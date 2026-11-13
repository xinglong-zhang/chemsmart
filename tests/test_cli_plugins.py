import importlib.metadata
import sys
from types import SimpleNamespace

import click
import pytest
from click.testing import CliRunner

from chemsmart.cli import plugins
from chemsmart.utils.cli import CtxObjArguments, MyGroup
from chemsmart.utils.mixins import RegistryMixin


@pytest.fixture(autouse=True)
def clear_cache():
    plugins.discover.cache_clear()
    yield
    plugins.discover.cache_clear()


@pytest.fixture
def shared_plugin(tmp_path, monkeypatch):
    """Real modules exercise Python's package cache and shared class identity."""
    name = "fixture_shared_plugin"
    package = tmp_path / name
    package.mkdir()
    (package / "__init__.py").write_text("")
    (package / "shared.py").write_text(
        "from chemsmart.utils.mixins import RegistryMixin\n"
        "class Shared(RegistryMixin):\n"
        "    pass\n"
    )
    (package / "good.py").write_text(
        "from . import shared\n"
        "import click\n"
        "good = click.Command('good')\n"
        "good.chemsmart_classes = (shared.Shared,)\n"
    )
    before = list(RegistryMixin._REGISTRY)
    monkeypatch.syspath_prepend(str(tmp_path))

    def discover(names, strict=False):
        entries = importlib.metadata.EntryPoints(
            importlib.metadata.EntryPoint(
                name=entry,
                value=f"{name}.{entry}:{entry}",
                group=plugins.JOB_COMMANDS,
            )
            for entry in names
        )
        monkeypatch.setattr(
            plugins.importlib.metadata, "entry_points", lambda: entries
        )
        plugins.discover.cache_clear()
        return plugins.discover(plugins.JOB_COMMANDS, strict=strict)

    yield name, package, discover
    RegistryMixin._REGISTRY[:] = before
    for module in list(sys.modules):
        if module == name or module.startswith(name + "."):
            sys.modules.pop(module)


@pytest.mark.parametrize("failure", ["import", "invalid", "declaration"])
@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("good_first", [False, True])
def test_shared_dependency_survives_failed_plugin(
    shared_plugin, failure, strict, good_first
):
    name, package, discover = shared_plugin
    bad_source = (
        "from . import shared\n"
        "from chemsmart.utils.mixins import RegistryMixin\n"
        "import click\n"
        "class Broken(RegistryMixin):\n"
        "    pass\n"
    )
    if failure == "import":
        bad_source += "raise ImportError('missing optional backend')\n"
    elif failure == "invalid":
        bad_source += "bad = click.Command('wrong-name')\n"
    else:
        bad_source += (
            "bad = click.Command('bad')\n"
            "bad.chemsmart_classes = (shared.Shared, object())\n"
        )
    (package / "bad.py").write_text(bad_source)
    if good_first:
        assert discover(["good"])[0].name == "good"
    before_bad = list(RegistryMixin._REGISTRY)
    if strict:
        with pytest.raises(plugins.PluginLoadError):
            discover(["bad"], strict=True)
    else:
        assert discover(["bad"]) == ()
    assert RegistryMixin._REGISTRY == before_bad
    cached_shared = sys.modules[f"{name}.shared"]
    assert getattr(sys.modules[name], "shared") is cached_shared

    # A successful plugin can activate a class imported by the failed plugin.
    for mode in (False, True, False):
        assert discover(["good"], strict=mode)[0].name == "good"
        assert sys.modules[f"{name}.shared"] is cached_shared
        assert RegistryMixin._REGISTRY.count(cached_shared.Shared) == 1


def test_deferred_registration_keeps_core_imports_and_resets_on_error():
    before = list(RegistryMixin._REGISTRY)
    try:
        with pytest.raises(SystemExit):
            with RegistryMixin.defer_plugin_registration():
                core = type(
                    "LazyCore",
                    (RegistryMixin,),
                    {"__module__": "chemsmart.test_lazy"},
                )
                external = type("Plugin", (RegistryMixin,), {})
                raise SystemExit(1)
        assert core in RegistryMixin._REGISTRY
        assert external not in RegistryMixin._REGISTRY
        after = type("NormalRegistration", (RegistryMixin,), {})
        assert after in RegistryMixin._REGISTRY
    finally:
        RegistryMixin._REGISTRY[:] = before


def test_existing_command_cannot_be_overwritten(monkeypatch):
    group = click.Group()
    original = click.Command("run")
    group.add_command(original)
    monkeypatch.setattr(
        plugins, "discover", lambda *args, **kwargs: (click.Command("run"),)
    )
    plugins.add_commands(plugins.COMMANDS, group, ())
    assert group.commands["run"] is original
    with pytest.raises(plugins.PluginLoadError, match="run"):
        plugins.add_commands(plugins.COMMANDS, group, (), strict=True)


def test_discovery_is_sorted_and_skips_broken(monkeypatch, caplog):
    class Entry:
        def __init__(self, name, broken=False):
            self.name = name
            self.value = name
            self.dist = SimpleNamespace(name="example")
            self.broken = broken

        def load(self):
            if self.broken:
                raise ImportError("missing dependency")
            return click.Command(self.name)

    class Entries:
        def select(self, **kwargs):
            return [Entry("z"), Entry("bad", True), Entry("a")]

    monkeypatch.setattr(
        plugins.importlib.metadata, "entry_points", lambda: Entries()
    )
    assert [c.name for c in plugins.discover(plugins.JOB_COMMANDS)] == [
        "a",
        "z",
    ]
    assert "missing dependency" in caplog.text
    with pytest.raises(plugins.PluginLoadError, match="missing dependency"):
        plugins.discover(plugins.JOB_COMMANDS, strict=True)


def test_duplicate_plugins_are_all_skipped(monkeypatch):
    entries = [
        SimpleNamespace(name="same", value=value, dist=None)
        for value in ("a", "b")
    ]
    monkeypatch.setattr(
        plugins.importlib.metadata,
        "entry_points",
        lambda: SimpleNamespace(select=lambda **kwargs: entries),
    )
    assert plugins.discover(plugins.JOB_COMMANDS) == ()
    with pytest.raises(plugins.PluginLoadError, match="multiple entry points"):
        plugins.discover(plugins.JOB_COMMANDS, strict=True)


def test_job_plugin_accepts_plain_click_command(monkeypatch):
    entry = SimpleNamespace(
        name="example",
        value="example:command",
        dist=None,
        load=lambda: click.Command("example"),
    )
    monkeypatch.setattr(
        plugins.importlib.metadata,
        "entry_points",
        lambda: SimpleNamespace(select=lambda **kwargs: [entry]),
    )
    assert plugins.discover(plugins.JOB_COMMANDS)[0].name == "example"
    assert (
        plugins.discover(plugins.JOB_COMMANDS, strict=True)[0].name
        == "example"
    )


def test_host_records_and_replays_nested_job_parameters(tmp_path):
    calls = []

    @click.group(cls=MyGroup)
    def sub():
        pass

    @click.group("example")
    @click.option("--project", "-p", type=click.Path(resolve_path=True))
    def example(project):
        calls.append(project)

    @example.group("example")
    def nested():
        pass

    @nested.command()
    @click.argument("filename")
    @click.option("--epochs", "steps", type=int)
    @click.option("--tag", multiple=True)
    @click.option("--pair", nargs=2, type=int)
    @click.option("--enabled/--disabled", default=True)
    @click.option("--fast", "mode", flag_value="fast")
    @click.option("--slow", "mode", flag_value="slow")
    @click.option("--truth", type=bool)
    @click.option("-v", count=True)
    def train(**values):
        calls.append(values)

    @sub.result_callback()
    @click.pass_context
    def arguments(ctx, result):
        return CtxObjArguments(
            ctx.obj["subcommand"], "sub"
        ).reconstruct_command_line()

    # run/sub share the command: repeated adaptation must not add wrappers.
    plugins._record_job_command(example)
    plugins._record_job_command(example)
    sub.add_command(example)
    argv = [
        "example",
        "-p",
        str(tmp_path / "project.yaml"),
        "example",
        "train",
        "--epochs",
        "3",
        "--tag",
        "",
        "--tag",
        "water",
        "--pair",
        "1",
        "2",
        "--disabled",
        "--slow",
        "--truth",
        "false",
        "-vv",
        "--",
        "-input.xyz",
    ]
    result = CliRunner().invoke(sub, argv, standalone_mode=False)
    assert result.exit_code == 0, result.output
    original = list(calls)
    calls.clear()
    replay = result.return_value
    assert replay[:2] == ["sub", "example"]
    assert "--epochs" in replay and "--steps" not in replay
    result = CliRunner().invoke(sub, replay[1:], standalone_mode=False)
    assert result.exit_code == 0, result.output
    assert calls == original


def test_host_records_lazy_plugin_commands():
    class LazyGroup(click.Group):
        def get_command(self, ctx, name):
            return leaf if name == "dynamic" else None

    @click.command("dynamic")
    @click.option("--value")
    @click.pass_context
    def leaf(ctx, value):
        return CtxObjArguments(
            ctx.obj["subcommand"]
        ).reconstruct_command_line()

    command = LazyGroup("lazy")
    plugins._record_job_command(command)
    result = CliRunner().invoke(
        command, ["dynamic", "--value", "x"], standalone_mode=False
    )
    assert result.exit_code == 0, result.output
    assert result.return_value == ["lazy", "dynamic", "--value", "x"]


@pytest.mark.parametrize("amount", [None, "4"])
@pytest.mark.parametrize("disable", [False, True])
def test_plugin_replay_preserves_callbacks_defaults_and_click_types(
    tmp_path, monkeypatch, amount, disable
):
    from datetime import datetime

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("PLUGIN_AMOUNT", raising=False)
    if amount is not None:
        monkeypatch.setenv("PLUGIN_AMOUNT", amount)
    (tmp_path / "input.txt").write_text("contents")
    calls = []

    def transform(ctx, param, value):
        ctx.params["derived"] = "from callback"
        return value * 2

    @click.command()
    @click.option(
        "--amount",
        type=int,
        callback=transform,
        default=lambda: 3,
        envvar="PLUGIN_AMOUNT",
    )
    @click.option("--date", type=click.DateTime(["%Y-%m-%d"]))
    @click.option(
        "--other-date",
        type=click.DateTime(["%Y-%m-%d"]),
        default=datetime(2026, 1, 1),
    )
    @click.option("--input", type=click.File())
    @click.option("--disable", is_flag=True, default=True)
    @click.option("--tag", is_flag=True, multiple=True)
    @click.option("--red", "colors", flag_value="red", multiple=True)
    @click.option("--blue", "colors", flag_value="blue", multiple=True)
    @click.pass_context
    def example(ctx, **values):
        values["input"] = values["input"].read()
        calls.append(values)
        return CtxObjArguments(
            ctx.obj["subcommand"]
        ).reconstruct_command_line()

    plugins._record_job_command(example)
    argv = [
        "--date",
        "2026-11-13",
        "--input",
        "input.txt",
        "--tag",
        "--tag",
        "--blue",
        "--red",
        "--blue",
    ]
    if disable:
        argv.append("--disable")
    result = CliRunner().invoke(example, argv, standalone_mode=False)
    assert result.exit_code == 0, result.output
    replay = result.return_value
    assert str(tmp_path / "input.txt") in replay
    result = CliRunner().invoke(example, replay[1:], standalone_mode=False)
    assert result.exit_code == 0, result.output
    assert calls[0] == calls[1]
    assert calls[0]["amount"] == (8 if amount else 6)
    assert calls[0]["disable"] is not disable
    assert calls[0]["tag"] == (True, True)
    assert calls[0]["colors"] == ("blue", "red", "blue")


def test_reserved_plugin_is_not_loaded(monkeypatch):
    loaded = []
    entry = SimpleNamespace(
        name="run",
        value="plugin:run",
        dist=None,
        load=lambda: loaded.append(True),
    )
    monkeypatch.setattr(
        plugins.importlib.metadata,
        "entry_points",
        lambda: SimpleNamespace(select=lambda **kwargs: [entry]),
    )
    group = click.Group()
    group.add_command(click.Command("run"))
    plugins.add_commands(plugins.COMMANDS, group, ())
    assert loaded == []
    with pytest.raises(plugins.PluginLoadError, match="conflicts"):
        plugins.add_commands(plugins.COMMANDS, group, (), strict=True)
    assert loaded == []


@pytest.mark.parametrize("failure", ["import", "invalid"])
def test_failed_plugin_rolls_back_registered_classes(monkeypatch, failure):
    from chemsmart.settings.executable import Executable
    from chemsmart.utils.mixins import RegistryMixin

    before = list(RegistryMixin._REGISTRY)

    def load():
        class PluginExecutable(Executable):
            PROGRAM = "rollback_test"

        if failure == "import":
            raise ImportError("missing dependency")
        return click.Command("wrong-name")

    entry = SimpleNamespace(
        name="example", value="example:command", dist=None, load=load
    )
    monkeypatch.setattr(
        plugins.importlib.metadata,
        "entry_points",
        lambda: SimpleNamespace(select=lambda **kwargs: [entry]),
    )
    assert plugins.discover(plugins.JOB_COMMANDS) == ()
    assert RegistryMixin._REGISTRY == before
    with pytest.raises(plugins.PluginLoadError):
        plugins.discover(plugins.JOB_COMMANDS, strict=True)
    assert RegistryMixin._REGISTRY == before


def test_submitter_resolves_plugin_executable(monkeypatch):
    from chemsmart.settings.executable import Executable
    from chemsmart.settings.submitters import Submitter

    expected = object()
    plugin = SimpleNamespace(
        PROGRAM="example", from_servername=lambda name: (name, expected)
    )
    monkeypatch.setattr(Executable, "subclasses", lambda: [plugin])
    submitter = SimpleNamespace(
        job=SimpleNamespace(PROGRAM="EXAMPLE"),
        server=SimpleNamespace(name="cluster"),
    )
    assert Submitter.executable.fget(submitter) == ("cluster", expected)


def test_duplicate_executable_rejected_by_runner_and_submitter(monkeypatch):
    from chemsmart.jobs.runner import _executable_class_for_program
    from chemsmart.settings.executable import Executable
    from chemsmart.settings.submitters import Submitter

    monkeypatch.setattr(
        Executable,
        "subclasses",
        lambda: [
            SimpleNamespace(PROGRAM="gaussian"),
            SimpleNamespace(PROGRAM="GAUSSIAN"),
        ],
    )
    submitter = SimpleNamespace(
        job=SimpleNamespace(PROGRAM="Gaussian"),
        server=SimpleNamespace(name="cluster"),
    )
    with pytest.raises(ValueError, match="Multiple executables"):
        _executable_class_for_program("gaussian")
    with pytest.raises(ValueError, match="Multiple executables"):
        Submitter.executable.fget(submitter)
