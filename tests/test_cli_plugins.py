import importlib.metadata
import sys
from types import SimpleNamespace

import click
import pytest

from chemsmart.cli import plugins
from chemsmart.utils.cli import MyCommand
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
        "from chemsmart.utils.cli import MyCommand\n"
        "good = MyCommand('good')\n"
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
        "from chemsmart.utils.cli import MyCommand\n"
        "class Broken(RegistryMixin):\n"
        "    pass\n"
    )
    if failure == "import":
        bad_source += "raise ImportError('missing optional backend')\n"
    elif failure == "invalid":
        bad_source += "bad = MyCommand('wrong-name')\n"
    else:
        bad_source += (
            "bad = MyCommand('bad')\n"
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
            return MyCommand(self.name)

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


def test_job_plugin_requires_recording_command(monkeypatch):
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
    assert plugins.discover(plugins.JOB_COMMANDS) == ()
    with pytest.raises(plugins.PluginLoadError, match="MyCommand"):
        plugins.discover(plugins.JOB_COMMANDS, strict=True)


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
