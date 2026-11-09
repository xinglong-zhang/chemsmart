from types import SimpleNamespace

import click
import pytest

from chemsmart.cli import plugins


@pytest.fixture(autouse=True)
def clear_cache():
    plugins.discover.cache_clear()
    yield
    plugins.discover.cache_clear()


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

    monkeypatch.setattr(plugins.importlib.metadata, "entry_points", lambda: Entries())
    assert [c.name for c in plugins.discover(plugins.JOB_COMMANDS)] == ["a", "z"]
    assert "missing dependency" in caplog.text
    with pytest.raises(plugins.PluginLoadError, match="missing dependency"):
        plugins.discover(plugins.JOB_COMMANDS, strict=True)


def test_duplicate_plugins_are_all_skipped(monkeypatch):
    entries = [SimpleNamespace(name="same", value=value, dist=None)
               for value in ("a", "b")]
    monkeypatch.setattr(
        plugins.importlib.metadata, "entry_points",
        lambda: SimpleNamespace(select=lambda **kwargs: entries),
    )
    assert plugins.discover(plugins.JOB_COMMANDS) == ()
    with pytest.raises(plugins.PluginLoadError, match="multiple entry points"):
        plugins.discover(plugins.JOB_COMMANDS, strict=True)


def test_reserved_plugin_is_not_loaded(monkeypatch):
    loaded = []
    entry = SimpleNamespace(
        name="run", value="plugin:run", dist=None,
        load=lambda: loaded.append(True),
    )
    monkeypatch.setattr(
        plugins.importlib.metadata, "entry_points",
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

    entry = SimpleNamespace(name="example", value="example:command", dist=None, load=load)
    monkeypatch.setattr(
        plugins.importlib.metadata, "entry_points",
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
        Executable, "subclasses",
        lambda: [SimpleNamespace(PROGRAM="gaussian"),
                 SimpleNamespace(PROGRAM="GAUSSIAN")],
    )
    submitter = SimpleNamespace(
        job=SimpleNamespace(PROGRAM="Gaussian"),
        server=SimpleNamespace(name="cluster"),
    )
    with pytest.raises(ValueError, match="Multiple executables"):
        _executable_class_for_program("gaussian")
    with pytest.raises(ValueError, match="Multiple executables"):
        Submitter.executable.fget(submitter)
