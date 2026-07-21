"""
Direct unit tests for ``YamlServerSettings``, ``ServerSettingsManager``,
and the scheduler-specific ``*Server`` subclasses in
``chemsmart.settings.server``.
"""

import os

import pytest

from chemsmart.settings.server import (
    LSFServer,
    PBSServer,
    ServerSettingsManager,
    SGE_Server,
    SLURMServer,
    YamlServerSettings,
)


class TestYamlServerSettings:
    def test_from_yaml_loads_server_section(self, server_yaml_file):
        settings = YamlServerSettings.from_yaml(server_yaml_file)
        assert settings.name == server_yaml_file
        assert settings.num_cores == 64
        assert settings.mem_gb == 375

    def test_repr_and_str(self, server_yaml_file):
        settings = YamlServerSettings.from_yaml(server_yaml_file)
        assert repr(settings) == f"YamlServerSettings(name={server_yaml_file})"
        assert str(settings) == f"YamlServerSettings: {server_yaml_file}"

    def test_eq_and_hash_by_name(self):
        a = YamlServerSettings(name="x")
        b = YamlServerSettings(name="x", NUM_CORES=99)
        c = YamlServerSettings(name="y")
        assert a == b
        assert a != c
        assert hash(a) == hash("x")

    def test_call_returns_self(self):
        settings = YamlServerSettings(name="x")
        assert settings() is settings

    def test_register_returns_self_without_side_effects(self):
        """Unlike ``Server.register()`` (see BUGS_FOUND.md #7),
        ``YamlServerSettings.register()`` is overridden to be a
        trivial, safe no-op that just returns ``self``."""
        settings = YamlServerSettings(name="x")
        assert settings.register() is settings


class TestServerSettingsManager:
    def test_rejects_none_filename(self):
        with pytest.raises(ValueError, match="filename is not specified"):
            ServerSettingsManager(filename=None)

    def test_stores_absolute_path(self, server_yaml_file):
        manager = ServerSettingsManager(filename=server_yaml_file)
        assert manager.filename == os.path.abspath(server_yaml_file)

    def test_create_loads_yaml_server_settings(self, server_yaml_file):
        manager = ServerSettingsManager(filename=server_yaml_file)
        settings = manager.create()
        assert isinstance(settings, YamlServerSettings)
        assert settings.num_cores == 64

    def test_create_raises_for_missing_file(self, tmp_path):
        manager = ServerSettingsManager(filename=str(tmp_path / "nope.yaml"))
        with pytest.raises(FileNotFoundError):
            manager.create()


class TestSchedulerSpecificServers:
    def test_pbs_server_construction(self):
        server = PBSServer(NUM_CORES=8)
        assert server.name == "PBS"
        assert server.SCHEDULER_TYPE == "PBS"
        assert server.num_cores == 8

    def test_lsf_server_construction(self):
        server = LSFServer(NUM_CORES=4)
        assert server.name == "LSF"
        assert server.SCHEDULER_TYPE == "LSF"

    def test_sge_server_construction(self):
        server = SGE_Server(NUM_CORES=2)
        assert server.name == "SGE"
        assert server.SCHEDULER_TYPE == "SGE"


class TestSLURMServerBug:
    """Documents a real bug: ``SLURMServer.__init__`` calls
    ``super().__init__(filename=..., **kwargs)``, but
    ``YamlServerSettings.__init__`` (and ``Server.__init__``) take a
    positional/``name`` argument, not ``filename`` — unlike
    ``PBSServer``/``LSFServer``/``SGE_Server``, which correctly pass
    ``self.NAME`` positionally. Constructing ``SLURMServer()`` directly
    always raises ``TypeError``. In practice this path isn't hit by
    ``Server.from_scheduler_type()`` (which goes through
    ``from_servername`` → ``ServerSettingsManager`` instead of calling
    the class directly), but the class is still public API and broken.
    See BUGS_FOUND.md."""

    def test_construction_raises_type_error(self):
        with pytest.raises(TypeError, match="missing 1 required positional"):
            SLURMServer(NUM_CORES=8)
