"""Shared paths and server/user-config fixtures."""

import os
from pathlib import Path

import pytest


@pytest.fixture()
def test_data_directory():
    """Absolute path to ``tests/data``."""
    tests_directory = os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))
    )
    return os.path.join(tests_directory, "data")


@pytest.fixture()
def utils_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "YAMLTests")


@pytest.fixture()
def server_yaml_file(utils_test_directory):
    return os.path.join(utils_test_directory, "server.yaml")


@pytest.fixture()
def gaussian_project_config_dir(tmp_path):
    """Minimal Gaussian project config under a temporary CHEMSMART config root."""
    config_root = tmp_path / "chemsmart_cfg"
    gaussian_cfg = config_root / "gaussian"
    gaussian_cfg.mkdir(parents=True)
    (gaussian_cfg / "test.yaml").write_text(
        "gas:\n  functional: B3LYP\n  basis: def2-SVP\n"
        "solv:\n  functional: B3LYP\n  basis: def2-SVP\n"
        "  solvent_model: smd\n  solvent_id: water\n"
    )
    return config_root


_LEGACY_SERVER_YAML = """\
SERVER:
    SCHEDULER: PBS
    QUEUE_NAME: normal
    NUM_HOURS: 24
    MEM_GB: 375
    NUM_CORES: 64
    NUM_GPUS: 0
    NUM_THREADS: 64
    SUBMIT_COMMAND: qsub
    SCRATCH_DIR: null
    USE_HOSTS: false
GAUSSIAN:
    EXEFOLDER: ~/programs/g16
    LOCAL_RUN: True
ORCA:
    EXEFOLDER: ~/programs/orca_6_0_0
    LOCAL_RUN: False
"""


@pytest.fixture()
def legacy_server_yaml(tmp_path):
    """Path to a temporary server YAML that has no XTB or CREST block.

    Simulates a server configuration written before xTB support was added to
    CHEMSMART, so that tests can verify the graceful fallback behaviour.
    """
    yaml_path = tmp_path / "legacy_server.yaml"
    yaml_path.write_text(_LEGACY_SERVER_YAML)
    return str(yaml_path)


@pytest.fixture()
def chemsmart_templates_config(mocker):
    """
    Point USER_CONFIG_DIR to the local templates directory.
    This avoids creating a mock directory
    and instead uses the provided templates.
    """
    # Locate templates: chemsmart/settings/templates/.chemsmart
    package_root = Path(__file__).resolve().parents[2]
    template_dir = (
        package_root / "chemsmart" / "settings" / "templates" / ".chemsmart"
    )

    if not template_dir.exists():
        raise FileNotFoundError(
            f"Template directory not found: {template_dir}"
        )

    # Patch the Class attribute
    mocker.patch(
        "chemsmart.settings.user.CHEMSMARTUserSettings.USER_CONFIG_DIR",
        str(template_dir),
    )

    # Patch the global instance in runner.py
    from chemsmart.settings.user import CHEMSMARTUserSettings

    new_settings = CHEMSMARTUserSettings()
    mocker.patch("chemsmart.jobs.runner.user_settings", new_settings)
    # Patch other module-level user_settings singletons used by the CLI path
    mocker.patch("chemsmart.settings.server.user_settings", new_settings)
    mocker.patch("chemsmart.settings.executable.user_settings", new_settings)

    return template_dir
