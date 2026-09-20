"""Program CLIs do not read chain YAML aliases."""

from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.user import CHEMSMARTUserSettings


@pytest.fixture()
def isolated_config_dir(tmp_path, monkeypatch):
    config_dir = tmp_path / "chemsmart_config"
    config_dir.mkdir()
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_dir))
    monkeypatch.setattr(
        CHEMSMARTUserSettings,
        "USER_CONFIG_DIR",
        str(config_dir),
    )
    return config_dir


def test_gaussian_dash_p_combined_does_not_load_chain_yaml(
    isolated_config_dir, single_molecule_xyz_file
):
    chain_dir = isolated_config_dir / "chain"
    chain_dir.mkdir()
    (chain_dir / "combined.yaml").write_text("gaussian: gas_solv\n")

    with patch(
        "chemsmart.settings.gaussian.GaussianProjectSettings.from_project",
        wraps=GaussianProjectSettings.from_project,
    ) as from_project:
        result = CliRunner().invoke(
            gaussian,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
            obj={"jobrunner": MagicMock()},
        )
    assert result.exit_code != 0
    from_project.assert_called_with("combined")
