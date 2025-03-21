import pytest
import logging
from unittest.mock import patch, MagicMock
from click.testing import CliRunner
from chemsmart.cli.gaussian.gaussian import gaussian

@pytest.fixture
def runner():
    return CliRunner()

@pytest.fixture
def mock_ctx():
    return {
        "settings": {"some_setting": "value"},
        "job_settings": MagicMock(),
        "keywords": MagicMock(),
        "molecules": ["TestMolecule"],
        "label": "test_label",
    }

@patch("chemsmart.jobs.gaussian.settings.GaussianQMMMJobSettings")
@patch("chemsmart.cli.gaussian.qmmm.get_setting_from_jobtype_for_gaussian")
@patch("chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob")
def test_qmmm_command(mock_job, mock_get_settings, mock_gaussian_settings, runner, mock_ctx):
    """Test the QM/MM CLI command execution"""

    # Mock external function calls
    mock_get_settings.return_value = MagicMock()
    mock_gaussian_settings.return_value = MagicMock()
    mock_job.return_value = MagicMock()

    # Define CLI parameters
    params = [
        "--functional-high", "B3LYP",
        "--basis-high", "6-31G(d)",
        "--force-field-high", "UFF",
        "--functional-medium", "PBE",
        "--basis-medium", "def2-SVP",
        "--force-field-medium", "GAFF",
        "--functional-low", "HF",
        "--basis-low", "STO-3G",
        "--force-field-low", "CHARMM",
        "--high-level-charges", "0",
        "--high-level-multiplicity", "1",
        "--low-level-charge", "0",
        "--high-level-atoms", "1,2,3",
        "--medium-level-atoms", "4,5,6",
        "--low-level-atoms", "7,8,9",
        "--bonded-atom", "1-4",
        "--scale-factor1", "1.0",
        "--scale-factor2", "1.0",
        "--scale-factor3", "1.0",
        "--num-atoms", "10",
        # "--jobtype", "optimization",
    ]

    # Run the CLI command
    with patch("chemsmart.cli.gaussian.qmmm.click.get_current_context") as mock_get_ctx:
        mock_get_ctx.return_value.obj = mock_ctx
        result = runner.invoke(gaussian, ["qmmm"] + params)
        print(result)

    # # Assertions
    # assert result.exit_code == 0, f"Command failed with output: {result.output}"
    # mock_get_settings.assert_called_once()
    # mock_job.assert_called_once()
    #
    # # Check if logging contains expected output
    # assert "ONIOM calculation of molecule" in result.output
    # assert "Running QM/MM job with settings" in result.output

    # mock_job.assert_called_once()

    # Get the actual arguments passed to GaussianQMMMJob
    _, kwargs = mock_job.call_args

    # Assertions to check if parameters were correctly assigned
    assert kwargs["settings"].functional_high == "B3LYP"
    # assert kwargs["settings"].basis_high == "6-31G"
    # assert kwargs["settings"].force_field_high == "UFF"
    # assert kwargs["settings"].functional_medium == "PBE"
    # assert kwargs["settings"].basis_medium == "DZVP"
    # assert kwargs["settings"].force_field_medium == "GAFF"
    # assert kwargs["settings"].functional_low == "PM6"
    # assert kwargs["settings"].basis_low == "STO-3G"
    # assert kwargs["settings"].force_field_low == "MMFF94"
    # assert kwargs["settings"].high_level_charges == "0"
    # assert kwargs["settings"].high_level_multiplicity == "1"
    # assert kwargs["settings"].low_level_charge == "-1"
    # assert kwargs["settings"].high_level_atoms == "1,2,3"
    # assert kwargs["settings"].medium_level_atoms == "4,5,6"
    # assert kwargs["settings"].low_level_atoms == "7,8,9"
    # assert kwargs["settings"].bonded_atom == "2,4"
    # assert kwargs["settings"].scale_factor1 == "1.2"
    # assert kwargs["settings"].scale_factor2 == "1.3"
    # assert kwargs["settings"].scale_factor3 == "1.4"
    # assert kwargs["settings"].num_atoms == "9"
