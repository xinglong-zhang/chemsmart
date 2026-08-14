"""
Tests for the ``nciplot`` CLI command.

Each test uses :class:`click.testing.CliRunner` to invoke the ``nciplot``
command and :mod:`unittest.mock` to intercept the job constructor so that
the call kwargs can be inspected without running an actual NCIPLOT job.
"""

import sys
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemsmart.cli.nciplot.nciplot import nciplot

# ``chemsmart.cli.nciplot.__init__`` does ``from .nciplot import nciplot``,
# which shadows the ``nciplot`` submodule attribute on the package with the
# click Command object. Fetch the real submodule from ``sys.modules``
# directly so we can patch ``NCIPLOTJob`` on it.
nciplot_module = sys.modules["chemsmart.cli.nciplot.nciplot"]


def run_nciplot_and_capture_kwargs(cli_args):
    runner = CliRunner()
    with patch.object(nciplot_module, "NCIPLOTJob") as mock_job_cls:
        mock_job_cls.return_value = MagicMock()
        result = runner.invoke(
            nciplot,
            cli_args,
            catch_exceptions=False,
        )
        call_kwargs = None
        if mock_job_cls.call_args is not None:
            call_kwargs = mock_job_cls.call_args.kwargs
    return result, call_kwargs


class TestNciplotCLI:
    def test_requires_file_or_pubchem(self):
        with pytest.raises(ValueError, match="Must provide at least one"):
            run_nciplot_and_capture_kwargs([])

    def test_single_file_job_creation_promolecular_label(
        self, single_molecule_xyz_file
    ):
        result, kwargs = run_nciplot_and_capture_kwargs(
            ["-f", single_molecule_xyz_file]
        )
        assert result.exit_code == 0, result.output
        assert kwargs is not None, "NCIPLOTJob was never instantiated"
        assert kwargs["label"].endswith("_promolecular")

    def test_custom_label_preserved(self, single_molecule_xyz_file):
        result, kwargs = run_nciplot_and_capture_kwargs(
            ["-f", single_molecule_xyz_file, "-l", "my_label"]
        )
        assert result.exit_code == 0, result.output
        assert kwargs["label"] == "my_label_promolecular"

    def test_multiple_files_label_joined(self, single_molecule_xyz_file):
        result, kwargs = run_nciplot_and_capture_kwargs(
            [
                "-f",
                single_molecule_xyz_file,
                "-f",
                single_molecule_xyz_file,
            ]
        )
        assert result.exit_code == 0, result.output
        assert "_and_" in kwargs["label"]

    def test_fragments_option_parsed(self, single_molecule_xyz_file):
        result, kwargs = run_nciplot_and_capture_kwargs(
            [
                "-f",
                single_molecule_xyz_file,
                "--fragments",
                "{1: [1, 2, 3], 2: [4, 5, 6]}",
            ]
        )
        assert result.exit_code == 0, result.output
        assert kwargs["settings"].fragments == {1: [1, 2, 3], 2: [4, 5, 6]}

    def test_invalid_fragments_option_raises(self, single_molecule_xyz_file):
        result, _ = run_nciplot_and_capture_kwargs(
            [
                "-f",
                single_molecule_xyz_file,
                "--fragments",
                "not-a-literal(",
            ]
        )
        assert result.exit_code != 0

    def test_pubchem_requires_label(self):
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_pubchem",
            return_value=MagicMock(),
        ):
            with pytest.raises(
                AssertionError, match="Label for file is required"
            ):
                run_nciplot_and_capture_kwargs(["-P", "water"])

    def test_pubchem_and_filenames_mutually_exclusive(
        self, single_molecule_xyz_file
    ):
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_pubchem",
            return_value=MagicMock(),
        ):
            with pytest.raises(AssertionError, match="Cannot provide both"):
                run_nciplot_and_capture_kwargs(
                    [
                        "-f",
                        single_molecule_xyz_file,
                        "-P",
                        "water",
                        "-l",
                        "water_label",
                    ]
                )
