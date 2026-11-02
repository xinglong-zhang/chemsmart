"""Tests for PySCF results reaching a CHEMSMART database."""

import os
import shutil

import h5py
import pytest
from click.testing import CliRunner

from chemsmart.cli.database.assemble import assemble
from chemsmart.database.assemble import PySCFAssembler
from chemsmart.database.database import Database
from chemsmart.database.utils import sha256_file


def _copy_output(source_folder, destination):
    """Copy a committed PySCF output, returning the path of its log."""
    os.makedirs(destination)
    label = os.path.basename(source_folder)
    for suffix in (".out", ".h5"):
        shutil.copy(
            os.path.join(source_folder, label + suffix),
            os.path.join(destination, "calc" + suffix),
        )
    return os.path.join(destination, "calc.out")


@pytest.mark.parametrize("extra", [[], ["-p", "pyscf"]])
def test_assemble_finds_pyscf_outputs_in_a_directory(
    pyscf_water_sp_outfolder, tmp_path, extra
):
    _copy_output(pyscf_water_sp_outfolder, tmp_path / "calculations")
    output = tmp_path / "pyscf.db"
    result = CliRunner().invoke(
        assemble,
        ["-d", str(tmp_path / "calculations"), "-o", str(output), *extra],
        obj={},
    )
    assert result.exit_code == 0, (result.output, result.exception)
    record = Database(str(output)).get_record(record_index=1)
    assert record["provenance"]["program"] == "PySCF"


@pytest.mark.parametrize(
    ("field", "value"),
    [("scf_tol", 1e-3), ("aux_basis", "def2-universal-jkfit")],
)
def test_settings_that_change_results_change_the_record_id(
    pyscf_water_sp_outfolder, tmp_path, field, value
):
    original = _copy_output(pyscf_water_sp_outfolder, tmp_path / "a")
    changed = _copy_output(pyscf_water_sp_outfolder, tmp_path / "b")
    with h5py.File(tmp_path / "b" / "calc.h5", "a") as handle:
        del handle[f"spec/{field}"]
        handle["spec"].create_dataset(field, data=value)

    assert (
        PySCFAssembler(filename=original).assemble().record_id
        != PySCFAssembler(filename=changed).assemble().record_id
    )


def test_record_identifies_the_results_file(pyscf_water_opt_outfolder):
    label = os.path.basename(pyscf_water_opt_outfolder)
    logfile = os.path.join(pyscf_water_opt_outfolder, f"{label}.out")
    resultsfile = os.path.join(pyscf_water_opt_outfolder, f"{label}.h5")
    record = PySCFAssembler(filename=logfile).assemble()
    assert record.provenance["source_file_hash"] == sha256_file(resultsfile)
    molecule = record.molecules[0]
    assert molecule["point_group"] == "C2V"
    assert molecule["dipole_moment_magnitude"] is not None
    assert molecule["is_optimized_structure"] is True
