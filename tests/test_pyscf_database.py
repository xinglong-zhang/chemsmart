"""Tests for PySCF results reaching a CHEMSMART database."""

import os

import h5py
import pytest
from click.testing import CliRunner

from chemsmart.cli.database.assemble import assemble
from chemsmart.database.assemble import PySCFAssembler, SingleFileAssembler
from chemsmart.database.database import Database
from chemsmart.database.utils import sha256_file
from chemsmart.io.pyscf.driver import APPLIED_SPEC_FIELDS
from chemsmart.io.pyscf.output import ROUTE_IDENTITY_EXCLUDED, PySCFOutput
from tests.helpers import (
    changed_value,
    copy_pyscf_output,
    replace_pyscf_spec_field,
)


@pytest.mark.parametrize("extra", [[], ["-p", "pyscf"]])
def test_assemble_finds_pyscf_outputs_in_a_directory(
    pyscf_water_sp_outfolder, tmp_path, extra
):
    copy_pyscf_output(pyscf_water_sp_outfolder, tmp_path / "calculations")
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
    original = copy_pyscf_output(pyscf_water_sp_outfolder, tmp_path / "a")
    changed = copy_pyscf_output(pyscf_water_sp_outfolder, tmp_path / "b")
    with h5py.File(tmp_path / "b" / "calc.h5", "a") as handle:
        del handle[f"spec/{field}"]
        handle["spec"].create_dataset(field, data=value)

    assert (
        PySCFAssembler(filename=original).assemble().record_id
        != PySCFAssembler(filename=changed).assemble().record_id
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        # The tolerance decides the point group, hence the rotational
        # symmetry number and the free energy.
        ("point_group_tol", 0.05),
        # A result storing fewer properties must not overwrite one storing
        # more: an optimisation without its Hessian has no thermochemistry,
        ("stages", ["scf", "opt"]),
        # and a run without forces has no gradient.
        ("forces", False),
        # A run that stopped at its step limit must not overwrite one that
        # had room to converge.
        ("opt_maxsteps", 5),
    ],
)
def test_settings_that_change_the_stored_results_change_the_record_id(
    pyscf_water_opt_distorted_outfolder, tmp_path, field, value
):
    original = copy_pyscf_output(
        pyscf_water_opt_distorted_outfolder, tmp_path / "a"
    )
    changed = copy_pyscf_output(
        pyscf_water_opt_distorted_outfolder, tmp_path / "b"
    )
    replace_pyscf_spec_field(changed, field, value)

    assert (
        PySCFAssembler(filename=original).assemble().record_id
        != PySCFAssembler(filename=changed).assemble().record_id
    )


@pytest.mark.parametrize("field", APPLIED_SPEC_FIELDS)
def test_only_identity_neutral_settings_keep_the_route(
    pyscf_water_opt_distorted_outfolder, tmp_path, field
):
    """Every applied setting reaches the identity unless it is excluded."""
    logfile = copy_pyscf_output(
        pyscf_water_opt_distorted_outfolder, tmp_path / "calc"
    )
    before = PySCFOutput(logfile)
    replace_pyscf_spec_field(
        logfile, field, changed_value(before.spec.get(field))
    )
    after = PySCFOutput(logfile)

    unchanged = after.route_string == before.route_string
    assert unchanged == (field in ROUTE_IDENTITY_EXCLUDED)


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


def test_assemble_keeps_converged_runs_and_skips_failed_ones(
    pyscf_water_sp_outfolder,
    pyscf_water_scf_failed_outfolder,
    pyscf_water_smd_failed_outfolder,
    tmp_path,
):
    root = tmp_path / "calculations"
    copy_pyscf_output(pyscf_water_sp_outfolder, root / "good")
    copy_pyscf_output(pyscf_water_scf_failed_outfolder, root / "scf_failed")
    copy_pyscf_output(pyscf_water_smd_failed_outfolder, root / "no_geometry")
    output = tmp_path / "pyscf.db"
    result = CliRunner().invoke(
        assemble, ["-d", str(root), "-p", "pyscf", "-o", str(output)], obj={}
    )
    assert result.exit_code == 0, (result.output, result.exception)
    assert Database(str(output)).count_records() == 1


def test_failed_run_with_a_geometry_can_be_included(
    pyscf_water_scf_failed_outfolder,
):
    label = os.path.basename(pyscf_water_scf_failed_outfolder)
    logfile = os.path.join(pyscf_water_scf_failed_outfolder, f"{label}.out")
    record = PySCFAssembler(filename=logfile, include_failed=True).assemble()
    assert record.provenance["normal_termination"] is False
    assert record.molecules[0]["is_optimized_structure"] is False


def test_failed_run_without_a_geometry_assembles_to_nothing(
    pyscf_water_smd_failed_outfolder,
):
    label = os.path.basename(pyscf_water_smd_failed_outfolder)
    logfile = os.path.join(pyscf_water_smd_failed_outfolder, f"{label}.out")
    assembler = SingleFileAssembler(logfile, include_failed=True)
    assert assembler.assemble_data is None


def test_open_shell_record_keeps_both_spin_channels(
    pyscf_oh_radical_outfolder,
):
    label = os.path.basename(pyscf_oh_radical_outfolder)
    logfile = os.path.join(pyscf_oh_radical_outfolder, f"{label}.out")
    record = SingleFileAssembler(logfile).assemble_data
    assert record.results["alpha_homo_energy"] is not None
    assert record.results["beta_homo_energy"] is not None
    # The same "restricted"/"unrestricted" meaning as Gaussian and ORCA.
    assert record.meta["spin"] == "unrestricted"


def _assemble_benzene(folder, tmp_path):
    output = tmp_path / "benzene.db"
    result = CliRunner().invoke(
        assemble, ["-d", folder, "-p", "pyscf", "-o", str(output)], obj={}
    )
    assert result.exit_code == 0, (result.output, result.exception)
    return Database(str(output)).get_record(record_index=1)


def test_benzene_record_has_its_energy_and_frequencies(
    pyscf_benzene_outfolder, tmp_path
):
    record = _assemble_benzene(pyscf_benzene_outfolder, tmp_path)
    molecule = record["molecules"][0]
    assert record["meta"]["jobtype"] == "opt"
    assert record["results"]["total_energy"] == pytest.approx(
        -232.08478438, abs=1e-7
    )
    assert record["results"]["zero_point_energy"] == pytest.approx(
        0.100499, abs=1e-5
    )
    assert molecule["is_optimized_structure"] is True
    assert molecule["num_vibrational_modes"] == 30


def test_benzene_record_has_its_true_symmetry(
    pyscf_benzene_outfolder, tmp_path
):
    record = _assemble_benzene(pyscf_benzene_outfolder, tmp_path)
    molecule = record["molecules"][0]
    assert molecule["point_group"] == "D6H"
    assert molecule["rotational_symmetry_number"] == 12
    assert record["results"]["gibbs_free_energy"] == pytest.approx(
        -232.009398, abs=1e-5
    )
