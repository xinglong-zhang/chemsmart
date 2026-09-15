import io
import os

import pytest

from chemsmart.jobs.xtb.opt import XTBOptJob
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.jobs.xtb.writer import XTBInputWriter


@pytest.mark.usefixtures("temporary_working_dir")
class TestXTBInputWriter:
    def test_write_distance_angle_dihedral_constraints(
        self, tmpdir, methane_molecule, xtb_jobrunner_no_scratch
    ):
        settings = XTBJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="opt",
            constraints=[[1, 2], [2, 1, 3], [2, 1, 3, 4]],
            force_constant=0.5,
        )
        job = XTBOptJob(
            molecule=methane_molecule,
            settings=settings,
            label="xtb_constrain",
            jobrunner=xtb_jobrunner_no_scratch,
        )
        writer = XTBInputWriter(job=job)
        writer.write(target_directory=tmpdir)

        inp_file = os.path.join(tmpdir, "xtb_constrain.inp")
        assert os.path.isfile(inp_file)
        content = open(inp_file).read()

        dist = methane_molecule.get_distance(1, 2)
        angle = methane_molecule.get_angle(2, 1, 3)
        dihedral = methane_molecule.get_dihedral(2, 1, 3, 4)

        assert content.startswith("$constrain\n")
        assert "   force constant=0.5\n" in content
        assert f"   distance: 1, 2, {dist:.4f}\n" in content
        assert f"   angle: 2, 1, 3, {angle:.4f}\n" in content
        assert f"   dihedral: 2, 1, 3, 4, {dihedral:.4f}\n" in content
        assert content.endswith("$end\n")

    def test_write_without_force_constant(
        self, tmpdir, water_molecule, xtb_jobrunner_no_scratch
    ):
        settings = XTBJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="opt",
            constraints=[[1, 2], [2, 1, 3]],
        )
        job = XTBOptJob(
            molecule=water_molecule,
            settings=settings,
            label="xtb_water_constrain",
            jobrunner=xtb_jobrunner_no_scratch,
        )
        XTBInputWriter(job=job).write(target_directory=tmpdir)

        content = open(os.path.join(tmpdir, "xtb_water_constrain.inp")).read()
        assert "force constant" not in content
        assert "distance: 1, 2," in content
        assert "angle: 2, 1, 3," in content

    def test_write_skips_inp_when_no_constraints(
        self, tmpdir, water_molecule, xtb_jobrunner_no_scratch
    ):
        settings = XTBJobSettings(charge=0, multiplicity=1, jobtype="opt")
        job = XTBOptJob(
            molecule=water_molecule,
            settings=settings,
            label="xtb_no_constrain",
            jobrunner=xtb_jobrunner_no_scratch,
        )
        XTBInputWriter(job=job).write(target_directory=tmpdir)
        assert not os.path.isfile(os.path.join(tmpdir, "xtb_no_constrain.inp"))

    def test_invalid_constraint_length_raises(
        self, water_molecule, xtb_jobrunner_no_scratch
    ):
        settings = XTBJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="opt",
            constraints=[[1]],
        )
        job = XTBOptJob(
            molecule=water_molecule,
            settings=settings,
            label="bad",
            jobrunner=xtb_jobrunner_no_scratch,
        )
        with pytest.raises(ValueError, match="2 \\(distance\\)"):
            XTBInputWriter(job=job).write_constraints(io.StringIO())

    def test_out_of_range_atom_index_raises(
        self, water_molecule, xtb_jobrunner_no_scratch
    ):
        settings = XTBJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="opt",
            constraints=[[1, 4]],
        )
        job = XTBOptJob(
            molecule=water_molecule,
            settings=settings,
            label="bad",
            jobrunner=xtb_jobrunner_no_scratch,
        )
        with pytest.raises(ValueError, match="out of range"):
            XTBInputWriter(job=job).write_constraints(io.StringIO())
