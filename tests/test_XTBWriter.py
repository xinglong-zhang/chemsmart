import io
import os
from unittest.mock import MagicMock

import pytest

from chemsmart.jobs.xtb.opt import XTBOptJob
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.jobs.xtb.writer import XTBInputWriter


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

    def test_write_creates_target_directory_when_missing(
        self, tmp_path, water_molecule, xtb_jobrunner_no_scratch
    ):
        settings = XTBJobSettings(
            charge=0, multiplicity=1, jobtype="opt", constraints=[[1, 2]]
        )
        job = XTBOptJob(
            molecule=water_molecule,
            settings=settings,
            label="xtb_water_constrain",
            jobrunner=xtb_jobrunner_no_scratch,
        )
        nested = tmp_path / "nested" / "target"
        assert not nested.exists()

        XTBInputWriter(job=job).write(target_directory=str(nested))

        assert nested.is_dir()
        assert (nested / "xtb_water_constrain.inp").is_file()

    def test_write_uses_job_folder_when_no_target_directory(
        self, water_molecule, monkeypatch
    ):
        settings = XTBJobSettings(
            charge=0, multiplicity=1, jobtype="opt", constraints=[[1, 2]]
        )
        job = MagicMock()
        job.settings = settings
        job.molecule = water_molecule
        job.label = "mylabel"

        captured = {}

        def fake_write_inp(self, folder):
            captured["folder"] = folder

        monkeypatch.setattr(XTBInputWriter, "_write_inp", fake_write_inp)
        XTBInputWriter(job=job).write()

        assert captured["folder"] == job.folder

    def test_write_constraints_no_op_on_empty_constraints(
        self, water_molecule, xtb_jobrunner_no_scratch
    ):
        settings = XTBJobSettings(
            charge=0, multiplicity=1, jobtype="opt", constraints=None
        )
        job = XTBOptJob(
            molecule=water_molecule,
            settings=settings,
            label="mylabel",
            jobrunner=xtb_jobrunner_no_scratch,
        )
        buf = io.StringIO()
        XTBInputWriter(job=job).write_constraints(buf)
        assert buf.getvalue() == ""

    def test_write_multiple_constraints_including_non_terminal_dihedral(
        self, tmpdir, methane_molecule, xtb_jobrunner_no_scratch
    ):
        """A 4-atom (dihedral) constraint followed by another constraint
        exercises the loop continuing past the dihedral branch."""
        settings = XTBJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="opt",
            constraints=[[2, 1, 3, 4], [1, 2]],
        )
        job = XTBOptJob(
            molecule=methane_molecule,
            settings=settings,
            label="xtb_multi_constrain",
            jobrunner=xtb_jobrunner_no_scratch,
        )
        XTBInputWriter(job=job).write(target_directory=tmpdir)

        content = open(os.path.join(tmpdir, "xtb_multi_constrain.inp")).read()
        dihedral = methane_molecule.get_dihedral(2, 1, 3, 4)
        dist = methane_molecule.get_distance(1, 2)
        assert f"dihedral: 2, 1, 3, 4, {dihedral:.4f}\n" in content
        assert f"distance: 1, 2, {dist:.4f}\n" in content

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
