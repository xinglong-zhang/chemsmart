import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.jobs.crest.writer import CRESTInputWriter


class TestCRESTInputWriter:
    def test_write_reuses_xtb_constraints(
        self, tmpdir, methane_molecule, crest_jobrunner_no_scratch
    ):
        settings = CRESTJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="conformers",
            constraints=[[1, 2], [2, 1, 3], [2, 1, 3, 4]],
            force_constant=0.5,
        )
        job = CRESTConformerSearchJob(
            molecule=methane_molecule,
            settings=settings,
            label="crest_constrain",
            jobrunner=crest_jobrunner_no_scratch,
        )
        CRESTInputWriter(job=job).write(target_directory=tmpdir)

        content = open(os.path.join(tmpdir, "constraints.inp")).read()
        dist = methane_molecule.get_distance(1, 2)
        angle = methane_molecule.get_angle(2, 1, 3)
        dihedral = methane_molecule.get_dihedral(2, 1, 3, 4)
        assert f"   distance: 1, 2, {dist:.4f}\n" in content
        assert f"   angle: 2, 1, 3, {angle:.4f}\n" in content
        assert f"   dihedral: 2, 1, 3, 4, {dihedral:.4f}\n" in content
        assert "   force constant=0.5\n" in content

    def test_write_creates_target_directory_when_missing(self, tmp_path):
        job = MagicMock()
        job.settings = CRESTJobSettings(constraints=[[1, 2]])
        job.folder = str(tmp_path / "unused")

        nested = tmp_path / "nested" / "target"
        assert not nested.exists()

        with patch(
            "chemsmart.jobs.crest.writer.XTBInputWriter"
        ) as mock_xtb_writer_cls:
            CRESTInputWriter(job=job).write(target_directory=str(nested))

        assert nested.is_dir()
        assert (nested / "constraints.inp").is_file()
        mock_xtb_writer_cls.assert_called_once_with(job=job)

    def test_write_uses_job_folder_when_no_target_directory(self, tmp_path):
        job = MagicMock()
        job.settings = CRESTJobSettings(constraints=[[1, 2]])
        job.folder = str(tmp_path)

        with patch("chemsmart.jobs.crest.writer.XTBInputWriter"):
            CRESTInputWriter(job=job).write()

        assert (tmp_path / "constraints.inp").is_file()

    def test_write_skips_inp_file_when_no_constraints(self, tmp_path):
        job = MagicMock()
        job.settings = CRESTJobSettings(constraints=None)
        job.folder = str(tmp_path)

        CRESTInputWriter(job=job).write(target_directory=str(tmp_path))

        assert not (tmp_path / "constraints.inp").exists()

    def test_write_toml_not_yet_implemented(self, tmp_path):
        job = MagicMock()
        job.settings = CRESTJobSettings()
        job.folder = str(tmp_path)

        with pytest.raises(NotImplementedError, match="not yet supported"):
            CRESTInputWriter(job=job)._write_toml(str(tmp_path))
