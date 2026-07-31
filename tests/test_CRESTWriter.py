import os

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
