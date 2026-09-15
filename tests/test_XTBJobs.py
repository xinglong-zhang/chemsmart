from pathlib import Path

from chemsmart.jobs.xtb.opt import XTBOptJob
from chemsmart.jobs.xtb.settings import XTBJobSettings


class TestXTBJobs:
    def test_opt_job_with_constraints(
        self, temporary_working_dir, water_molecule, xtb_jobrunner_no_scratch
    ):
        settings = XTBJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="opt",
            freq=False,
            solvent_model="alpb",
            solvent_id="water",
            constraints=[[1, 2]],
            force_constant=0.5,
            gfn_version="gfn2",
        )
        job = XTBOptJob(
            molecule=water_molecule,
            settings=settings,
            label="water_opt",
            jobrunner=xtb_jobrunner_no_scratch,
        )
        assert isinstance(job, XTBOptJob)
        assert xtb_jobrunner_no_scratch.run(job) == 0

        folder = Path(job.folder)
        assert (folder / f"{job.label}.xyz").is_file()
        assert (folder / f"{job.label}.inp").is_file()
        assert (folder / f"{job.label}.out").is_file()

        inp = (folder / f"{job.label}.inp").read_text()
        assert "$constrain" in inp
        assert "distance: 1, 2" in inp
        assert "force constant=0.5" in inp

        out = (folder / f"{job.label}.out").read_text()
        assert "--opt" in out
        assert f"--input {job.label}.inp" in out
        assert "--gfn 2" in out
        assert "--alpb water" in out
