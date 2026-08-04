import pytest

from chemsmart.cli.gromacs.factory import (
    GromacsJobFactory as CliGromacsJobFactory,
)
from chemsmart.jobs.gromacs.factory import (
    GromacsJobFactory as CoreGromacsJobFactory,
)
from chemsmart.jobs.gromacs.job import GromacsMDJob, GromacsNPTJob
from chemsmart.settings.gromacs import GromacsProjectSettings


@pytest.mark.parametrize(
    ("factory_cls", "job_type", "expected_job_cls"),
    [
        (CoreGromacsJobFactory, "npt", GromacsNPTJob),
        (CoreGromacsJobFactory, "md", GromacsMDJob),
        (CliGromacsJobFactory, "npt", GromacsNPTJob),
        (CliGromacsJobFactory, "md", GromacsMDJob),
    ],
)
def test_gromacs_factory_creates_npt_and_md_jobs(
    tmp_path,
    factory_cls,
    job_type,
    expected_job_cls,
):
    settings = GromacsProjectSettings.from_dict(
        {
            "project_name": f"prepared_{job_type}",
            "workflow": "prepared",
            "job_type": job_type,
            "structure_file": tmp_path / "input.gro",
            "top_file": tmp_path / "topol.top",
        }
    )

    job = factory_cls.create_from_settings(
        settings=settings,
        molecule=None,
        jobrunner=None,
    )

    assert isinstance(job, expected_job_cls)
    assert job.TYPE == expected_job_cls.TYPE
    assert job.label == f"prepared_{job_type}"
    assert job.workflow == "prepared"


@pytest.mark.parametrize(
    "factory_cls",
    [
        CoreGromacsJobFactory,
        CliGromacsJobFactory,
    ],
)
def test_gromacs_factory_rejects_unsupported_job_type(
    tmp_path,
    factory_cls,
):
    settings = GromacsProjectSettings.from_dict(
        {
            "project_name": "unsupported",
            "workflow": "prepared",
            "job_type": "invalid",
            "structure_file": tmp_path / "input.gro",
            "top_file": tmp_path / "topol.top",
        }
    )

    with pytest.raises(ValueError, match="Unsupported GROMACS job type"):
        factory_cls.create_from_settings(
            settings=settings,
            molecule=None,
            jobrunner=None,
        )