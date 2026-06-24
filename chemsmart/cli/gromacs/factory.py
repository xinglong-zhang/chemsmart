from chemsmart.jobs.gromacs.job import GromacsEMJob
from chemsmart.settings.gromacs import GromacsProjectSettings


class GromacsJobFactory:
    """
    Factory for creating GROMACS jobs from project YAML or settings
    """

    JOB_TYPE_MAP = {
        "em": GromacsEMJob,
    }

    @classmethod
    def create_from_project_yaml(
        cls, project_yaml, molecule=None, label=None, jobrunner=None, **kwargs
    ):
        settings = GromacsProjectSettings.from_yaml(project_yaml)
        return cls.create_from_settings(
            settings,
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def create_from_settings(
        cls, settings, molecule=None, label=None, jobrunner=None, **kwargs
    ):
        settings.validate()
        job_cls = cls.JOB_TYPE_MAP.get(settings.job_type)
        if job_cls is None:
            raise ValueError(
                f"Unsupported GROMACS job type: {settings.job_type}"
            )
        return job_cls.from_project_settings(
            settings=settings,
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
