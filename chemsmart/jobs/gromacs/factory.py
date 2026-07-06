from chemsmart.jobs.gromacs.job import GromacsEMJob, GromacsNVTJob
from chemsmart.settings.gromacs import GromacsProjectSettings


class GromacsJobFactory:
    """
    Create concrete GROMACS job objects from settings or project YAML files.
    """

    JOB_TYPE_MAP = {
        "em": GromacsEMJob,
        "nvt": GromacsNVTJob,
    }

    @classmethod
    def create_from_project_yaml(
        cls,
        project_yaml,
        molecule=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Create a GROMACS job from a project YAML file.
        """
        settings = GromacsProjectSettings.from_yaml(project_yaml)

        return cls.create_from_settings(
            settings=settings,
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def create_from_settings(
        cls,
        settings,
        molecule=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Create a GROMACS job from GromacsProjectSettings.
        """
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
