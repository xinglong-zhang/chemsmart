from pathlib import Path

from chemsmart.jobs.job import Job


class GromacsJob(Job):
    PROGRAM = "gromacs"
    TYPE = "gromacs"

    def __init__(
        self,
        molecule,
        label,
        jobrunner,
        mdp_file=None,
        structure_file=None,
        top_file=None,
        tpr_file=None,
        itp_files=None,
        index_file=None,
        workflow="prepared",
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

        self.mdp_file = Path(mdp_file) if mdp_file else None
        self.structure_file = Path(structure_file) if structure_file else None
        self.top_file = Path(top_file) if top_file else None

        self._use_default_tpr_file = tpr_file is None
        self.tpr_file = (
            Path(tpr_file)
            if tpr_file is not None
            else Path(self.folder) / f"{self.label}.tpr"
        )
        self.index_file = Path(index_file) if index_file else None
        self.itp_files = [Path(f) for f in itp_files] if itp_files else []
        self.workflow = workflow

    @classmethod
    def from_project_settings(
        cls,
        settings,
        molecule=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Create a GROMACS job from GromacsProjectSettings.

        Project settings store reusable project-level information, while the
        job object stores the concrete inputs needed by the runner.
        """
        job_kwargs = settings.to_job_kwargs()

        if label is None:
            label = settings.project_name or "gromacs_job"

        job_kwargs.update(kwargs)

        return cls(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **job_kwargs,
        )

    def _run(self, **kwargs):
        self.jobrunner.run(self, **kwargs)

    def has_tpr(self):
        return self.tpr_file is not None and self.tpr_file.exists()

    def has_topology(self):
        return self.top_file is not None and self.top_file.exists()

    def has_required_prepared_inputs(self):
        """
        Check whether the job has the minimum user-provided files
        required to assemble a TPR file.
        """
        required_files = [
            self.mdp_file,
            self.structure_file,
            self.top_file,
        ]

        return all(
            path is not None and Path(path).exists()
            for path in required_files
        )

    def set_folder(self, folder):
        """
        Set the job folder and update the default TPR path if it was not
         explicitly provided by the user.
        """
        super().set_folder(folder)

        if self._use_default_tpr_file:
             self.tpr_file = Path(self.folder) / f"{self.label}.tpr"

class GromacsEMJob(GromacsJob):
    """
    Energy minimization job for GROMACS.
    """

    TYPE = "gmxem"

    def __init__(
        self,
        molecule=None,
        label="gromacs_em",
        jobrunner=None,
        mdp_file=None,
        structure_file=None,
        top_file=None,
        tpr_file=None,
        itp_files=None,
        index_file=None,
        workflow="prepared",
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            mdp_file=mdp_file,
            structure_file=structure_file,
            top_file=top_file,
            tpr_file=tpr_file,
            itp_files=itp_files,
            index_file=index_file,
            workflow=workflow,
            **kwargs,
        )
