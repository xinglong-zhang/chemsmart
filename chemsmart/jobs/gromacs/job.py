from __future__ import annotations

"""
GROMACS job definitions.
"""

from pathlib import Path

from chemsmart.jobs.job import Job


class GromacsJob(Job):
    """
    Base GROMACS job.

    This class stores concrete input files and workflow options used by the
    runner. Reusable project-level settings should be stored in
    GromacsProjectSettings.
    """

    PROGRAM = "gromacs"
    TYPE = "gromacs"

    def __init__(
        self,
        molecule,
        label,
        jobrunner,
        mdp_file=None,
        ions_mdp_file=None,
        structure_file=None,
        input_pdb=None,
        top_file=None,
        tpr_file=None,
        itp_files=None,
        index_file=None,
        workflow="prepared",
        processed_structure_file=None,
        boxed_structure_file=None,
        solvated_structure_file=None,
        ions_tpr_file=None,
        ionized_structure_file=None,
        force_field=None,
        water_model=None,
        timestep=None,
        temperature=None,
        pressure=None,
        thermostat=None,
        barostat=None,
        constraints=None,
        constraint_algorithm=None,
        nsteps=None,
        emtol=None,
        emstep=None,
        tau_t=None,
        tc_grps=None,
        tau_p=None,
        compressibility=None,
        box_type="cubic",
        box_distance=1.0,
        solvent_file=None,
        positive_ion="NA",
        negative_ion="CL",
        neutral=True,
        genion_group="SOL",
        grompp_maxwarn=None,
        mdrun_threads=None,
        mdrun_ntmpi=None,
        mdrun_ntomp=None,
        mdrun_extra_args=None,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

        self.mdp_file = Path(mdp_file) if mdp_file else None
        self.ions_mdp_file = Path(ions_mdp_file) if ions_mdp_file else None

        self.structure_file = Path(structure_file) if structure_file else None
        self.input_pdb = Path(input_pdb) if input_pdb else None

        self.top_file = Path(top_file) if top_file else None
        self.index_file = Path(index_file) if index_file else None
        self.itp_files = [Path(f) for f in itp_files] if itp_files else []

        self.workflow = workflow

        self._use_default_tpr_file = tpr_file is None
        self.tpr_file = (
            Path(tpr_file)
            if tpr_file is not None
            else Path(self.folder) / f"{self.label}.tpr"
        )

        self.processed_structure_file = (
            Path(processed_structure_file)
            if processed_structure_file
            else None
        )
        self.boxed_structure_file = (
            Path(boxed_structure_file)
            if boxed_structure_file
            else None
        )
        self.solvated_structure_file = (
            Path(solvated_structure_file)
            if solvated_structure_file
            else None
        )
        self.ions_tpr_file = Path(ions_tpr_file) if ions_tpr_file else None
        self.ionized_structure_file = (
            Path(ionized_structure_file)
            if ionized_structure_file
            else None
        )

        self.force_field = force_field
        self.water_model = water_model

        self.timestep = timestep
        self.temperature = temperature
        self.pressure = pressure
        self.thermostat = thermostat
        self.barostat = barostat
        self.constraints = constraints
        self.constraint_algorithm = constraint_algorithm

        self.nsteps = nsteps
        self.emtol = emtol
        self.emstep = emstep
        self.tau_t = tau_t
        self.tc_grps = tc_grps
        self.tau_p = tau_p
        self.compressibility = compressibility

        self.box_type = box_type
        self.box_distance = box_distance
        self.solvent_file = Path(solvent_file) if solvent_file else None

        self.positive_ion = positive_ion
        self.negative_ion = negative_ion
        self.neutral = neutral
        self.genion_group = genion_group

        self.grompp_maxwarn = grompp_maxwarn
        self.mdrun_threads = mdrun_threads
        self.mdrun_ntmpi = mdrun_ntmpi
        self.mdrun_ntomp = mdrun_ntomp
        self.mdrun_extra_args = list(mdrun_extra_args or [])

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
        """
        Return True if the job already has an existing TPR file.
        """
        return self.tpr_file is not None and self.tpr_file.exists()

    def has_topology(self):
        """
        Return True if the job already has an existing topology file.
        """
        return self.top_file is not None and self.top_file.exists()

    def has_required_prepared_inputs(self):
        """
        Check whether the job has the minimum files for prepared workflow.
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

    def has_required_full_setup_inputs(self):
        """
        Check whether the job has the minimum inputs for full setup workflow.
        """
        input_file = self.input_pdb or self.structure_file

        return (
            input_file is not None
            and Path(input_file).exists()
            and self.mdp_file is not None
            and Path(self.mdp_file).exists()
            and self.force_field is not None
            and self.water_model is not None
        )

    def set_folder(self, folder):
        """
        Set the job folder and update default TPR path if needed.
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
        ions_mdp_file=None,
        structure_file=None,
        input_pdb=None,
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
            ions_mdp_file=ions_mdp_file,
            structure_file=structure_file,
            input_pdb=input_pdb,
            top_file=top_file,
            tpr_file=tpr_file,
            itp_files=itp_files,
            index_file=index_file,
            workflow=workflow,
            **kwargs,
        )


class GromacsNVTJob(GromacsJob):
    """
    NVT equilibration job for GROMACS.
    """

    TYPE = "gmxnvt"

    def __init__(
        self,
        molecule=None,
        label="gromacs_nvt",
        jobrunner=None,
        mdp_file=None,
        structure_file=None,
        input_pdb=None,
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
            input_pdb=input_pdb,
            top_file=top_file,
            tpr_file=tpr_file,
            itp_files=itp_files,
            index_file=index_file,
            workflow=workflow,
            **kwargs,
        )


class GromacsNPTJob(GromacsJob):
    """
    NPT equilibration job for GROMACS.
    """

    TYPE = "gmxnpt"

    def __init__(
        self,
        molecule=None,
        label="gromacs_npt",
        jobrunner=None,
        mdp_file=None,
        structure_file=None,
        input_pdb=None,
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
            input_pdb=input_pdb,
            top_file=top_file,
            tpr_file=tpr_file,
            itp_files=itp_files,
            index_file=index_file,
            workflow=workflow,
            **kwargs,
        )
