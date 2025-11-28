"""
PyMOL base job and helpers for molecular visualization.

Defines the core `PyMOLJob` used by higher-level visualization jobs.
It provides common file paths, backup utilities, and factory methods
to create visualization jobs from files or PubChem records.
"""

import logging
import os
import re

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class PyMOLJob(Job):
    """
    Base class for PyMOL molecular visualization jobs.

    Provides core functionality for creating and managing PyMOL
    visualization jobs including structure loading, styling options,
    and output management for molecular visualization tasks.

    Attributes:
        PROGRAM (str): Program identifier ('PyMOL').
        molecule: Molecule or list[Molecule] to visualize.
        label (str): Job identifier used for file naming and outputs.
        jobrunner (JobRunner): Execution backend responsible for running PyMOL.
        pymol_script (str | None): Optional custom PyMOL script path.
        style: Visualization style configuration.
        trace: Path/trace configuration for trajectories.
        vdw: Van der Waals representation settings.
        quiet_mode (bool): Whether to run PyMOL in quiet mode.
        command_line_only (bool): Whether to use CLI-only PyMOL execution.
        coordinates: Optional coordinates for labeling.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    PROGRAM = "PyMOL"

    def __init__(
        self,
        molecule=None,
        label=None,
        jobrunner=None,
        pymol_script=None,
        style=None,
        trace=None,
        vdw=None,
        quiet_mode=True,
        command_line_only=True,
        coordinates=None,
        isosurface_value=None,
        transparency_value=None,
        surface_quality=None,
        antialias_value=None,
        ray_trace_mode=None,
        label_offset=None,
        **kwargs,
    ):
        """
        Initialize a PyMOL visualization job.

        Sets up a PyMOL job with molecular structure, styling options,
        and execution parameters for generating molecular visualizations.

        Args:
            molecule: Molecule or list[Molecule] to visualize.
            label: Job identifier string (default: None).
            jobrunner: Runner for executing the job (default: None).
            pymol_script: Custom PyMOL script path (default: None).
            style: Visualization style settings (default: None).
            trace: Whether to trace molecular paths (default: None).
            vdw: Van der Waals representation settings (default: None).
            quiet_mode: Run PyMOL in quiet mode (default: True).
            command_line_only: Use command line interface only (default: True).
            coordinates: Coordinates for atom labeling (default: None).
            isosurface_value: Isosurface value for NCI visualization (optional).
            transparency_value: Transparency value for surfaces (optional).
            surface_quality: Surface quality setting (optional).
            antialias_value: Antialiasing level for rendering (optional).
            ray_trace_mode: Ray tracing mode for rendering (optional).
            label_offset: Offset for pymol labels (default: None).
            **kwargs: Additional arguments passed to parent Job class.
        """
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )
        self.pymol_script = pymol_script
        self.style = style
        self.trace = trace
        self.vdw = vdw
        self.quiet_mode = quiet_mode
        self.command_line_only = command_line_only
        self.coordinates = coordinates  # coordinates for labelling

        # Set defaults for pml files
        if isosurface_value is None:
            isosurface_value = 0.05
        if transparency_value is None:
            transparency_value = 0.2
        if surface_quality is None:
            surface_quality = 3
        if antialias_value is None:
            antialias_value = 3
        if ray_trace_mode is None:
            ray_trace_mode = 1

        self.isosurface_value = isosurface_value
        self.transparency_value = transparency_value
        self.surface_quality = surface_quality
        self.antialias_value = antialias_value
        self.ray_trace_mode = ray_trace_mode
        self.label_offset = label_offset

        logger.debug(f"Isosurface value: {self.isosurface_value}")
        logger.debug(f"Transparency value: {self.transparency_value}")
        logger.debug(f"Surface quality: {self.surface_quality}")
        logger.debug(f"Antialias value: {self.antialias_value}")
        logger.debug(f"Ray trace mode: {self.ray_trace_mode}")

    @property
    def job_basename(self):
        """
        Get the base name for job-related files.
        Returns:
            str: Base name derived from _get_job_basename method.
        """
        return self._get_job_basename()

    def _get_job_basename(self):
        """
        Internal method to derive the job base name.

        Returns:
            str: Base name derived from the job label.
        """
        return self.label

    @property
    def inputfile(self):
        """
        Get the path to the input XYZ file for PyMOL.

        Returns:
            str: Absolute path to the input XYZ coordinate file.
        """
        inputfile = self.label + ".xyz"
        return os.path.join(self.folder, inputfile)

    @property
    def logfile(self):
        """
        Get the path to the PyMOL log file.

        Returns:
            str: Absolute path to the job log file.
        """
        logfile = "log." + self.job_basename
        return os.path.join(self.folder, logfile)

    @property
    def outputfile(self):
        """
        Get the path to the PyMOL session file output.

        Returns:
            str: Absolute path to the output PSE session file.
        """
        outputfile = self.job_basename + ".pse"
        return os.path.join(self.folder, outputfile)

    @property
    def errfile(self):
        """
        Get the path to the error file for job diagnostics.

        Returns:
            str: Absolute path to the error log file.
        """
        errfile = self.job_basename + ".err"
        return os.path.join(self.folder, errfile)

    def _backup_files(self, backup_chk=False, **kwargs):
        """
        Create backup copies of job-related files.

        Backs up input, output, and optionally checkpoint files to
        a designated backup folder for data preservation.

        Args:
            backup_chk: Whether to backup checkpoint files (default: False).
            **kwargs: Additional arguments for backup operations.
        """
        folder = self._create_backup_folder_name()
        self.backup_file(self.inputfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)
        if backup_chk:
            self.backup_file(self.chkfile, folder=folder, **kwargs)

    def _output(self):
        """
        Get the absolute path to the job output file if it exists.

        Returns:
            str or None: Absolute path to output file, None if not found.
        """
        if not os.path.exists(self.outputfile):
            return None
        return os.path.abspath(self.outputfile)

    def _job_is_complete(self):
        """
        Check if the PyMOL job has completed successfully.

        Determines job completion by checking for the existence of
        the expected output PSE session file.

        Returns:
            bool: True if job completed, False otherwise.
        """
        return os.path.exists(self.outputfile)

    def _run(self):
        """
        Execute the PyMOL job using the configured job runner.

        Delegates job execution to the assigned jobrunner instance
        which handles the actual PyMOL process management.
        """
        self.jobrunner.run(self)

    @classmethod
    def from_filename(
        cls,
        filename,
        pymol_script=None,
        index="-1",
        label=None,
        jobrunner=None,
        style=None,
        vdw=None,
        quiet_mode=True,
        command_line_only=True,
        **kwargs,
    ):
        """
        Create a PyMOL job from a molecular structure file.

        Reads molecular structures from various file formats and creates
        a PyMOL visualization job with specified parameters and styling.

        Args:
            filename: Path to the molecular structure file.
            pymol_script: Custom PyMOL script path (default: None).
            index: Molecule index selection string (default: "-1").
            label: Job identifier (default: derived from filename).
            jobrunner: Job execution runner (default: auto-created).
            style: Visualization style settings (default: None).
            vdw: Van der Waals representation (default: None).
            quiet_mode: Run PyMOL quietly (default: True).
            command_line_only: Use CLI interface only (default: True).
            **kwargs: Additional arguments for job configuration.

        Returns:
            PyMOLJob: Configured PyMOL visualization job instance.
        """
        # get all molecule in a file and give the result as a list
        filename = re.sub(r" ", r"\\", filename)
        logger.info(f"Reading molecules from file: {filename}.")
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True, **kwargs
        )
        logger.info(f"Num of molecules read: {len(molecules)}.")

        if label is None:
            # by default, if no label is given and the job is read in
            # from a file, the label is set to the file basename
            label = os.path.basename(filename).split(".")[0]

        logger.info(f"Num of molecules read: {len(molecules)}.")
        molecules = molecules[string2index_1based(index)]
        logger.info(f"Num of molecules to use: {len(molecules)}.")

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecules,
                    label=label,
                    pymol_script=pymol_script,
                    style=style,
                    vdw=vdw,
                    quiet_mode=quiet_mode,
                    command_line_only=command_line_only,
                    jobrunner=None,
                    **kwargs,
                ),
                server=kwargs.get("server"),
                scratch=kwargs.get("scratch"),
                fake=kwargs.get("fake", False),
                **kwargs,
            )

        return cls(
            molecule=molecules,
            label=label,
            jobrunner=jobrunner,
            pymol_script=pymol_script,
            style=style,
            vdw=vdw,
            quiet_mode=quiet_mode,
            command_line_only=command_line_only,
            **kwargs,
        )

    @classmethod
    def from_pubchem(cls, identifier, label=None, jobrunner=None, **kwargs):
        """
        Create a PyMOL job from a PubChem molecular database entry.

        Downloads molecular structure data from PubChem database and
        creates a PyMOL visualization job for the specified compound.

        Args:
            identifier: PubChem identifier (CID, name, etc.).
            label: Job identifier string (default: None).
            jobrunner: Job execution runner (default: auto-created).
            **kwargs: Additional arguments for job configuration.

        Returns:
            PyMOLJob: Configured PyMOL job with PubChem structure data.
        """
        molecules = Molecule.from_pubchem(identifier=identifier)

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecules,
                    label=label,
                    jobrunner=None,
                    **kwargs,
                ),
                server=kwargs.get("server"),
                scratch=kwargs.get("scratch"),
                fake=kwargs.get("fake", False),
                **kwargs,
            )

        return cls(
            molecule=molecules,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
