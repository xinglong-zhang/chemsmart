"""PySCF job classes.

Unlike Gaussian and ORCA, PySCF is a Python library rather than an external
binary. ChemSmart therefore *generates* the calculation script and controls
both ends of the interface, which changes the file contract:

``label.py``   the generated, standalone-runnable driver script
``label.out``  PySCF's own log (human debugging; nothing parses it)
``label.h5``   the structured results -- the sole programmatic contract
``label.err``  stderr from the subprocess
"""

import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class PySCFJob(Job):
    """Base PySCF job.

    Attributes:
        PROGRAM (str): Program identifier ('PySCF').
        molecule (Molecule): Molecular structure used for the calculation.
        settings (PySCFJobSettings): Configuration options for the job.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    PROGRAM = "PySCF"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        if not isinstance(settings, PySCFJobSettings):
            raise ValueError(
                f"Settings must be instance of {PySCFJobSettings} for "
                f"{self}, but is {settings} instead!"
            )

        if not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is "
                f"{molecule} instead!"
            )

        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = settings.copy()

        if label is None:
            label = molecule.get_chemical_formula(empirical=True)
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[PySCFJobSettings]:
        return PySCFJobSettings

    @property
    def inputfile(self):
        """Path to the generated driver script."""
        return os.path.join(self.folder, self.label + ".py")

    @property
    def outputfile(self):
        """Path to PySCF's own log file."""
        return os.path.join(self.folder, self.label + ".out")

    @property
    def resultsfile(self):
        """Path to the structured results file.

        This is the only file read programmatically. ``label.out`` exists for
        humans, and is deliberately never parsed.
        """
        return os.path.join(self.folder, self.label + ".h5")

    @property
    def errfile(self):
        return os.path.join(self.folder, self.label + ".err")

    def _backup_files(self, **kwargs):
        """Back up the script, log and results together.

        All three are needed to reproduce or audit a run, so backing up a
        subset would leave an unusable record.
        """
        folder = self._create_backup_folder_name()
        self.backup_file(self.inputfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)
        self.backup_file(self.resultsfile, folder=folder, **kwargs)

    def _output(self):
        """Return a ``PySCFOutput``, or None when results are absent.

        Keyed on the results file rather than the log: a run that produced a
        log but no ``.h5`` did not finish, and must not read as complete.
        """
        if not os.path.exists(self.resultsfile):
            logger.debug(f"Results file not found: {self.resultsfile}")
            return None

        try:
            from chemsmart.io.pyscf.output import PySCFOutput

            return PySCFOutput(self.outputfile)
        except (OSError, KeyError, ValueError) as e:
            logger.error(f"Error creating PySCFOutput object: {e}")
            return None

    def _run(self, **kwargs):
        logger.info(f"Running PySCFJob {self} with jobrunner {self.jobrunner}")
        self.jobrunner.run(self, **kwargs)

    @classmethod
    def from_filename(
        cls,
        filename,
        settings=None,
        index="-1",
        label=None,
        jobrunner=None,
        keywords=("charge", "multiplicity"),
        **kwargs,
    ):
        """Create a PySCF job from a molecular structure file."""
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        logger.info(f"Number of molecules read: {len(molecules)}")
        molecules = molecules[string2index_1based(index)]

        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecules,
                    settings=settings,
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
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def from_pubchem(
        cls, identifier, settings=None, label=None, jobrunner=None, **kwargs
    ):
        """Create a PySCF job from a PubChem identifier."""
        molecule = Molecule.from_pubchem(identifier=identifier)

        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecule,
                    settings=settings,
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
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )


class PySCFGeneralJob(PySCFJob):
    """General PySCF job.

    Exists so that specialised subclasses can create and run a plain PySCF
    job internally without recursing into their own type.
    """

    TYPE = "pyscfjob"

    def __init__(self, molecule, settings=None, label=None, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
