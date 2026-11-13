"""Job that calculates the Sterimol parameters of a substituent."""

import logging
import os

from chemsmart.jobs.job import Job

logger = logging.getLogger(__name__)


class SterimolJob(Job):
    """
    Calculate L, B1 and B5 of a substituent from one structure.

    The calculation runs in Python, with no external program.

    Attributes:
        PROGRAM (str): Program identifier ('Sterimol').
        TYPE (str): Job type identifier ('sterimol').
        molecule (Molecule): Structure to measure.
        atom1 (int): 1-based parent atom of the axis.
        atom2 (int): 1-based first substituent atom of the axis.
        substituent_atoms (list[int]): 1-based atoms to measure.
        filename (str): File that the structure was read from.
        index (str): Structure selection within the file.
        selection (str): How the atoms were chosen, for the report.
        output (str | None): Report path. Defaults to a name from the label.
    """

    PROGRAM = "Sterimol"
    TYPE = "sterimol"

    def __init__(
        self,
        molecule,
        atom1,
        atom2,
        substituent_atoms,
        filename,
        index="-1",
        selection="explicit",
        output=None,
        label=None,
        jobrunner=None,
        skip_completed=False,
        **kwargs,
    ):
        """
        Initialize the job.

        Completed jobs rerun by default: the calculation is cheap, and a
        report from the same file and axis may come from different atoms.
        """
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )
        self.atom1 = atom1
        self.atom2 = atom2
        self.substituent_atoms = substituent_atoms
        self.filename = filename
        self.index = index
        self.selection = selection
        self.output = output

    @property
    def outputfile(self):
        """Absolute path of the JSON report."""
        if self.output is not None:
            return os.path.abspath(self.output)
        return os.path.join(self.folder, f"{self.label}_sterimol.json")

    def _job_is_complete(self):
        """Return True if the report exists."""
        return os.path.exists(self.outputfile)

    def _run(self, **kwargs):
        """Run the job with its jobrunner."""
        self.jobrunner.run(self, **kwargs)
