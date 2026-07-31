"""
CREST input file writer for conformational search jobs.

Generates constraint input files (constraints.inp, xcontrol syntax) for use
with crest --cinp, and CREST 3.0 TOML input files for use with crest --input.
"""

import logging
import os

from chemsmart.jobs.writer import InputWriter
from chemsmart.jobs.xtb.writer import XTBInputWriter

logger = logging.getLogger(__name__)


class CRESTInputWriter(InputWriter):
    """
    Writer class for generating CREST input files from job settings.

    Currently supports writing the xTB-syntax constraints file via
    XTBInputWriter.write_constraints (distance, angle, dihedral).
    Future versions will add TOML input file generation for CREST 3.0.
    """

    def write(self, **kwargs):
        """
        Write CREST input files.

        Args:
            **kwargs: Additional arguments passed to _write method.
        """
        self._write(**kwargs)

    def _write(self, target_directory=None):
        if target_directory is not None:
            if not os.path.exists(target_directory):
                logger.debug(f"Creating target directory: {target_directory}")
                os.makedirs(target_directory)
            folder = target_directory
        else:
            folder = self.job.folder

        if self.settings.constraints:
            self._write_inp(folder)

    def _write_inp(self, folder):
        """Write constraints.inp (xTB xcontrol syntax) for --cinp."""
        job_inputfile = os.path.join(folder, "constraints.inp")
        logger.info(f"Writing CREST constraints file to: {job_inputfile}")

        with open(job_inputfile, "w") as f:
            self.write_constraints(f)

        logger.info(
            f"Wrote CREST constraints file with "
            f"{len(self.settings.constraints)} constraint(s)"
        )

    def write_constraints(self, f):
        """Write the $constrain block using XTBInputWriter."""
        XTBInputWriter(job=self.job).write_constraints(f)

    def _write_toml(self, folder):
        """Write CREST 3.0 TOML input file for --input.

        Not yet implemented. Reserved for future CREST 3.0 support.
        """
        raise NotImplementedError(
            "CREST 3.0 TOML input file generation is not yet supported."
        )
