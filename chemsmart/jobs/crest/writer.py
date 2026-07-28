"""
CREST input file writer for conformational search jobs.

Generates constraint input files (constraints.inp, xcontrol syntax) for use
with crest --cinp, and CREST 3.0 TOML input files for use with crest --input.
"""

import logging
import os

from chemsmart.jobs.writer import InputWriter

logger = logging.getLogger(__name__)


class CRESTInputWriter(InputWriter):
    """
    Writer class for generating CREST input files from job settings.

    Currently supports writing the xTB-syntax constraints file.
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
            self._write_constraints(f)

        logger.info(
            f"Wrote CREST constraints file with "
            f"{len(self.settings.constraints)} distance constraints"
        )

    def _write_constraints(self, f):
        """Write the $constrain block to an open file handle.

        Raises:
            ValueError: If atom_pairs contains invalid indices.
        """
        atom_pairs = self.settings.constraints
        force_constant = self.settings.force_constant

        n_atoms = len(self.job.molecule)
        for pair in atom_pairs:
            if len(pair) != 2:
                raise ValueError(
                    f"Each constraint must be a pair of atom indices, "
                    f"got {pair}."
                )
            for idx in pair:
                if idx < 1 or idx > n_atoms:
                    raise ValueError(
                        f"Atom index {idx} out of range [1, {n_atoms}]."
                    )

        f.write("$constrain\n")
        if force_constant is not None:
            f.write(f"  force constant={force_constant}\n")

        for pair in atom_pairs:
            i, j = pair[0], pair[1]
            distance = self.job.molecule.get_distance(i - 1, j - 1)
            f.write(f"  distance: {i}, {j}, {distance:.4f}\n")

        f.write("$end\n")

    def _write_toml(self, folder):
        """Write CREST 3.0 TOML input file for --input.

        Not yet implemented. Reserved for future CREST 3.0 support.
        """
        raise NotImplementedError(
            "CREST 3.0 TOML input file generation is not yet supported."
        )
