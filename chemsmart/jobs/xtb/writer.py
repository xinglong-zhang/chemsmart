"""
xTB input file writer.

Generates xcontrol detailed input files for use with xtb --input.
Currently supports the $constrain block (distance, angle, dihedral).
"""

import logging
import os

from chemsmart.jobs.writer import InputWriter

logger = logging.getLogger(__name__)


class XTBInputWriter(InputWriter):
    """
    Writer class for generating xTB xcontrol input files from job settings.

    Currently supports writing the $constrain block. Additional xcontrol
    groups can be added in future versions.
    """

    def write(self, **kwargs):
        """
        Write xTB input files.

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
        """Write {label}.inp (xTB xcontrol syntax) for --input."""
        job_inputfile = os.path.join(folder, f"{self.job.label}.inp")
        logger.info(f"Writing xTB input file to: {job_inputfile}")

        with open(job_inputfile, "w") as f:
            self._write_all(f)

        logger.info(
            f"Wrote xTB input file with "
            f"{len(self.settings.constraints)} constraint(s)"
        )

    def _write_all(self, f):
        """Write all currently supported xcontrol sections."""
        self.write_constraints(f)

    def write_constraints(self, f):
        """Write the $constrain block to an open file handle.

        Each entry in settings.constraints is a list of 1-based atom indices:
        - length 2: distance constraint
        - length 3: angle constraint
        - length 4: dihedral constraint

        Constraint values are computed from the current molecular geometry.

        Raises:
            ValueError: If a constraint has invalid length or atom indices.
        """
        constraints = self.settings.constraints
        if not constraints:
            return

        force_constant = self.settings.force_constant
        molecule = self.job.molecule
        n_atoms = len(molecule)

        for atoms in constraints:
            if len(atoms) not in (2, 3, 4):
                raise ValueError(
                    f"Each constraint must specify 2 (distance), 3 (angle), "
                    f"or 4 (dihedral) atom indices, got {atoms}."
                )
            for idx in atoms:
                if idx < 1 or idx > n_atoms:
                    raise ValueError(
                        f"Atom index {idx} out of range [1, {n_atoms}]."
                    )

        f.write("$constrain\n")
        if force_constant is not None:
            f.write(f"   force constant={force_constant}\n")

        for atoms in constraints:
            if len(atoms) == 2:
                idx1, idx2 = atoms
                distance = molecule.get_distance(idx1, idx2)
                f.write(f"   distance: {idx1}, {idx2}, {distance:.4f}\n")
            elif len(atoms) == 3:
                idx1, idx2, idx3 = atoms
                angle = molecule.get_angle(idx1, idx2, idx3)
                f.write(f"   angle: {idx1}, {idx2}, {idx3}, {angle:.4f}\n")
            elif len(atoms) == 4:
                idx1, idx2, idx3, idx4 = atoms
                dihedral = molecule.get_dihedral(idx1, idx2, idx3, idx4)
                f.write(
                    f"   dihedral: {idx1}, {idx2}, {idx3}, {idx4}, "
                    f"{dihedral:.4f}\n"
                )

        f.write("$end\n")
