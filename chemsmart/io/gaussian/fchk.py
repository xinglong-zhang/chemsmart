"""Gaussian formatted checkpoint files and Cartesian force constants."""

import re
from functools import cached_property

import numpy as np

from chemsmart.utils.mixins import FileMixin


class GaussianFormattedCheckpoint(FileMixin):
    """Read Gaussian .fchk files, including unweighted Cartesian Hessians."""

    def __init__(self, filename):
        self.filename = filename

    @cached_property
    def num_atoms(self):
        """Read the atom count stored in the formatted checkpoint."""
        for line in self.contents:
            match = re.fullmatch(r"Number of atoms\s+I\s+(\d+)", line)
            if match and int(match[1]) > 0:
                return int(match[1])
        raise RuntimeError(
            f"Missing or invalid atom count in {self.filepath}."
        )

    def read_fchk_hessian(self):
        """Read the packed lower triangle in Hartree/Bohr², without rounding."""
        path = self.filepath
        n_atoms = self.num_atoms
        size = 3 * n_atoms
        expected = size * (size + 1) // 2
        values = None
        atom_count = None
        with open(path, encoding="utf-8") as stream:
            for line in stream:
                match = re.match(r"Number of atoms\s+I\s+(\d+)\s*$", line)
                if match:
                    atom_count = int(match[1])
                match = re.match(
                    r"Cartesian Force Constants\s+R\s+N=\s*(\d+)\s*$", line
                )
                if not match:
                    continue
                if int(match[1]) != expected:
                    raise RuntimeError(f"Hessian size mismatch in {path}.")
                values = []
                while len(values) < expected:
                    fields = next(stream, "").split()
                    if not fields:
                        raise RuntimeError(f"Incomplete Hessian in {path}.")
                    try:
                        values.extend(
                            float(v.replace("D", "E")) for v in fields
                        )
                    except ValueError as error:
                        raise RuntimeError(
                            f"Invalid Hessian data in {path}."
                        ) from error
                if len(values) != expected or not np.isfinite(values).all():
                    raise RuntimeError(f"Invalid Hessian data in {path}.")
        if atom_count != n_atoms or values is None:
            raise RuntimeError(
                f"Missing Hessian or atom count mismatch in {path}."
            )
        matrix = np.zeros((size, size))
        indices = np.tril_indices(size)
        matrix[indices] = values
        matrix[(indices[1], indices[0])] = values
        return matrix
