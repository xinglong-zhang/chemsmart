"""Read full Cartesian force constants from ORCA .hess files."""

import numpy as np

from chemsmart.utils.mixins import FileMixin


class ORCAHessianFile(FileMixin):
    """Read Cartesian force constants from an ORCA .hess file."""

    def __init__(self, filename):
        self.filename = filename

    def read_hessian(self, n_atoms):
        """Return the complete symmetric Hessian in Hartree/Bohr squared.

        ORCA writes numbered column blocks, each followed by every numbered
        matrix row. Missing, repeated, nonfinite or asymmetric data is rejected.
        """
        lines = self.contents
        try:
            start = lines.index("$hessian") + 1
            n = int(lines[start])
        except (ValueError, IndexError) as error:
            raise ValueError(
                "Missing or invalid ORCA $hessian section"
            ) from error
        if n != 3 * n_atoms:
            raise ValueError(
                "ORCA Hessian dimension does not match atom count"
            )
        matrix = np.full((n, n), np.nan)
        filled = np.zeros((n, n), dtype=bool)
        cursor = start + 1
        try:
            while cursor < len(lines):
                line = lines[cursor]
                cursor += 1
                if not line:
                    continue
                if line.startswith("$"):
                    break
                columns = [int(value) for value in line.split()]
                if not columns or len(set(columns)) != len(columns):
                    raise ValueError("Invalid Hessian columns")
                if min(columns) < 0 or max(columns) >= n:
                    raise ValueError("Hessian column out of range")
                for row in range(n):
                    values = lines[cursor].split()
                    cursor += 1
                    if (
                        len(values) != len(columns) + 1
                        or int(values[0]) != row
                    ):
                        raise ValueError("Invalid Hessian row")
                    if filled[row, columns].any():
                        raise ValueError("Repeated Hessian entries")
                    matrix[row, columns] = [
                        float(value.replace("D", "E").replace("d", "e"))
                        for value in values[1:]
                    ]
                    filled[row, columns] = True
        except (ValueError, IndexError) as error:
            raise ValueError("Malformed ORCA Hessian matrix") from error
        if not np.isfinite(matrix).all():
            raise ValueError("Incomplete or nonfinite ORCA Hessian")
        if not np.allclose(matrix, matrix.T, rtol=1e-7, atol=1e-9):
            raise ValueError("ORCA Hessian is not symmetric")
        return (matrix + matrix.T) / 2
