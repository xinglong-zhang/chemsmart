"""Read unweighted Cartesian force constants from Gaussian files."""

import re

import numpy as np


def read_fchk_hessian(path, n_atoms):
    """Read the packed lower triangle in Hartree/Bohr², without rounding."""
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
                    values.extend(float(v.replace("D", "E")) for v in fields)
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


def read_cartesian_hessian(path, n_atoms):
    """Read the last complete lower-triangular Hessian (Hartree/Bohr²).

    Requires Gaussian #p output. Reject incomplete matrices rather than
    silently treating missing force constants as zeros.
    """
    size = 3 * n_atoms
    result = None
    matrix = None
    columns = []
    with open(path, encoding="utf-8") as stream:
        for line in stream:
            if "Force constants in Cartesian coordinates:" in line:
                matrix = np.full((size, size), np.nan)
                columns = []
                continue
            if matrix is None:
                continue
            fields = line.split()
            if not fields:
                continue
            if all(field.isdigit() for field in fields):
                columns = [int(field) - 1 for field in fields]
                continue
            if not fields[0].isdigit() or not columns:
                continue
            row = int(fields[0]) - 1
            try:
                values = [
                    float(value.replace("D", "E")) for value in fields[1:]
                ]
            except ValueError:
                continue
            if not 0 <= row < size:
                continue
            for col, value in zip(columns, values):
                if 0 <= col <= row:
                    matrix[row, col] = matrix[col, row] = value
            if np.isfinite(matrix).all():
                result = matrix.copy()
                matrix = None
    if matrix is not None or result is None:
        raise RuntimeError(
            f"Missing or incomplete Cartesian Hessian in {path}."
        )
    return result
