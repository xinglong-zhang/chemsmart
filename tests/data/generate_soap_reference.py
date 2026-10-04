#!/usr/bin/env python
"""Generate ``soap_reference.npz`` golden vectors from DScribe 2.1.2.

Requires DScribe 2.1.2 and ASE. Run from the repo root with the chemsmart
conda environment (or any env that has those packages)::

    python tests/data/generate_soap_reference.py

Existing keys are overwritten with freshly generated values. Regenerate only
when intentionally changing SOAP conventions.
"""

from __future__ import annotations

import importlib.metadata as md
from pathlib import Path

import numpy as np
from ase import Atoms
from dscribe.descriptors import SOAP

OUT_PATH = Path(__file__).resolve().parent / "soap_reference.npz"

WATER_SYM = ["O", "H", "H"]
WATER_POS = np.array(
    [
        [0.00000, 0.00000, 0.11779],
        [0.00000, 0.75545, -0.47116],
        [0.00000, -0.75545, -0.47116],
    ],
    dtype=float,
)

METHANOL_SYM = ["C", "O", "H", "H", "H", "H"]
METHANOL_POS = np.array(
    [
        [-0.6626, 0.0000, 0.0000],
        [0.6626, 0.0000, 0.0000],
        [-1.0280, 1.0280, 0.0000],
        [-1.0280, -0.5140, 0.8900],
        [-1.0280, -0.5140, -0.8900],
        [1.0280, 0.0000, 0.0000],
    ],
    dtype=float,
)

# Water + He spectator at 8.0 Å from O (inside r_cut=6 + ~3.72·σ padding shell).
PADDING_SYM = ["O", "H", "H", "He"]
PADDING_POS = np.vstack([WATER_POS, [[8.0, 0.0, 0.11779]]])

SULFUR_SYM = ["S", "H", "H"]
SULFUR_POS = WATER_POS.copy()


def _dscribe_soap(symbols, positions, **kwargs) -> np.ndarray:
    atoms = Atoms(symbols=symbols, positions=positions)
    soap = SOAP(rbf="gto", average="off", periodic=False, **kwargs)
    return np.asarray(soap.create(atoms), dtype=np.float64)


def main() -> None:
    arrays = {
        "water_sigma_0_5": _dscribe_soap(
            WATER_SYM,
            WATER_POS,
            species=["H", "O"],
            r_cut=6.0,
            n_max=4,
            l_max=2,
            sigma=0.5,
        ),
        "methanol_centers_2_1": _dscribe_soap(
            METHANOL_SYM,
            METHANOL_POS,
            species=["C", "H", "O"],
            r_cut=6.0,
            n_max=4,
            l_max=2,
            sigma=1.0,
        )[
            [1, 0], :
        ],  # 1-based centers [2, 1] -> 0-based [1, 0]
        "water_species_HCON": _dscribe_soap(
            WATER_SYM,
            WATER_POS,
            species=["H", "C", "O", "N"],
            r_cut=6.0,
            n_max=4,
            l_max=2,
            sigma=1.0,
        ),
        "sulfur_species_HOS": _dscribe_soap(
            SULFUR_SYM,
            SULFUR_POS,
            species=["H", "O", "S"],
            r_cut=6.0,
            n_max=4,
            l_max=2,
            sigma=1.0,
        ),
        "water_defaults": _dscribe_soap(
            WATER_SYM,
            WATER_POS,
            species=["H", "O"],
            r_cut=6.0,
            n_max=8,
            l_max=6,
            sigma=1.0,
        ),
        "water_lmax4_sigma1": _dscribe_soap(
            WATER_SYM,
            WATER_POS,
            species=["H", "O"],
            r_cut=6.0,
            n_max=4,
            l_max=4,
            sigma=1.0,
        ),
        "methanol_sigma1": _dscribe_soap(
            METHANOL_SYM,
            METHANOL_POS,
            species=["C", "H", "O"],
            r_cut=6.0,
            n_max=4,
            l_max=2,
            sigma=1.0,
        ),
        "padding_shell_He": _dscribe_soap(
            PADDING_SYM,
            PADDING_POS,
            species=["H", "He", "O"],
            r_cut=6.0,
            n_max=4,
            l_max=2,
            sigma=1.0,
        ),
    }

    np.savez(OUT_PATH, **arrays)
    print(f"Wrote {OUT_PATH}")
    for key, arr in sorted(arrays.items()):
        print(f"  {key}: {arr.shape}")
    print(
        "versions:",
        f"dscribe={md.version('dscribe')}",
        f"ase={md.version('ase')}",
        f"numpy={np.__version__}",
    )


if __name__ == "__main__":
    main()
