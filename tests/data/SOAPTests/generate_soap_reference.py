#!/usr/bin/env python
"""Generate ``soap_reference.npz`` golden vectors from DScribe 2.1.2.

Requires DScribe 2.1.2 and ASE. Run from the repo root with the chemsmart
conda environment (or any env that has those packages)::

    python tests/data/SOAPTests/generate_soap_reference.py

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
        [0.0000, 0.0000, 0.1173],
        [0.0000, 0.7572, -0.4692],
        [0.0000, -0.7572, -0.4692],
    ],
    dtype=float,
)

METHANOL_SYM = ["C", "O", "H", "H", "H", "H"]
METHANOL_POS = np.array(
    [
        [0.000000, 0.000000, 0.000000],  # C
        [1.430000, 0.000000, 0.000000],  # O
        [1.109545, 0.904936, 0.000000],  # H (hydroxyl)
        [-0.363849, 1.027479, 0.000000],  # H
        [-0.363849, -0.513740, 0.889823],  # H
        [-0.363849, -0.513740, -0.889823],  # H
    ],
    dtype=float,
)

# Water + He spectator at 8.0 Å from O (inside r_cut=6 + ~3.72·σ padding shell).
PADDING_SYM = ["O", "H", "H", "He"]
PADDING_POS = np.vstack([WATER_POS, [[8.0, 0.0, 0.11779]]])

SULFUR_SYM = ["S", "H", "H"]
SULFUR_POS = WATER_POS.copy()

METHANE_SYM = ["C", "H", "H", "H", "H"]
METHANE_POS = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.09],
        [0.0, 1.03, -0.36],
        [0.89, -0.51, -0.36],
        [-0.89, -0.51, -0.36],
    ],
    dtype=float,
)

# Every atom lies on the z axis (the poles of the spherical harmonics).
HCN_SYM = ["H", "C", "N"]
HCN_POS = np.array(
    [[0.0, 0.0, -1.06], [0.0, 0.0, 0.0], [0.0, 0.0, 1.15]], dtype=float
)

# Asymmetric water, distinct from WATER_POS, used for the high-l check.
BENT_WATER_SYM = ["O", "H", "H"]
BENT_WATER_POS = np.array(
    [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]], dtype=float
)

# 17-atom CH3-CH(CH3)-CH2-CHFCl: four species and no symmetry.
HALO_CHAIN_SYM = ["C"] * 5 + ["F", "Cl"] + ["H"] * 10
HALO_CHAIN_POS = np.array(
    [
        [-1.0144, 1.0305, 0.9144],
        [-0.5694, 0.4225, -0.4231],
        [-1.2529, -0.9297, -0.664],
        [0.9715, 0.3516, -0.5595],
        [1.6781, -0.5734, 0.4436],
        [1.4033, -1.899, 0.1621],
        [3.441, -0.3363, 0.3383],
        [-0.4859, 1.99, 1.098],
        [-2.1056, 1.2374, 0.8889],
        [-0.8143, 0.3382, 1.7586],
        [-0.9238, 1.1101, -1.2233],
        [-2.3501, -0.7845, -0.7606],
        [-0.8811, -1.3889, -1.6046],
        [-1.071, -1.6266, 0.1808],
        [1.3772, 1.3807, -0.4425],
        [1.2285, 0.0227, -1.5904],
        [1.369, -0.3454, 1.4833],
    ],
    dtype=float,
)

ARGON_SYM = ["Ar"]
ARGON_POS = np.array([[0.0, 0.0, 0.0]], dtype=float)


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
        "methane_rcut5_lmax3": _dscribe_soap(
            METHANE_SYM,
            METHANE_POS,
            species=["C", "H"],
            r_cut=5.0,
            n_max=4,
            l_max=3,
            sigma=1.0,
        ),
        "hcn_linear_sigma_0_5": _dscribe_soap(
            HCN_SYM,
            HCN_POS,
            species=["C", "H", "N"],
            r_cut=5.0,
            n_max=4,
            l_max=3,
            sigma=0.5,
        ),
        "halo_chain_4species": _dscribe_soap(
            HALO_CHAIN_SYM,
            HALO_CHAIN_POS,
            species=["C", "Cl", "F", "H"],
            r_cut=5.0,
            n_max=4,
            l_max=3,
            sigma=1.0,
        ),
        "water_lmax9": _dscribe_soap(
            BENT_WATER_SYM,
            BENT_WATER_POS,
            species=["H", "O"],
            r_cut=4.0,
            n_max=3,
            l_max=9,
            sigma=0.7,
        ),
        "methane_lmax0_rcut3": _dscribe_soap(
            METHANE_SYM,
            METHANE_POS,
            species=["C", "H"],
            r_cut=3.0,
            n_max=2,
            l_max=0,
            sigma=0.3,
        ),
        "argon_single_atom": _dscribe_soap(
            ARGON_SYM,
            ARGON_POS,
            species=["Ar"],
            r_cut=5.0,
            n_max=4,
            l_max=3,
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
