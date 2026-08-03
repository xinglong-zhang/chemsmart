"""SOAP descriptor calculations from molecular geometry.

Smooth Overlap of Atomic Positions (SOAP) descriptors are computed from
nuclear coordinates and elemental species via the optional DScribe package.
No bonding / connectivity information is required.

Standard geometric SOAP does **not** encode charge, spin multiplicity,
oxidation state, electronic configuration, or molecular chirality.
"""

from __future__ import annotations

import logging
from typing import Sequence

import numpy as np

from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)

_VALID_AGGREGATIONS = frozenset({None, "mean", "sum"})
# DScribe's default GTO radial basis requires r_cut > 1 Å.
_MIN_RCUT_ANGSTROM = 1.0


def _require_dscribe():
    """Import DScribe lazily and raise a clear error if missing."""
    try:
        from dscribe.descriptors import SOAP
    except ImportError as exc:  # pragma: no cover - exercised via mock
        raise ImportError(
            "SOAP calculations require DScribe. "
            "Install it with: pip install 'chemsmart[soap]'"
        ) from exc
    return SOAP


def _has_active_pbc(molecule: Molecule) -> bool:
    """Return True if the molecule has any active periodic directions.

    ``Molecule.pbc_conditions`` stores per-axis flags (truthy = periodic).
    Finite-system SOAP rejects any active periodicity.
    """
    conditions = molecule.pbc_conditions
    if conditions is None:
        return False
    return any(bool(c) for c in conditions)


def _normalize_species(species: Sequence[str] | None, symbols: list[str]):
    """Return the SOAP elemental species basis.

    When *species* is ``None``, returns a sorted unique list inferred from
    *symbols*. When *species* is provided, returns that list unchanged
    (order and duplicates preserved) after verifying it covers all atoms.
    """
    if species is None:
        return sorted(set(symbols))
    if len(species) == 0:
        raise ValueError("species must be a non-empty sequence of symbols.")
    normalized = [str(s) for s in species]
    missing = sorted(set(symbols) - set(normalized))
    if missing:
        raise ValueError(
            f"Molecule contains species not in the SOAP basis: {missing}. "
            "Pass an explicit shared species list that covers all atoms."
        )
    return normalized


def _normalize_centers(
    centers: Sequence[int] | None, num_atoms: int
) -> list[int] | None:
    """Convert optional 1-based center indices to 0-based DScribe indices."""
    if centers is None:
        return None
    if len(centers) == 0:
        raise ValueError("centers must be a non-empty sequence when provided.")

    zero_based: list[int] = []
    for raw in centers:
        if isinstance(raw, bool) or not isinstance(raw, (int, np.integer)):
            raise TypeError(
                "centers must be 1-based integer atom indices; "
                f"got {raw!r} ({type(raw).__name__})."
            )
        idx = int(raw)
        if idx < 1 or idx > num_atoms:
            raise ValueError(
                f"centers index {idx} is out of range for a molecule with "
                f"{num_atoms} atoms (expected 1..{num_atoms})."
            )
        zero_based.append(idx - 1)
    return zero_based


def _validate_hyperparameters(
    r_cut: float, n_max: int, l_max: int, sigma: float
) -> None:
    """Validate SOAP hyperparameters for DScribe's default GTO basis."""
    if not np.isfinite(r_cut) or r_cut <= _MIN_RCUT_ANGSTROM:
        raise ValueError(
            f"r_cut must be a finite number greater than "
            f"{_MIN_RCUT_ANGSTROM} Å (DScribe GTO basis requirement); "
            f"got {r_cut}."
        )
    if isinstance(n_max, bool) or not isinstance(n_max, (int, np.integer)):
        raise TypeError(
            f"n_max must be an integer; got {type(n_max).__name__}."
        )
    if int(n_max) < 1:
        raise ValueError(f"n_max must be >= 1; got {n_max}.")
    if isinstance(l_max, bool) or not isinstance(l_max, (int, np.integer)):
        raise TypeError(
            f"l_max must be an integer; got {type(l_max).__name__}."
        )
    if int(l_max) < 0:
        raise ValueError(f"l_max must be >= 0; got {l_max}.")
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError(
            f"sigma must be a positive finite number; got {sigma}."
        )


def _validate_positions(positions: np.ndarray, num_atoms: int) -> np.ndarray:
    """Validate and return a float64 copy of atomic positions."""
    arr = np.asarray(positions, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"positions must have shape (N, 3); got {arr.shape}.")
    if arr.shape[0] != num_atoms:
        raise ValueError(
            f"Number of positions ({arr.shape[0]}) does not match "
            f"number of symbols ({num_atoms})."
        )
    if not np.isfinite(arr).all():
        raise ValueError("positions must contain only finite values.")
    return arr


def _live_symbols(molecule: Molecule) -> list[str]:
    """Return current chemical symbols without using cached properties."""
    # Prefer ``symbols`` over the ``chemical_symbols`` cached_property so
    # SOAP always reflects the molecule's current elemental composition.
    if molecule.symbols is None:
        raise ValueError("Cannot compute SOAP for a molecule with no atoms.")
    return [str(s) for s in list(molecule.symbols)]


def calculate_soap(
    molecule: Molecule,
    *,
    r_cut: float = 6.0,
    n_max: int = 8,
    l_max: int = 6,
    sigma: float = 1.0,
    species: Sequence[str] | None = None,
    centers: Sequence[int] | None = None,
    aggregation: str | None = None,
) -> np.ndarray:
    """Compute SOAP descriptors for a finite (non-periodic) molecule.

    Args:
        molecule: ``Molecule`` instance with symbols and positions (Å).
        r_cut: Cutoff radius in Å for the local atomic environment.
            Must be greater than 1 Å for DScribe's default GTO basis.
        n_max: Number of radial basis functions.
        l_max: Maximum angular momentum (spherical harmonics degree).
        sigma: Gaussian width in Å (standard deviation) used to expand
            densities.
        species: Explicit elemental species basis for the SOAP feature space.
            When ``None``, a sorted unique list of the molecule's own elements
            is used (convenient for a single structure). For comparable
            descriptors across a dataset, pass a shared species list that
            covers every element present in the dataset.
        centers: Optional 1-based atom indices on which to evaluate SOAP.
            When ``None``, every atom is used as a center. Order is preserved.
        aggregation: ``None`` returns local (per-center) vectors with shape
            ``(n_centers, n_features)``. ``"mean"`` returns the arithmetic
            mean of local power spectra over the selected centers (DScribe
            outer averaging). ``"sum"`` returns the extensive sum over the
            selected centers. Mean is preferred for size-comparable
            fingerprints; sum scales with the number of centers.

    Returns:
        Dense ``float64`` NumPy array of SOAP features.

    Raises:
        ImportError: If DScribe is not installed.
        ValueError: For invalid geometry, hyperparameters, centers, species,
            aggregation, active periodic boundary conditions, or unrecognized
            elemental symbols.
        TypeError: For invalid center / hyperparameter types.
    """
    if not isinstance(molecule, Molecule):
        raise TypeError(
            f"molecule must be a Molecule instance; got {type(molecule)}."
        )
    if _has_active_pbc(molecule):
        raise ValueError(
            "Periodic SOAP is not supported. The molecule has active PBC "
            "conditions. Use a finite (non-periodic) molecular geometry."
        )
    if aggregation not in _VALID_AGGREGATIONS:
        raise ValueError(
            "aggregation must be one of "
            f"{sorted(_VALID_AGGREGATIONS - {None})!r} "
            f"or None; got {aggregation!r}."
        )

    _validate_hyperparameters(r_cut, n_max, l_max, sigma)

    symbols = _live_symbols(molecule)
    if len(symbols) == 0:
        raise ValueError("Cannot compute SOAP for a molecule with no atoms.")
    positions = _validate_positions(molecule.positions, len(symbols))
    species_list = _normalize_species(species, symbols)
    center_indices = _normalize_centers(centers, len(symbols))

    SOAP = _require_dscribe()
    from ase import Atoms

    try:
        atoms = Atoms(symbols=symbols, positions=positions)
    except KeyError as exc:
        raise ValueError(
            f"Unrecognized elemental symbol for ASE/SOAP: {exc.args[0]!r}. "
            "Use standard chemical element symbols (e.g. 'H' not 'D')."
        ) from exc

    soap = SOAP(
        species=species_list,
        periodic=False,
        r_cut=float(r_cut),
        n_max=int(n_max),
        l_max=int(l_max),
        sigma=float(sigma),
        average="off",
        sparse=False,
        dtype="float64",
    )

    logger.debug(
        "Computing SOAP for %d atoms (species=%s, r_cut=%s, n_max=%s, "
        "l_max=%s, sigma=%s, centers=%s, aggregation=%s)",
        len(symbols),
        species_list,
        r_cut,
        n_max,
        l_max,
        sigma,
        centers,
        aggregation,
    )

    features = soap.create(atoms, centers=center_indices)
    features = np.asarray(features, dtype=np.float64)

    # Defensive: older DScribe versions returned 1-D for a single center.
    if features.ndim == 1:  # pragma: no cover
        features = features.reshape(1, -1)
    if features.ndim != 2:  # pragma: no cover
        raise RuntimeError(
            f"Unexpected SOAP output shape {features.shape}; expected 2-D."
        )

    if aggregation == "mean":
        return np.mean(features, axis=0)
    if aggregation == "sum":
        return np.sum(features, axis=0)
    return features
