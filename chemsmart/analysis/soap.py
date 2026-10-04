"""SOAP descriptor calculations from molecular geometry.

Smooth Overlap of Atomic Positions (SOAP) descriptors are computed from
nuclear coordinates and elemental species using a built-in NumPy/SciPy
implementation. No bonding / connectivity information is required.

The GTO SOAP formulation matches DScribe 2.1.2 with ``rbf="gto"``,
``average="off"``, and default ``compression="off"``.

Standard geometric SOAP does **not** encode charge, spin multiplicity,
oxidation state, electronic configuration, or molecular chirality.
"""

from __future__ import annotations

import logging
from typing import Sequence

import numpy as np
from scipy.linalg import sqrtm
from scipy.special import factorial, gamma, lpmv

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.periodictable import PeriodicTable

logger = logging.getLogger(__name__)

_VALID_AGGREGATIONS = frozenset({None, "mean", "sum"})
# GTO radial basis requires r_cut > 1 Å (same constraint as DScribe).
_MIN_RCUT_ANGSTROM = 1.0
# Match DScribe 2.1.2's hard upper bound on spherical-harmonic degree.
_MAX_LMAX = 20
_GTO_DECAY_THRESHOLD = 1e-3
# Guard divisions that are analytically non-zero but can underflow numerically.
_EPS = 1e-15
# Per-l angular prefactors in DScribe's analytic GTO integral (soapGTO.cpp).
_RADIAL_ANGULAR_SCALE = {
    0: np.pi / 2.0,
    1: np.pi * np.sqrt(3.0) / 2.0,
}
_PERIODIC_TABLE = PeriodicTable()


def _has_active_pbc(molecule: Molecule) -> bool:
    """Return True if the molecule has any active periodic directions.

    ``Molecule.pbc_conditions`` stores per-axis flags (truthy = periodic).
    Finite-system SOAP rejects any active periodicity.
    """
    conditions = molecule.pbc_conditions
    if conditions is None:
        return False
    return any(bool(c) for c in conditions)


def _has_translation_vectors(molecule: Molecule) -> bool:
    """Return True if the molecule carries a non-empty cell / TV list.

    Some constructors (e.g. ``Molecule.from_coordinate_block_text``) set
    ``translation_vectors`` without setting ``pbc_conditions``. Finite-system
    SOAP rejects those geometries as well.
    """
    tvs = molecule.translation_vectors
    if tvs is None:
        return False
    try:
        return len(tvs) > 0
    except TypeError:
        return True


def _unrecognized_symbol_error(symbol: str) -> ValueError:
    """SOAP error for a symbol that is not a standard chemical element."""
    return ValueError(
        f"Unrecognized elemental symbol for SOAP: {symbol!r}. "
        "Use standard chemical element symbols (e.g. 'H' not 'D')."
    )


def _canonicalize_symbol(symbol: str) -> str:
    """Return a canonical element symbol (``fe`` → ``Fe``).

    Uses ``PeriodicTable.to_element`` for the same capitalization rules as
    the rest of CHEMSMART. The stripped input must be that canonical symbol
    ignoring case, so prefix matches such as ``Xx`` → ``X`` are rejected.
    """
    raw = str(symbol).strip()
    try:
        canonical = _PERIODIC_TABLE.to_element(raw)
    except ValueError as exc:
        raise _unrecognized_symbol_error(symbol) from exc
    if (
        canonical not in _PERIODIC_TABLE.PERIODIC_TABLE
        or raw.casefold() != canonical.casefold()
    ):
        raise _unrecognized_symbol_error(symbol)
    return canonical


def _normalize_species(species: Sequence[str] | None, symbols: list[str]):
    """Return the SOAP elemental species basis.

    When *species* is ``None``, returns a sorted unique list inferred from
    *symbols*. When *species* is provided, returns that list after verifying
    each entry is a recognized element and that the list covers all atoms.
    Downstream channel indexing deduplicates by atomic number (DScribe
    convention); user list order does not control feature layout.
    """
    if species is None:
        return sorted(set(symbols))
    if len(species) == 0:
        raise ValueError("species must be a non-empty sequence of symbols.")
    normalized = [_canonicalize_symbol(s) for s in species]
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
    """Convert optional 1-based center indices to 0-based center indices."""
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
    """Validate SOAP hyperparameters for the default GTO basis."""
    if not np.isfinite(r_cut) or r_cut <= _MIN_RCUT_ANGSTROM:
        raise ValueError(
            f"r_cut must be a finite number greater than "
            f"{_MIN_RCUT_ANGSTROM} Å (GTO basis requirement); "
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
    if int(l_max) > _MAX_LMAX:
        raise ValueError(
            f"l_max must be <= {_MAX_LMAX} (DScribe GTO parity limit); "
            f"got {l_max}."
        )
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
    """Return current chemical symbols from the live ``symbols`` attribute.

    Symbols are case-normalized (``fe`` → ``Fe``) so channel indexing
    matches ``PeriodicTable.to_atomic_number``.
    """
    # Prefer ``symbols`` over ``chemical_symbols`` so SOAP always reflects
    # the molecule's current elemental composition after in-place mutation.
    if molecule.symbols is None:
        raise ValueError("Cannot compute SOAP for a molecule with no atoms.")
    return [_canonicalize_symbol(s) for s in list(molecule.symbols)]


def _dscribe_species_index_map(species_list: Sequence[str]) -> dict[str, int]:
    """Map element symbols to SOAP channel indices (DScribe atomic-number order).

    DScribe's ``get_atomic_numbers`` deduplicates and sorts by atomic number,
    not by the order symbols appear in the user-supplied species list.
    """
    atomic_numbers = sorted(
        {_PERIODIC_TABLE.to_atomic_number(sym) for sym in species_list}
    )
    return {
        _PERIODIC_TABLE.to_symbol(Z): idx
        for idx, Z in enumerate(atomic_numbers)
    }


def _gto_basis(
    r_cut: float, n_max: int, l_max: int
) -> tuple[np.ndarray, np.ndarray]:
    """Return GTO radial ``alphas`` and Löwdin ``betas`` (DScribe convention)."""
    a = np.linspace(1.0, float(r_cut), int(n_max))
    alphas_full = np.zeros((l_max + 1, n_max), dtype=np.float64)
    betas_full = np.zeros((l_max + 1, n_max, n_max), dtype=np.float64)

    for ang_l in range(l_max + 1):
        # a is linspace(1, r_cut, n_max) with r_cut > 1, so a >= 1 analytically;
        # still floor denominators against underflow / accidental zeros.
        a_pow = np.maximum(np.power(a, ang_l), _EPS)
        a_sq = np.maximum(a * a, _EPS)
        alphas = -np.log(_GTO_DECAY_THRESHOLD / a_pow) / a_sq
        m = np.maximum(alphas[:, None] + alphas[None, :], _EPS)
        overlap = 0.5 * gamma(ang_l + 1.5) * m ** (-(ang_l + 1.5))
        betas = sqrtm(np.linalg.inv(overlap))
        if np.iscomplexobj(betas):
            raise ValueError(
                "Could not calculate normalization factors for the radial "
                "basis in the domain of real numbers. Lowering the number of "
                "radial basis functions (n_max) or increasing the radial "
                "cutoff (r_cut) is advised."
            )
        alphas_full[ang_l, :] = alphas
        betas_full[ang_l, :, :] = np.real(betas)

    return alphas_full, betas_full


def _solid_harmonics_polynomial(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, l_max: int
) -> np.ndarray:
    """Real solid harmonics ``R_l^m(x, y, z)`` for ``l >= 2``.

    Returns an array of shape ``(n_points, (l_max + 1) ** 2)``. Entries for
    ``l < 2`` are left zero; callers handle ``l = 0`` and ``l = 1`` directly.
    """
    n_pts = x.shape[0]
    n_harm = (l_max + 1) ** 2
    result = np.zeros((n_pts, n_harm), dtype=np.float64)
    if l_max < 2 or n_pts == 0:
        return result

    r2 = x * x + y * y + z * z
    r = np.sqrt(r2)
    # Floor r before dividing so coincident atoms (r=0) and tiny separations
    # cannot produce 0/0 or huge cos_theta; clip for SciPy lpmv domain.
    safe_r = np.maximum(r, _EPS)
    cos_theta = np.clip(z / safe_r, -1.0, 1.0)
    cos_theta = np.where(r > 0.0, cos_theta, 0.0)
    phi = np.arctan2(y, x)

    for ang_l in range(2, l_max + 1):
        rl = np.power(r, ang_l)
        for im, m in enumerate(range(-ang_l, ang_l + 1)):
            idx = ang_l * ang_l + im
            abs_m = abs(m)
            # factorial(ang_l + abs_m) >= 1 for ang_l >= 2, abs_m >= 0.
            denom = max(float(factorial(ang_l + abs_m, exact=True)), _EPS)
            norm = np.sqrt(
                (2 * ang_l + 1)
                / (4.0 * np.pi)
                * float(factorial(ang_l - abs_m, exact=True))
                / denom
            )
            plm = lpmv(abs_m, ang_l, cos_theta)
            if m < 0:
                trig = np.sin(abs_m * phi)
            elif m > 0:
                trig = np.cos(m * phi)
            else:
                trig = 1.0
            if m != 0:
                trig *= np.sqrt(2.0)
            result[:, idx] = rl * norm * plm * trig

    return result


def _power_spectrum_prefactor(ang_l: int) -> float:
    """Wigner-D normalization prefactor used by DScribe's ``getPD``."""
    base = np.pi * np.sqrt(8.0 / max(2 * ang_l + 1, 1))
    if ang_l > 1:
        return base * np.pi**3
    return base


def _power_spectrum(
    cnnd: np.ndarray, n_species: int, n_max: int, l_max: int
) -> np.ndarray:
    """Flatten SOAP power spectrum in DScribe ``compression='off'`` layout."""
    n_features = (
        (n_species * n_max) * (n_species * n_max + 1) // 2 * (l_max + 1)
    )
    features = np.zeros(n_features, dtype=np.float64)
    shift = 0
    for j in range(n_species):
        for jd in range(j, n_species):
            for ang_l in range(l_max + 1):
                prel = _power_spectrum_prefactor(ang_l)
                m_start = ang_l * ang_l
                m_end = (ang_l + 1) ** 2
                # (n_max, 2l+1) @ (n_max, 2l+1).T -> (n_max, n_max) dots
                dots = cnnd[j, :, m_start:m_end] @ cnnd[jd, :, m_start:m_end].T
                if j == jd:
                    # Upper triangle including diagonal, row-major (k, kd>=k)
                    for k in range(n_max):
                        for kd in range(k, n_max):
                            features[shift] = prel * dots[k, kd]
                            shift += 1
                else:
                    features[shift : shift + n_max * n_max] = (
                        prel * dots.ravel()
                    )
                    shift += n_max * n_max
    return features


def _soap_gto_features(
    positions: np.ndarray,
    species_indices: np.ndarray,
    *,
    n_species: int,
    r_cut: float,
    n_max: int,
    l_max: int,
    sigma: float,
    center_indices: Sequence[int],
) -> np.ndarray:
    """Compute per-center SOAP vectors (DScribe 2.1.2 GTO parity)."""
    n_centers = len(center_indices)
    n_features = (
        (n_species * n_max) * (n_species * n_max + 1) // 2 * (l_max + 1)
    )
    result = np.zeros((n_centers, n_features), dtype=np.float64)

    alphas, betas = _gto_basis(r_cut, n_max, l_max)
    two_sigma_sq = 2.0 * sigma**2
    inv_eta = 2.0 * sigma**2
    inv_eta_3_2 = np.sqrt(inv_eta**3)

    a_oa = np.zeros((l_max + 1, n_max), dtype=np.float64)
    b_oa = np.zeros((l_max + 1, n_max, n_max), dtype=np.float64)
    for ang_l in range(l_max + 1):
        for k in range(n_max):
            denom = max(1.0 + two_sigma_sq * float(alphas[ang_l, k]), _EPS)
            one_over = 1.0 / denom
            a_oa[ang_l, k] = -alphas[ang_l, k] * one_over
            scale = np.sqrt(one_over) * (one_over ** (ang_l + 1))
            b_oa[ang_l, :, k] = inv_eta_3_2 * betas[ang_l, :, k] * scale

    # DScribe 2.1.2 includes neighbors out to r_cut + cutoff_padding, where
    # cutoff_padding = sigma * sqrt(-2 ln threshold). Analytic GTO integrals
    # have no hard r_cut filter; the padded neighbor list is the cutoff.
    cutoff_padding = float(sigma) * np.sqrt(
        -2.0 * np.log(_GTO_DECAY_THRESHOLD)
    )
    neighbor_cut_sq = (float(r_cut) + cutoff_padding) ** 2

    for ci, center_atom in enumerate(center_indices):
        center = positions[center_atom]
        deltas = positions - center
        r2_all = np.sum(deltas * deltas, axis=1)
        neighbor_idx = np.flatnonzero(r2_all <= neighbor_cut_sq)
        # The center itself has r=0, so this is empty only for empty geometries
        # that are already rejected upstream; keep the guard for safety.
        if neighbor_idx.size == 0:  # pragma: no cover
            continue

        x = deltas[neighbor_idx, 0]
        y = deltas[neighbor_idx, 1]
        z = deltas[neighbor_idx, 2]
        r2 = r2_all[neighbor_idx]
        neigh_species = species_indices[neighbor_idx]

        harmonics = _solid_harmonics_polynomial(x, y, z, l_max)
        cnnd = np.zeros((n_species, n_max, (l_max + 1) ** 2), dtype=np.float64)

        for type_j in range(n_species):
            mask = neigh_species == type_j
            if not np.any(mask):
                continue

            xj = x[mask]
            yj = y[mask]
            zj = z[mask]
            r2j = r2[mask]
            harm_j = harmonics[mask]

            for ang_l in range(l_max + 1):
                k_l = _RADIAL_ANGULAR_SCALE.get(ang_l, 1.0)
                # pre_exp: (n_max, n_neighbors) — one radial Gaussian per k
                pre_exp = k_l * np.exp(a_oa[ang_l, :, None] * r2j[None, :])

                if ang_l == 0:
                    # sums: (n_max, 1)
                    sums = pre_exp.sum(axis=1, keepdims=True)
                    m_start, m_end = 0, 1
                elif ang_l == 1:
                    # DScribe l=1 m-order: z, x, y -> indices 1, 2, 3
                    sums = np.column_stack(
                        (
                            pre_exp @ zj,
                            pre_exp @ xj,
                            pre_exp @ yj,
                        )
                    )
                    m_start, m_end = 1, 4
                else:
                    m_start = ang_l * ang_l
                    m_end = (ang_l + 1) ** 2
                    # (n_max, n_neigh) @ (n_neigh, 2l+1) -> (n_max, 2l+1)
                    sums = pre_exp @ harm_j[:, m_start:m_end]

                # b_oa[l]: (n, k); sums: (k, m) -> cnnd n-block += b @ sums
                cnnd[type_j, :, m_start:m_end] += b_oa[ang_l] @ sums

        result[ci] = _power_spectrum(cnnd, n_species, n_max, l_max)

    return result


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
            Must be greater than 1 Å for the default GTO basis. Matching
            DScribe 2.1.2, neighbors within
            ``r_cut + sigma * sqrt(-2 * ln(1e-3))`` also contribute
            (Gaussian density tails).
        n_max: Number of radial basis functions.
        l_max: Maximum angular momentum (spherical harmonics degree).
            Must satisfy ``0 <= l_max <= 20`` (DScribe GTO parity limit).
        sigma: Gaussian width in Å (standard deviation) used to expand
            densities.
        species: Explicit elemental species basis for the SOAP feature space.
            When ``None``, a sorted unique list of the molecule's own elements
            is used (convenient for a single structure). For comparable
            descriptors across a dataset, pass a shared species list that
            covers every element present in the dataset. Symbols are
            case-normalized (``fe`` → ``Fe``).
        centers: Optional 1-based atom indices on which to evaluate SOAP.
            When ``None``, every atom is used as a center. Order is preserved
            and duplicate indices are kept (they overweight ``"mean"`` /
            ``"sum"`` aggregations).
        aggregation: ``None`` returns local (per-center) vectors with shape
            ``(n_centers, n_features)``. ``"mean"`` / ``"sum"`` are computed
            post-hoc over local power spectra (outer-average / extensive-sum
            equivalent). Mean is preferred for size-comparable fingerprints;
            sum scales with the number of centers.

    Returns:
        Dense ``float64`` NumPy array of SOAP features.

    Raises:
        ValueError: For invalid geometry, hyperparameters, centers, species,
            aggregation, active periodic boundary conditions, non-empty
            translation vectors, or unrecognized elemental symbols.
        TypeError: For invalid center / hyperparameter types.
    """
    if not isinstance(molecule, Molecule):
        raise TypeError(
            f"molecule must be a Molecule instance; got {type(molecule)}."
        )
    if _has_active_pbc(molecule) or _has_translation_vectors(molecule):
        raise ValueError(
            "Periodic SOAP is not supported. The molecule has active PBC "
            "conditions and/or translation vectors (cell). Use a finite "
            "(non-periodic) molecular geometry."
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
    if center_indices is None:
        center_indices = list(range(len(symbols)))

    species_to_index = _dscribe_species_index_map(species_list)
    n_species = len(species_to_index)
    species_indices = np.array(
        [species_to_index[sym] for sym in symbols], dtype=np.int64
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

    features = _soap_gto_features(
        positions,
        species_indices,
        n_species=n_species,
        r_cut=float(r_cut),
        n_max=int(n_max),
        l_max=int(l_max),
        sigma=float(sigma),
        center_indices=center_indices,
    )

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
