"""The one place this host decides which atoms are adjacent.

A bond is not a distance, and this module does not pretend otherwise. What
it owns is a *named, declared convention* for turning coordinates into an
adjacency relation, together with the margin by which every pair cleared
or missed that convention -- so a reader can tell a structural fact from a
threshold artifact without re-deriving the threshold.

Why it exists
-------------
Five perceivers answered "which atoms are bonded" with four calibrations.
The analysis plane shrank every X-H tolerance to 0.05 A on top of ASE
covalent radii, which put the bond/no-bond line at 1.120 A for C-H,
1.020 A for O-H and **0.740 A for H-H**. Measured consequences, all on
real structures:

- H2 at its experimental 0.7414 A had no perceived bond;
- a converged B3LYP/def2-SVP formaldehyde (C-H 1.1215 A) had *neither*
  C-H bond, while the same molecule at def2-TZVP (1.1078 A) had both, so
  the host's molecular graph changed with the basis set (live,
  ``sm1-formaldehyde``, 2026-09-11);
- a strongly hydrogen-bonded O-H at 1.030 A had no bond, in a laboratory
  that runs aqueous pKa cycles;
- SiH4 at its experimental geometry was reported as five separated
  pieces.

The cause was structural rather than a tuning accident. Every X-H
equilibrium bond is *longer* than the sum of covalent radii, and hydrogen
is the worst case: H2 exceeds it by 0.121 A. Hydrogen's radius
under-describes its bonds more than any other element's, so X-H needs the
**widest** tolerance and the code gave it the narrowest.

The form, and why it has two regimes
------------------------------------
An additive buffer is scale-inconsistent: 0.05 A is 8.1% of the H-H radius
sum and 3.3% of C-C, so a fixed addend is tightest exactly where relative
bond-length variation is largest. A pure multiplicative factor has the
opposite failure -- 30% of a large radius sum is a large absolute slack --
and it was measured against this repository's own corpus before being
rejected: the admissible single factor is only ``(1.1958, 1.2473)``, its
lower bound set by H2 and its upper bound by a **C...Ti contact at
2.944 A** in ``tests/data/StructuresTests/conformers``. A single factor of
1.30 calls that contact a bond.

So the declared convention takes whichever of the two is tighter::

    cutoff(A, B) = min(f * (r_A + r_B), (r_A + r_B) + cap)

with ``f = 1.30`` and ``cap = 0.45`` A. The factor governs pairs whose
radii sum below ``cap / (f - 1) = 1.5`` A -- which is every pair involving
hydrogen -- and the cap governs the rest. Verified against the corpus: it
drops **no** pair the conservative additive-0.3 convention calls bonded,
and adds exactly one, an O...C at 1.793 A inside a transition state, which
is a genuine partial bond and now arrives with its margin.

What stays outside
------------------
Bond *order*, aromaticity and valence saturation are not here and are not
derived from this cutoff. A distance rule cannot see where electrons are,
and a number derived from the cutoff moves whenever the cutoff moves --
which is how one earlier repair traded a false connectivity for a false
bond order without anything noticing. Those belong to the scientist, or to
a program that computes them.

Nor does this module decide the cases that have no distance answer.
[FHF]-, B-H-B bridges, agostic interactions and every proton-transfer
saddle sit near the line by their nature. For those the margin is the
answer: the host reports what its convention said *and how narrowly*, and
the reader draws the chemical conclusion.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

from chemsmart.io.molecules import get_covalent_radius

#: The convention's identity, carried on every receipt derived from it, so
#: a stored adjacency says which rule produced it and a policy change is
#: visible in provenance rather than silently retrofitted onto old
#: records.
BOND_PERCEPTION_POLICY_ID = "chemsmart.bond.radii-two-regime.v1"

#: Multiplicative regime: scale-consistent, and it governs every pair
#: whose radii sum is small -- which is every pair involving hydrogen.
BOND_PERCEPTION_FACTOR = 1.30

#: Additive ceiling: caps the absolute slack a large radius sum would
#: otherwise receive. Set from the measured corpus bound (a C...Ti contact
#: at 2.944 A must not be a bond) rather than chosen.
BOND_PERCEPTION_CAP_ANGSTROM = 0.45

#: Beyond this, atoms are not examined at all. Purely an arithmetic bound:
#: the largest cutoff any declared pair can have is far below it.
_MAX_EXAMINED_ANGSTROM = 6.0


@dataclass(frozen=True)
class PerceivedPairV1:
    """One pair, its distance, the cutoff applied, and the margin.

    ``margin`` is ``cutoff - distance``: positive for a perceived
    adjacency and negative for a rejected one, so its magnitude says how
    far the decision was from flipping. A reader that ignores it gets the
    same boolean as before; a reader that reads it can tell a structural
    fact from a threshold artifact.
    """

    schema_version: str
    first_index: int
    second_index: int
    first_symbol: str
    second_symbol: str
    distance_angstrom: float
    cutoff_angstrom: float
    margin_angstrom: float
    adjacent: bool
    policy_id: str

    @property
    def relative_margin(self) -> float:
        """The margin as a fraction of the cutoff, for scale-free reading."""

        return (
            self.margin_angstrom / self.cutoff_angstrom
            if self.cutoff_angstrom
            else 0.0
        )


def bond_cutoff(first: str, second: str) -> float:
    """The declared adjacency cutoff for one element pair, in angstrom."""

    total = get_covalent_radius(str(first)) + get_covalent_radius(str(second))
    return min(
        BOND_PERCEPTION_FACTOR * total,
        total + BOND_PERCEPTION_CAP_ANGSTROM,
    )


def binding_regime(first: str, second: str) -> str:
    """Which of the two regimes decided this pair's cutoff."""

    total = get_covalent_radius(str(first)) + get_covalent_radius(str(second))
    factor_cut = BOND_PERCEPTION_FACTOR * total
    return (
        "factor"
        if factor_cut <= total + BOND_PERCEPTION_CAP_ANGSTROM
        else "cap"
    )


def perceive_pairs(
    symbols: Sequence[str],
    positions: Any,
    *,
    include_rejected: bool = False,
) -> tuple[PerceivedPairV1, ...]:
    """Every examined pair, with its distance, cutoff and margin.

    With ``include_rejected`` the result also carries pairs the convention
    rejected, which is what a reader needs to see a near miss. Pairs
    beyond an arithmetic screening distance are not examined at all and
    never appear: absence there means "not examined", not "rejected".
    """

    labels = [str(item) for item in symbols]
    coordinates = np.asarray(positions, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(
            "adjacency perception needs an (n, 3) coordinate array; got "
            f"{coordinates.shape}"
        )
    if len(labels) != coordinates.shape[0]:
        raise ValueError(
            f"{len(labels)} symbol(s) against {coordinates.shape[0]} "
            "coordinate row(s)"
        )
    found: list[PerceivedPairV1] = []
    for i, j in itertools.combinations(range(len(labels)), 2):
        distance = float(np.linalg.norm(coordinates[i] - coordinates[j]))
        if distance > _MAX_EXAMINED_ANGSTROM:
            continue
        cutoff = bond_cutoff(labels[i], labels[j])
        adjacent = distance < cutoff
        if not adjacent and not include_rejected:
            continue
        found.append(
            PerceivedPairV1(
                schema_version="chemsmart.perceived-pair.v1",
                first_index=i,
                second_index=j,
                first_symbol=labels[i],
                second_symbol=labels[j],
                distance_angstrom=distance,
                cutoff_angstrom=cutoff,
                margin_angstrom=cutoff - distance,
                adjacent=adjacent,
                policy_id=BOND_PERCEPTION_POLICY_ID,
            )
        )
    return tuple(found)


def adjacency_matrix(
    symbols: Sequence[str], positions: Any
) -> list[list[int]]:
    """The binary adjacency matrix in source atom order."""

    size = len(list(symbols))
    matrix = [[0 for _ in range(size)] for _ in range(size)]
    for pair in perceive_pairs(symbols, positions):
        matrix[pair.first_index][pair.second_index] = 1
        matrix[pair.second_index][pair.first_index] = 1
    return matrix


def adjacency_graph(symbols: Sequence[str], positions: Any) -> Any:
    """A zero-based ``networkx`` graph of the perceived adjacency."""

    import networkx as nx

    graph = nx.Graph()
    graph.add_nodes_from(range(len(list(symbols))))
    for pair in perceive_pairs(symbols, positions):
        graph.add_edge(pair.first_index, pair.second_index)
    return graph


def molecule_adjacency_matrix(molecule: Any) -> list[list[int]]:
    """Convenience over a ``Molecule``-shaped object."""

    return adjacency_matrix(molecule.chemical_symbols, molecule.positions)


def molecule_adjacency_graph(molecule: Any) -> Any:
    """Convenience over a ``Molecule``-shaped object."""

    return adjacency_graph(molecule.chemical_symbols, molecule.positions)


def separated_piece_count(molecule: Any) -> int:
    """How many connected pieces the perceived adjacency has."""

    import networkx as nx

    return int(
        nx.number_connected_components(molecule_adjacency_graph(molecule))
    )


def margin_observations(molecule: Any, *, nearest: int = 4) -> tuple[str, ...]:
    """The pairs closest to the line, as lines a reader can act on.

    The host says what its convention decided and how narrowly. It does
    not say whether the chemistry agrees -- that is the reader's call, and
    the whole reason the margin travels.
    """

    pairs = perceive_pairs(
        molecule.chemical_symbols, molecule.positions, include_rejected=True
    )
    if not pairs:
        return ()
    ranked = sorted(pairs, key=lambda item: abs(item.margin_angstrom))
    lines = [
        f"adjacency policy {BOND_PERCEPTION_POLICY_ID}: "
        f"cutoff = min({BOND_PERCEPTION_FACTOR:g} x sum of covalent radii, "
        f"sum + {BOND_PERCEPTION_CAP_ANGSTROM:g} A)"
    ]
    for pair in ranked[: max(0, int(nearest))]:
        verdict = "adjacent" if pair.adjacent else "not adjacent"
        lines.append(
            f"  {pair.first_symbol}{pair.first_index + 1}-"
            f"{pair.second_symbol}{pair.second_index + 1} "
            f"{pair.distance_angstrom:.4f} A vs cutoff "
            f"{pair.cutoff_angstrom:.4f} A: {verdict} by "
            f"{abs(pair.margin_angstrom):.4f} A "
            f"({abs(pair.relative_margin) * 100:.1f}%)"
        )
    return tuple(lines)


def declared_policy() -> dict[str, Any]:
    """The convention as a record, for the policy registry and receipts."""

    return {
        "policy_id": BOND_PERCEPTION_POLICY_ID,
        "question": "which atoms are adjacent",
        "form": "min(factor * (r_a + r_b), (r_a + r_b) + cap)",
        "factor": BOND_PERCEPTION_FACTOR,
        "cap_angstrom": BOND_PERCEPTION_CAP_ANGSTROM,
        "radii_source": "ase.data.covalent_radii",
        "regime_crossover_radii_sum_angstrom": (
            BOND_PERCEPTION_CAP_ANGSTROM / (BOND_PERCEPTION_FACTOR - 1.0)
        ),
        "delivers": "binary adjacency plus per-pair distance, cutoff, margin",
        "excludes": "bond order, aromaticity, valence saturation",
    }


__all__ = [
    "BOND_PERCEPTION_CAP_ANGSTROM",
    "BOND_PERCEPTION_FACTOR",
    "BOND_PERCEPTION_POLICY_ID",
    "PerceivedPairV1",
    "adjacency_graph",
    "adjacency_matrix",
    "binding_regime",
    "bond_cutoff",
    "declared_policy",
    "margin_observations",
    "molecule_adjacency_graph",
    "molecule_adjacency_matrix",
    "perceive_pairs",
    "separated_piece_count",
]
