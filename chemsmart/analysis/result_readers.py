"""One result reader per program, behind a single semantic vocabulary.

Typed quantity extraction was written against PySCF because ChemSmart controls
its structured HDF5 schema, and its contracts were deliberately left
program-neutral "so that Gaussian, ORCA, and xTB readers can be registered
later without changing the expression evaluator".  Until they are, every
observable from those programs has to be read by eye out of a log, which is
exactly the channel the project-YAML hub exists to close: a number a model
typed is not a number ChemSmart measured.

Registering them here needs no new parser.  Each program's reader already
exposes the same few accessors under the same names -- ``energies``,
``gibbs_free_energy``, ``vibrational_frequencies``, ``molecule`` -- so a
selector maps to an attribute rather than to a per-program branch ladder.  A
program is supported exactly when it appears in ``RESULT_READERS``; nothing
else in the extraction path needs to know which programs exist.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from chemsmart.io.molecules.perception import (
    BOND_PERCEPTION_POLICY_ID,
    molecule_adjacency_matrix,
    perceive_pairs,
)

__all__ = [
    "DECLARED_SELECTORS",
    "PRINTED_THERMOCHEMISTRY_CONVENTIONS",
    "RESULT_READERS",
    "MissingQuantityError",
    "ResultReaderV1",
    "SELECTOR_UNITS",
    "reader_for",
    "registered_reader_jobtype_selectors",
    "registered_reader_programs",
    "registered_reader_selectors",
    "atom_resolved_selector_metadata",
    "merge_selector_declarations",
]


class MissingQuantityError(LookupError):
    """The run did not produce this quantity.

    Distinct from an unsupported selector: a single point genuinely has no
    Gibbs energy, and reporting that honestly is what lets a caller decide to
    run a frequency job rather than treat a parse failure as a bug.
    """


#: The canonical unit each semantic selector is reported in.  A reader that
#: cannot produce a quantity in its declared unit must fail rather than guess:
#: a silently mis-scaled energy is worse than an absent one.
SELECTOR_UNITS = {
    "absorption_wavelengths": "nm",
    "energy": "Eh",
    "energies": "Eh",
    "entropy_times_temperature": "Eh",
    "excitation_energies": "eV",
    "singlet_excitation_energies": "eV",
    "triplet_excitation_energies": "eV",
    "excited_state_indices": "",
    "excited_state_manifold_roots": "",
    "excited_state_multiplicities": "",
    "excited_state_labels": "",
    "excited_state_spin_square": "",
    "gibbs_free_energy": "Eh",
    "oscillator_strengths": "",
    "singlet_oscillator_strengths": "",
    "triplet_oscillator_strengths": "",
    # Per-root transition dipole vectors; PySCF returns atomic units and
    # the driver converts with the CODATA factor it records.
    "transition_dipole_moments": "Debye",
    # 1 per converged root, 0 otherwise; the followed root of an
    # excited-surface optimisation as a 1-based index.
    "excited_state_converged": "",
    "excited_state_followed_root": "",
    # The coupled-cluster components beside the final method's whole
    # correlation energy.
    "ccsd_correlation_energy": "Eh",
    "triples_correction": "Eh",
    "vibrational_frequencies": "cm^-1",
    "ir_intensities": "km/mol",
    "vibrational_mode_atom_participation": "",
    "vibrational_mode_degeneracy_group": "",
    "scan_energies": "Eh",
    # A scanned coordinate is a distance for a bond and an angle for a torsion,
    # so a single declared unit here would be false for half of all scans.  The
    # values carry the kind the run actually drove, which the scan coordinate
    # record states.
    "scan_coordinate_values": "",
    # The position of each point in the surface, so a session can name the
    # point it chose rather than describe it.
    "scan_point_indices": "",
    "scan_steps_reached": "1",
    "scan_steps_planned": "1",
    "vpt2_harmonic_frequencies": "cm^-1",
    "vpt2_fundamental_frequencies": "cm^-1",
    "vpt2_zero_point_rovibrational_energy": "cm^-1",
    "positions": "Angstrom",
    "reached_positions": "Angstrom",
    "supplied_positions": "Angstrom",
    "connectivity": "",
    "symbols": "",
    "charge": "",
    "multiplicity": "",
    "scf_energy": "Eh",
    "reference_energy": "Eh",
    "correlation_energy": "Eh",
    "dispersion_energy": "Eh",
    "auxiliary_basis": "",
    "auxiliary_basis_role": "",
    "dipole_moment": "Debye",
    "dipole_moment_magnitude": "Debye",
    "homo": "eV",
    "lumo": "eV",
    "gap": "eV",
    # An unrestricted reference has two frontier pairs, and which one answers
    # the question is a scientific choice rather than a parser detail, so the
    # channels stay separately nameable instead of being collapsed.
    "alpha_homo": "eV",
    "alpha_lumo": "eV",
    "beta_homo": "eV",
    "beta_lumo": "eV",
    "spin_square": "",
    "spin_square_after_annihilation": "",
    "spin_square_target": "",
    "spin_square_deviation": "",
    "effective_multiplicity": "",
    "wavefunction_stability_verdict": "",
    "wavefunction_stability_history": "",
    "trajectory_frame_count": "",
    # The energy at every frame of a path, and the spectrum at its start.
    "trajectory_energies": "Eh",
    "trajectory_start_frequencies": "cm^-1",
    "trajectory_start_positions": "Angstrom",
    "trajectory_end_positions": "Angstrom",
    "trajectory_start_connectivity": "",
    "trajectory_end_connectivity": "",
    "trajectory_connectivity_changed": "",
    "irc_direction": "",
    "solvation_model": "",
    "solvation_electrostatic_energy": "Eh",
    "solvation_nonelectrostatic_energy": "Eh",
    "solvation_cavity_surface_area": "angstrom^2",
    "mulliken_atomic_charges": "e",
    "hirshfeld_atomic_charges": "e",
    "loewdin_atomic_charges": "e",
    "xtb_scc_atomic_charges": "e",
    "mulliken_atomic_spin_populations": "1",
    "loewdin_atomic_spin_populations": "1",
    "functional": "",
    "method": "",
    "ab_initio": "",
    "basis": "",
    "converged": "1",
    "irc_converged": "1",
    "solvent": "",
    "surface_id": "",
    "wiberg_bond_orders": "1",
}


# A selector is the semantic identity of an atom-resolved population.  This
# small registry is deliberately descriptive: it tells the Agent what the
# host extracted and how vector indices align, but never turns a Mulliken,
# Loewdin, Hirshfeld, or xTB SCC population into a common "atomic charge".
_ATOM_RESOLVED_SELECTOR_METADATA: Mapping[str, Mapping[str, str]] = {
    "mulliken_atomic_charges": {
        "semantic_quantity": "atomic_partial_charge",
        "population_scheme": "Mulliken",
        "atom_order": "zero-based molecular atom order",
    },
    "loewdin_atomic_charges": {
        "semantic_quantity": "atomic_partial_charge",
        "population_scheme": "Loewdin",
        "atom_order": "zero-based molecular atom order",
    },
    "hirshfeld_atomic_charges": {
        "semantic_quantity": "atomic_partial_charge",
        "population_scheme": "Hirshfeld",
        "atom_order": "zero-based molecular atom order",
    },
    "xtb_scc_atomic_charges": {
        "semantic_quantity": "atomic_partial_charge",
        "population_scheme": "xTB self-consistent-charge population",
        "atom_order": "zero-based molecular atom order",
    },
    "mulliken_atomic_spin_populations": {
        "semantic_quantity": "atomic_spin_population",
        "population_scheme": "Mulliken",
        "atom_order": "zero-based molecular atom order",
    },
    "loewdin_atomic_spin_populations": {
        "semantic_quantity": "atomic_spin_population",
        "population_scheme": "Loewdin",
        "atom_order": "zero-based molecular atom order",
    },
    "wiberg_bond_orders": {
        "semantic_quantity": "bond_order",
        "population_scheme": "Wiberg",
        "atom_order": "zero-based molecular atom order",
        "data_shape": "rows of [atom_i, atom_j, wiberg_bond_order]",
        "sparsity": (
            "the native xTB sidecar is thresholded; omitted pairs have no "
            "reported Wiberg value and are not zero"
        ),
    },
}


def atom_resolved_selector_metadata(selector: str) -> dict[str, str]:
    """Return declared population semantics, or an explicit empty mapping."""

    return dict(_ATOM_RESOLVED_SELECTOR_METADATA.get(selector, {}))


def _last_energy(output: Any) -> float:
    """Return the converged energy, preferring an explicit final value."""

    final = getattr(output, "final_energy", None)
    if final is not None:
        return float(final)
    energies = getattr(output, "energies", None)
    if not energies:
        raise ValueError("result exposes no energy")
    return float(energies[-1])


def _positions(output: Any) -> list[list[float]]:
    molecule = output.molecule
    return [[float(value) for value in row] for row in molecule.positions]


#: Frequencies closer than this are treated as one degenerate set.  Two
#: programs print frequencies to two decimals, so a tighter window would
#: split a genuinely degenerate pair on rounding alone; a much looser one
#: would merge distinct modes of a floppy molecule.
MODE_DEGENERACY_TOLERANCE_CM1 = 1.0


def _mode_displacement_matrices(output: Any) -> list[Any]:
    """The program's own per-mode displacement matrices, index-aligned.

    Every registered reader reaches this through ``vibrational_modes``, and
    the pairing with ``vibrational_frequencies`` is the invariant that makes
    a mode index mean one thing: mode k is the mode whose frequency is
    frequency k.  A result whose two lists disagree is refused rather than
    silently mispaired.
    """

    import numpy as np

    # PySCF hands the modes back as one numpy array (nmode, natm, 3) and
    # the log readers as a list of matrices; ``or []`` on an array is an
    # ambiguous truth value, and that single expression hid the whole
    # quantity on every real PySCF Hessian behind a swallowed ValueError.
    modes = getattr(output, "vibrational_modes", None)
    modes = [] if modes is None else list(modes)
    frequencies = getattr(output, "vibrational_frequencies", None)
    frequencies = [] if frequencies is None else list(frequencies)
    if not modes:
        raise MissingQuantityError(
            "result establishes no vibrational normal modes"
        )
    if len(modes) != len(frequencies):
        raise MissingQuantityError(
            f"result reports {len(frequencies)} vibrational frequencies but "
            f"{len(modes)} normal modes, so a mode index would not name the "
            "mode its frequency names"
        )
    matrices = []
    for mode in modes:
        matrix = np.asarray(mode, dtype=float)
        if matrix.ndim != 2 or matrix.shape[1] != 3:
            raise MissingQuantityError(
                "a vibrational normal mode must carry three displacement "
                f"components per atom; got shape {matrix.shape}"
            )
        matrices.append(matrix)
    return matrices


def _vibrational_mode_atom_participation(output: Any) -> list[list[float]]:
    """Per-atom share of each mode's squared displacement.

    ``s_i = ||x_i||^2 / sum_j ||x_j||^2`` for atom ``i`` of one mode, so each
    row sums to one and the table reads "how much of mode k is atom i".

    This is the one quantity the four supported programs agree on.  They do
    not agree on the vector: ORCA, Gaussian and xTB print Cartesian
    displacements normalised to unit norm, while PySCF returns the same
    physical displacement scaled by ``1/sqrt(reduced mass)`` in amu^-1/2 --
    a per-mode scalar.  Dividing by the row's own total removes exactly that
    scalar, and with it the arbitrary eigenvector sign and the program's
    choice of coordinate frame, none of which survive into a ratio of
    squared magnitudes.  What it does not remove is each program's atomic
    mass table (Gaussian uses most-abundant-isotope masses where the others
    use standard atomic weights), which perturbs the eigenvector itself by
    roughly a tenth of a percent for C/H/O and about one percent for heavy
    halogens.

    The number is an observation.  Reading "atoms 6, 7 and 8 carry 97% of
    this imaginary mode, and they are the three hydrogens of a methyl group"
    is chemistry, and it belongs to the scientist, not to this function.
    """

    participation = []
    for matrix in _mode_displacement_matrices(output):
        squared = (matrix**2).sum(axis=1)
        total = float(squared.sum())
        if total <= 0.0:
            raise MissingQuantityError(
                "a vibrational normal mode has zero displacement, so no "
                "atom carries a share of it"
            )
        participation.append([float(value) for value in squared / total])
    return participation


def _vibrational_mode_degeneracy_group(output: Any) -> list[int]:
    """Which modes share a frequency, as 1-based group labels.

    Within a degenerate set the individual eigenvectors are an arbitrary
    basis: the two bending modes of a linear triatomic can be printed as
    pure-x and pure-y, or as any rotation of that pair, and the program's
    choice carries no physics.  Per-mode participation is therefore
    ill-posed inside such a set -- only the sum over the set is meaningful.

    Modes whose frequencies lie within
    :data:`MODE_DEGENERACY_TOLERANCE_CM1` of the group's first member share
    a label, so a reader can see that mode k has company before drawing a
    conclusion about which atoms move in it.  Singleton labels are the
    ordinary case; the grouping states a fact and judges nothing.
    """

    frequencies = getattr(output, "vibrational_frequencies", None)
    frequencies = [] if frequencies is None else list(frequencies)
    if not frequencies:
        raise MissingQuantityError(
            "result establishes no vibrational frequencies"
        )
    labels: list[int] = []
    group = 0
    anchor: float | None = None
    for frequency in (float(item) for item in frequencies):
        if (
            anchor is None
            or abs(frequency - anchor) > MODE_DEGENERACY_TOLERANCE_CM1
        ):
            group += 1
            anchor = frequency
        labels.append(group)
    return labels


def _symbols(output: Any) -> list[str]:
    return [str(item) for item in output.molecule.chemical_symbols]


_ATOM_LABEL = re.compile(r"^([A-Za-z]{1,3})(\d+)$")


def _orca_spin_populations(output: Any, *, quantity: str) -> list[float]:
    """Per-atom spin populations of an open-shell ORCA result, in order.

    Refused on a closed shell, where the block has one column and there
    is no spin to partition; checked against the sum ORCA prints beneath
    the block, so a dropped or duplicated atom cannot pass as a
    population. That sum is 2 Ms of the determinant the run ended on:
    2S = multiplicity - 1 for the state the coordinate line names, and
    2 Ms of the state ORCA says it converged to where a FlipSpin or
    BrokenSym request flipped it. Checking every result against the
    coordinate line refused both routes to ino2's antiferromagnetic
    Ni(II)2 state (R11 truth-2, CUHK 2157086).
    """

    labelled = getattr(output, quantity)
    if labelled is None:
        raise MissingQuantityError(
            f"{quantity}: this result printed no spin populations -- the "
            "population block has one column, as for a closed shell"
        )
    values = _per_atom_vector(
        labelled, _orca_symbols(output), quantity=quantity
    )
    try:
        multiplicity = int(output.multiplicity)
    except (TypeError, ValueError):
        multiplicity = None
    final_ms = getattr(output, "broken_symmetry_ms", None)
    expected = basis = None
    if final_ms is not None:
        expected = 2.0 * float(final_ms)
        basis = (
            f"2 Ms for the Ms ORCA converged its broken-symmetry request "
            f"to ({float(final_ms):.1f}) is {expected:.1f}"
        )
    elif multiplicity is not None:
        expected = float(multiplicity - 1)
        basis = f"2S for multiplicity {multiplicity} is {expected:.1f}"
    if expected is not None:
        total = sum(values)
        if abs(total - expected) > 0.05:
            raise MissingQuantityError(
                f"{quantity} sums to {total:.3f} where {basis}; the vector "
                "is not the complete molecule in order"
            )
    return values


def _per_atom_vector(
    labelled: Any, symbols: list[str], *, quantity: str
) -> list[float]:
    """Order an atom-labelled per-atom mapping into molecular order.

    Every program's per-atom parsers hand back a mapping keyed by an atom
    label, and the labelling schemes do not agree: ORCA and Gaussian number
    atoms globally, so carbon dioxide's carbon is ``"C3"``, while xTB counts
    within each element and calls the same atom ``"C1"``.  Reading one
    scheme as the other returns a different atom without any error.

    The typed layer cannot repair the order later either -- freezing a
    mapping sorts it by label, so ``O1, C2, H3`` becomes ``C2, H3, O1``, and
    per-atom values silently stop lining up with the geometry they describe.

    So the label is resolved here, against the molecule's own symbols, and
    every atom must be named exactly once by a label whose element matches
    the atom at that position.  A per-element scheme fails that check rather
    than passing quietly, which is the point: this refuses to guess.
    """

    if not isinstance(labelled, Mapping):
        values = [float(item) for item in labelled]
        if len(values) != len(symbols):
            raise MissingQuantityError(
                f"{quantity} carries {len(values)} values for "
                f"{len(symbols)} atoms"
            )
        return values
    if len(labelled) != len(symbols):
        raise MissingQuantityError(
            f"{quantity} carries {len(labelled)} values for "
            f"{len(symbols)} atoms"
        )
    ordered: list[float | None] = [None] * len(symbols)
    for label, value in labelled.items():
        match = _ATOM_LABEL.match(str(label).strip())
        if match is None:
            raise MissingQuantityError(
                f"{quantity} carries an unreadable atom label {label!r}"
            )
        element, index_text = match.groups()
        index = int(index_text)
        if not 1 <= index <= len(symbols):
            raise MissingQuantityError(
                f"{quantity} names atom {index} of {len(symbols)}"
            )
        expected = symbols[index - 1]
        if expected.strip().lower() != element.strip().lower():
            raise MissingQuantityError(
                f"{quantity} labels atom {index} {element!r} while the "
                f"geometry has {expected!r}; the atom-label scheme does "
                "not match this molecule's atom order"
            )
        if ordered[index - 1] is not None:
            raise MissingQuantityError(
                f"{quantity} names atom {index} more than once"
            )
        ordered[index - 1] = float(value)
    if any(item is None for item in ordered):
        raise MissingQuantityError(f"{quantity} does not name every atom")
    return [float(item) for item in ordered]


def _required_record_values(
    records: list[dict[str, Any]], key: str, description: str
) -> list[Any]:
    if not records or any(record.get(key) is None for record in records):
        raise MissingQuantityError(
            f"result does not establish {description} for every excited state"
        )
    return [record[key] for record in records]


def _excited_spin_squares(records, response_method, manifold):
    """Each root's <S^2>, where the quantity has one meaning.

    A spin-adapted root's <S^2> is exact (0 or 2) and a Tamm-Dancoff root
    of an open-shell reference is a CIS-like wavefunction whose <S^2> every
    program computes alike -- Gaussian and ORCA agree to the printed digit
    on the allyl radical (0.756/0.756 for D1, CUHK Slurm 2150194).  A full
    TD-DFT root of an open-shell reference is not a wavefunction, and the
    two programs print two different approximations for it (0.713 in
    Gaussian, 0.801 in ORCA for the same D1): under one name that would be
    two quantities, so none is served there.
    """

    if response_method == "tddft" and manifold == "unrestricted":
        raise MissingQuantityError(
            "a full TD-DFT root of an open-shell reference has no unique "
            "<S^2>: Gaussian and ORCA print two different approximations "
            "for the same root (allyl D1: 0.713 vs 0.801); a Tamm-Dancoff "
            "run (response_method: tda) gives the CIS-like <S^2> both "
            "programs agree on"
        )
    return [
        float(item)
        for item in _required_record_values(
            records, "spin_square", "a printed excited-state <S^2>"
        )
    ]


def _frontier_orbital(base: str, offset: float) -> str:
    step = int(round(offset))
    return base if step == 0 else f"{base}{step:+d}"


def _dominant_excitation_values(records, what: str) -> list[Any]:
    """Each root's largest single excitation, in one program-neutral form.

    Every reader records ``(occupied_offset, virtual_offset, weight,
    channel)`` from its own program's print: the occupied orbital counted
    from the HOMO and the virtual from the LUMO of its own spin, and the
    spin channel of an unrestricted root.  ``labels`` spells it
    ``HOMO-1 -> LUMO`` (``beta HOMO -> LUMO`` on an open shell) and
    ``weights`` gives the program's own weight of that excitation (2c^2 of
    a spin-adapted coefficient, c^2 of an unrestricted one, ORCA's printed
    weight).  A root's character is what it is made of, so it is how a
    root is recognised in another program's list, where its position may
    differ.
    """

    values = []
    for record in records:
        dominant = record.get("dominant_excitation")
        if dominant is None:
            raise MissingQuantityError(
                "result does not print the excitations of every root"
            )
        occupied, virtual, weight, channel = dominant
        if what == "weights":
            values.append(float(weight))
            continue
        spin = {1: "alpha ", -1: "beta "}.get(int(round(channel)), "")
        values.append(
            f"{spin}{_frontier_orbital('HOMO', occupied)} -> "
            f"{_frontier_orbital('LUMO', virtual)}"
        )
    return values


def _orca_served_records(output: Any) -> list[dict[str, Any]]:
    """The roots of the manifold the request named, ranked by energy.

    ORCA has no triplet-only solve, so a ``triplet`` request is written as
    ``Triplets true`` and its output also holds the singlet block nobody
    asked for.  The aggregate selectors serve the requested manifold -- as
    Gaussian's and PySCF's triplet runs do, ranked 1..n -- and the
    singlets stay readable by name (``singlet_*``).
    """

    records = list(output.excited_state_records or ())
    if getattr(output, "state_manifold", None) == "triplet":
        triplets = [item for item in records if item["multiplicity"] == 3]
        records = [
            {**item, "state_index": rank}
            for rank, item in enumerate(triplets, start=1)
        ]
    return records


def _orca_excitation_energies(output: Any) -> list[float]:
    records = _orca_served_records(output)
    if records:
        return [float(item["energy_eV"]) for item in records]
    # ORCA 6 spectrum-only fragments remain useful for the legacy aggregate
    # selector, but they cannot support multiplicity-specific selectors.
    return [float(item) for item in output.excitation_energies_eV]


def _orca_absorption_rows(output: Any, records, key: str) -> list[float]:
    """The absorption-table value of each record, in the records' order.

    Each root is paired with the row whose ``N-MA`` label is ORCA's own
    label for it (the ``STATE`` number within its block and the printed
    multiplicity), never by position: ORCA prints the ``STATE`` blocks
    singlets first and the absorption table in energy order, and pairing
    by position served acrolein's bright S2 (6.536 eV) beside f = 0 and
    its T2 (3.196 eV) beside f = 0.381 (CUHK Slurm 2150076).
    """

    rows = list(output.electronic_absorption_transition_records or ())
    by_label = {
        (row["manifold_root"], row["multiplicity"]): row for row in rows
    }
    values = []
    for record in records:
        label = (record["orca_state"], record["orca_multiplicity"])
        row = by_label.get(label)
        if row is None:
            raise MissingQuantityError(
                "ORCA printed no electric-dipole absorption row for its "
                f"state {label[0]}-{label[1]}A"
            )
        values.append(float(row[key]))
    return values


def _orca_absorption_per_state(output: Any, key: str) -> list[float]:
    """One absorption-table value per state, in ``excitation_energies`` order."""

    records = _orca_served_records(output)
    if not records:
        # A spectrum-only fragment: the table is the only record, in the
        # order the energies fall back to as well.
        rows = list(output.electronic_absorption_transition_records or ())
        return [float(row[key]) for row in rows]
    return _orca_absorption_rows(output, records, key)


def _orca_manifold_values(output: Any, multiplicity: int, key: str):
    """One spin block's values in manifold-root order (S_k or T_k)."""

    records = [
        item
        for item in (output.excited_state_records or ())
        if item["multiplicity"] == multiplicity
    ]
    if key == "energy_eV":
        return [float(item["energy_eV"]) for item in records]
    return _orca_absorption_rows(output, records, key)


def _last_spin_square(output: Any, key: str | None = None) -> float:
    """Return the last printed ``<S^2>``, or say why the run has none.

    A spin-restricted closed-shell calculation never prints an expectation
    value, because its wavefunction is an eigenfunction of ``S^2`` by
    construction.  That is a property of the calculation rather than a gap in
    the parser, so name the reason and point at the state evidence the result
    does carry instead of surfacing an index error.
    """

    history = getattr(output, "spin_square_history", None) or ()
    if not history:
        raise MissingQuantityError(
            "this result prints no <S^2> expectation value; a spin-restricted "
            "closed-shell calculation is an eigenfunction of S^2 by "
            "construction, so read 'multiplicity' for its electronic state"
        )
    entry = history[-1]
    return float(entry if key is None else entry[key])


def _spin_square_target(output: Any) -> float:
    multiplicity = getattr(output, "multiplicity", None)
    if not isinstance(multiplicity, int) or multiplicity <= 0:
        raise MissingQuantityError(
            "result does not establish a positive integer multiplicity"
        )
    return (float(multiplicity) ** 2 - 1.0) / 4.0


def _effective_multiplicity(spin_square: float) -> float:
    if spin_square < 0.0:
        raise MissingQuantityError("negative <S^2> cannot define multiplicity")
    return float((1.0 + 4.0 * spin_square) ** 0.5)


#: The determinant a result ran, in one vocabulary for every program:
#: PySCF's own family names (``reference_family``), which Gaussian's
#: ``SCF Done: E(UB3LYP)`` label and ORCA's ``HFTyp`` line are read into.
SCF_REFERENCE_WORDS = ("rhf", "rks", "rohf", "roks", "uhf", "uks")
#: The references that are spin eigenfunctions by construction: a
#: restricted closed shell and a restricted open shell cannot break spin
#: symmetry, whatever guess started them.
SPIN_EIGENFUNCTION_REFERENCES = frozenset({"rhf", "rks", "rohf", "roks"})
#: How far <S**2> must exceed the bound state's S(S+1) before the host
#: calls an unrestricted solution's spin symmetry broken.  A collapsed
#: (spin-symmetric) solution sits at the target to the precision the SCF
#: converged -- 1e-10 to 1e-15 in PySCF and 0.0000 as Gaussian prints it --
#: while the broken-symmetry singlets measured sit at 0.70-1.04 and the
#: weakest real break measured, UHF's pi instability in planar ethylene,
#: at 0.0335 (R10 Q18 O0, CUHK Slurm 2153330).  The number always rides
#: beside the word, which is the host's reading of it and never a verdict
#: on the chemistry: a broken solution may be a real diradical or an
#: artefact of the method, and saying which is the session's.
SPIN_SYMMETRY_BROKEN_THRESHOLD = 0.01


def spin_symmetry_word(reference: str, spin_square: float, target: float):
    """``broken`` or ``unbroken`` for an unrestricted singlet's <S**2>.

    None for a reference that is a spin eigenfunction by construction
    (whose <S**2> says nothing about symmetry breaking), for one no reader
    stated, and for an open shell: there <S**2> above S(S+1) is spin
    contamination of a state that has unpaired spin anyway, which the
    measurement itself reports, not the open-shell singlet this word is
    about.
    """

    word = str(reference or "").strip().lower()
    if word not in SCF_REFERENCE_WORDS or word in (
        SPIN_EIGENFUNCTION_REFERENCES
    ):
        return None
    if abs(float(target)) > 1e-12:
        return None
    deviation = float(spin_square) - float(target)
    return (
        "broken"
        if deviation >= SPIN_SYMMETRY_BROKEN_THRESHOLD
        else ("unbroken")
    )


def spin_symmetry_record(reader: Any, output: Any) -> dict[str, Any] | None:
    """What one result says about its own spin symmetry, or None.

    The reference that ran and whether the program's own record shows the
    broken-symmetry request applied -- both from the result's level, as
    every other organ reads them -- and <S**2> against the bound state's
    S(S+1) through the reader's own selectors, with the host's reading of
    it.  These are the facts a session needs to tell a broken-symmetry
    diradical from a solution that collapsed to the spin-symmetric one.
    None where the level states no reference.
    """

    try:
        level = reader.level_for_output(output) or {}
    except Exception:  # noqa: BLE001 - a reader that cannot say says nothing
        return None
    reference = level.get("reference")
    if reference not in SCF_REFERENCE_WORDS:
        return None
    record: dict[str, Any] = {
        "reference": str(reference),
        "broken_symmetry_requested": bool(level.get("broken_symmetry")),
    }
    followed = getattr(output, "broken_symmetry_record", None)
    if isinstance(followed, Mapping):
        # The program's own stability answer about the restricted solution
        # it started from, where its translation asked (PySCF).
        record["followed_instability"] = {
            name: followed[name]
            for name in (
                "external_stable",
                "external_lowest_eigenvalue",
                "internal_stable",
                "restricted_energy",
                "unrestricted_energy",
                "eigenvalue_unit",
            )
            if name in followed
        }
    try:
        spin_square = float(reader.read(output, "spin_square")[0])
        target = float(reader.read(output, "spin_square_target")[0])
    except Exception:  # noqa: BLE001 - a restricted run prints no <S**2>
        return record
    record["spin_square"] = spin_square
    record["spin_square_target"] = target
    word = spin_symmetry_word(reference, spin_square, target)
    if word is not None:
        record["spin_symmetry"] = word
        record["threshold"] = SPIN_SYMMETRY_BROKEN_THRESHOLD
    return record


def _connectivity_matrix(molecule: Any) -> list[list[int]]:
    """Return geometry-perceived covalent connectivity in source atom order.

    This is intentionally a binary connectivity observation rather than an
    asserted electronic bond order.  IRC endpoints are not fully optimized,
    and covalent-radius perception is appropriate for deciding whether two
    path ends have different molecular graphs while leaving chemical
    interpretation to the scientist.

    Perception is delegated to the one declared convention
    (:mod:`chemsmart.io.molecules.perception`), so this plane and the
    execution plane cannot answer the same question differently. They did:
    this function passed ``adjust_H=True``, which shrank every X-H
    tolerance to 0.05 A and *ignored the buffer argument above it*, while
    the execution plane passed ``adjust_H=False``. The bond/no-bond line
    for C-H therefore sat at 1.120 A here and 1.370 A there, and a
    converged B3LYP/def2-SVP formaldehyde (C-H 1.1215 A) was delivered
    with **neither** C-H bond while the same molecule at def2-TZVP
    (1.1078 A) had both -- a molecular graph that changed with the basis
    set, in a claim (live, ``sm1-formaldehyde``, 2026-09-11).

    The docstring this replaces asserted the opposite of what the code
    did, which is why it read as safe for a year.
    """

    return molecule_adjacency_matrix(molecule)


def _irc_structures(output: Any) -> list[Any]:
    jobtype = str(getattr(output, "jobtype", "") or "").casefold()
    if jobtype not in {"irc", "ircf", "ircr"}:
        raise MissingQuantityError(
            "trajectory selectors require an IRC program result"
        )
    structures = list(getattr(output, "all_structures", ()) or ())
    if len(structures) < 2:
        raise MissingQuantityError(
            "IRC output does not contain at least two parsed trajectory frames"
        )
    symbols = tuple(structures[0].chemical_symbols)
    if any(tuple(item.chemical_symbols) != symbols for item in structures[1:]):
        raise MissingQuantityError(
            "IRC trajectory changes atom identity or order"
        )
    return structures


def _irc_direction(output: Any) -> str:
    jobtype = str(getattr(output, "jobtype", "") or "").casefold()
    directions = {"ircf": "forward", "ircr": "reverse", "irc": "combined"}
    if jobtype not in directions:
        raise MissingQuantityError(
            "result does not establish an IRC direction"
        )
    return directions[jobtype]


def _orca_irc_direction(output: Any) -> str:
    """Return the direction explicitly echoed by an ORCA IRC input block."""

    if str(getattr(output, "jobtype", "") or "").casefold() != "irc":
        raise MissingQuantityError("result is not an ORCA IRC calculation")
    observed = getattr(output, "irc_direction", None)
    if observed is None:
        # Preserve the result-reader protocol for light-weight parser fixtures
        # while production ORCAInput/ORCAOutput objects use the shared typed
        # ``irc_direction`` property above.
        pattern = re.compile(
            r"\bDirection\s+(both|forward|backward|down)\b",
            re.IGNORECASE,
        )
        matches = [
            match.group(1).casefold()
            for line in getattr(output, "contents", ())
            if (match := pattern.search(str(line))) is not None
        ]
        observed = matches[-1] if matches else None
    if observed not in {"backward", "both", "down", "forward"}:
        raise MissingQuantityError(
            "ORCA output does not explicitly establish the IRC direction"
        )
    return observed


def _trajectory_connectivity_changed(output: Any) -> int:
    structures = _irc_structures(output)
    return int(
        _connectivity_matrix(structures[0])
        != _connectivity_matrix(structures[-1])
    )


def _xyz_trajectory_structures(output: Any) -> list[Any]:
    """Return a multi-frame XYZ trajectory without assigning program meaning.

    ORCA writes IRC paths to XYZ sidecars, which already enter Runtime V2 as
    ``geometry_xyz`` artifacts.  The file itself establishes ordered frames,
    but not whether a path is forward, reverse, or the combined IRC.  Generic
    trajectory selectors therefore expose only observations available from
    the registered bytes; ``irc_direction`` remains a program-log selector.
    """

    structures = list(output.get_molecules(index=":", return_list=True) or ())
    if len(structures) < 2:
        raise MissingQuantityError(
            "XYZ artifact does not contain at least two trajectory frames"
        )
    symbols = tuple(structures[0].chemical_symbols)
    if any(tuple(item.chemical_symbols) != symbols for item in structures[1:]):
        raise MissingQuantityError(
            "XYZ trajectory changes atom identity or order"
        )
    return structures


def _xyz_trajectory_connectivity_changed(output: Any) -> int:
    structures = _xyz_trajectory_structures(output)
    return int(
        _connectivity_matrix(structures[0])
        != _connectivity_matrix(structures[-1])
    )


#: The molecular states one completed result can carry.
#:
#: A single ORCA output holds several: the geometry it was handed, the
#: geometry the optimiser reached, the geometry the Hessian was computed
#: at, and every frame in between. Quantities read from different ones
#: are not interchangeable -- on a live unconverged OptTS the supplied
#: and reached structures sat 1.23 A apart with energies 182.2 kcal/mol
#: apart -- so a selector says which it reads and a consumer asks for
#: the role it needs.
#:
#: ``stateless`` is the honest answer for a value no structure changes:
#: the method name, the solvation model, the atom symbols.
STRUCTURAL_STATES = (
    "as_supplied",
    "as_reached",
    "thermochemistry_reference",
    "trajectory_endpoint",
    "scan_point",
    "stateless",
)

#: Which electronic state or method a selector's value belongs to.
#:
#: A structural state identifies a geometry, not a density.  One PySCF
#: artifact of an excited-root optimisation carries the followed root's
#: total energy beside the ground-state reference's dipole, populations and
#: orbital energies, all at one geometry; a correlated result carries a
#: CCSD total beside an HF dipole.  Every hash, unit and geometry check
#: passes while a session reads "the S1 dipole" off a ground-state number,
#: so the axis is declared like the structural one and rendered beside it.
#:
#: ``reference`` -- the SCF the job built on (mean-field properties);
#: ``excited_root`` -- a root of the response manifold, an index at this
#: geometry and never a state identity; ``correlated`` -- the correlated
#: method's own component; ``computed_surface`` -- the surface the job
#: computed on, resolved per artifact (the reference for HF/DFT and for a
#: spectrum, the followed root for an excited-root optimisation, the
#: correlated method otherwise); ``stateless`` -- a value no electronic
#: state changes, and the honest answer for a selector nobody declared.
ELECTRONIC_PROVENANCES = (
    "reference",
    "excited_root",
    "correlated",
    "computed_surface",
    "stateless",
    # A result whose recorded surface names a method family this release
    # has never audited. It is not "reference": that word would be a true
    # sentence about the wrong density, which is exactly how an excited
    # or correlated number came to be read as a ground-state one.
    "unknown",
)

#: Method words that mark a post-SCF surface when a program's ``ab_initio``
#: selector names one; the resolver reads the selector, never a route.
_CORRELATED_METHOD_MARKERS = ("mp2", "ccsd", "cepa", "cisd", "nevpt", "caspt")

#: The fields of a surface identity, in the order the token spells them.
#: A reader that cannot say what a field was leaves it absent, and an
#: absent field is never equal to another absent field: two results
#: agree about a frozen core when both recorded one, and about nothing
#: when neither did.
SURFACE_IDENTITY_FIELDS = (
    "method_family",
    "reference",
    "functional_applied",
    "basis",
    "charge",
    "multiplicity",
    "solvent_model",
    "solvent",
    "dispersion",
    "density_fitting",
    "frozen_core",
    "excited_root",
    "state_manifold",
    "constraints",
)
#: What an absent field spells in the token. ``-`` is a field the reader
#: knows does not apply -- no solvent, no followed root -- and ``?`` is a
#: field it cannot determine. Two ``?`` are never agreement: a frozen
#: core nobody recorded on either side is two unknowns, not a match.
_SURFACE_FIELD_ABSENT = "-"
_SURFACE_FIELD_UNKNOWN = "?"
#: The word a reader writes into a surface field it cannot determine.
SURFACE_UNKNOWN = "unknown"


def surface_token(surface: Mapping[str, Any] | None) -> str:
    """Render a surface identity as one comparable, readable token.

    Two results are on one surface when their tokens are equal *and*
    neither carries an unknown field: the token is what a human reads on
    the level line and what a receipt carries, and the join that decides
    whether a Hessian characterises a geometry asks the reader rather
    than the string.
    """

    if not surface:
        return ""
    parts = []
    for name in SURFACE_IDENTITY_FIELDS:
        value = surface.get(name)
        if value == SURFACE_UNKNOWN:
            parts.append(_SURFACE_FIELD_UNKNOWN)
        elif value is None or value == [] or value == ():
            parts.append(_SURFACE_FIELD_ABSENT)
        elif isinstance(value, bool):
            parts.append("yes" if value else "no")
        elif isinstance(value, (list, tuple)):
            parts.append("+".join(str(item) for item in value))
        else:
            parts.append(str(value).strip().lower())
    return ":".join(parts)


def surfaces_agree(first, second):
    """Whether two results are on one electronic surface.

    ``True`` when every field agrees, ``False`` when one differs, and
    ``None`` when the question cannot be answered -- a result that
    recorded no surface, or a field either reader could not determine.
    A caller must not read ``None`` as agreement: that is precisely the
    reading under which a ground-state Hessian characterised an excited
    minimum, and the caller says "on another surface" or "not
    comparable" rather than staying silent.
    """

    if not first or not second:
        return None
    for name in SURFACE_IDENTITY_FIELDS:
        left = first.get(name)
        right = second.get(name)
        if SURFACE_UNKNOWN in (left, right):
            return None
        if left != right:
            return False
    return True


def surface_from_accessors(reader, output):
    """A surface identity assembled from what one reader can answer.

    For a program whose result records no identity of its own, this is
    what the host can honestly say: the level fields its reader already
    parses, and the word ``unknown`` for every field it cannot -- never
    a default, because a default would make two different surfaces
    compare equal.
    """

    def _read(selector):
        accessor = reader.accessors.get(selector)
        if accessor is None:
            return SURFACE_UNKNOWN
        try:
            value = accessor(output)
        except MissingQuantityError:
            return None
        except Exception:  # noqa: BLE001 - an unreadable field is unknown
            return SURFACE_UNKNOWN
        if isinstance(value, str):
            value = value.strip().lower() or None
        return value

    ab_initio = _read("ab_initio")
    functional = _read("functional")
    if isinstance(ab_initio, str) and any(
        marker in ab_initio for marker in _CORRELATED_METHOD_MARKERS
    ):
        method_family = ab_initio
    elif ab_initio == "hf":
        method_family = "hf"
    elif isinstance(functional, str) and functional:
        method_family = "dft"
    else:
        method_family = SURFACE_UNKNOWN
    charge = _read("charge")
    multiplicity = _read("multiplicity")
    return {
        "basis": _read("basis"),
        "charge": None if charge is None else charge,
        "constraints": [],
        # Not parsed from a log by this release; an unknown field makes
        # the surface incomparable rather than falsely equal.
        "density_fitting": SURFACE_UNKNOWN,
        "dispersion": SURFACE_UNKNOWN,
        "excited_root": None,
        "frozen_core": SURFACE_UNKNOWN,
        "functional_applied": functional,
        "method_family": method_family,
        "multiplicity": None if multiplicity is None else multiplicity,
        "reference": SURFACE_UNKNOWN,
        "solvent": _read("solvent"),
        "solvent_model": _read("solvation_model"),
        "state_manifold": None,
    }


def _resolve_computed_surface(reader, output, selector, word):
    """Resolve ``computed_surface`` to the surface *this* artifact computed on.

    Read from the result's own words: the followed root of an excited-root
    optimisation, the correlated method a stage recorded, or the
    ``ab_initio`` selector a log reader serves.  Any other declared word
    passes through unchanged.
    """

    if word != "computed_surface":
        return word
    # The surface the result recorded, where it recorded one: a method
    # family this resolver has never heard of used to fall through every
    # branch below and answer "reference", which is a true word about the
    # wrong density.
    family = (
        str(
            (getattr(output, "surface", None) or {}).get("method_family") or ""
        )
        .strip()
        .lower()
    )
    if family:
        if getattr(output, "excited_state_followed_root", None) or family in {
            "tda",
            "tddft",
            "eom_ccsd",
        }:
            return "excited_root"
        if any(marker in family for marker in _CORRELATED_METHOD_MARKERS):
            return "correlated"
        if family in {"hf", "dft"}:
            return "reference"
        return "unknown"
    if getattr(output, "excited_state_followed_root", None):
        return "excited_root"
    if getattr(output, "correlated_method", None):
        return "correlated"
    accessor = reader.accessors.get("ab_initio")
    if accessor is not None:
        try:
            value = str(accessor(output) or "").strip().lower()
        except Exception:  # noqa: BLE001 - absence resolves to the reference
            value = ""
        if any(marker in value for marker in _CORRELATED_METHOD_MARKERS):
            return "correlated"
    return "reference"


#: The two questions a reference-stability analysis asks, whatever
#: program asks them.  A program names its own orbital-rotation space
#: inside each answer -- PySCF says ``RHF/RKS -> UHF/UKS`` for a
#: restricted reference and ``UHF/UKS -> GHF/GKS`` for an unrestricted
#: one, Gaussian names none -- so the space rides the answer and is
#: never assumed from the question.  ``considered_perturbations`` is
#: Gaussian's own hedge, kept as its own word rather than folded into
#: either of the other two.  ``real_to_complex`` is the question PySCF
#: solves inside its external analysis and does not return: a record
#: that heard PySCF's answer to it carries it as a third question.
REFERENCE_STABILITY_QUESTIONS = (
    "internal",
    "external",
    "real_to_complex",
    "considered_perturbations",
)


def _stability_answer(
    question, *, rotation_space=None, reason=None, lowest_eigenvalue=None
):
    """One question's answer, in the shape every reader writes.

    ``lowest_eigenvalue`` is the number the verdict was drawn from, in the
    program's own normalisation and unit, where the program printed one:
    an anomaly is recorded with the numbers that tripped it.
    """

    answer = {"question": str(question)}
    if rotation_space:
        answer["rotation_space"] = str(rotation_space)
    if reason:
        answer["reason"] = str(reason)
    if lowest_eigenvalue is not None:
        answer["lowest_eigenvalue"] = float(lowest_eigenvalue)
    return answer


def _pyscf_stationarity_gradient(output: Any) -> float | None:
    """max|g| at the geometry a PySCF Hessian was taken on, in Eh/Bohr.

    ``results/forces`` is the negative gradient the ``hess`` stage recorded
    at the one structure it differentiated, and ``forces_unit`` declares
    the unit; both are checked here so a future contract that stores
    another unit is an absence rather than a number read as Eh/Bohr.  A
    stage that recorded no gradient -- a single point, an optimisation --
    answers nothing, which is what "the host cannot say" means.
    """

    import numpy as np

    if str(getattr(output, "forces_unit", "") or "") != "Eh/Bohr":
        return None
    forces = getattr(output, "forces", None)
    if forces is None:
        return None
    values = np.asarray(forces, dtype=float)
    if not values.size or not bool(np.isfinite(values).all()):
        return None
    return float(np.max(np.abs(values)))


@dataclass(frozen=True)
class CartesianHessianV1:
    """The Cartesian Hessian at the structure a result's modes belong to.

    What a projected harmonic analysis needs and a program's printed
    spectrum does not keep: the second-derivative matrix itself
    (``hessian``, 3N x 3N in Eh/Bohr^2), the structure it was taken at in
    the same frame (``positions_bohr``), and the atomic masses the
    program's own printed spectrum was computed with, so the host can show
    that it read the matrix the printed frequencies came from before it
    removes anything.  ``gradient`` (N x 3, Eh/Bohr) is present only where
    the result records the gradient at that same structure.  ``source``
    names where the matrix was read and ``source_sha256`` the digest of a
    sidecar file when it was not the result artifact itself.
    """

    positions_bohr: Any
    hessian: Any
    masses_amu: Any
    symbols: tuple[str, ...]
    source: str
    mass_convention: str
    source_sha256: str = ""
    gradient: Any = None


def _file_sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _orca_cartesian_hessian(output: Any) -> CartesianHessianV1 | None:
    """ORCA's ``<basename>.hess``, the sidecar every ORCA Freq run writes.

    ORCA prints its normal modes to six decimals and never the matrix; the
    ``.hess`` beside the output holds the matrix to eleven significant
    figures with the structure (Bohr) and the masses its analysis used
    (isotope-averaged: 15.999 for O, 1.008 for H).  A sidecar left by
    another run is not this result's: it is served only when its atoms
    are the result's in order, and the host then checks that it
    reproduces the printed spectrum.  No gradient: the ``.engrad`` beside
    an ``Opt Freq`` pairs the final coordinates and energy with the
    gradient of the cycle before (R10 Q21 g2-hooh, measured 1.6e-4 A
    apart), which belongs to no structure a Hessian was taken at.
    """

    import numpy as np

    filename = getattr(output, "filename", None)
    if not filename:
        return None
    path = Path(str(filename)).with_suffix(".hess")
    if not path.is_file():
        return None
    lines = path.read_text(errors="replace").splitlines()
    starts = {
        line.strip()[1:]: index
        for index, line in enumerate(lines)
        if line.strip().startswith("$")
    }
    if "hessian" not in starts or "atoms" not in starts:
        return None
    row = starts["hessian"] + 1
    size = int(lines[row].split()[0])
    hessian = np.full((size, size), np.nan)
    row += 1
    while row < len(lines) and not lines[row].strip().startswith("$"):
        header = lines[row].split()
        if not header:
            row += 1
            continue
        columns = [int(value) for value in header]
        for offset in range(size):
            parts = lines[row + 1 + offset].split()
            index = int(parts[0])
            for column, value in zip(columns, parts[1:]):
                hessian[index, column] = float(value)
        row += 1 + size
    row = starts["atoms"] + 1
    count = int(lines[row].split()[0])
    symbols, masses, positions = [], [], []
    for offset in range(count):
        parts = lines[row + 1 + offset].split()
        symbols.append(parts[0])
        masses.append(float(parts[1]))
        positions.append([float(value) for value in parts[2:5]])
    if 3 * count != size or not bool(np.isfinite(hessian).all()):
        return None
    try:
        expected = [str(item) for item in _orca_symbols(output)]
    except Exception:  # noqa: BLE001 - a reader that cannot say says nothing
        return None
    if expected != symbols:
        return None
    return CartesianHessianV1(
        positions_bohr=np.asarray(positions, dtype=float),
        hessian=hessian,
        masses_amu=np.asarray(masses, dtype=float),
        symbols=tuple(symbols),
        source=f"the ORCA sidecar {path.name}",
        mass_convention="the masses ORCA's .hess records",
        source_sha256=_file_sha256(path),
    )


def _gaussian_archive_sections(output: Any) -> list[str] | None:
    """The last archive entry of a Gaussian log, split at its ``\\\\``.

    Gaussian wraps the entry at a fixed width wherever the characters fall,
    the closing ``\\\\@`` included, so the lines are joined before the end
    is looked for: searched in the raw text, an entry wrapped between its
    last backslash and the ``@`` read as no entry, and the Hessian of 3 of
    122 archive-bearing Gaussian logs in this repository and R10 Q30's
    oracle (a frozen-atom optimisation; ethane held at 5 deg) was absent.
    """

    text = getattr(output, "content_lines_string", None)
    if not text:
        return None
    start = text.rfind(" 1\\1\\")
    if start < 0:
        return None
    joined = "".join(
        line[1:] if line.startswith(" ") else line
        for line in text[start:].splitlines()
    )
    end = joined.find("\\\\@")
    if end < 0:
        return None
    return joined[: end + 3].split("\\\\")


def _gaussian_cartesian_hessian(output: Any) -> CartesianHessianV1 | None:
    """The force constants a Gaussian frequency job writes into its archive.

    A ``Freq`` archive entry ends with the lower triangle of the Cartesian
    Hessian (Eh/Bohr^2) and then the gradient (Eh/Bohr, the negative of the
    printed forces), both in the frame of the geometry the same entry
    carries -- the input orientation, not the standard one the log's
    tables use.  The masses are the ones the log says its analysis used
    ("Atom N has atomic number Z and mass M").
    """

    import numpy as np

    sections = _gaussian_archive_sections(output)
    if not sections:
        return None
    try:
        marker = next(
            index
            for index, section in enumerate(sections)
            if "NImag=" in section
        )
        geometry = sections[3].split("\\")[1:]
        symbols = [
            re.sub(r"\(.*\)", "", row.split(",")[0]) for row in geometry
        ]
        positions = np.array(
            [
                [float(value) for value in row.split(",")[-3:]]
                for row in geometry
            ]
        )
        lower = [float(value) for value in sections[marker + 1].split(",")]
        gradient = np.array(
            [float(value) for value in sections[marker + 2].split(",")]
        )
    except (StopIteration, IndexError, ValueError):
        return None
    size = 3 * len(symbols)
    if len(lower) != size * (size + 1) // 2 or gradient.size != size:
        return None
    hessian = np.zeros((size, size))
    hessian[np.tril_indices(size)] = lower
    hessian = hessian + np.tril(hessian, -1).T
    masses = {}
    for line in getattr(output, "contents", ()) or ():
        found = re.match(
            r"\s*Atom\s+(\d+)\s+has atomic number\s+\d+\s+and mass\s+"
            r"([0-9.]+)",
            line,
        )
        if found:
            masses[int(found.group(1))] = float(found.group(2))
    if sorted(masses) != list(range(1, len(symbols) + 1)):
        return None
    return CartesianHessianV1(
        positions_bohr=positions / 0.529177210903,
        hessian=hessian,
        masses_amu=np.array(
            [masses[index + 1] for index in range(len(symbols))]
        ),
        symbols=tuple(symbols),
        source="the log's own archive entry",
        mass_convention="the masses the Gaussian log states it used",
        gradient=gradient.reshape(-1, 3),
    )


def _pyscf_cartesian_hessian(output: Any) -> CartesianHessianV1 | None:
    """``results/hessian`` of a PySCF ``hess`` stage, with its own gradient.

    The Hessian is stored per atom pair as (N, N, 3, 3) in Eh/Bohr^2 at
    ``results/positions``, and ``results/forces`` is the negative gradient
    at that same structure.  PySCF's analysis used isotope-averaged
    masses (``atom_mass_list(isotope_avg=True)``, the artifact's
    ``mass_convention``); the standard atomic weights ase tabulates are
    those numbers for the light elements, and the host's reproduction of
    the printed spectrum is what shows it for the atoms at hand.
    """

    import numpy as np
    from ase.data import atomic_masses, atomic_numbers

    results = getattr(output, "results", None) or {}
    stored = results.get("hessian") if hasattr(results, "get") else None
    if stored is None:
        return None
    if str(
        (getattr(output, "result_units", None) or {}).get(
            "results/hessian", "Eh/Bohr^2"
        )
    ) not in {"Eh/Bohr^2"}:
        return None
    blocks = np.asarray(stored, dtype=float)
    if blocks.ndim != 4:
        return None
    count = blocks.shape[0]
    symbols = [str(item) for item in (output.chemical_symbols or ())]
    if len(symbols) != count:
        return None
    positions = np.asarray(output.positions, dtype=float) / 0.529177210903
    forces = getattr(output, "forces", None)
    gradient = None
    if forces is not None and str(getattr(output, "forces_unit", "")) == (
        "Eh/Bohr"
    ):
        gradient = -np.asarray(forces, dtype=float).reshape(count, 3)
    return CartesianHessianV1(
        positions_bohr=positions,
        hessian=blocks.transpose(0, 2, 1, 3).reshape(3 * count, 3 * count),
        masses_amu=np.array(
            [
                float(atomic_masses[atomic_numbers[symbol]])
                for symbol in symbols
            ]
        ),
        symbols=tuple(symbols),
        source="the result's results/hessian",
        mass_convention="isotope-averaged standard atomic weights",
        gradient=gradient,
    )


@dataclass(frozen=True)
class ConvergenceCheckV1:
    """A program's own convergence check at the structure a result ends on.

    What the program printed when it judged whether the structure its
    modes (or its reached geometry) belong to is stationary: each
    criterion as ``(name, value, threshold, met)`` in the program's own
    coordinates and units, which rows measure the gradient
    (``force_rows``), the program's verdict on the whole check
    (``converged``), and the energy change its own quadratic model
    predicts where it prints one.  The criteria are the program's:
    Gaussian and ORCA judge forces on redundant internal coordinates, xTB
    the norm of the Cartesian gradient, geomeTRIC each Cartesian
    component.  One structure passes all of them only by accident -- a
    Gaussian-converged bromochloromethane has internal forces of at most
    2.9e-4 and a largest Cartesian component of 4.9e-4 -- so a host
    criterion applied to a structure another program judged is not that
    program's criterion (R10 Q33 census: 3 Gaussian and 9 ORCA converged
    minima above geomeTRIC's 4.5e-4).
    """

    program: str
    criterion: str
    source: str
    rows: tuple[tuple[str, float, float, bool], ...]
    force_rows: tuple[str, ...]
    #: The program's own answer, by its own rule: ORCA's and xTB's printed
    #: verdict on the search (ORCA converges a held H2O2 with its RMS
    #: gradient row above tolerance, by rules its table does not print),
    #: Gaussian's force rows at the structure (its frequency step's
    #: displacement rows are the Newton step along soft modes with the
    #: exact Hessian: they fail at 1e-8 Eh of predicted change on five
    #: archived minima its optimiser converged).
    stationary: bool
    converged: bool | None = None
    predicted_energy_change_eh: float | None = None
    #: The check was taken at the structure the result's modes were
    #: computed at (Gaussian's frequency step), not only at the end of a
    #: search: it then speaks for a result that searched for nothing.
    at_modes: bool = False

    @property
    def forces_met(self) -> bool:
        """Whether every gradient row of the check is within its threshold."""

        rows = [row for row in self.rows if row[0] in self.force_rows]
        return bool(rows) and all(row[3] for row in rows)

    def verdict(self) -> str:
        """The program's printed verdict, or '' when it printed none."""

        if self.converged is None:
            return ""
        return (
            "the program's verdict: converged"
            if self.converged
            else "the program's verdict: not converged"
        )

    def words(self, force_rows: bool = True) -> str:
        """``maximum force 1.62e-04 (threshold 4.5e-04)``, joined."""

        return ", ".join(
            f"{name.lower()} {value:.3g} (threshold {threshold:.3g})"
            for name, value, threshold, _met in self.rows
            if (name in self.force_rows) == force_rows
        )


def _gaussian_job_blocks(output: Any) -> list[list[str]]:
    """A Gaussian log's lines, one list per job step (Link1)."""

    blocks, current = [], []
    for line in getattr(output, "contents", None) or ():
        current.append(line)
        if line.startswith("Normal termination of Gaussian") or (
            line.startswith("Error termination")
        ):
            blocks.append(current)
            current = []
    if current:
        blocks.append(current)
    return blocks


_GAUSSIAN_CHECK_ROW = re.compile(
    r"^(Maximum Force|RMS\s+Force|Maximum Displacement|RMS\s+Displacement)"
    r"\s+([-+0-9.DE]+)\s+([-+0-9.DE]+)\s+(YES|NO)\b"
)


def _gaussian_float(text: str) -> float:
    return float(str(text).replace("D", "E"))


def _gaussian_convergence_check(output: Any) -> ConvergenceCheckV1 | None:
    """Gaussian's own convergence table at the structure its modes belong to.

    Gaussian judges a structure with its Berny criteria on redundant
    internal coordinates -- maximum and RMS force (Hartree/Bohr or
    Hartree/radian), maximum and RMS displacement of the Newton step --
    and prints the table wherever it computes forces: every optimisation
    step, and once more in a frequency job, with the exact Hessian, at the
    very structure the modes are computed at.  That last table is read:
    the frequency step's when the log has modes (an ``opt freq`` or a
    ``freq`` job), else the optimisation's last.  A frozen coordinate is
    left out of Gaussian's force check, so on a constrained optimisation
    the check is on the surface it held.
    """

    blocks = _gaussian_job_blocks(output)
    with_modes = [
        block
        for block in blocks
        if any(line.startswith("Harmonic frequencies") for line in block)
    ]
    with_check = [
        block
        for block in blocks
        if any("Converged?" in line and "Threshold" in line for line in block)
        and any(line.startswith("Berny optimization") for line in block)
    ]
    chosen = with_modes[-1:] or with_check[-1:]
    if not chosen:
        return None
    block = chosen[0]
    starts = [
        index
        for index, line in enumerate(block)
        if "Converged?" in line and "Threshold" in line
    ]
    if not starts:
        return None
    start = starts[-1]
    rows = []
    for line in block[start + 1 : start + 5]:
        found = _GAUSSIAN_CHECK_ROW.match(line)
        if found:
            rows.append(
                (
                    re.sub(r"\s+", " ", found.group(1)),
                    _gaussian_float(found.group(2)),
                    _gaussian_float(found.group(3)),
                    found.group(4) == "YES",
                )
            )
    if len(rows) != 4:
        return None
    predicted = None
    converged = False
    for line in block[start + 5 : start + 10]:
        if line.startswith("Predicted change in Energy="):
            predicted = _gaussian_float(line.split("=", 1)[1].split()[0])
        if line.startswith("Optimization completed"):
            converged = True
    return ConvergenceCheckV1(
        program="gaussian",
        criterion=(
            "Gaussian's Berny criteria on redundant internal coordinates "
            "(forces in Hartree/Bohr or Hartree/radian)"
        ),
        source=(
            "its frequency step's check with the exact Hessian"
            if with_modes
            else "its optimisation's last check"
        ),
        rows=tuple(rows),
        force_rows=("Maximum Force", "RMS Force"),
        stationary=all(row[3] for row in rows[:2]),
        converged=converged,
        predicted_energy_change_eh=predicted,
        at_modes=bool(with_modes),
    )


def _orca_convergence_check(output: Any) -> ConvergenceCheckV1 | None:
    """ORCA's last geometry-convergence table and its own verdict.

    ORCA checks the energy change, the RMS and maximum gradient and the RMS
    and maximum step of its optimisation coordinates (redundant internal
    coordinates unless told otherwise) against the tolerances the table
    prints, and says whether the search converged in its own words.  The
    gradient rows are at the structure whose gradient it last computed; the
    step rows bound how far the final structure is from it.
    """

    lines = list(getattr(output, "contents", None) or ())
    starts = [
        index
        for index, line in enumerate(lines)
        if "|Geometry convergence|" in line
    ]
    if not starts:
        return None
    rows = []
    for line in lines[starts[-1] + 3 : starts[-1] + 12]:
        parts = line.split()
        if len(parts) >= 5 and parts[-1] in {"YES", "NO"}:
            try:
                rows.append(
                    (
                        " ".join(parts[:-3]),
                        float(parts[-3]),
                        float(parts[-2]),
                        parts[-1] == "YES",
                    )
                )
            except ValueError:
                continue
        elif rows and line.startswith("---"):
            break
    if not rows:
        return None
    converged = None
    for line in lines:
        if "THE OPTIMIZATION HAS CONVERGED" in line:
            converged = True
        elif "The optimization did not converge" in line:
            converged = False
    gradient_rows = [
        row for row in rows if row[0] in {"RMS gradient", "MAX gradient"}
    ]
    return ConvergenceCheckV1(
        program="orca",
        criterion=(
            "ORCA's geometry-convergence tolerances on its optimisation "
            "coordinates (gradient in Eh/Bohr)"
        ),
        source="its last geometry-convergence check",
        rows=tuple(rows),
        force_rows=("RMS gradient", "MAX gradient"),
        stationary=(
            bool(converged)
            if converged is not None
            else bool(gradient_rows) and all(row[3] for row in gradient_rows)
        ),
        converged=converged,
    )


def _xtb_convergence_check(output: Any) -> ConvergenceCheckV1 | None:
    """xTB's own criterion for the optimisation level it ran.

    xTB's optimiser (ANCopt) converges when the energy change and the norm
    of the Cartesian gradient fall below the thresholds of the level asked
    for (normal: 5e-6 Eh and 1e-3 Eh/Bohr; loose: 5e-5 and 4e-3); the
    setup block prints both, and the summary prints the gradient norm at
    the structure the run reached.  A run that optimised nothing (a
    ``--hess`` at a handed geometry) prints no level and has no check.
    """

    try:
        threshold = output.gradient_convergence
        level = output.optimization_level
        norm = output.gradient_norm
        energy = output.energy_convergence
        converged = output.geometry_optimization_converged
    except AttributeError:
        return None
    if threshold is None or norm is None:
        return None
    return ConvergenceCheckV1(
        program="xtb",
        criterion=(
            f"xTB's optimisation level {level}: gradient norm at most "
            f"{float(threshold):g} Eh/Bohr"
            + (
                f" and energy change at most {float(energy):g} Eh"
                if energy is not None
                else ""
            )
        ),
        source="the gradient norm at the structure it reached",
        rows=(
            (
                "gradient norm",
                float(norm),
                float(threshold),
                float(norm) <= float(threshold),
            ),
        ),
        force_rows=("gradient norm",),
        stationary=bool(converged),
        converged=bool(converged),
    )


_GAUSSIAN_PARAMETER_ROW = re.compile(
    r"^!\s*(\S+)\s+([RADL])\(([-0-9,]+)\)\s+(-?[0-9.]+)\s+(.*)$"
)


def _gaussian_frozen_parameters(
    output: Any,
) -> tuple[list[dict[str, Any]], tuple[int, ...]] | None:
    """What Gaussian froze, from its own parameter table; None if none.

    Gaussian prints every coordinate its optimiser works in, and whether
    it is frozen, in the "Initial Parameters" table of each optimising or
    frequency step -- whatever put the freeze there: a ModRedundant
    section, a Cartesian ``-1`` flag in the geometry, or a checkpoint read
    with ``geom=check`` that still carries a previous job's constraints.
    The last table is the one that holds for the structure the log ends
    on.  Returns the frozen internal coordinates as the ModRedundant rows
    would name them (``kind``, one-based ``atoms``, the ``value`` in
    Angstrom or degrees, ``label``), and the zero-based atoms whose three
    Cartesian coordinates it froze.
    """

    lines = list(getattr(output, "contents", None) or ())
    starts = [
        index
        for index, line in enumerate(lines)
        if "Initial Parameters" in line
    ]
    if not starts:
        return None
    kinds = {2: "bond", 3: "angle", 4: "dihedral"}
    internal, cartesian = [], {}
    for line in lines[starts[-1] + 1 :]:
        if line.startswith("GradGradGrad") or line.startswith("Trust Radius"):
            break
        found = _GAUSSIAN_PARAMETER_ROW.match(line)
        if not found or "frozen" not in found.group(5):
            continue
        atoms = [int(value) for value in found.group(3).split(",")]
        if len(atoms) == 2 and atoms[1] < 0:
            cartesian.setdefault(atoms[0] - 1, set()).add(atoms[1])
            continue
        kind = kinds.get(len(atoms)) if found.group(2) != "L" else None
        if kind is None or min(atoms) < 1:
            continue
        internal.append(
            {
                "kind": kind,
                "atoms": tuple(atoms),
                "value": float(found.group(4)),
                "label": f"{found.group(1)} {found.group(2)}"
                f"({found.group(3)}) frozen",
            }
        )
    frozen_atoms = tuple(
        sorted(atom for atom, axes in cartesian.items() if len(axes) == 3)
    )
    return internal, frozen_atoms


def _gaussian_frozen_atoms(output: Any) -> tuple[int, ...] | None:
    parameters = _gaussian_frozen_parameters(output)
    return None if parameters is None else parameters[1]


def _orca_frozen_atoms(output: Any) -> tuple[int, ...] | None:
    """The zero-based atoms ORCA held fixed in space (``{C n C}``).

    ORCA reports a Cartesian constraint as "Will constrain atom n
    coordinate k" for each of an atom's three coordinates.  It also holds
    only the out-of-plane coordinate of every atom when a dihedral is
    constrained at an exactly planar value, which holds that dihedral and
    no atom, so only atoms held in all three coordinates are counted.
    """

    held: dict[int, set[int]] = {}
    for line in getattr(output, "contents", None) or ():
        found = re.match(r"Will constrain atom (\d+) coordinate (\d+)", line)
        if found:
            held.setdefault(int(found.group(1)), set()).add(
                int(found.group(2))
            )
    return tuple(sorted(atom for atom, axes in held.items() if len(axes) == 3))


@dataclass(frozen=True)
class TorsionalScanV1:
    """A relaxed scan of one dihedral, as the program that drove it ran it.

    ``atoms`` are the driven dihedral's one-based atoms in the program's
    own order; ``values_deg`` the dihedral at each converged point and
    ``energies_eh`` its energy (the program's total, dispersion included
    where the method has it).  ``measured`` says whether each value was
    measured in the structure the program wrote for that point
    (``positions_angstrom``) or is the target the program was asked to
    hold, used only where no structure was kept.  ``planned_points`` is
    how many points the scan declared, ``source`` where the rows were read.
    """

    atoms: tuple[int, ...]
    values_deg: tuple[float, ...]
    energies_eh: tuple[float, ...]
    symbols: tuple[str, ...]
    source: str
    measured: bool
    planned_points: int | None = None
    positions_angstrom: tuple[Any, ...] | None = None


def _dihedral_degrees(positions: Any, atoms: Sequence[int]) -> float:
    """One-based ``atoms``' dihedral in degrees, (-180, 180]."""

    from chemsmart.analysis.thermochemistry import internal_coordinate_value

    return math.degrees(
        internal_coordinate_value(positions, [int(i) - 1 for i in atoms])
    )


def _read_xyz_rows(path: Path) -> tuple[list[str], Any] | None:
    import numpy as np

    try:
        lines = Path(path).read_text(errors="replace").splitlines()
        count = int(lines[0].split()[0])
        rows = [line.split() for line in lines[2 : 2 + count]]
        return (
            [row[0] for row in rows],
            np.array([[float(value) for value in row[1:4]] for row in rows]),
        )
    except (OSError, ValueError, IndexError):
        return None


def _orca_torsional_scan(output: Any) -> TorsionalScanV1 | None:
    """ORCA's relaxed surface scan of one dihedral.

    The surface is ORCA's own 'Actual Energy' table (``scan_profile``) and
    each point's structure the ``<stem>.NNN.xyz`` ORCA writes for it
    (``scan_point_records``); the dihedral is measured in that structure.
    """

    coordinate = getattr(output, "scan_coordinate", None)
    if not coordinate or coordinate.get("kind") != "dihedral":
        return None
    records = tuple(getattr(output, "scan_point_records", ()) or ())
    if not records:
        return None
    symbols = [str(item) for item in _orca_symbols(output)]
    atoms = tuple(int(item) for item in coordinate["atoms"])
    values, energies, structures = [], [], []
    measured = True
    for record in records:
        energies.append(float(record["energy"]))
        read = (
            _read_xyz_rows(Path(record["geometry_file"]))
            if record.get("geometry_file")
            else None
        )
        if read is None or read[0] != symbols:
            measured = False
            structures.append(None)
            values.append(float(record["coordinate"]))
            continue
        structures.append(read[1])
        values.append(_dihedral_degrees(read[1], atoms))
    return TorsionalScanV1(
        atoms=atoms,
        values_deg=tuple(values),
        energies_eh=tuple(energies),
        symbols=tuple(symbols),
        source="ORCA's relaxed-surface-scan table and its per-point .xyz files",
        measured=measured,
        planned_points=int(coordinate.get("points") or 0) or None,
        positions_angstrom=(
            tuple(structures)
            if all(s is not None for s in structures)
            else None
        ),
    )


def _gaussian_torsional_scan(output: Any) -> TorsionalScanV1 | None:
    """Gaussian's relaxed scan of one dihedral (``opt=modredundant`` ``S``).

    Gaussian prints no surface table; the points are the structures the log
    marks as converged scan points, each with its own energy, and the
    dihedral is measured in each (``scan_profile``).
    """

    import numpy as np

    coordinate = getattr(output, "scan_coordinate", None)
    if not coordinate or coordinate.get("kind") != "dihedral":
        return None
    profile = getattr(output, "scan_profile", None)
    structures = getattr(output, "all_structures", None)
    if not profile or not structures or len(profile) != len(structures):
        return None
    symbols = [str(item) for item in structures[0].chemical_symbols]
    atoms = tuple(int(item) for item in coordinate["atoms"])
    positions = [
        np.asarray(item.positions, dtype=float) for item in structures
    ]
    return TorsionalScanV1(
        atoms=atoms,
        values_deg=tuple(_dihedral_degrees(x, atoms) for x in positions),
        energies_eh=tuple(float(row["energy"]) for row in profile),
        symbols=tuple(symbols),
        source="the structures Gaussian's log marks as converged scan points",
        measured=True,
        planned_points=int(coordinate.get("points") or 0) or None,
        positions_angstrom=tuple(positions),
    )


def _pyscf_reference_diagnostics(output: Any) -> Mapping[str, Any] | None:
    """What PySCF's own analysis said about the reference this result
    stands on, or ``None`` when nothing was recorded.

    ``None`` covers three different absences on purpose -- an artifact
    written under a contract older than v7, a run nobody asked, and an
    analysis that raised -- because none of them is a stable reference
    and no consumer should be able to tell them apart by accident.  The
    third is already a ``property_failures`` entry the validator reports
    on its own.
    """

    record = _pyscf_stability_record(output)
    if not isinstance(record, Mapping):
        return None
    analyses = record.get("analyses")
    if not isinstance(analyses, Mapping):
        return None
    unstable, stable, unavailable = [], [], []
    for question in REFERENCE_STABILITY_QUESTIONS:
        entry = analyses.get(question)
        if not isinstance(entry, Mapping):
            continue
        if entry.get("unavailable"):
            unavailable.append(
                _stability_answer(
                    question, reason=str(entry.get("unavailable"))
                )
            )
            continue
        answered = entry.get("stable")
        if answered is None:
            continue
        target = stable if bool(answered) else unstable
        eigenvalues = entry.get("lowest_eigenvalues") or ()
        target.append(
            _stability_answer(
                question,
                rotation_space=entry.get("rotation_space"),
                lowest_eigenvalue=(
                    eigenvalues[0]
                    if eigenvalues
                    and record.get("eigenvalue_unit")
                    == _pyscf_stability_eigenvalue_unit()
                    else None
                ),
            )
        )
    not_determined = [
        _stability_answer(
            str(entry.get("question") or name),
            rotation_space=entry.get("rotation_space"),
            reason=(
                str(entry.get("reason") or "")
                + _printed_pointer(output, str(name))
            ),
        )
        for name, entry in sorted((record.get("not_determined") or {}).items())
        if isinstance(entry, Mapping)
    ]
    if not (unstable or stable or unavailable or not_determined):
        return None
    return {
        "analysis": "scf_stability",
        "applies_to": str(record.get("applies_to") or "reference"),
        "source": str(record.get("source") or ""),
        "reference_class": str(record.get("reference_class") or ""),
        "reference_family": str(record.get("reference_family") or ""),
        # The orbitals the analysis ran on.  PySCF answers for an SCF
        # that never converged, about orbitals that are not stationary,
        # and a reader of the answer is owed that fact beside it.
        "reference_converged": record.get("scf_converged"),
        "unstable": tuple(unstable),
        "stable": tuple(stable),
        "unavailable": tuple(unavailable),
        "not_determined": tuple(not_determined),
    }


#: Gaussian's own verdict words, and what each says about the
#: wavefunction the run ended on.  A linked job can find an instability,
#: reoptimise the wavefunction and end stable, so the *last* verdict is
#: the state of the reference and the earlier ones are its history: the
#: archived ``dna_link_sp`` log finds an internal instability and then
#: reports stability, and reading anything but the last entry would call
#: that reference unstable.
_GAUSSIAN_STABILITY_VERDICTS = {
    "internal_instability": ("unstable", "internal"),
    "external_instability": ("unstable", "external"),
    "stable_under_considered_perturbations": (
        "stable",
        "considered_perturbations",
    ),
}


def _gaussian_reference_diagnostics(output: Any) -> Mapping[str, Any] | None:
    """Gaussian's stability verdict for the wavefunction it ended on.

    The rotation space is the one Gaussian's own sentence names, and only
    where it names one: an unstable answer says "internal" or the larger
    space the reference fell into ("RHF -> UHF"), and a stable answer --
    "stable under the perturbations considered" -- names none, so none is
    recorded and none is invented.  The lowest eigenvalue of the stability
    matrix printed with the verdict rides the answer, in hartree.
    """

    history = tuple(
        str(item)
        for item in (
            getattr(output, "wavefunction_stability_history", None) or ()
        )
    )
    if not history:
        return None
    standing, question = _GAUSSIAN_STABILITY_VERDICTS.get(
        history[-1], (None, None)
    )
    if standing is None:
        return None
    records = list(
        getattr(output, "wavefunction_stability_records", None) or ()
    )
    last = records[-1] if records else {}
    eigenvalues = last.get("eigenvalues") or ()
    answer = (
        _stability_answer(
            question,
            rotation_space=last.get("rotation_space"),
            lowest_eigenvalue=(
                eigenvalues[0]["eigenvalue"] if eigenvalues else None
            ),
        ),
    )
    return {
        "analysis": "wavefunction_stability",
        "applies_to": "reference",
        "source": "gaussian stability analysis",
        "verdict": history[-1],
        "history": history,
        "unstable": answer if standing == "unstable" else (),
        "stable": answer if standing == "stable" else (),
        "unavailable": (),
        "not_determined": (),
    }


def _electronic_provenance_table(accessors, declared):
    """Keep the declared provenance rows this reader actually implements."""

    return tuple(
        sorted(
            (selector, word)
            for selector, word in declared
            if selector in accessors
        )
    )


@dataclass(frozen=True)
class ResultReaderV1:
    """How one program's result file answers the shared selector vocabulary."""

    program: str
    #: Host artifact kind this reader consumes; extraction refuses any other.
    artifact_kind: str
    parser_id: str
    #: Constructs the parser from an already host-verified path.
    open_output: Callable[[Path], Any]
    #: Selector name to a callable reading it from the parser object.
    accessors: dict[str, Callable[[Any], Any]]
    #: Optional native file that supplied one structural selector.  Most
    #: readers obtain coordinates from the bound result artifact itself; xTB
    #: reaches its final optimisation frame through ``xtbopt.*`` beside the
    #: main log.  The shared handoff seals that sidecar's digest when present.
    geometry_source_path_for_selector: (
        Callable[[Any, str], Path | None] | None
    ) = None
    #: Optional native sidecars a selector reads.  Extraction, unlike the
    #: geometry handoff, owns the receipt that carries these bytes, so it
    #: verifies and digests them beside the primary result artifact.
    native_evidence_paths_for_selector: (
        Callable[[Any, str], tuple[Path, ...]] | None
    ) = None
    #: Program-native units that differ from the shared selector display unit.
    source_units: Mapping[str, str] = field(default_factory=dict)
    #: Selectors this parser can extract from a program job type when the
    #: chosen method/settings emit them. Missing coverage means unknown, never
    #: that the job produces no quantities.
    jobtype_selectors: tuple[tuple[str, tuple[str, ...]], ...] = ()
    #: Where the word a plan names a stage with is not the word a finished
    #: result of that stage answers to, as ``((stage, (jobtype, ...)), ...)``.
    #: A stage absent here is its own result word.  See
    #: ``result_jobtypes_for_stage``.
    stage_result_jobtypes: tuple[tuple[str, tuple[str, ...]], ...] = ()
    #: Which structural state each selector reads, as
    #: ``((selector, state), ...)``. A selector absent from this tuple is
    #: ``stateless``.
    #:
    #: One ORCA result carries several molecular states, and two
    #: selectors of the same result were serving different ones under
    #: one receipt. Measured on a live unconverged ``OptTS``
    #: (po3-r18 ``ts-ester-c4``, 2026-09-11): ``positions`` returns
    #: ``thermochemistry_molecule``, the geometry the Hessian was
    #: computed at, while ``energy`` returns the last of a hundred
    #: printed ``FINAL SINGLE POINT ENERGY`` lines. The two describe
    #: structures **1.231588 A** apart whose energies differ by
    #: **182.2 kcal/mol**, and nothing on the receipt said which
    #: structure either belonged to. ``build_reached_geometry`` then
    #: asked for ``positions`` and wrote "the structure orca reached"
    #: onto the *input* coordinates, so the recovery route the repair
    #: menu offers for a non-converged run returned the seed.
    #:
    #: Which state a parser reads is a fact about the parser, not a
    #: chemical judgement. Declaring it lets a consumer ask for the role
    #: it needs and be refused rather than silently served another.
    selector_structural_states: tuple[tuple[str, str], ...] = ()
    #: Selectors this reader *introduces* to the shared vocabulary, as
    #: ``((selector, unit, dimension), ...)``.  A selector several programs
    #: already serve stays in the shared tables below; one that only this
    #: program's parser can answer is declared here, beside the accessor
    #: that reads it, so adding it is an edit to this reader and to nothing
    #: else.  ``merge_selector_declarations`` folds these into the shared
    #: tables and refuses a unit or dimension another declaration disagrees
    #: with: one selector name has one physical meaning on every program.
    selector_declarations: tuple[tuple[str, str, str], ...] = ()
    #: Atom-resolved metadata for a selector this reader introduces, as
    #: ``((selector, ((key, value), ...)), ...)`` -- the same record
    #: ``atom_resolved_selector_metadata`` answers for the shared ones.
    atom_resolved_declarations: tuple[
        tuple[str, tuple[tuple[str, str], ...]], ...
    ] = ()
    #: Native program outputs must prove normal termination.  A standalone
    #: geometry artifact is data rather than an engine run, so that format can
    #: opt out while retaining the same typed quantity path.
    requires_normal_termination: bool = True
    #: A reader may hold a stricter admission for *quantities* than for
    #: opening.  PySCF opens any receipt-bound artifact -- a failed run is
    #: inspectable, its last geometry bindable -- while every extracted
    #: number and every free energy also demands the green receipt.  Log
    #: readers have no such split: their normal-termination gate is the
    #: whole of it.
    admit_for_analysis: Callable[[Path], Any] | None = None
    #: Which electronic state or method each selector's value belongs to,
    #: as ``((selector, word), ...)`` over ``ELECTRONIC_PROVENANCES``; a
    #: selector absent from the tuple is ``stateless``.  Declared beside
    #: the structural state because the two answer different questions: a
    #: geometry identity says nothing about whose density a dipole is.
    selector_electronic_provenance: tuple[tuple[str, str], ...] = ()
    #: Resolves a declared word against one opened result; the shared
    #: resolver turns ``computed_surface`` into the concrete word this
    #: artifact's own record supports.  None keeps every word as declared.
    resolve_electronic_provenance: Callable[..., str] | None = None
    #: The level of theory one opened result computed at, as its own
    #: record names it -- method, basis, solvent, the frozen-core count a
    #: correlated stage applied, the response an excited stage ran on and
    #: the root it followed -- for the inspection reply.  None renders no
    #: level.  A level states what the program applied in the vocabulary
    #: every program's level uses (a functional is the literal its applied
    #: form names), so an operation combining two results' numbers can be
    #: told when their levels differ; it is compared, never refused.
    resolve_level: Callable[[Any], Mapping[str, Any]] | None = None
    #: The electronic surface one opened result is on, as a mapping over
    #: ``SURFACE_IDENTITY_FIELDS``.  Unlike the level, a surface is *for*
    #: comparing: two results are one surface when every field agrees and
    #: neither reader wrote ``unknown``.  None means this reader cannot
    #: say, and the organs that ask treat that as "not comparable" rather
    #: than as agreement.
    resolve_surface: Callable[[Any], Mapping[str, Any] | None] | None = None
    #: What this result's own run recorded about the *reference* every
    #: number above it stands on -- today, whether a stability analysis
    #: found the converged orbitals to be a minimum in orbital-rotation
    #: space.  A diagnostic about the wavefunction rather than about a
    #: quantity, so it is not a selector: it has no unit, enters no
    #: arithmetic, and flattening two named questions into one verdict
    #: string is exactly the shape the recording round rejected.  It is
    #: the one path such a diagnostic takes into the host's sensors, so
    #: the next one is a reader function rather than another hand-thread.
    #: ``None`` means this reader cannot say, and no organ reads that as
    #: a stable reference.
    resolve_reference_diagnostics: (
        Callable[[Any], Mapping[str, Any] | None] | None
    ) = None
    #: The largest absolute Cartesian gradient component, in Eh/Bohr, at
    #: *the structure this result's spectrum belongs to* -- the one number
    #: that says whether the point a Hessian was taken at is stationary at
    #: all.  It is not a selector: it enters no arithmetic and answers a
    #: question about a geometry rather than carrying a quantity from it.
    #:
    #: A reader declares it only where the result records a gradient bound
    #: to that one structure.  ORCA's and Gaussian's ``forces`` are a list
    #: of the gradient at *every* optimisation step, so a maximum over them
    #: belongs to no single geometry and those readers answer nothing here
    #: rather than a number read from the wrong structure.  ``None`` is an
    #: absence, never zero, and no organ reads it as stationarity.
    resolve_stationarity_gradient: Callable[[Any], float | None] | None = None
    #: The Cartesian Hessian at the structure this result's spectrum
    #: belongs to (:class:`CartesianHessianV1`), for the analysis that
    #: removes a held coordinate from it.  Like the gradient above it is
    #: not a selector: a matrix enters no expression, it is what a free
    #: energy along a held coordinate is computed from.  ``None`` means
    #: this reader cannot serve one for this result.
    resolve_cartesian_hessian: (
        Callable[[Any], CartesianHessianV1 | None] | None
    ) = None
    #: The relaxed scan of one dihedral this result ran
    #: (:class:`TorsionalScanV1`), for the hindered-rotor treatment that
    #: takes a torsion's potential from it.  Not a selector: it is a table
    #: of structures and energies a partition function is built from.
    #: ``None`` means this reader cannot serve one for this result.
    resolve_torsional_scan: Callable[[Any], TorsionalScanV1 | None] | None = (
        None
    )
    #: The program's own convergence check at the structure this result's
    #: modes, or its reached geometry, belong to (:class:`ConvergenceCheckV1`):
    #: the criteria the program judged that structure by, in its own
    #: coordinates, with the numbers it printed.  Not a selector: it says
    #: whether a structure is stationary by the criterion it was converged
    #: with, which a host criterion on another coordinate system is not.
    #: ``None`` means the program printed no such check for this result.
    resolve_convergence_check: (
        Callable[[Any], ConvergenceCheckV1 | None] | None
    ) = None
    #: The zero-based atoms this result held fixed in space in all three
    #: Cartesian coordinates, from the program's own record.  Not a
    #: selector: it says what the structure was not relaxed along.  Gaussian
    #: prints no force on such an atom (its archive gradient reads zero
    #: there), so no gradient this host can read shows the structure
    #: stationary.  ``None`` means this reader cannot say; ``()`` that the
    #: program froze none.
    resolve_frozen_atoms: Callable[[Any], tuple[int, ...] | None] | None = None

    def convergence_check_for_output(
        self, output: Any
    ) -> ConvergenceCheckV1 | None:
        """The program's own check at this result's structure, or None."""

        if self.resolve_convergence_check is None:
            return None
        try:
            return self.resolve_convergence_check(output)
        except (
            Exception
        ):  # noqa: BLE001 - a reader that cannot say says nothing
            return None

    def frozen_atoms_for_output(self, output: Any) -> tuple[int, ...] | None:
        """The atoms this result froze in space, or None when unknowable."""

        if self.resolve_frozen_atoms is None:
            return None
        try:
            return self.resolve_frozen_atoms(output)
        except (
            Exception
        ):  # noqa: BLE001 - a reader that cannot say says nothing
            return None

    def torsional_scan_for_output(self, output: Any) -> TorsionalScanV1 | None:
        """The dihedral scan this result ran, or None when unservable."""

        if self.resolve_torsional_scan is None:
            return None
        try:
            return self.resolve_torsional_scan(output)
        except (
            Exception
        ):  # noqa: BLE001 - a reader that cannot say says nothing
            return None

    def cartesian_hessian_for_output(
        self, output: Any
    ) -> CartesianHessianV1 | None:
        """The Hessian at this result's structure, or None when unservable."""

        if self.resolve_cartesian_hessian is None:
            return None
        try:
            return self.resolve_cartesian_hessian(output)
        except (
            Exception
        ):  # noqa: BLE001 - a reader that cannot say says nothing
            return None

    def stationarity_gradient_for_output(self, output: Any) -> float | None:
        """max|g| (Eh/Bohr) where this result's spectrum lives, or None.

        One function, because two organs ask it: the run sensor that
        raises ``stationary_point.gradient_above_optimizer_criterion`` and
        the characterisation that refuses to call a non-stationary point a
        stationary point of any order.  They had each reached into the
        reader's attributes with their own unit guard.
        """

        if self.resolve_stationarity_gradient is None:
            return None
        try:
            value = self.resolve_stationarity_gradient(output)
        except (
            Exception
        ):  # noqa: BLE001 - a reader that cannot say says nothing
            return None
        if value is None:
            return None
        value = float(value)
        return value if value == value else None

    def surface_for_output(self, output: Any) -> Mapping[str, Any] | None:
        """The surface this result is on, or None when unknowable."""

        if self.resolve_surface is None:
            return None
        return self.resolve_surface(output)

    def reference_diagnostics_for_output(
        self, output: Any
    ) -> Mapping[str, Any] | None:
        """What this result recorded about its own reference, or None."""

        if self.resolve_reference_diagnostics is None:
            return None
        return self.resolve_reference_diagnostics(output)

    def spin_symmetry_for_output(self, output: Any) -> dict[str, Any] | None:
        """Which reference ran and whether its spin symmetry broke, or None.

        One function for every program (:func:`spin_symmetry_record`), read
        from the level and the <S**2> selectors this reader already serves,
        so the sensor that raises a collapsed broken-symmetry request and a
        session reading the same result see one answer.
        """

        return spin_symmetry_record(self, output)

    def __post_init__(self) -> None:
        jobtypes = tuple(item[0] for item in self.jobtype_selectors)
        if jobtypes != tuple(sorted(set(jobtypes))):
            raise ValueError(
                "jobtype selector declarations must be sorted and unique"
            )
        for jobtype, selectors in self.jobtype_selectors:
            if not jobtype or jobtype != jobtype.strip().lower():
                raise ValueError("jobtype selector keys must be normalized")
            if selectors != tuple(sorted(set(selectors))):
                raise ValueError(
                    f"{jobtype} selector declaration must be sorted and unique"
                )
            undeclared = set(selectors) - self.selectors
            if undeclared:
                raise ValueError(
                    f"{jobtype} selector declaration is not implemented: "
                    f"{sorted(undeclared)}"
                )
        named = tuple(item[0] for item in self.selector_structural_states)
        if named != tuple(sorted(set(named))):
            raise ValueError(
                "selector structural states must be sorted and unique"
            )
        for selector, state in self.selector_structural_states:
            if selector not in self.selectors:
                raise ValueError(
                    f"structural state declared for {selector!r}, which "
                    "this reader does not implement"
                )
            if state not in STRUCTURAL_STATES:
                raise ValueError(
                    f"{selector}: {state!r} is not one of "
                    f"{STRUCTURAL_STATES}"
                )
        provenance_names = tuple(
            item[0] for item in self.selector_electronic_provenance
        )
        if provenance_names != tuple(sorted(set(provenance_names))):
            raise ValueError(
                "selector electronic provenance must be sorted and unique"
            )
        for selector, word in self.selector_electronic_provenance:
            if selector not in self.selectors:
                raise ValueError(
                    f"electronic provenance declared for {selector!r}, "
                    "which this reader does not implement"
                )
            if word not in ELECTRONIC_PROVENANCES:
                raise ValueError(
                    f"{selector}: {word!r} is not one of "
                    f"{ELECTRONIC_PROVENANCES}"
                )

    @property
    def selectors(self) -> frozenset[str]:
        return frozenset(self.accessors)

    def structural_state(self, selector: str) -> str:
        """Which molecular state this selector's value belongs to."""

        for name, state in self.selector_structural_states:
            if name == selector:
                return state
        return "stateless"

    def electronic_provenance(self, selector: str) -> str:
        """Which electronic state or method this selector's value belongs to,
        as declared; ``computed_surface`` stays symbolic here."""

        for name, word in self.selector_electronic_provenance:
            if name == selector:
                return word
        return "stateless"

    def electronic_provenance_for_output(
        self, output: Any, selector: str
    ) -> str:
        """The declared provenance resolved against one opened result."""

        word = self.electronic_provenance(selector)
        if self.resolve_electronic_provenance is None:
            return word
        return str(
            self.resolve_electronic_provenance(self, output, selector, word)
        )

    def level_for_output(self, output: Any) -> dict[str, Any]:
        """The level this result computed at, from its own record."""

        if self.resolve_level is None:
            return {}
        return dict(self.resolve_level(output))

    def geometry_source_path_for_output(
        self, output: Any, selector: str
    ) -> Path | None:
        """Return the exact native geometry file behind a selector, if any."""

        if self.geometry_source_path_for_selector is None:
            return None
        path = self.geometry_source_path_for_selector(output, selector)
        return Path(path) if path is not None else None

    def native_evidence_paths_for_output(
        self, output: Any, selector: str
    ) -> tuple[Path, ...]:
        """Return every native sidecar that this selector actually reads."""

        if self.native_evidence_paths_for_selector is None:
            return ()
        return tuple(
            Path(path)
            for path in self.native_evidence_paths_for_selector(
                output, selector
            )
        )

    def selectors_in_state(self, state: str) -> tuple[str, ...]:
        """Every selector this reader serves for one structural state."""

        return tuple(
            sorted(
                name
                for name in self.accessors
                if self.structural_state(name) == state
            )
        )

    def selectors_in_state_for_output(
        self, output: Any, state: str
    ) -> tuple[str, ...]:
        """Selectors in one structural state that *this result* declares.

        ``selectors_in_state`` answers for the program; a jobtype
        declaration is the semantic claim about what a value means for
        the job that actually ran, and the two questions had two
        answers.  A host organ asking for a role must get the
        intersection, or the role is served by a selector nobody audited
        for this jobtype: ORCA's IRC log prints only its starting
        structure, so no selector there can honestly answer
        ``as_reached``, and the recovery route that asked the program-
        level question alone would have been handed the transition state
        labelled as the structure the path reached -- which is the defect
        the IRC selector restriction exists to prevent, arriving through
        a different door.

        A reader that declares no jobtype coverage at all (the xyz
        reader over registered geometry artifacts) is ungated here for
        the same reason it is ungated at extraction: its values carry no
        jobtype semantics to misread.
        """

        named = self.selectors_in_state(state)
        if not self.jobtype_selectors:
            return named
        jobtype = str(getattr(output, "jobtype", "") or "").strip().lower()
        declared = self.selectors_for_jobtype(jobtype)
        if declared is None:
            return ()
        return tuple(name for name in named if name in declared)

    def selectors_for_jobtype(self, jobtype: str) -> tuple[str, ...] | None:
        """Return exact declared coverage, or ``None`` when it is unknown."""

        normalized = str(jobtype).strip().lower()
        return next(
            (
                selectors
                for declared_jobtype, selectors in self.jobtype_selectors
                if declared_jobtype == normalized
            ),
            None,
        )

    def result_jobtypes_for_stage(self, stage: str) -> tuple[str, ...]:
        """The result words one planned ChemSmart stage's results answer to.

        A plan names a stage with the word the public CLI takes, and a
        reader is keyed on the word a finished log says it is.  They are
        usually the same word and this answers ``(stage,)``.  They are not
        always: ChemSmart runs one Gaussian ``irc`` job as two native
        one-direction inputs, and each log's route says ``ircf`` or
        ``ircr``, so the stage word appears in no Gaussian log ever
        written by this hub.

        Declaring the correspondence here rather than teaching each
        consumer the stage word keeps it one fact.  It was three: the
        plan-time selector gate and the geometry-producer predicate both
        assumed stage-is-key and so found a Gaussian IRC unreadable and
        unusable, while the preview verifier spelled the branch pair out
        by hand.
        """

        normalized = str(stage).strip().lower()
        declared = next(
            (
                jobtypes
                for declared_stage, jobtypes in self.stage_result_jobtypes
                if declared_stage == normalized
            ),
            None,
        )
        if declared is not None:
            return declared
        if self.selectors_for_jobtype(normalized) is None:
            return ()
        return (normalized,)

    def selectors_for_stage(self, stage: str) -> tuple[str, ...] | None:
        """Coverage every result of one planned stage declares.

        The intersection, not the union: a stage that can produce more
        than one result promises only what all of them carry, because a
        plan built before any of them exists cannot know which one an
        extraction will be handed.  For Gaussian's two IRC branches the
        two declarations are identical, so the intersection costs the
        stage nothing and states the rule for the next stage that splits.
        """

        jobtypes = self.result_jobtypes_for_stage(stage)
        if not jobtypes:
            return None
        declared = [self.selectors_for_jobtype(item) for item in jobtypes]
        if any(item is None for item in declared):
            return None
        shared = set(declared[0] or ())
        for item in declared[1:]:
            shared &= set(item or ())
        return tuple(sorted(shared))

    def available_selectors(self, output: Any) -> tuple[str, ...]:
        """Return the selectors this one opened result actually resolves.

        A program's selector set is what its parser implements.  What any
        single result carries is narrower, and the difference is a property of
        the job that produced it rather than of the program: a single point has
        no Hessian, a spin-restricted run prints no ``<S^2>``, and a job run
        without a frequency step has no thermochemistry.  ``jobtype_selectors``
        can only ever describe a job type, so it cannot answer this question
        for a specific artifact; probing the accessors can, exactly.

        Reporting the inventory keeps a reader from discovering the shape of a
        result one failed extraction at a time, and keeps the host descriptive:
        it states what this artifact holds without choosing which quantity the
        scientific question needs.
        """

        available = []
        for selector in sorted(self.accessors):
            try:
                self.read(output, selector)
            except Exception:  # noqa: BLE001 - absence is the answer here
                continue
            available.append(selector)
        return tuple(available)

    def read(self, output: Any, selector: str) -> tuple[Any, str]:
        """Return ``(value, unit)`` for ``selector``.

        A selector this reader does not implement is refused by name, and a
        quantity this particular run did not produce -- a single point has no
        Gibbs energy and no optimised geometry -- is refused as absent.  Both
        are stated gaps rather than a silent ``None`` or a parser traceback
        leaking out of the extraction path.
        """

        if selector not in self.accessors:
            raise ValueError(
                f"{self.program} result reader does not provide {selector!r}; "
                f"it provides {sorted(self.selectors)}"
            )
        from chemsmart.analysis.result_quantities import (
            QuantityExtractionError,
        )

        try:
            value = self.accessors[selector](output)
        except (MissingQuantityError, QuantityExtractionError):
            # Absence and divergence are different answers.  A block this run
            # never wrote is an absent quantity; a stored unit that is not the
            # one we read the dataset as means the writer's contract and this
            # reader have diverged, and reporting that as "no value" would
            # hide a defect behind the ordinary meaning of a missing block.
            raise
        except Exception as exc:
            # The parsers raise IndexError/TypeError when a block the run never
            # wrote is requested. That is an absent quantity, not a defect.
            raise MissingQuantityError(
                f"{self.program} result contains no {selector!r} value "
                f"({type(exc).__name__})"
            ) from exc
        if value is None or (
            isinstance(value, (list, tuple, Mapping, set, frozenset))
            and not value
        ):
            raise MissingQuantityError(
                f"{self.program} result contains no {selector!r} value"
            )
        if isinstance(value, Mapping):
            # Every program's per-atom parsers key by an atom label, and the
            # labels do not agree: ORCA and Gaussian use a global 1-based
            # index while xTB uses a per-element counter, so CO2 with atom
            # order O,O,C is {"O1","O2","C1"} there and {"O1","O2","C3"}
            # elsewhere -- the same key names a different atom. The typed
            # layer cannot repair that later either: freezing a mapping sorts
            # it by label, which reorders per-atom data out of molecular
            # order, and the resulting object is accepted as a "matrix" whose
            # cells are half strings and only fails at first arithmetic.
            #
            # So a mapping never leaves an accessor. Normalise it to a
            # positional vector in molecular order and pair it with symbols,
            # which is what connectivity already does. Refusing here keeps
            # that a typed contract failure rather than an uncaught TypeError
            # deep in numpy.
            from chemsmart.analysis import result_quantities as rq

            raise rq.QuantityExtractionError(
                f"{self.program} accessor for {selector!r} returned an "
                "atom-labelled mapping; per-atom quantities must be "
                "normalised to a positional vector in molecular order "
                "before they leave the accessor, because atom-label "
                "schemes differ between programs"
            )
        return value, self.source_units.get(selector, SELECTOR_UNITS[selector])


def _text_output_accessors(
    *,
    thermochemistry: bool = True,
    mode_composition: bool = False,
) -> dict[str, Callable[[Any], Any]]:
    """Accessors shared by the log-parsing programs.

    ``mode_composition`` is opt-in rather than shared because an accessor is
    what ``registered_reader_selectors`` shows the model.  Gaussian's own
    displacement block comes in variants this reader cannot yet tell apart
    -- ``freq=HPModes`` prints a second, higher-precision block, and
    ``freq=raman`` moves the row header -- and we never run Gaussian, so the
    variant is the user's choice and not an observable of ours.  Listing the
    quantity for a reader that may be looking at the wrong block would
    advertise support we cannot stand behind, so Gaussian is left out until
    the reader can detect what it is reading.
    """

    accessors: dict[str, Callable[[Any], Any]] = {
        "energy": _last_energy,
        "energies": lambda output: [float(item) for item in output.energies],
        "vibrational_frequencies": lambda output: [
            float(item) for item in output.vibrational_frequencies
        ],
        "positions": _positions,
        "connectivity": lambda output: _connectivity_matrix(output.molecule),
        "symbols": _symbols,
        "charge": lambda output: int(output.charge),
        "multiplicity": lambda output: int(output.multiplicity),
    }
    if mode_composition:
        accessors["vibrational_mode_atom_participation"] = (
            _vibrational_mode_atom_participation
        )
        accessors["vibrational_mode_degeneracy_group"] = (
            _vibrational_mode_degeneracy_group
        )
    if thermochemistry:
        accessors["gibbs_free_energy"] = lambda output: float(
            output.gibbs_free_energy
        )
    return accessors


def _orca_output(path: Path) -> Any:
    from chemsmart.io.orca.output import ORCAOutput

    return ORCAOutput(filename=str(path))


def _orca_total_energy(output: Any) -> float:
    """Return the final total ORCA energy, including post-SCF correlation.

    The last explicit ORCA final-energy record is the program's total result
    and remains correct for DFT, empirical-dispersion, and post-HF jobs.
    """

    values: list[float] = []
    for line in getattr(output, "contents", ()):
        if "FINAL SINGLE POINT ENERGY" not in str(line).upper():
            continue
        try:
            values.append(float(str(line).split()[-1]))
        except (TypeError, ValueError):
            continue
    if values:
        return values[-1] - _orca_fixed_geometry_root_shift(output)
    return _last_energy(output)


def _orca_fixed_geometry_root_shift(output: Any) -> float:
    """What ORCA added to a spectrum's final energy, which ``energy`` removes.

    A fixed-geometry %tddft run prints ``E(tot) = E(SCF) + DE(CIS)`` of its
    ``IRoot`` (root 1 unless set) as ``FINAL SINGLE POINT ENERGY``, even
    with ``Follow IRoot ... off``: on the archived water TDA that is
    -76.080796713 Eh beside a reference of -76.358315131
    (orca_differential/water_td_b3lypg.out).  ``energy`` on a spectrum is
    the surface the job computed on, the reference, as PySCF's ``td``
    answers it; so the printed excitation is taken back off.  A moving
    job on a root -- an excited-state optimisation -- keeps ORCA's total,
    because that root is the surface it walked.
    """

    if getattr(output, "jobtype", None) != "td":
        return 0.0
    shift = 0.0
    pattern = re.compile(r"^\s*DE\(CIS\)\s*=\s*(-?\d+\.\d+)\s*Eh")
    for line in getattr(output, "contents", ()):
        match = pattern.match(str(line))
        if match:
            shift = float(match.group(1))
    return shift


def _orca_scf_energy(output: Any) -> float:
    # ``ORCAOutput.final_scf_energy`` currently treats an empty optimization
    # slice as an optimization result and can therefore return ``None`` for a
    # perfectly valid correlated single point.  Read the last explicit TOTAL
    # SCF ENERGY value first; for post-HF output that is the reference energy
    # paired with the final correlated total below it.
    lines = tuple(str(line) for line in getattr(output, "contents", ()))
    number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
    cbs_pattern = re.compile(
        rf"^\s*Extrapolated CBS SCF energy\b.*?:\s*(?P<value>{number})\b",
        re.IGNORECASE,
    )
    cbs_values = [
        float(match.group("value"))
        for line in lines
        if (match := cbs_pattern.match(line)) is not None
    ]
    if cbs_values:
        return cbs_values[-1]

    values: list[float] = []
    for line in lines:
        text = str(line)
        if "Total Energy" not in text or ":" not in text:
            continue
        fields = text.split()
        try:
            colon = fields.index(":")
            values.append(float(fields[colon + 1]))
        except (ValueError, IndexError):
            continue
    if values:
        return values[-1]
    value = getattr(output, "final_scf_energy", None)
    if value is None:
        raise MissingQuantityError("ORCA result contains no final SCF energy")
    return float(value)


def _orca_correlation_energy(output: Any) -> float:
    """Return ORCA's final native post-SCF correlation energy.

    DFT-D's separately printed D3/D4 correction is not an electronic
    correlation energy.  A DFT result therefore cannot satisfy this selector
    merely because its total differs from its SCF component.  Parse the native
    MP2/CC result record itself rather than inferring correlation from a method
    token or subtracting independently rounded total-energy components.

    ORCA prints progressively more final records for correlated workflows.  A
    later CBS extrapolation supersedes the preceding finite-basis result; a
    later final CC correlation energy (including a printed triples correction)
    supersedes corrected CCSD; and corrected CCSD supersedes its preceding MP2
    record.  Compound outputs may then start another job, so chronology across
    all supported native records—not record-class priority—identifies the last
    completed correlation result.
    """

    number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
    record_patterns = (
        re.compile(
            rf"^\s*Extrapolated CBS correlation energy\b.*?:\s*"
            rf"(?P<value>{number})\b",
            re.IGNORECASE,
        ),
        re.compile(
            rf"^\s*Final correlation energy\s+(?:\.\.\.\s+)?"
            rf"(?P<value>{number})(?:\s+Eh)?\s*$",
            re.IGNORECASE,
        ),
        re.compile(
            rf"^\s*E\(CORR\)\(corrected\)\s+(?:\.\.\.\s+)?"
            rf"(?P<value>{number})\b",
            re.IGNORECASE,
        ),
        re.compile(
            rf"^\s*(?:RI-)?(?:DLPNO-)?MP2 correlation energy\s*"
            rf"(?:\.\.\.|:)\s*(?P<value>{number})(?:\s+Eh)?\s*$",
            re.IGNORECASE,
        ),
    )
    lines = tuple(str(line) for line in getattr(output, "contents", ()))
    records: list[tuple[int, float]] = []
    for line_index, line in enumerate(lines):
        for pattern in record_patterns:
            match = pattern.match(line)
            if match is not None:
                records.append((line_index, float(match.group("value"))))
                break
    if not records:
        raise MissingQuantityError(
            "ORCA result contains no explicit final post-SCF correlation "
            "energy record"
        )
    return records[-1][1]


def _orca_mayer_atoms_checked(output: Any) -> list[tuple[Any, ...]]:
    """The last Mayer table, in molecular order, checked against symbols."""

    from chemsmart.analysis import result_quantities as rq

    rows = list(getattr(output, "mayer_atom_rows", None) or ())
    if not rows:
        raise MissingQuantityError(
            "this ORCA result printed no Mayer population analysis"
        )
    symbols = _orca_symbols(output)
    if [row[0] for row in rows] != list(range(len(symbols))) or [
        str(row[1]).capitalize() for row in rows
    ] != [str(symbol).capitalize() for symbol in symbols]:
        raise rq.QuantityExtractionError(
            "the Mayer table's atoms are not this molecule's atoms in "
            "order; a per-atom vector is not reordered after the fact"
        )
    return rows


def _orca_mayer_bond_orders(output: Any) -> list[list[Any]]:
    """Mayer bond orders ORCA printed, as ``[atom_i, atom_j, order]`` rows.

    Printed by default beneath every ORCA population analysis (238
    archived outputs) and served by no reader: the bond order of a partial
    bond at a saddle, of a delocalised radical, of a metal-ligand bond.
    Zero-based molecular atom order, i < j.  Sparse as printed: ORCA lists
    only orders above 0.1, so an omitted pair has no value here, not zero.
    The last block belongs to the final density.
    """

    from chemsmart.analysis import result_quantities as rq

    n_atoms = len(_orca_mayer_atoms_checked(output))
    pairs = list(getattr(output, "mayer_bond_order_rows", None) or ())
    if not pairs:
        raise MissingQuantityError(
            "this ORCA result printed no Mayer bond order above 0.1"
        )
    rows: list[list[Any]] = []
    seen: set[tuple[int, int]] = set()
    for atom_i, atom_j, order in pairs:
        i0, j0 = sorted((int(atom_i), int(atom_j)))
        if i0 == j0 or i0 < 0 or j0 >= n_atoms or (i0, j0) in seen:
            raise rq.QuantityExtractionError(
                f"ORCA's Mayer bond-order list names an impossible pair "
                f"({atom_i}, {atom_j}) for {n_atoms} atoms"
            )
        seen.add((i0, j0))
        rows.append([i0, j0, float(order)])
    rows.sort(key=lambda row: (row[0], row[1]))
    return rows


def _orca_mayer_free_valence(output: Any) -> list[float]:
    """Mayer's free valence per atom (FA), in molecular order.

    Zero for a closed-shell atom whose valence is all bonding; on an open
    shell it measures the unpaired population an atom carries (the Fe(II)
    quintet archived here: 3.66 on iron).  Served beside the bond orders
    it completes: total valence = bonded valence + free valence.
    """

    return [float(row[-1]) for row in _orca_mayer_atoms_checked(output)]


def _orca_t1_diagnostic(output: Any) -> float:
    """The T1 diagnostic of the last coupled-cluster calculation printed.

    The standard single-reference check of a CCSD(T) number (Lee and
    Taylor's T1: the singles-amplitude norm over the square root of twice
    the number of correlated electrons).  ORCA prints it for every
    canonical and DLPNO coupled-cluster calculation; 16 archived outputs
    carry it and no reader served it, so a delivered CCSD(T) energy could
    not be questioned from its own result.  A basis-set extrapolation
    prints one per basis, and the last is the largest basis's.
    """

    values = list(getattr(output, "t1_diagnostics", None) or ())
    if not values:
        raise MissingQuantityError(
            "this ORCA result printed no T1 diagnostic: no coupled-cluster "
            "calculation ran"
        )
    return float(values[-1])


def _orca_dispersion_energy(output: Any) -> float:
    """Return the final explicit empirical dispersion correction."""

    value = getattr(output, "final_dispersion_energy", None)
    if value is None:
        raise MissingQuantityError(
            "ORCA result contains no empirical dispersion correction"
        )
    return float(value)


def _orca_auxiliary_basis(output: Any) -> str:
    value = getattr(output, "aux_basis", None)
    if not value:
        raise MissingQuantityError(
            "ORCA final route contains no explicit auxiliary basis"
        )
    return str(value)


def _orca_auxiliary_basis_role(output: Any) -> str:
    value = output.route_object.auxiliary_basis_role
    if not value or value == "missing":
        raise MissingQuantityError(
            "ORCA final route contains no explicit auxiliary-basis role"
        )
    return str(value)


def _orca_solvation_context(output: Any) -> tuple[str, str | None]:
    """Read solvent treatment from the final echoed ORCA job route."""

    routes = tuple(getattr(output, "_input_route_strings", ()))
    if not routes:
        raise MissingQuantityError(
            "ORCA result contains no echoed final route for solvent evidence"
        )
    final_route = str(routes[-1]).lower()
    # Simple-input keywords are whitespace-delimited route tokens.  Match a
    # complete token so a different model such as ``CPCM-X(water)`` cannot be
    # silently shortened to CPCM.
    # The leading boundary admits the route's own "!" prefix: "!SMD(water)"
    # with no space after the bang is the common ORCA idiom, and requiring
    # preceding whitespace made exactly that input report as gas phase --
    # the unrecognized-treatment guard below never fired either, because
    # "smd" contains none of its substrings.  Found on review, not in a run:
    # hub-written inputs always space the bang; user-supplied results need
    # not.
    match = re.search(
        r"(?<![a-z0-9_-])(?P<model>smd|cpcmc|cpcm|cosmors)"
        r"(?:\((?P<solvent>[^)]+)\))?(?!\S)",
        final_route,
    )
    if match is None:
        if re.search(
            r"(?i)(?:"
            r"\b\S*pcm\S*|\b\S*solv\S*|\b\S*cosmo\S*|\bsmd\b|"
            r"\b(?:alpb|ddcosmo|cpcmx|gbsa|tmcosmo)(?:\([^)]*\))?"
            r")",
            final_route,
        ):
            raise MissingQuantityError(
                "ORCA final route uses an unrecognized solvent treatment"
            )
        return "gas_phase", None

    model = match.group("model")
    solvent = match.group("solvent")
    if solvent:
        return model, solvent.strip().lower()

    lines = tuple(str(line) for line in getattr(output, "contents", ()))
    markers = tuple(getattr(output, "_orca_job_markers", ()))
    start = markers[-1][0] if markers else 0
    pattern = re.compile(
        r"^\s*Solvent:\s*(?:\.\.\.\s*)?(?P<name>\S.*?)\s*$",
        re.I,
    )
    values = []
    for line in lines[start:]:
        native_match = pattern.match(line)
        if native_match is not None:
            values.append(native_match.group("name").strip().lower())
    return model, values[-1] if values else None


def _orca_solvation_model(output: Any) -> str:
    return _orca_solvation_context(output)[0]


def _chemical_core_table() -> tuple[int, ...]:
    table = [0] * 119
    for first, last, orbitals in (
        (5, 12, 1),
        (13, 30, 5),
        (31, 38, 9),
        (39, 48, 14),
        (49, 70, 18),
        (71, 80, 23),
        (81, 103, 34),
        (104, 112, 50),
        (113, 118, 55),
    ):
        for number in range(first, last + 1):
            table[number] = orbitals
    return tuple(table)


#: The chemical core, in orbitals, of each element by atomic number: the
#: rule PySCF 2.14 applies for ``frozen_core: auto`` (``chemcore_atm`` in
#: pyscf/data/elements.py, read on the CUHK deployment 2026-09-24), which
#: ORCA 6.1.1's default frozen core matched on every molecule measured --
#: water, CH2O, CH3, CH3Cl, HBr, ZnH2, HI, CH3I, I (CUHK Slurm 2149277,
#: 2149487, 2151772, 2151773).  Gaussian's default matches it there too,
#: except zinc (9 orbitals against 5).  A core potential's orbitals leave
#: the count, as PySCF's own rule removes them.
CHEMICAL_CORE_ORBITALS = _chemical_core_table()

#: The frozen-core word of an operand whose molecule has no core to
#: freeze (H, He, Li, Be): every convention freezes nothing there.
FROZEN_CORE_NO_CORE = "no_core"


def chemical_core_orbitals(
    symbols: Any, ecp_core_electrons: Mapping[str, int] | None = None
) -> int:
    """Orbitals the chemical-core rule freezes in one molecule."""

    from ase.data import atomic_numbers

    cores = ecp_core_electrons or {}
    total = 0
    for symbol in symbols:
        orbitals = CHEMICAL_CORE_ORBITALS[int(atomic_numbers[str(symbol)])]
        replaced = int(cores.get(str(symbol), 0)) // 2
        total += 0 if replaced > orbitals else orbitals - replaced
    return total


def frozen_core_conventions(
    frozen: int,
    *,
    symbols: Any,
    ecp_core_electrons: Mapping[str, int] | None,
    program_rule: str | None,
) -> tuple[str, ...]:
    """The conventions one frozen-core count is consistent with.

    A count is a fact about one molecule -- HI freezes 4 orbitals where H
    freezes none under one rule -- so two results of one level compare
    by rule, never by count: ``chemical_core`` when the count is the
    molecule's chemical core, ``all_electrons`` when it is zero, and the
    rule the program applied by default (``orca_default``,
    ``gaussian_default``, ``pyscf_default``, ``pyscf_auto``) when it was
    not told one.  A molecule with no core is consistent with every rule.
    """

    chemical = chemical_core_orbitals(symbols, ecp_core_electrons)
    if chemical == 0:
        return (FROZEN_CORE_NO_CORE,)
    words = set()
    if int(frozen) == chemical:
        words.add("chemical_core")
    if int(frozen) == 0:
        words.add("all_electrons")
    if program_rule:
        words.add(program_rule)
    return tuple(sorted(words))


def _add_frozen_core_conventions(
    level: dict[str, Any], symbols: Any, program_rule: str | None
) -> None:
    """State, beside a correlated level's frozen count, the rules it fits.

    A molecule with no core is consistent with every rule whether or not
    the program printed a count for it (ORCA prints none for a hydrogen
    atom); otherwise the stated count is classified, and a level that
    states no count states no convention.
    """

    method = str(level.get("ab_initio") or "").lower()
    cores = level.get("ecp_core_electrons")
    if (
        not any(marker in method for marker in _CORRELATED_METHOD_MARKERS)
        or cores is None
        or not symbols
    ):
        return
    if chemical_core_orbitals(symbols, cores) == 0:
        level["frozen_core_conventions"] = (FROZEN_CORE_NO_CORE,)
        return
    frozen = level.get("frozen_core")
    if frozen is None:
        return
    level["frozen_core_conventions"] = frozen_core_conventions(
        int(frozen),
        symbols=symbols,
        ecp_core_electrons=cores,
        program_rule=program_rule,
    )


def _molecule_symbols(output: Any) -> list[str]:
    try:
        return [str(symbol) for symbol in output.molecule.chemical_symbols]
    except Exception:  # noqa: BLE001 - no structure, no symbols
        return []


def _normalized_dispersion(value: Any) -> str:
    """One word for an empirical dispersion, whatever program spelled it."""

    word = str(value or "").strip().lower()
    if word.startswith("empiricaldispersion="):
        word = word.split("=", 1)[1]
    return {
        "": "none",
        "gd3bj": "d3bj",
        "gd3": "d3zero",
        "d3": "d3zero",
        "gd2": "d2",
    }.get(word, word)


def _orca_scf_reference(output: Any) -> str | None:
    """The determinant ORCA ran, from its own SCF settings block.

    ORCA prints ``Hartree-Fock type HFTyp .... RHF|UHF|ROHF`` and, above
    it, ``Density Functional Method .... DFT(GTOs)`` or ``Ab initio
    Hamiltonian Method .... Hartree-Fock(GTOs)``; the last block printed
    is the SCF the result's numbers come from.  None where neither was
    printed.
    """

    family = None
    kohn_sham = None
    for line in getattr(output, "contents", ()) or ():
        text = str(line)
        match = re.search(r"Hartree-Fock type\s+HFTyp\s+\.+\s+(\S+)", text)
        if match:
            family = match.group(1).strip().casefold()
            continue
        if re.search(r"Density Functional\s+Method\s+\.+", text):
            kohn_sham = True
        elif re.search(r"Ab initio Hamiltonian\s+Method\s+\.+", text):
            kohn_sham = False
    if family not in {"rhf", "uhf", "rohf"} or kohn_sham is None:
        return None
    return family[:-2] + "ks" if kohn_sham else family


def _orca_level(output: Any) -> dict[str, Any]:
    """The Hamiltonian an ORCA result computed with, from its own record.

    The applied functional (the literal its printed form names), the
    correlated method, the basis, the dispersion and the continuum; the
    frozen core a correlated stage applied, in orbitals, from ORCA's own
    ``NCore``/``chemical core (N el)`` line.  Numerics (grid, RI, COSX)
    are not a level.
    """

    level: dict[str, Any] = {}
    try:
        level["functional"] = _orca_functional(output)
    except Exception:  # noqa: BLE001 - a post-HF run applies no functional
        pass
    ab_initio = getattr(output, "ab_initio", None)
    if ab_initio:
        level["ab_initio"] = str(ab_initio).lower()
    basis = getattr(output, "basis", None)
    if basis:
        level["basis"] = str(basis).lower()
    if "functional" in level:
        level["dispersion"] = _normalized_dispersion(
            getattr(output, "dispersion", None)
        )
    try:
        model, solvent = _orca_solvation_context(output)
    except MissingQuantityError:
        model, solvent = None, None
    if model not in (None, "gas_phase"):
        level["solvent_model"] = model
        if solvent:
            level["solvent"] = solvent
    if ab_initio:
        electrons = None
        for line in getattr(output, "contents", ()):
            match = re.search(
                r"Freezing NCore=(\d+)|chemical core \((\d+) el\)", str(line)
            )
            if match:
                electrons = int(match.group(1) or match.group(2))
        if electrons is not None:
            level["frozen_core"] = electrons // 2
    # What the basis was built of: ORCA's angular form is a fact about
    # ORCA, and its core potentials are named where it assigns them.
    from chemsmart.jobs.orca.settings import ORCA_BASIS_FUNCTIONS

    if "basis" in level:
        level["basis_functions"] = ORCA_BASIS_FUNCTIONS
    try:
        cores = output.ecp_core_electrons
    except Exception:  # noqa: BLE001 - an unreadable record states nothing
        cores = None
    if cores is not None:
        level["ecp_core_electrons"] = dict(cores)
    # ORCA freezes its chemical core unless its input says otherwise.
    told = any(
        str(line).lstrip().startswith("|")
        and "frozencore" in str(line).lower()
        for line in getattr(output, "contents", ())
    )
    _add_frozen_core_conventions(
        level, _molecule_symbols(output), None if told else "orca_default"
    )
    # The response an excited stage ran on, in the words Gaussian's and
    # PySCF's levels use: ORCA's own header says which approximation it
    # applied and how many roots of each block it determined, and the
    # manifold is the word the request named (ORCA spells ``triplet`` and
    # ``singlet_triplet`` alike).  TDA and full TD-DFT are two calculations
    # of the same roots, and ORCA's level had said neither.
    applied = getattr(output, "excited_state_applied", None)
    if isinstance(applied, Mapping) and applied.get("response_method"):
        level["response_method"] = applied["response_method"]
        manifold = getattr(output, "state_manifold", None)
        if manifold:
            level["state_manifold"] = manifold
        if applied.get("nstates") is not None:
            level["nstates"] = int(applied["nstates"])
    try:
        requested = bool(getattr(output, "broken_symmetry", False))
    except Exception:  # noqa: BLE001 - an unreadable echo requests nothing
        requested = False
    _add_reference_identity(level, _orca_scf_reference(output), requested)
    return level


def _gaussian_level(output: Any) -> dict[str, Any]:
    """The Hamiltonian a Gaussian result computed with, from its own record.

    The applied functional (the name ``SCF Done`` states), the correlated
    method whose total ``energy`` is, the basis, the empirical dispersion,
    the continuum, and the frozen orbitals Gaussian printed (``NFC=``) for
    a correlated method.
    """

    level: dict[str, Any] = {}
    try:
        functional = _gaussian_functional(output)
    except Exception:  # noqa: BLE001 - a post-HF run applies no functional
        functional = None
    if functional:
        level["functional"] = re.sub(
            r"-d[234](?:bj|zero)?$", "", str(functional).lower()
        )
    correlated = getattr(output, "correlated_method", None)
    if correlated and correlated != "double_hybrid":
        level["ab_initio"] = str(correlated)
    basis = getattr(output, "basis", None)
    if basis:
        level["basis"] = str(basis).lower()
    if "functional" in level:
        level["dispersion"] = _normalized_dispersion(
            getattr(output, "dispersion", None)
        )
    model = getattr(output, "solvent_model", None)
    if model:
        level["solvent_model"] = str(model).lower()
        solvent = getattr(output, "solvent_id", None)
        if solvent:
            level["solvent"] = str(solvent).lower()
    if correlated:
        frozen = None
        for line in getattr(output, "contents", ()):
            match = re.search(r"\bNFC=\s*(\d+)", str(line))
            if match:
                frozen = int(match.group(1))
        if frozen is not None:
            level["frozen_core"] = frozen
    # What the basis was built of, in Gaussian's own statement: the
    # angular form it printed beside the basis, and the electrons its
    # core potentials replaced.
    form = getattr(output, "basis_angular_form", None)
    if form and "basis" in level:
        level["basis_functions"] = form
    try:
        cores = output.ecp_core_electrons
    except Exception:  # noqa: BLE001 - an unreadable record states nothing
        cores = None
    if cores is not None:
        level["ecp_core_electrons"] = dict(cores)
    # Gaussian's FC is its default; a route that names another rule
    # (Full, a window, a named freeze) chose its own.
    route_words = set(
        re.split(
            r"[^a-z0-9]+", str(getattr(output, "route_string", "")).lower()
        )
    )
    told = bool(
        route_words
        & {
            "full",
            "window",
            "readwindow",
            "rw",
            "freezeg2",
            "freezeg3",
            "freezeg4",
            "freezenoblegascore",
            "freezeinnernoblegascore",
        }
    )
    if correlated and correlated != "double_hybrid":
        _add_frozen_core_conventions(
            level,
            _molecule_symbols(output),
            None if told else "gaussian_default",
        )
    # The response an excited stage ran on, in the words PySCF's level
    # uses: full TD-DFT and TDA are different calculations of the same
    # roots, and the route is where Gaussian says which ran.
    request = getattr(output, "excited_state_request", None)
    if isinstance(request, Mapping):
        for name in ("response_method", "state_manifold", "nstates"):
            if request.get(name) is not None:
                level[name] = request[name]
        root = getattr(output, "excited_state_followed_root", None)
        if root is not None:
            level["excited_state_root"] = int(root)
    try:
        requested = bool(getattr(output, "broken_symmetry", False))
    except Exception:  # noqa: BLE001 - an unreadable route requests nothing
        requested = False
    _add_reference_identity(level, _gaussian_scf_reference(output), requested)
    return level


def _add_reference_identity(level, reference, broken_symmetry):
    """State, on a level, the determinant that ran and the request behind it.

    ``reference`` is the family the program itself printed (a word of
    ``SCF_REFERENCE_WORDS``); ``broken_symmetry`` is true where the
    program's own record shows the broken-symmetry request applied.  Both
    are shown and never compared: a closed-shell molecule and a radical in
    one reaction energy run different determinants by necessity, so the
    reference is not a ``LEVEL_IDENTITY_FIELDS`` field, and an operation is
    never told its operands differ in level because their states do.
    """

    if reference in SCF_REFERENCE_WORDS:
        level["reference"] = reference
    if broken_symmetry is True:
        level["broken_symmetry"] = True


def _orca_solvent(output: Any) -> str:
    model, solvent = _orca_solvation_context(output)
    if model == "gas_phase":
        raise MissingQuantityError("ORCA final job is explicitly gas phase")
    if not solvent:
        raise MissingQuantityError(
            "ORCA final solvated job contains no explicit solvent name"
        )
    return solvent


def _orca_positions(output: Any) -> list[list[float]]:
    molecule = output.thermochemistry_molecule
    return [[float(value) for value in row] for row in molecule.positions]


def _orca_reached_positions(output: Any) -> list[list[float]]:
    """The structure the run actually reached, not the one it was handed.

    ``positions`` reads ``thermochemistry_molecule`` -- the geometry the
    Hessian was computed at, which for an ``OptTS Freq`` is step 0,
    because ChemSmart's ORCA ts writer computes the Hessian first. On a
    live unconverged saddle search the two sat **1.231588 A** apart, and
    ``build_reached_geometry`` -- the recovery route the repair menu
    offers for exactly that ending -- asked for ``positions`` and wrote
    "the structure orca reached" onto the seed coordinates. A session
    caught it by measuring the bytes and recorded the route as spent.

    ``output.molecule`` is the last printed structure, which is what
    "reached" means for an optimiser that ran and did not converge.
    """

    molecule = output.molecule
    if isinstance(molecule, (list, tuple)):
        if not molecule:
            raise MissingQuantityError(
                "this ORCA result prints no structure to have reached"
            )
        molecule = molecule[-1]
    positions = getattr(molecule, "positions", None)
    if positions is None:
        raise MissingQuantityError(
            "this ORCA result prints no reached structure"
        )
    return [[float(value) for value in row] for row in positions]


def _orca_symbols(output: Any) -> list[str]:
    return [
        str(item) for item in output.thermochemistry_molecule.chemical_symbols
    ]


def _orca_hirshfeld_charges(output: Any) -> list[float] | None:
    """Return the Hirshfeld charges as a positional vector, or ``None``.

    Parsing once matters -- the property rescans the whole file -- and
    absence must reach ``read`` as a value rather than as an exception, so
    that a run which never asked for the analysis is refused in those words
    instead of by naming a Python type at a scientist.
    """

    charges = output.hirshfeld_charges
    if charges is None:
        return None
    return _per_atom_vector(
        charges,
        _orca_symbols(output),
        quantity="hirshfeld_atomic_charges",
    )


def _orca_hirshfeld_spins(output: Any) -> list[float] | None:
    """The Hirshfeld spin populations, positionally, or ``None``.

    The second column of the block the ``Hirshfeld`` directive prints; on
    a doublet it closes on one unpaired electron (the methyl radical,
    CUHK 2152359: C 0.839, each H 0.054, total 1.000000).
    """

    spins = output.hirshfeld_spin_densities
    if spins is None:
        return None
    return _per_atom_vector(
        spins,
        _orca_symbols(output),
        quantity="hirshfeld_atomic_spin_populations",
    )


#: Hirshfeld's spin partition, declared by each program that prints it: the
#: same record the Mulliken and Loewdin spin populations carry.
_HIRSHFELD_SPIN_DECLARATION = (
    ("hirshfeld_atomic_spin_populations", "1", "DIMENSIONLESS"),
)
_HIRSHFELD_SPIN_ATOM_DECLARATION = (
    (
        "hirshfeld_atomic_spin_populations",
        tuple(
            {
                "semantic_quantity": "atomic_spin_population",
                "population_scheme": "Hirshfeld",
                "atom_order": "zero-based molecular atom order",
            }.items()
        ),
    ),
)


def _orca_channel_eigenvalues(
    output: Any, channel: str
) -> tuple[list[float], list[float]]:
    """Return one spin channel's occupied and virtual energies, both in eV.

    ORCA prints either a single restricted block or a spin-up/spin-down pair,
    and ``ORCAOutput`` already normalises both to eV.  For a restricted
    reference the two channels are the same doubly-occupied set, which is why
    the alpha accessors answer it unchanged.
    """

    if channel == "beta":
        occupied = output.beta_occ_eigenvalues
        virtual = output.beta_virtual_eigenvalues
    else:
        occupied = output.alpha_occ_eigenvalues
        virtual = output.alpha_virtual_eigenvalues
    return (
        [float(item) for item in occupied or ()],
        [float(item) for item in virtual or ()],
    )


def _orca_frontier_pair(output: Any) -> tuple[float, float]:
    """Return (HOMO, LUMO) in eV as extrema over every occupied channel.

    For an unrestricted reference the frontier orbitals need not share a spin
    channel -- the highest occupied level can be alpha while the lowest
    virtual is beta -- so the extremum over both channels is the definition
    that survives that case.  A question about one channel specifically is
    asked through the spin-resolved selectors instead.
    """

    channels = (
        ("alpha", "beta")
        if bool(getattr(output, "is_unrestricted", False))
        else ("alpha",)
    )
    occupied: list[float] = []
    virtual: list[float] = []
    for channel in channels:
        channel_occupied, channel_virtual = _orca_channel_eigenvalues(
            output, channel
        )
        occupied.extend(channel_occupied)
        virtual.extend(channel_virtual)
    if not occupied or not virtual:
        raise MissingQuantityError(
            "ORCA result establishes no occupied and virtual orbital pair; "
            "orbital energies are printed alongside the population analysis"
        )
    return max(occupied), min(virtual)


def _orca_frontier_gap(output: Any) -> float:
    homo, lumo = _orca_frontier_pair(output)
    return lumo - homo


def _orca_channel_frontier(output: Any, channel: str, occupied: bool) -> float:
    values = _orca_channel_eigenvalues(output, channel)[0 if occupied else 1]
    if not values:
        raise MissingQuantityError(
            f"ORCA result establishes no {'occupied' if occupied else 'virtual'}"
            f" {channel} orbital"
        )
    return max(values) if occupied else min(values)


def _route_functional(output: Any) -> str:
    """The DFT functional named in the echoed route; a semantic identity.

    Post-HF results record no functional -- that identity is ``ab_initio``.
    Refusing here, rather than returning an empty string, keeps "which
    method produced this number" a fact the host can compare across a
    series instead of a blank that equals every other blank.
    """

    value = getattr(output, "functional", None)
    if not value:
        raise MissingQuantityError(
            "this result records no DFT functional in its route; a post-HF "
            "method identity is read by 'ab_initio'"
        )
    return str(value)


def _orca_functional(output: Any) -> str:
    """The functional an ORCA run applied, in the ChemSmart vocabulary.

    The route parser answers the literal the route's keyword means
    (ORCA's ``B3LYP`` is ``b3lyp5``); the Hamiltonian block ORCA printed
    says which local correlation it actually ran, and where the writer's
    table names the label that form prints, the two must agree -- a run
    whose printed ``LDAOpt`` is not its route's form is not reported as
    either functional.
    """

    from chemsmart.jobs.orca.settings import orca_functional_lda_label

    value = _route_functional(output)
    expected = orca_functional_lda_label(value)
    printed = getattr(output, "lda_correlation", None)
    if expected and printed and printed.upper() != expected.upper():
        raise MissingQuantityError(
            f"the route names {value}, whose ORCA form prints LDAOpt "
            f"{expected}, and this run printed {printed}: the functional "
            "ORCA applied is not the one its route names"
        )
    return value


def _gaussian_functional(output: Any) -> str:
    """The functional a Gaussian run applied, in the ChemSmart vocabulary.

    The route word says what was asked; ``SCF Done:  E(R<name>)`` says
    what Gaussian ran.  They differ when Gaussian completed a route word
    to another keyword -- ``pbe0`` ran as ``RPBE0DH`` (CUHK Slurm 2149277)
    -- and then the applied name is the answer, so a result never reports
    the functional its route asked for in place of the one it computed.
    """

    from chemsmart.io.gaussian import GAUSSIAN_ALL_FUNCTIONALS
    from chemsmart.jobs.gaussian.settings import (
        GAUSSIAN_FUNCTIONAL_NATIVE,
        gaussian_functional_literal,
    )

    value = _route_functional(output)
    label = None
    for line in reversed(getattr(output, "contents", ()) or ()):
        match = re.search(r"SCF Done:\s+E\(([^)]+)\)", str(line))
        if match:
            label = match.group(1).strip().casefold()
            break
    if label is None or ":" in value:
        return value
    if label in {"rhf", "uhf", "rohf"}:
        raise MissingQuantityError(
            "this result ran a Hartree-Fock reference and applied no "
            "functional; the method identity is read by 'ab_initio'"
        )
    # Gaussian labels a pure combination with a hyphen (``RB-LYP``,
    # ``RB-P86``) that its keyword does not carry, so names are compared
    # without hyphens on both sides.
    label = label.replace("-", "")
    asked = re.sub(r"-d[234](?:bj|zero)?$", "", value.casefold())
    native = GAUSSIAN_FUNCTIONAL_NATIVE.get(asked, asked).casefold()
    native = native.replace("-", "")
    if label in {f"{prefix}{native}" for prefix in ("r", "u", "ro")}:
        return value
    keywords = {
        word.replace("-", ""): word for word in GAUSSIAN_ALL_FUNCTIONALS
    }
    for prefix in ("ro", "r", "u"):
        applied = label[len(prefix) :]
        if label.startswith(prefix) and applied in keywords:
            return gaussian_functional_literal(keywords[applied])
    return gaussian_functional_literal(
        label[1:] if label[:1] in "ru" else label
    )


#: Gaussian's semiempirical labels in ``SCF Done: E(...)``: NDDO
#: references, restricted or unrestricted like Hartree-Fock's.
_GAUSSIAN_SEMIEMPIRICAL_LABELS = (
    "am1",
    "cndo",
    "indo",
    "mndo",
    "pddg",
    "pm3",
    "pm6",
    "pm7",
    "zindo",
)


def _gaussian_scf_reference(output: Any) -> str | None:
    """The determinant Gaussian ran, from its own ``SCF Done`` label.

    ``E(RB3LYP)``, ``E(UB3LYP)``, ``E(ROB3LYP)``, ``E(RHF)``, ``E(UHF)``,
    ``E(ROHF)``: the prefix is what Gaussian ran, whatever the route asked
    for -- a bare method on a singlet runs restricted and on an open shell
    unrestricted, and a restricted route with ``guess=mix`` stays restricted
    (R10 Q18 O0, CUHK Slurm 2153330).  The functional reader strips this
    prefix to name the functional; it is kept here.  None where the log
    printed no ``SCF Done`` line.
    """

    label = None
    for line in reversed(getattr(output, "contents", ()) or ()):
        match = re.search(r"SCF Done:\s+E\(([^)]+)\)", str(line))
        if match:
            label = match.group(1).strip().casefold()
            break
    if not label:
        return None
    for prefix in ("ro", "u", "r"):
        if label.startswith(prefix):
            rest = label[len(prefix) :]
            hartree_fock = (
                rest == "hf" or rest in _GAUSSIAN_SEMIEMPIRICAL_LABELS
            )
            return f"{prefix}{'hf' if hartree_fock else 'ks'}"
    return None


def _route_ab_initio(output: Any) -> str:
    value = getattr(output, "ab_initio", None)
    if not value:
        raise MissingQuantityError(
            "this result records no post-HF method in its route; a DFT "
            "functional identity is read by 'functional'"
        )
    return str(value)


def _route_basis(output: Any) -> str:
    value = getattr(output, "basis", None)
    if not value:
        raise MissingQuantityError(
            "this result records no basis set in its echoed route"
        )
    return str(value)


def _irc_run_converged(output: Any) -> int:
    """1 only when every IRC branch printed its convergence marker.

    The first Agent-executed IRC stopped at ORCA's default iteration limit
    and validated with nothing typed able to see it; this selector was
    added only after both phrasings -- converged and exhausted -- were
    observed in real artifacts.
    """

    value = getattr(output, "irc_converged", None)
    if value is None:
        raise MissingQuantityError(
            "this result records no IRC convergence marker"
        )
    return int(bool(value))


def _optimization_converged(output: Any) -> int:
    """1 when the program printed its optimization-converged marker.

    ``False`` is an observation -- the program printed its own
    exhaustion marker ("did not converge ... maximum number of
    optimization cycles") -- and is served as 0. ``None`` means the log
    carries no convergence marker of either kind, a non-optimizing run,
    and extraction refuses rather than manufacturing a 0 that would
    read as an observed failure.
    """

    value = getattr(output, "converged", None)
    if value is None:
        raise MissingQuantityError(
            "this result records no optimization convergence marker"
        )
    return int(bool(value))


def _scan_steps_reached(output: Any) -> int:
    """How many relaxed-scan steps this run actually started."""
    value = getattr(output, "scan_step_count", None)
    if not value:
        raise MissingQuantityError(
            "this result announces no relaxed-scan steps"
        )
    return int(value)


def _scan_steps_planned(output: Any) -> int:
    """How many steps the scan declared before its first point."""
    coordinate = getattr(output, "scan_coordinate", None)
    if not coordinate:
        raise MissingQuantityError(
            "this result declares no driven scan coordinate"
        )
    return int(coordinate["points"])


#: How far two hartree values printed at different precisions may sit
#: apart and still be the same number.  ORCA prints the IRC path table to
#: six decimals and the geometry sidecars to twelve, so the pairing check
#: below compares them at the coarser one with room to spare; a wrong
#: pairing is millihartrees out, never microhartrees.
_ORCA_IRC_ENERGY_PAIRING_TOLERANCE = 1e-5


def _orca_irc_path_records(output: Any) -> tuple[Any, ...]:
    """ORCA's printed IRC path, refused when the log carries none."""

    records = tuple(getattr(output, "irc_path_records", ()) or ())
    if len(records) < 2:
        raise MissingQuantityError(
            "this ORCA result prints no IRC PATH SUMMARY table, so it "
            "establishes no reaction path"
        )
    return records


def _orca_irc_branch(output: Any) -> Any:
    """The one branch this IRC walked, or a refusal naming the route.

    ORCA writes a branch's final structure beside the log and a
    ``direction both`` run writes two of them.  Which end of a reaction
    coordinate a structure belongs to is then a scientific fact the host
    cannot settle, and serving either as "the structure the run reached"
    is how one geometry comes to be delivered under another's name.  One
    direction per irc node keeps the answer unambiguous, which is also
    the shape the qualified Agent route already plans.
    """

    records = tuple(getattr(output, "irc_endpoint_records", ()) or ())
    if not records:
        raise MissingQuantityError(
            "this ORCA IRC left no branch endpoint structure beside its "
            "output, so it establishes no path geometry"
        )
    if len(records) > 1:
        directions = ", ".join(str(item["direction"]) for item in records)
        raise MissingQuantityError(
            f"this ORCA IRC walked {len(records)} branches ({directions}), "
            "so no single structure is the one it reached; plan one "
            "direction per irc node and each branch's endpoint is its own"
        )
    return records[0]


def _orca_irc_start_molecule(output: Any) -> Any:
    """The transition state the branch was handed, checked against the path.

    An ORCA IRC log prints one structure and it is the starting point.
    That premise is verified rather than trusted: the printed structure
    must be the only one in the log and its energy must be the first row
    of the path table, which is the row ORCA itself marks ``<= TS``.
    """

    _orca_irc_branch(output)
    records = _orca_irc_path_records(output)
    structures = list(getattr(output, "all_structures", ()) or ())
    if len(structures) != 1:
        raise MissingQuantityError(
            "this ORCA IRC log prints "
            f"{len(structures)} structures, so which one the branch "
            "started from is not established"
        )
    printed_energy = getattr(output, "final_energy", None)
    if printed_energy is None or not math.isfinite(float(printed_energy)):
        raise MissingQuantityError(
            "this ORCA IRC records no energy for its printed structure"
        )
    if (
        abs(float(printed_energy) - float(records[0]["energy"]))
        > _ORCA_IRC_ENERGY_PAIRING_TOLERANCE
    ):
        raise MissingQuantityError(
            "this ORCA IRC's printed structure does not carry the path's "
            "first energy, so the log's structure is not the point the "
            "path starts from"
        )
    return structures[0]


def _orca_irc_end_molecule(output: Any) -> Any:
    """The structure the branch reached, read from ORCA's own endpoint file.

    The endpoint is the whole product of an IRC and it lives in a sidecar
    rather than in the log.  It is read through ChemSmart's XYZ parser and
    bound to the log by two checks: the atom identities and order must be
    the log's, and the energy the sidecar carries must be the last row of
    the path table.  A sidecar that belongs to another run fails both.
    """

    record = _orca_irc_branch(output)
    records = _orca_irc_path_records(output)
    path = Path(str(record["geometry_file"]))
    from chemsmart.io.xyz.xyzfile import XYZFile

    try:
        molecule = XYZFile(str(path)).get_molecules(index="-1")
        comment = str(XYZFile(str(path)).get_comments(index="-1") or "")
    except (OSError, IndexError, TypeError, ValueError) as error:
        raise MissingQuantityError(
            f"ORCA's IRC endpoint file {path.name} is not readable"
        ) from error
    if molecule is None:
        raise MissingQuantityError(
            f"ORCA's IRC endpoint file {path.name} holds no structure"
        )
    expected = tuple(
        str(item) for item in output.thermochemistry_molecule.chemical_symbols
    )
    if tuple(str(item) for item in molecule.chemical_symbols) != expected:
        raise MissingQuantityError(
            "ORCA's IRC endpoint file changes atom identity or atom order "
            "against the log it sits beside"
        )
    energy = None
    for pattern in _XYZ_HARTREE_PATTERNS:
        match = pattern.search(comment.strip())
        if match is not None:
            energy = float(match.group(1).replace("D", "E").replace("d", "e"))
            break
    if energy is None:
        raise MissingQuantityError(
            "ORCA's IRC endpoint file records no energy, so it cannot be "
            "bound to the path it is supposed to end"
        )
    if (
        abs(energy - float(records[-1]["energy"]))
        > _ORCA_IRC_ENERGY_PAIRING_TOLERANCE
    ):
        raise MissingQuantityError(
            "ORCA's IRC endpoint file does not carry the path's last "
            "energy, so it is not this path's endpoint"
        )
    return molecule


def _orca_constraint_records(output: Any) -> tuple[Any, ...]:
    """The internal coordinates ORCA held, checked against the molecule.

    The parser reads ORCA's own constraint table, whose definitions carry
    an element symbol beside each index.  Those labels are resolved
    against the molecule's own symbols rather than trusted, exactly as the
    per-atom population vectors are: a label that does not match is a
    reader that has drifted from the structure, and guessing is what the
    populations lesson cost this repository once already.
    """

    records = tuple(getattr(output, "constrained_coordinate_records", ()))
    if not records:
        # True of a plain optimisation and equally true of one that froze
        # Cartesian positions: ORCA constrains those atom by atom and
        # prints no internal coordinate for them, so this family has
        # nothing to say about such a run rather than nothing being held.
        raise MissingQuantityError(
            "this ORCA result held no internal coordinate; this family "
            "answers bonds, angles and dihedrals, and a Cartesian atom "
            "freeze is not one of them"
        )
    symbols = _orca_symbols(output)
    for record in records:
        for index, symbol in zip(record["atoms"], record["symbols"]):
            if not 1 <= int(index) <= len(symbols):
                raise MissingQuantityError(
                    f"constraint {record['label']} names atom {index}, "
                    f"outside this molecule's {len(symbols)} atoms"
                )
            if symbols[int(index) - 1] != str(symbol):
                raise MissingQuantityError(
                    f"constraint {record['label']} calls atom {index} "
                    f"{symbol!r} while the molecule has "
                    f"{symbols[int(index) - 1]!r} there"
                )
    return records


#: How far a held coordinate may sit from the value ORCA declared before
#: the reader stops calling it held.  ORCA holds a constraint to its own
#: optimiser tolerance: across this repository's archived constrained
#: results the largest disagreement between the declared value and the
#: returned structure is 5.5e-7 Angstrom and 1.1e-5 degrees.  These bounds
#: are three orders of magnitude above that and far below any change a
#: chemist would call a different structure, so they separate "the reader
#: and the geometry disagree" from optimiser noise and from chemistry.
_ORCA_CONSTRAINT_TOLERANCE = {"bond": 1e-3, "angle": 1e-2, "dihedral": 1e-2}

#: What a held internal coordinate is measured in, for every reader that
#: answers a constrained optimisation, so ORCA's and Gaussian's answers
#: are one declaration rather than two copies of it.
_CONSTRAINED_COORDINATE_DECLARATIONS = (
    ("constrained_angle_atoms", "1", "DIMENSIONLESS"),
    ("constrained_bond_angles", "degree", "ANGLE"),
    ("constrained_bond_atoms", "1", "DIMENSIONLESS"),
    ("constrained_bond_lengths", "Angstrom", "LENGTH"),
    ("constrained_coordinate_count", "1", "DIMENSIONLESS"),
    ("constrained_dihedral_angles", "degree", "ANGLE"),
    ("constrained_dihedral_atoms", "1", "DIMENSIONLESS"),
)
#: What each excited root is made of, served by all three td readers: the
#: root's largest single excitation in frontier-orbital words (``HOMO-1 ->
#: LUMO``, ``beta HOMO -> LUMO``) and the program's own weight of it.
_EXCITED_CHARACTER_DECLARATIONS = (
    ("excited_state_dominant_excitations", "", "DIMENSIONLESS"),
    ("excited_state_dominant_weights", "1", "DIMENSIONLESS"),
)
_CONSTRAINED_COORDINATE_ATOM_DECLARATIONS = (
    (
        "constrained_bond_atoms",
        (
            ("semantic_quantity", "constrained_internal_coordinate"),
            ("atom_order", "zero-based molecular atom order"),
            ("data_shape", "rows of [atom_i, atom_j]"),
        ),
    ),
    (
        "constrained_angle_atoms",
        (
            ("semantic_quantity", "constrained_internal_coordinate"),
            ("atom_order", "zero-based molecular atom order"),
            (
                "data_shape",
                "rows of [atom_i, atom_j, atom_k], vertex in the middle",
            ),
        ),
    ),
    (
        "constrained_dihedral_atoms",
        (
            ("semantic_quantity", "constrained_internal_coordinate"),
            ("atom_order", "zero-based molecular atom order"),
            (
                "data_shape",
                "rows of [atom_i, atom_j, atom_k, atom_l], about the j-k bond",
            ),
        ),
    ),
)


def _orca_held_coordinates(output: Any, kind: str) -> list[dict[str, Any]]:
    """The held coordinates of one kind, measured in the reached structure.

    ORCA declares the value it imposed in its constraint table; this
    measures the same coordinate in the structure ORCA handed back, which
    is the value the geometry a later stage consumes actually has.  The
    two must agree, and disagreement is reported with both numbers rather
    than resolved: a constrained optimisation whose constraint did not
    hold is not the calculation that was approved, and that is a finding,
    not a parse error.

    One kind per selector because a bond is a length and an angle is an
    angle.  Mixing them would force a dimensionless vector and lose the
    unit, and it would also make the atom rows ragged, which the typed
    layer cannot carry.
    """

    records = [
        record
        for record in _orca_constraint_records(output)
        if record["kind"] == kind
    ]
    if not records:
        raise MissingQuantityError(
            f"this ORCA result held no {kind}; it held "
            + ", ".join(
                sorted(
                    {
                        str(item["kind"])
                        for item in _orca_constraint_records(output)
                    }
                )
            )
        )
    molecule = output.molecule
    if isinstance(molecule, (list, tuple)):
        molecule = molecule[-1] if molecule else None
    if molecule is None:
        raise MissingQuantityError(
            "this ORCA result prints no reached structure to measure a "
            "held coordinate in"
        )
    measure = {
        "bond": molecule.get_distance,
        "angle": molecule.get_angle,
        "dihedral": molecule.get_dihedral,
    }[kind]
    held = []
    for record in records:
        # ``Molecule.get_distance`` and its siblings number atoms from
        # one, as everything the Agent *writes* does. What the extraction
        # plane *delivers* is zero-based, because it indexes the vectors
        # this same plane delivers -- symbols, positions, populations --
        # so the two bases are converted here and never mixed.
        reached = float(measure(*record["atoms"]))
        declared = float(record["value"])
        if abs(reached - declared) > _ORCA_CONSTRAINT_TOLERANCE[kind]:
            raise MissingQuantityError(
                f"ORCA declared {record['label']} held at {declared} but "
                f"the structure it returned has {reached}; the constraint "
                "this result was run under is not the one its geometry "
                "carries, and the geometry itself is readable as "
                "reached_positions"
            )
        held.append(
            {
                "atoms": tuple(int(index) - 1 for index in record["atoms"]),
                "label": record["label"],
                "value": reached,
            }
        )
    return held


def _orca_geometry_source_path(output: Any, selector: str) -> Path | None:
    """Name the ORCA sidecar behind the structure an IRC branch reached."""

    if selector != "trajectory_end_positions":
        return None
    try:
        record = _orca_irc_branch(output)
    except MissingQuantityError:
        return None
    return Path(str(record["geometry_file"]))


def _orca_native_evidence_paths(
    output: Any, selector: str
) -> tuple[Path, ...]:
    """Return the ORCA sidecar whose bytes a selector directly consumes."""

    if selector not in {
        "trajectory_connectivity_changed",
        "trajectory_end_connectivity",
        "trajectory_end_positions",
    }:
        return ()
    try:
        record = _orca_irc_branch(output)
    except MissingQuantityError:
        return ()
    return (Path(str(record["geometry_file"])),)


def _orca_printed_hessian(
    accessor: Callable[[Any], Any],
) -> Callable[[Any], Any]:
    """Refuse a frequency quantity ORCA did not print, saying why.

    A ScanTS with Freq prints the frequency table for its first Hessian
    only; the output says so (``unprinted_frequency_table_reason``) rather
    than handing over another Hessian's table (R10 Q31, CUHK 2154022).
    """

    def read(output: Any) -> Any:
        reason = getattr(output, "unprinted_frequency_table_reason", None)
        if reason:
            raise MissingQuantityError(reason)
        return accessor(output)

    return read


def _orca_accessors() -> dict[str, Callable[[Any], Any]]:
    accessors = _text_output_accessors(mode_composition=True)
    for selector in (
        "vibrational_frequencies",
        "vibrational_mode_atom_participation",
        "vibrational_mode_degeneracy_group",
    ):
        accessors[selector] = _orca_printed_hessian(accessors[selector])
    accessors.update(
        {
            # A relaxed scan is a surface, so it reaches the typed layer as
            # two parallel vectors rather than one opaque record list: the
            # existing operations then compose against it directly -- the
            # height of a torsional barrier is the spread of the energies --
            # instead of waiting on a bespoke profile type.
            "scan_coordinate_values": lambda output: [
                float(point["coordinate"]) for point in output.scan_profile
            ],
            "scan_energies": lambda output: [
                float(point["energy"]) for point in output.scan_profile
            ],
            "scan_point_indices": lambda output: [
                float(record["index"]) for record in output.scan_point_records
            ],
            # A constrained optimisation is defined by what it held. Each
            # kind answers under its own name so the value keeps its unit
            # and the atom rows keep one width, with the atoms beside the
            # values in the same order. The count is the whole of them and
            # is the number a claim is rendered from, which is how a
            # session states that the program applied the constraint it
            # was approved for instead of describing it.
            "constrained_coordinate_count": lambda output: float(
                len(_orca_constraint_records(output))
            ),
            "constrained_bond_atoms": lambda output: [
                [float(index) for index in held["atoms"]]
                for held in _orca_held_coordinates(output, "bond")
            ],
            "constrained_bond_lengths": lambda output: [
                held["value"]
                for held in _orca_held_coordinates(output, "bond")
            ],
            "constrained_angle_atoms": lambda output: [
                [float(index) for index in held["atoms"]]
                for held in _orca_held_coordinates(output, "angle")
            ],
            "constrained_bond_angles": lambda output: [
                held["value"]
                for held in _orca_held_coordinates(output, "angle")
            ],
            "constrained_dihedral_atoms": lambda output: [
                [float(index) for index in held["atoms"]]
                for held in _orca_held_coordinates(output, "dihedral")
            ],
            "constrained_dihedral_angles": lambda output: [
                held["value"]
                for held in _orca_held_coordinates(output, "dihedral")
            ],
            # Reached versus planned is the whole diagnosis when a scan
            # dies partway. The parser has always known both -- the step
            # announcements it counted, and the step total ORCA stated
            # before the first point -- and threw them away: a truncated
            # scan reached the typed layer only as shorter vectors, with
            # nothing to say what the target was. A live scan died at
            # step 2 of 12 and no tool could state either number.
            "scan_steps_reached": _scan_steps_reached,
            "scan_steps_planned": _scan_steps_planned,
            "absorption_wavelengths": lambda output: (
                _orca_absorption_per_state(output, "wavelength_nm")
            ),
            "excitation_energies": _orca_excitation_energies,
            "oscillator_strengths": lambda output: (
                _orca_absorption_per_state(output, "oscillator_strength")
            ),
            "excited_state_indices": lambda output: [
                int(item["state_index"])
                for item in _orca_served_records(output)
            ],
            "excited_state_manifold_roots": lambda output: [
                int(item["manifold_root"])
                for item in _orca_served_records(output)
            ],
            "excited_state_multiplicities": lambda output: [
                int(item)
                for item in _required_record_values(
                    _orca_served_records(output),
                    "multiplicity",
                    "a spin multiplicity (an unrestricted manifold has none)",
                )
            ],
            "excited_state_spin_square": lambda output: (
                _excited_spin_squares(
                    _orca_served_records(output),
                    (output.excited_state_applied or {}).get(
                        "response_method"
                    ),
                    output.state_manifold,
                )
            ),
            "excited_state_dominant_excitations": lambda output: (
                _dominant_excitation_values(
                    _orca_served_records(output), "labels"
                )
            ),
            "excited_state_dominant_weights": lambda output: (
                _dominant_excitation_values(
                    _orca_served_records(output), "weights"
                )
            ),
            "singlet_excitation_energies": lambda output: (
                _orca_manifold_values(output, 1, "energy_eV")
            ),
            "triplet_excitation_energies": lambda output: (
                _orca_manifold_values(output, 3, "energy_eV")
            ),
            "singlet_oscillator_strengths": lambda output: (
                _orca_manifold_values(output, 1, "oscillator_strength")
            ),
            "triplet_oscillator_strengths": lambda output: (
                _orca_manifold_values(output, 3, "oscillator_strength")
            ),
            "energy": _orca_total_energy,
            "energies": lambda output: [
                float(item) - _orca_fixed_geometry_root_shift(output)
                for item in output.energies
            ],
            "entropy_times_temperature": lambda output: float(
                output.entropy_times_temperature
            ),
            "positions": _orca_positions,
            "reached_positions": _orca_reached_positions,
            "connectivity": lambda output: _connectivity_matrix(
                output.thermochemistry_molecule
            ),
            "symbols": _orca_symbols,
            "charge": lambda output: int(output.thermochemistry_charge),
            "multiplicity": lambda output: int(
                output.thermochemistry_multiplicity
            ),
            "scf_energy": _orca_scf_energy,
            "reference_energy": _orca_scf_energy,
            "correlation_energy": _orca_correlation_energy,
            "t1_diagnostic": _orca_t1_diagnostic,
            "mayer_bond_orders": _orca_mayer_bond_orders,
            "mayer_free_valence": _orca_mayer_free_valence,
            "dispersion_energy": _orca_dispersion_energy,
            "auxiliary_basis": _orca_auxiliary_basis,
            "auxiliary_basis_role": _orca_auxiliary_basis_role,
            "solvation_model": _orca_solvation_model,
            "solvation_electrostatic_energy": (
                lambda output: output.solvation_electrostatic_energy
            ),
            "solvation_nonelectrostatic_energy": (
                lambda output: output.solvation_nonelectrostatic_energy
            ),
            "solvation_cavity_surface_area": (
                lambda output: output.solvation_cavity_surface_area
            ),
            "mulliken_atomic_charges": lambda output: _per_atom_vector(
                output.mulliken_atomic_charges,
                _orca_symbols(output),
                quantity="mulliken_atomic_charges",
            ),
            "loewdin_atomic_charges": lambda output: _per_atom_vector(
                output.loewdin_atomic_charges,
                _orca_symbols(output),
                quantity="loewdin_atomic_charges",
            ),
            # Spin populations sit in the second column of the same block
            # and were discarded for years; the sum is 2S by construction,
            # which is the checksum that a per-atom vector is complete.
            "mulliken_atomic_spin_populations": lambda output: (
                _orca_spin_populations(
                    output, quantity="mulliken_atomic_spin_populations"
                )
            ),
            "loewdin_atomic_spin_populations": lambda output: (
                _orca_spin_populations(
                    output, quantity="loewdin_atomic_spin_populations"
                )
            ),
            # A partition of the density into atomic basins rather than
            # over basis functions, so it does not carry Mulliken's basis
            # sensitivity -- which is why the condensed-Fukui literature
            # asks for it.  ORCA prints the block only when the route says
            # so, and a run that did not ask has no Hirshfeld analysis at
            # all rather than a failed one.
            "hirshfeld_atomic_charges": _orca_hirshfeld_charges,
            "hirshfeld_atomic_spin_populations": _orca_hirshfeld_spins,
            "functional": _orca_functional,
            "ab_initio": _route_ab_initio,
            "basis": _route_basis,
            "converged": _optimization_converged,
            "irc_converged": _irc_run_converged,
            "solvent": _orca_solvent,
            "dipole_moment": lambda output: [
                float(item)
                for item in output.dipole_moment_in_debye.reshape(-1)
            ],
            "dipole_moment_magnitude": lambda output: float(
                output.dipole_moment_magnitude_in_debye
            ),
            # Frontier orbital energies were already parsed into eV by
            # ORCAOutput and were reachable through no selector, so a session
            # that ran ORCA for a HOMO-LUMO gap executed the calculation and
            # could bind nothing from it.  Nothing new is parsed here.
            "homo": lambda output: _orca_frontier_pair(output)[0],
            "lumo": lambda output: _orca_frontier_pair(output)[1],
            "gap": _orca_frontier_gap,
            "alpha_homo": lambda output: _orca_channel_frontier(
                output, "alpha", True
            ),
            "alpha_lumo": lambda output: _orca_channel_frontier(
                output, "alpha", False
            ),
            "beta_homo": lambda output: _orca_channel_frontier(
                output, "beta", True
            ),
            "beta_lumo": lambda output: _orca_channel_frontier(
                output, "beta", False
            ),
            "vpt2_harmonic_frequencies": lambda output: [
                float(item) for item in output.vpt2_harmonic_frequencies
            ],
            "vpt2_fundamental_frequencies": lambda output: [
                float(item) for item in output.vpt2_fundamental_frequencies
            ],
            "vpt2_zero_point_rovibrational_energy": lambda output: float(
                output.vpt2_zero_point_rovibrational_energy
            ),
            "spin_square": _last_spin_square,
            "spin_square_target": _spin_square_target,
            "spin_square_deviation": lambda output: float(
                _last_spin_square(output) - _spin_square_target(output)
            ),
            "effective_multiplicity": lambda output: _effective_multiplicity(
                _last_spin_square(output)
            ),
            # An ORCA IRC writes its path to the log's PATH SUMMARY table
            # and its endpoint to a sidecar; the log body prints only the
            # saddle it was handed.  Reading the trajectory family out of
            # ``all_structures`` -- as the log-parsing readers do for
            # programs that print every frame -- would return that one
            # structure as both ends of the path, which is the loss this
            # jobtype's selector declarations were narrowed to prevent.
            "trajectory_frame_count": lambda output: len(
                _orca_irc_path_records(output)
            ),
            "trajectory_energies": lambda output: [
                float(record["energy"])
                for record in _orca_irc_path_records(output)
            ],
            "trajectory_start_positions": lambda output: [
                [float(value) for value in row]
                for row in _orca_irc_start_molecule(output).positions
            ],
            "trajectory_end_positions": lambda output: [
                [float(value) for value in row]
                for row in _orca_irc_end_molecule(output).positions
            ],
            "trajectory_start_connectivity": lambda output: (
                _connectivity_matrix(_orca_irc_start_molecule(output))
            ),
            "trajectory_end_connectivity": lambda output: _connectivity_matrix(
                _orca_irc_end_molecule(output)
            ),
            "trajectory_connectivity_changed": lambda output: int(
                _connectivity_matrix(_orca_irc_start_molecule(output))
                != _connectivity_matrix(_orca_irc_end_molecule(output))
            ),
            "irc_direction": _orca_irc_direction,
        }
    )
    return accessors


def _gaussian_output(path: Path) -> Any:
    from chemsmart.io.gaussian.output import Gaussian16Output

    return Gaussian16Output(filename=str(path))


def _gaussian_spin_square_after_annihilation(output: Any) -> float:
    """Return Gaussian's printed post-annihilation spin diagnostic.

    Restricted calculations and some incomplete unrestricted blocks record
    ``None`` for the post-annihilation value.  That is an absent scientific
    quantity, not a value that should reach ``float(None)`` and be reported as
    a parser-type failure.
    """

    value = output.spin_square_history[-1]["after_annihilation"]
    if value is None:
        raise MissingQuantityError(
            "Gaussian result does not print <S^2> after annihilation"
        )
    return float(value)


#: What one Gaussian IRC branch log answers.  ``energies`` is declared here
#: and nowhere else among the path job types because for an IRC it is the
#: path's own profile -- one printed energy per accepted point, index for
#: index with the frames -- rather than an optimiser's trace.  No frequency
#: table is printed by an IRC step, so nothing vibrational or
#: thermochemical is declared.
_GAUSSIAN_IRC_BRANCH_SELECTORS = (
    "ab_initio",
    "basis",
    "charge",
    "connectivity",
    "dipole_moment",
    "dipole_moment_magnitude",
    "effective_multiplicity",
    "energies",
    "energy",
    "functional",
    "gap",
    "homo",
    "irc_direction",
    "lumo",
    "multiplicity",
    "positions",
    "reached_positions",
    "scf_energy",
    "spin_square",
    "spin_square_after_annihilation",
    "spin_square_deviation",
    "spin_square_target",
    "symbols",
    "trajectory_connectivity_changed",
    "trajectory_end_connectivity",
    "trajectory_end_positions",
    "trajectory_frame_count",
    "trajectory_start_connectivity",
    "trajectory_start_positions",
    "wavefunction_stability_history",
    "wavefunction_stability_lowest_eigenvalue",
    "wavefunction_stability_rotation_space",
    "wavefunction_stability_verdict",
)


def _gaussian_reached_positions(output: Any) -> list[list[float]]:
    """The structure a Gaussian run's optimiser or path walk stopped on.

    ``positions`` answers with ``molecule``, which is the last parsed
    orientation whatever produced it -- for a single point that is the
    geometry the job was handed, and calling it "reached" would hand a
    recovery route its own seed.  This answers the *role*: a structure
    exists here only where a stage moved the geometry and the log records
    that it stopped on the frame it printed.

    Which job types those are is read back from this reader's own
    ``jobtype_selectors``, which is where the claim is made, rather than
    from a second list beside it: a probe that calls the accessor without
    the job-type gate (``available_selectors``) then answers exactly what
    the declaration promises, and the two cannot drift apart.

    Convergence is deliberately not asked.  "Reached" means the last
    structure the run printed, which is exactly what ORCA and PySCF mean
    by it, and the route that consumes it -- the geometry lift the repair
    menu offers -- exists *for* the run that stopped without converging:
    refusing there sends a session back to its own seed and throws away
    every step the optimiser took.  Measured on the archived
    error-terminated triplet optimisation in this repository, that is 72
    atoms carried 0.3777 A from where they started, over seven complete
    frames.  Whether the point is stationary is a different question, and
    ``converged``, the spectrum and the validity verdict answer it; the
    producer edge inside an approval keeps its own convergence test,
    because that one feeds a calculation a human approved.

    What is asked is that a frame exists.  ``molecule`` falls back to the
    input coordinate block when a log printed no orientation at all, so
    reading it would hand back the seed under this name -- the precise
    defect ORCA's accessor records -- and the parser's own frame list is
    read instead.  Every abnormally terminated structure-moving Gaussian
    log in this repository parses fewer complete frames than it printed
    orientations, because the truncation that drops a half-written block
    happens before this point.
    """

    jobtype = str(getattr(output, "jobtype", "") or "").strip().lower()
    declared = reader_for("gaussian").selectors_for_jobtype(jobtype) or ()
    if "reached_positions" not in declared:
        raise MissingQuantityError(
            f"a gaussian {jobtype or 'unknown'} result reaches no structure "
            "beyond the one it was handed; bind the supplied geometry or "
            "the producing optimisation's result instead"
        )
    frames = list(getattr(output, "all_structures", ()) or ())
    if not frames:
        raise MissingQuantityError(
            "this gaussian result printed no complete structure, so the "
            "only geometry it carries is the one it was supplied; bind "
            "that instead"
        )
    return [[float(value) for value in row] for row in frames[-1].positions]


def _gaussian_held_rows(output: Any) -> list[dict[str, Any]]:
    """The internal coordinates a Gaussian run froze, from its own record.

    The rows its echoed ModRedundant section froze, and every other
    internal coordinate its parameter table marks frozen: a ``ts`` search
    that read its geometry from a constrained optimisation's checkpoint
    (``geom=check``) keeps that job's frozen bonds without echoing a
    ModRedundant section, and was read as holding nothing (the archived
    ``Pd_insertion_ts_r.log``: two Pd-C bonds frozen, a Cartesian gradient
    of 0.012 Eh/Bohr, called stationary; R10 Q33).  A row the table adds
    carries no value: it was frozen where the run's first structure has it.
    """

    held = list(getattr(output, "held_internal_coordinates", None) or ())
    parameters = _gaussian_frozen_parameters(output)
    if parameters is None:
        return held

    def canonical(atoms: Sequence[int]) -> tuple[int, ...]:
        atoms = tuple(int(index) for index in atoms)
        return min(atoms, tuple(reversed(atoms)))

    named = {canonical(row["atoms"]) for row in held}
    for row in parameters[0]:
        if canonical(row["atoms"]) not in named:
            named.add(canonical(row["atoms"]))
            held.append({**row, "value": None})
    return held


def _gaussian_held_count(output: Any) -> float:
    held = _gaussian_held_rows(output)
    if not held:
        raise MissingQuantityError(
            "this gaussian result held no internal coordinate; this family "
            "answers the bonds, angles and dihedrals a run froze, and a "
            "Cartesian atom freeze is not one of them"
        )
    return float(len(held))


def _gaussian_held_coordinates(output: Any, kind: str) -> list[dict[str, Any]]:
    """The coordinates of one kind a Gaussian run held, where it ended.

    What ORCA's reader answers for a constrained optimisation, from
    Gaussian's own record (``_gaussian_held_rows``), each measured in the
    structure the run returned.  Gaussian states no
    held value unless the row carried one -- it freezes a coordinate where
    the input geometry has it -- so the value it was held at is the first
    printed structure's, or the row's own.  The two must agree within the
    tolerance ORCA's reader uses; disagreement is reported with both
    numbers, because a constraint that did not hold is a finding, not a
    parse error.  Atoms are delivered zero-based like every other atom
    index this plane serves.
    """

    all_held = _gaussian_held_rows(output)
    rows = [row for row in all_held if row["kind"] == kind]
    if not rows:
        kinds = sorted({str(row["kind"]) for row in all_held})
        raise MissingQuantityError(
            f"this gaussian result held no {kind}"
            + (f"; it held {', '.join(kinds)}" if kinds else "")
        )
    frames = list(getattr(output, "all_structures", ()) or ())
    if not frames:
        raise MissingQuantityError(
            "this gaussian result prints no reached structure to measure a "
            "held coordinate in"
        )
    method = {
        "bond": "get_distance",
        "angle": "get_angle",
        "dihedral": "get_dihedral",
    }[kind]
    held = []
    for row in rows:
        reached = float(getattr(frames[-1], method)(*row["atoms"]))
        declared = (
            float(row["value"])
            if row["value"] is not None
            else float(getattr(frames[0], method)(*row["atoms"]))
        )
        offset = reached - declared
        if kind == "dihedral":
            offset = ((offset + 180.0) % 360.0) - 180.0
        if abs(offset) > _ORCA_CONSTRAINT_TOLERANCE[kind]:
            raise MissingQuantityError(
                f"gaussian held {row['label']} at {declared} but the "
                f"structure it returned has {reached}; the constraint this "
                "result was run under is not the one its geometry carries, "
                "and the geometry itself is readable as reached_positions"
            )
        held.append(
            {
                "atoms": tuple(int(index) - 1 for index in row["atoms"]),
                "label": row["label"],
                "value": reached,
            }
        )
    return held


def _gaussian_population(
    attribute: str, *, quantity: str
) -> Callable[[Any], list[float]]:
    """Read one Gaussian population block into molecular atom order.

    Gaussian labels its population rows by global atom index, which is the
    scheme ``_per_atom_vector`` resolves against the molecule's own
    symbols; a row set that names another molecule's atoms is refused
    rather than reordered.  A run that never asked for the analysis has no
    block, and the parser answers that with ``None`` for Mulliken and by
    indexing an empty list for Hirshfeld -- both are an absent quantity
    here, not a parser failure.
    """

    def _read(output: Any) -> list[float]:
        try:
            labelled = getattr(output, attribute)
        except IndexError as error:
            raise MissingQuantityError(
                f"this gaussian result prints no {quantity} block; the run "
                "did not request that population analysis"
            ) from error
        if labelled is None:
            raise MissingQuantityError(
                f"this gaussian result prints no {quantity} block; the run "
                "did not request that population analysis"
            )
        return _per_atom_vector(labelled, _symbols(output), quantity=quantity)

    return _read


def _gaussian_ir_intensities(output: Any) -> list[float]:
    """Per-normal-mode IR absorption intensities, paired with frequencies.

    Both lists come from the same frequency job step, so a disagreement in
    length can only mean the parser read two different tables; refuse
    rather than deliver an intensity under another mode's index.
    """

    intensities = [float(item) for item in output.ir_intensities or ()]
    if not intensities:
        raise MissingQuantityError(
            "this gaussian result records no IR intensities"
        )
    frequencies = [
        float(item) for item in output.vibrational_frequencies or ()
    ]
    if len(intensities) != len(frequencies):
        raise MissingQuantityError(
            f"this gaussian result reports {len(frequencies)} vibrational "
            f"frequencies and {len(intensities)} IR intensities, so an "
            "intensity index would not name the mode its frequency names"
        )
    return intensities


def _gaussian_scan_profile(output: Any) -> list[Mapping[str, Any]]:
    """The relaxed-scan surface, or an honest absence."""

    profile = getattr(output, "scan_profile", None)
    if not profile:
        raise MissingQuantityError(
            "this gaussian result establishes no relaxed-scan surface; a "
            "constrained optimisation drives no coordinate, and a section "
            "driving several has no single coordinate value per point"
        )
    return list(profile)


def _gaussian_energies(output: Any) -> list[float]:
    """Totals on the surface the route asks for, one per geometry.

    The parser names which printed line that is (``energy_source``); a
    post-HF method whose total it does not read answers nothing rather
    than the lower level Gaussian printed on the way to it.
    """

    values = [float(item) for item in output.energies]
    if values:
        return values
    if getattr(output, "energy_source", None) == "unrecognized_post_hf":
        raise MissingQuantityError(
            "this Gaussian route's method prints a total this reader does "
            "not read; the lower levels printed on the way to it are not "
            "its energy and are not served as one"
        )
    raise MissingQuantityError("this Gaussian result printed no energy")


def _gaussian_last_stability_record(output: Any) -> Mapping[str, Any]:
    """The stability analysis the run ended on, or why there is none."""

    records = list(
        getattr(output, "wavefunction_stability_records", None) or ()
    )
    if not records:
        raise MissingQuantityError(
            "this Gaussian result printed no stability analysis (the route "
            "asked for no Stable)"
        )
    return records[-1]


def _gaussian_stability_lowest_eigenvalue(output: Any) -> float:
    """The lowest stability-matrix eigenvalue of the last analysis, in Eh.

    Printed as ``Eigenvector 1: <label> Eigenvalue= X`` beside every
    verdict and served by no reader: the number the verdict is drawn from.
    For a restricted reference it is the lowest root over the singlet
    (internal) and triplet (RHF -> UHF) blocks Gaussian tests; on singlet
    O2 at RB3LYP/def2-SVP it is -0.0926178 Eh, the triplet root, equal to
    PySCF's RHF/RKS -> UHF/UKS eigenvalue at the same level to 1e-6 Eh
    (CUHK 2152098, 2151881) -- while PySCF's internal eigenvalue is four
    times Gaussian's singlet root, a normalisation of PySCF's own.
    """

    record = _gaussian_last_stability_record(output)
    eigenvalues = record.get("eigenvalues") or ()
    if not eigenvalues:
        raise MissingQuantityError(
            "the last stability verdict this log prints has no eigenvalue "
            "printed before it"
        )
    return float(eigenvalues[0]["eigenvalue"])


def _gaussian_stability_rotation_space(output: Any) -> str:
    """The space Gaussian's last verdict names, where it names one."""

    record = _gaussian_last_stability_record(output)
    space = record.get("rotation_space")
    if not space:
        raise MissingQuantityError(
            "Gaussian's last verdict names no rotation space: 'stable under "
            "the perturbations considered' does not say which were"
        )
    return str(space)


def _gaussian_molecular_volume(output: Any) -> float:
    """The molecule's volume as Gaussian's ``volume`` keyword measured it.

    The volume inside the 0.001 e/bohr^3 density contour, a Monte Carlo
    estimate per molecule, in bohr^3 (Gaussian prints it beside the molar
    figure in cm^3/mol, which is N_A times it).  Asked for in R10 Q6
    (ar04/ar10) and verified unreachable because no reader served it.
    Printed only when the route asks for ``volume``; the SMD solvent
    parameter printed under the same words is not read.
    """

    values = list(getattr(output, "molecular_volumes_bohr3", None) or ())
    if not values:
        raise MissingQuantityError(
            "this Gaussian result printed no molecular volume: the route did "
            "not ask for volume"
        )
    return float(values[-1])


def _gaussian_solvation_model(output: Any) -> str:
    """The continuum the route applied, in the route's own word, or
    ``gas_phase`` when the route asks for none.

    Read from the route Gaussian echoes, as the level record already does,
    so a solvation term is read beside the model that gives it its
    meaning.  A route that asks for a continuum this reader cannot name is
    an absence, never gas phase.
    """

    model = getattr(output, "solvent_model", None)
    if model:
        return str(model).strip().lower()
    route = str(getattr(output, "route_string", "") or "").lower()
    if "scrf" in route:
        raise MissingQuantityError(
            "this Gaussian route asks for a continuum (scrf) this reader "
            "does not name"
        )
    return "gas_phase"


def _gaussian_solvent(output: Any) -> str:
    if _gaussian_solvation_model(output) == "gas_phase":
        raise MissingQuantityError("this Gaussian run is gas phase")
    solvent = getattr(output, "solvent_id", None)
    if not solvent:
        raise MissingQuantityError(
            "this Gaussian solvated route names no solvent"
        )
    return str(solvent).strip().lower()


def _gaussian_smd_cds_energy(output: Any) -> float:
    """SMD's non-electrostatic term, as Gaussian printed it, in kcal/mol.

    The same quantity ORCA and PySCF already serve under this name: the
    cavity-dispersion-solvent-structure part of the solvation free energy
    the model put into the total.  The solvation charter topic recorded
    "No archived Gaussian log carries the printed terms"; fourteen archived
    Gaussian SMD logs print this one ("SMD-CDS (non-electrostatic) energy
    (kcal/mol) = ...  (included in total energy above)").  Gaussian prints
    no separate electrostatic term, so that one stays undeclared.  The last
    print belongs to the final structure; a run with no SMD continuum
    printed none, and the absence says so.
    """

    values = list(getattr(output, "smd_cds_energies_kcal_per_mol", None) or ())
    if not values:
        raise MissingQuantityError(
            "this Gaussian result printed no SMD-CDS term: it ran in the gas "
            "phase or in a continuum model without one"
        )
    return float(values[-1])


def _gaussian_electronic_spatial_extent(output: Any) -> float:
    """<R**2> of the SCF density at the structure the run ended on, bohr^2.

    Printed by every Gaussian population analysis and served by no reader
    until now, although it is one of the few numbers that say how diffuse
    a density is: an anion or a Rydberg-like state that the basis cannot
    hold shows up here before it shows up anywhere else, and Q10's LG1
    session (CUHK 2151662) could only be told that the log printed it.
    The last print belongs to the final density.  The value depends on the
    origin, so it is served only where Gaussian computed in its standard
    orientation (origin at the centre of nuclear charge), and only for the
    SCF density the population analysis names.
    """

    records = list(getattr(output, "electronic_spatial_extents", None) or ())
    if not records:
        raise MissingQuantityError(
            "this Gaussian result printed no <R**2>: no population analysis "
            "ran"
        )
    last = records[-1]
    density = str(last.get("density") or "")
    if not density.casefold().startswith("scf"):
        raise MissingQuantityError(
            "the last <R**2> this log prints belongs to the "
            f"{density or 'unnamed'} density, not the SCF reference this "
            "selector serves"
        )
    if not getattr(output, "standard_orientations", None):
        raise MissingQuantityError(
            "this Gaussian run computed in the input orientation "
            "(reorientation suppressed), so its <R**2> is about the "
            "coordinates' own origin; the spatial extent depends on the "
            "origin and is served only about the centre of nuclear charge"
        )
    return float(last["value"])


def _gaussian_scf_energy(output: Any) -> float:
    """The SCF reference total at the last geometry (``SCF Done``).

    For an SCF route it is ``energy``; under a correlated method, a double
    hybrid or a TD optimisation it is the reference beneath the surface
    ``energy`` names, as ``scf_energy`` is on ORCA and PySCF.
    """

    values = list(getattr(output, "scf_energies", None) or ())
    if not values:
        raise MissingQuantityError("this Gaussian result printed no SCF Done")
    return float(values[-1])


#: Gaussian's word for the electronic-provenance axis over the selectors
#: its reader implements: ``energy`` is the total of the surface the route
#: computed on -- the correlated method, a double hybrid's total, the
#: followed root of a TD optimisation -- resolved per artifact, and the
#: density-derived properties beside it are the SCF reference's.
_GAUSSIAN_ELECTRONIC_PROVENANCE_DECLARED = (
    ("absorption_wavelengths", "excited_root"),
    ("dipole_moment", "reference"),
    ("dipole_moment_magnitude", "reference"),
    ("effective_multiplicity", "reference"),
    ("electronic_spatial_extent", "reference"),
    ("molecular_volume", "reference"),
    ("energies", "computed_surface"),
    ("energy", "computed_surface"),
    ("excitation_energies", "excited_root"),
    ("excited_state_dominant_excitations", "excited_root"),
    ("excited_state_dominant_weights", "excited_root"),
    ("excited_state_indices", "excited_root"),
    ("excited_state_labels", "excited_root"),
    ("excited_state_manifold_roots", "excited_root"),
    ("excited_state_multiplicities", "excited_root"),
    ("excited_state_spin_square", "excited_root"),
    ("gap", "reference"),
    ("hirshfeld_atomic_charges", "reference"),
    ("hirshfeld_atomic_spin_populations", "reference"),
    ("homo", "reference"),
    ("lumo", "reference"),
    ("mulliken_atomic_charges", "reference"),
    ("mulliken_atomic_spin_populations", "reference"),
    ("oscillator_strengths", "excited_root"),
    ("scf_energy", "reference"),
    ("singlet_excitation_energies", "excited_root"),
    ("singlet_oscillator_strengths", "excited_root"),
    # SMD's non-electrostatic term is part of the SCF total it printed
    # beneath, the reference's own energy.
    ("solvation_nonelectrostatic_energy", "reference"),
    ("spin_square", "reference"),
    ("spin_square_after_annihilation", "reference"),
    ("spin_square_deviation", "reference"),
    ("spin_square_target", "reference"),
    ("triplet_excitation_energies", "excited_root"),
    ("triplet_oscillator_strengths", "excited_root"),
)


def _gaussian_accessors() -> dict[str, Callable[[Any], Any]]:
    accessors = _text_output_accessors()
    accessors.update(
        {
            "energy": lambda output: _gaussian_energies(output)[-1],
            "energies": _gaussian_energies,
            "scf_energy": _gaussian_scf_energy,
            "reached_positions": _gaussian_reached_positions,
            # Whether the optimiser said it converged, and what a
            # constrained optimisation held -- the answers ORCA's modred
            # gives, from Gaussian's own record.
            "converged": _optimization_converged,
            "constrained_coordinate_count": _gaussian_held_count,
            "constrained_bond_atoms": lambda output: [
                [float(index) for index in held["atoms"]]
                for held in _gaussian_held_coordinates(output, "bond")
            ],
            "constrained_bond_lengths": lambda output: [
                held["value"]
                for held in _gaussian_held_coordinates(output, "bond")
            ],
            "constrained_angle_atoms": lambda output: [
                [float(index) for index in held["atoms"]]
                for held in _gaussian_held_coordinates(output, "angle")
            ],
            "constrained_bond_angles": lambda output: [
                held["value"]
                for held in _gaussian_held_coordinates(output, "angle")
            ],
            "constrained_dihedral_atoms": lambda output: [
                [float(index) for index in held["atoms"]]
                for held in _gaussian_held_coordinates(output, "dihedral")
            ],
            "constrained_dihedral_angles": lambda output: [
                held["value"]
                for held in _gaussian_held_coordinates(output, "dihedral")
            ],
            # The surface reaches the typed layer as two parallel vectors,
            # exactly as ORCA's does, so the existing operations compose
            # against it: the height of a torsional barrier is the spread
            # of the energies, on either program.
            "scan_coordinate_values": lambda output: [
                float(point["coordinate"])
                for point in _gaussian_scan_profile(output)
            ],
            "scan_energies": lambda output: [
                float(point["energy"])
                for point in _gaussian_scan_profile(output)
            ],
            "scan_point_indices": lambda output: [
                float(point["index"])
                for point in _gaussian_scan_profile(output)
            ],
            "scan_steps_reached": _scan_steps_reached,
            "scan_steps_planned": _scan_steps_planned,
            "ir_intensities": _gaussian_ir_intensities,
            "mulliken_atomic_charges": _gaussian_population(
                "mulliken_atomic_charges",
                quantity="mulliken_atomic_charges",
            ),
            "mulliken_atomic_spin_populations": _gaussian_population(
                "mulliken_spin_densities",
                quantity="mulliken_atomic_spin_populations",
            ),
            "hirshfeld_atomic_charges": _gaussian_population(
                "hirshfeld_charges",
                quantity="hirshfeld_atomic_charges",
            ),
            # The S-H column of the same block, parsed and never served.
            "hirshfeld_atomic_spin_populations": _gaussian_population(
                "hirshfeld_spin_densities",
                quantity="hirshfeld_atomic_spin_populations",
            ),
            "absorption_wavelengths": lambda output: [
                float(item) for item in output.absorptions_in_nm
            ],
            "excitation_energies": lambda output: [
                float(item["energy_eV"])
                for item in output.excited_state_records
            ],
            "oscillator_strengths": lambda output: [
                float(item) for item in output.oscillatory_strengths
            ],
            "excited_state_indices": lambda output: [
                int(item["state_index"])
                for item in output.excited_state_records
            ],
            "excited_state_manifold_roots": lambda output: [
                int(item)
                for item in _required_record_values(
                    output.excited_state_records,
                    "manifold_root",
                    "a manifold-local root index",
                )
            ],
            "excited_state_multiplicities": lambda output: [
                int(item)
                for item in _required_record_values(
                    output.excited_state_records,
                    "multiplicity",
                    "a source-labelled spin multiplicity",
                )
            ],
            "excited_state_labels": lambda output: [
                str(item["state_label"])
                for item in output.excited_state_records
            ],
            "excited_state_spin_square": lambda output: (
                _excited_spin_squares(
                    output.excited_state_records,
                    (output.excited_state_request or {}).get(
                        "response_method"
                    ),
                    (output.excited_state_request or {}).get("state_manifold"),
                )
            ),
            "excited_state_dominant_excitations": lambda output: (
                _dominant_excitation_values(
                    output.excited_state_records, "labels"
                )
            ),
            "excited_state_dominant_weights": lambda output: (
                _dominant_excitation_values(
                    output.excited_state_records, "weights"
                )
            ),
            "singlet_excitation_energies": lambda output: [
                float(item["energy_eV"])
                for item in output.excited_state_records
                if item["multiplicity"] == 1
            ],
            "triplet_excitation_energies": lambda output: [
                float(item["energy_eV"])
                for item in output.excited_state_records
                if item["multiplicity"] == 3
            ],
            "singlet_oscillator_strengths": lambda output: [
                float(item["oscillator_strength"])
                for item in output.excited_state_records
                if item["multiplicity"] == 1
                and item["oscillator_strength"] is not None
            ],
            "triplet_oscillator_strengths": lambda output: [
                float(item["oscillator_strength"])
                for item in output.excited_state_records
                if item["multiplicity"] == 3
                and item["oscillator_strength"] is not None
            ],
            "dipole_moment": lambda output: [
                float(item) for item in output.all_dipole_moments[-1]
            ],
            "dipole_moment_magnitude": lambda output: float(
                output.all_dipole_moment_magnitudes[-1]
            ),
            "electronic_spatial_extent": _gaussian_electronic_spatial_extent,
            "molecular_volume": _gaussian_molecular_volume,
            "solvation_model": _gaussian_solvation_model,
            "solvent": _gaussian_solvent,
            "solvation_nonelectrostatic_energy": _gaussian_smd_cds_energy,
            "spin_square": lambda output: _last_spin_square(
                output, "before_annihilation"
            ),
            "spin_square_after_annihilation": (
                _gaussian_spin_square_after_annihilation
            ),
            "spin_square_target": _spin_square_target,
            "spin_square_deviation": lambda output: float(
                _last_spin_square(output, "before_annihilation")
                - _spin_square_target(output)
            ),
            "effective_multiplicity": lambda output: _effective_multiplicity(
                _last_spin_square(output, "before_annihilation")
            ),
            "wavefunction_stability_verdict": lambda output: str(
                output.wavefunction_stability_history[-1]
            ),
            "wavefunction_stability_history": lambda output: [
                str(item) for item in output.wavefunction_stability_history
            ],
            "wavefunction_stability_lowest_eigenvalue": (
                _gaussian_stability_lowest_eigenvalue
            ),
            "wavefunction_stability_rotation_space": (
                _gaussian_stability_rotation_space
            ),
            "trajectory_frame_count": lambda output: len(
                _irc_structures(output)
            ),
            "trajectory_start_positions": lambda output: [
                [float(value) for value in row]
                for row in _irc_structures(output)[0].positions
            ],
            "trajectory_end_positions": lambda output: [
                [float(value) for value in row]
                for row in _irc_structures(output)[-1].positions
            ],
            "trajectory_start_connectivity": lambda output: (
                _connectivity_matrix(_irc_structures(output)[0])
            ),
            "trajectory_end_connectivity": lambda output: _connectivity_matrix(
                _irc_structures(output)[-1]
            ),
            "trajectory_connectivity_changed": (
                _trajectory_connectivity_changed
            ),
            "irc_direction": _irc_direction,
            "functional": _gaussian_functional,
            "homo": _gaussian_frontier("homo_energy"),
            "lumo": _gaussian_frontier("lumo_energy"),
            "gap": _gaussian_frontier("fmo_gap"),
            "ab_initio": _route_ab_initio,
            "basis": _route_basis,
        }
    )
    return accessors


def _xtb_output(path: Path | str) -> Any:
    from chemsmart.io.xtb.output import XTBOutput

    # Unlike the Gaussian and ORCA readers, XTBOutput represents the complete
    # calculation directory because xTB can distribute one result across the
    # main log, geometry, Hessian, and vibrational-spectrum files.  The host
    # has already resolved and verified the exact main output artifact; its
    # parent is therefore the corresponding calculation directory.
    #
    # The host's organs do not agree on the type they hand a reader: the mode
    # displacement passes a Path and the stationary-point characterisation
    # passes ``str(artifact.path)``.  Every other reader takes either, because
    # each only wraps the value in its own parser; this one navigates it, so a
    # str reached ``.parent`` and a live session asking what its saddle was
    # got ``AttributeError: 'str' object has no attribute 'parent'`` instead
    # of an answer (goal r8x-xtb-saddle-escape, CUHK 2142404).  Normalised
    # here, where the navigation happens.
    return XTBOutput(folder=str(Path(path).parent))


def _xtb_state_integer(attribute: str) -> Callable[[Any], int]:
    def _read(output: Any) -> int:
        value = getattr(output, attribute, None)
        if value is None:
            raise MissingQuantityError(
                f"this xtb result records no {attribute}"
            )
        return int(value)

    return _read


#: What a program's own printed free energy is, in the host's terms, for the
#: readers that serve one under ``gibbs_free_energy``.  Each sentence was
#: measured, not recalled: oracle O1 (R10 Q5, CUHK Slurm 2149853/2149909)
#: re-derived every printed value from the same program's frequencies.
#: Gaussian 16's equals the host's RRHO derivation with Gaussian's own
#: symmetry number to <= 5e-7 Eh over ten molecules -- and that number was
#: 1 for a C3v NH3. xTB 6.7.1's equals the host's Grimme derivation with a
#: 50 cm-1 cutoff and alpha 4 to within 0.005 kcal/mol for nine closed-shell
#: molecules, and sits RT ln 2 = 0.41 kcal/mol above it for doublet NO2: it
#: carries no electronic spin-degeneracy entropy.  ORCA's printed value
#: (Grimme, 100 cm-1) is not served at all; its reader says why.
PRINTED_THERMOCHEMISTRY_CONVENTIONS: Mapping[str, str] = {
    "gaussian": (
        "Gaussian's printed free energy: harmonic (RRHO), ideal gas at the "
        "temperature and pressure the route set (298.15 K and 1 atm unless "
        "it said otherwise), Gaussian's own rotational symmetry number, "
        "imaginary modes left out"
    ),
    "xtb": (
        "xTB's printed free energy: modified RRHO (Grimme's free-rotor "
        "interpolation below xTB's rotor cutoff, 50 cm-1 by default), "
        "298.15 K unless set, xTB's own rotational symmetry number (its "
        "symmetry search tolerates 0.1 A), and no electronic "
        "spin-degeneracy entropy -- an open shell sits RT ln(multiplicity) "
        "above the host's derivation"
    ),
}


def _xtb_gibbs(output: Any) -> float:
    """xTB prints G only when a Hessian ran; None is an absent quantity."""

    value = getattr(output, "gibbs_free_energy", None)
    if value is None:
        raise MissingQuantityError(
            "this xtb result records no thermochemistry; free energy "
            "requires a Hessian calculation"
        )
    return float(value)


def _xtb_reached_positions(output: Any) -> list[list[float]]:
    """Return only a converged optimisation's own final ``xtbopt`` frame.

    ``XTBOutput.molecule`` is deliberately broad: a Hessian may reconstruct
    a frame from a Gaussian-format sidecar and a single point from its input.
    Those are inspectable structures, not proof that the stage reached a new
    geometry.  The recovery consumer asks for this stricter role.
    """

    if str(getattr(output, "jobtype", "") or "").casefold() != "opt":
        raise MissingQuantityError(
            "reached positions require an xTB optimisation result"
        )
    molecule = getattr(output, "optimized_structure", None)
    if molecule is None:
        raise MissingQuantityError(
            "this xTB optimisation did not normally converge with a readable "
            "xtbopt final geometry"
        )
    return [[float(value) for value in row] for row in molecule.positions]


def _xtb_scc_atomic_charges(output: Any) -> list[float]:
    """Read xTB's positional SCC population without relabelling its scheme."""

    charges_file = getattr(output, "charges_file", None)
    values = getattr(charges_file, "partial_charges", None)
    symbols = _symbols(output)
    if values is None:
        raise MissingQuantityError(
            "this xTB result wrote no SCC atomic-population charges sidecar"
        )
    charges = [float(value) for value in values]
    if len(charges) != len(symbols):
        raise MissingQuantityError(
            "xTB SCC atomic-population sidecar has "
            f"{len(charges)} values for {len(symbols)} atoms"
        )
    total = getattr(output, "charge", None)
    if total is not None and abs(sum(charges) - float(total)) > 1e-3:
        raise MissingQuantityError(
            "xTB SCC atomic-population charges sum to "
            f"{sum(charges):.6f} e while the result's total charge is "
            f"{float(total):.6f} e"
        )
    return charges


def _xtb_wiberg_bond_orders(output: Any) -> list[list[Any]]:
    """Return sparse Wiberg records in canonical zero-based atom order."""

    wbo_file = getattr(output, "wbo_file", None)
    if wbo_file is None:
        raise MissingQuantityError(
            "this xTB result wrote no Wiberg bond order (wbo) sidecar"
        )
    symbols = _symbols(output)
    n_atoms = len(symbols)
    try:
        pairs = getattr(wbo_file, "bond_orders", None)
    except ValueError as exc:
        from chemsmart.analysis import result_quantities as rq

        raise rq.QuantityExtractionError(
            f"xTB WBO sidecar is malformed: {exc}"
        ) from exc
    if not pairs:
        raise MissingQuantityError(
            "this xTB result wrote an empty WBO sidecar (no native WBO "
            "pair records)"
        )
    rows: list[list[Any]] = []
    seen: set[tuple[int, int]] = set()
    for atom_i, atom_j, order in pairs:
        if atom_i < 1 or atom_i > n_atoms or atom_j < 1 or atom_j > n_atoms:
            from chemsmart.analysis import result_quantities as rq

            raise rq.QuantityExtractionError(
                f"xTB WBO atom index out of bounds: ({atom_i}, {atom_j}) "
                f"for molecule with {n_atoms} atoms"
            )
        i0, j0 = atom_i - 1, atom_j - 1
        if i0 == j0:
            from chemsmart.analysis import result_quantities as rq

            raise rq.QuantityExtractionError(
                f"xTB WBO sidecar names a self-pair: ({atom_i}, {atom_j})"
            )
        if i0 > j0:
            i0, j0 = j0, i0
        pair = (i0, j0)
        if pair in seen:
            from chemsmart.analysis import result_quantities as rq

            raise rq.QuantityExtractionError(
                "xTB WBO sidecar repeats unordered atom pair " f"({i0}, {j0})"
            )
        numeric_order = float(order)
        if not math.isfinite(numeric_order) or numeric_order < 0.0:
            from chemsmart.analysis import result_quantities as rq

            raise rq.QuantityExtractionError(
                "xTB WBO sidecar contains a non-finite or negative bond "
                f"order for pair ({i0}, {j0})"
            )
        seen.add(pair)
        rows.append([i0, j0, numeric_order])
    rows.sort(key=lambda r: (r[0], r[1]))
    return rows


def _xtb_dispersion_energy(output: Any) -> float:
    """Return the parsed dispersion energy contribution from the summary block."""

    main_out = getattr(output, "main_out", None)
    value = getattr(main_out, "dispersion_energy", None) if main_out else None
    if value is None:
        raise MissingQuantityError(
            "this xTB result records no dispersion energy in its summary block"
        )
    return float(value)


#: One solvation free-energy term, as this vocabulary names it and as the
#: xTB parser names it.  ALPB and GBSA print the total and its four parts
#: in one SUMMARY block: Gsolv = Gelec + Gsasa + Ghb + Gshift, exactly, to
#: the twelve decimals xTB prints.  Gelec is the *whole* electrostatic part
#: of Gsolv, which is what the shared name means on every program, so it
#: answers to it.  The other three are not: the non-electrostatic part of
#: an ALPB solvation free energy is Gsasa **and** Ghb **and** Gshift, so
#: calling Gsasa ``solvation_nonelectrostatic_energy`` would be a name that
#: is false about the model.  They keep the scheme in the name, as the xTB
#: population does.
_XTB_SOLVATION_TERMS: tuple[tuple[str, str], ...] = (
    ("solvation_free_energy", "solvation_energy_gsolv"),
    ("solvation_electrostatic_energy", "electronic_solvation_energy_gelec"),
    ("xtb_solvation_sasa_energy", "surface_area_solvation_energy_gsasa"),
    (
        "xtb_solvation_hydrogen_bond_energy",
        "hydrogen_bonding_solvation_energy_ghb",
    ),
    ("xtb_solvation_shift_energy", "empirical_shift_correction_gshift"),
)

#: The parts whose sum is the total, in the order xTB prints them.
_XTB_SOLVATION_COMPONENTS: tuple[str, ...] = tuple(
    selector for selector, _attribute in _XTB_SOLVATION_TERMS[1:]
)

#: Print precision is twelve decimals, so four rounded parts and a rounded
#: total can disagree by a few times 1e-13.  Anything above this is not
#: rounding: it is terms read from different SUMMARY blocks, which is the
#: one failure a single-term accessor cannot see.
_XTB_SOLVATION_CLOSURE_TOLERANCE = 1e-9


def _xtb_solvation_context(output: Any) -> tuple[str, str | None]:
    """The solvent treatment this result's own setup block recorded.

    A gas-phase result answers ``gas_phase`` rather than nothing, the same
    word ORCA's reader uses, because "this ran without solvent" is an
    answer and not a missing value.
    """

    if not bool(getattr(output, "solvent_on", False)):
        return "gas_phase", None
    model = getattr(output, "solvent_model", None)
    if model in (None, ""):
        raise MissingQuantityError(
            "this xTB result records solvation but names no solvation model"
        )
    solvent = getattr(output, "solvent_id", None)
    if solvent in (None, ""):
        return str(model).strip().lower(), None
    return str(model).strip().lower(), str(solvent).strip().lower()


def _xtb_solvation_terms(output: Any) -> dict[str, float]:
    """Every solvation free-energy term this result printed, checked as one.

    The parser resolves each term by scanning the last SUMMARY block for
    its own keyword, so nothing in a single read establishes that five
    numbers came from one block.  The model's own identity does:
    ``Gsolv`` is the sum of its four parts.  A result whose parts do not
    add to its total is evidence the block was misread, not a fact about
    the solvent, so it is refused rather than served.
    """

    from chemsmart.analysis import result_quantities as rq

    model, _solvent = _xtb_solvation_context(output)
    values: dict[str, float] = {}
    for selector, attribute in _XTB_SOLVATION_TERMS:
        raw = getattr(output, attribute, None)
        if raw is None:
            continue
        number = float(raw)
        if not math.isfinite(number):
            raise rq.QuantityExtractionError(
                f"xTB energy summary prints a non-finite {selector}: {number}"
            )
        values[selector] = number
    if model == "gas_phase":
        if values:
            raise rq.QuantityExtractionError(
                "this xTB result's setup block reports no solvation while "
                "its energy summary prints "
                f"{', '.join(sorted(values))}; the two disagree about "
                "whether a solvent was applied"
            )
        return values
    total = values.get("solvation_free_energy")
    parts = tuple(
        values.get(selector) for selector in _XTB_SOLVATION_COMPONENTS
    )
    if total is not None and all(part is not None for part in parts):
        residual = total - math.fsum(
            part for part in parts if part is not None
        )
        if abs(residual) > _XTB_SOLVATION_CLOSURE_TOLERANCE:
            raise rq.QuantityExtractionError(
                "xTB solvation terms do not close: Gsolv "
                f"{total:.12f} Eh differs from the sum of its parts by "
                f"{residual:.3e} Eh"
            )
    return values


def _xtb_solvation_model(output: Any) -> str:
    """The solvation model this result applied, or ``gas_phase``."""

    return _xtb_solvation_context(output)[0]


def _xtb_solvent(output: Any) -> str:
    """The solvent this result applied; a gas-phase run has none."""

    model, solvent = _xtb_solvation_context(output)
    if model == "gas_phase":
        raise MissingQuantityError("this xTB result ran in the gas phase")
    if not solvent:
        raise MissingQuantityError(
            "this solvated xTB result names no solvent in its setup block"
        )
    return solvent


def _xtb_solvation_energy(selector: str) -> Callable[[Any], float]:
    """Read one solvation free-energy term of the applied model."""

    def _read(output: Any) -> float:
        values = _xtb_solvation_terms(output)
        if selector in values:
            return values[selector]
        model, _solvent = _xtb_solvation_context(output)
        if model == "gas_phase":
            raise MissingQuantityError(
                "this xTB result ran in the gas phase, so it has no "
                f"{selector}"
            )
        raise MissingQuantityError(
            f"this {model} xTB result printed no {selector} term in its "
            "energy summary"
        )

    return _read


def _xtb_level(output: Any) -> dict[str, Any]:
    """Return the applied xTB Hamiltonian and solvent from result bytes."""

    level: dict[str, Any] = {}
    hamiltonian = getattr(output, "hamiltonian", None)
    if hamiltonian not in (None, ""):
        level["method"] = str(hamiltonian)
    if getattr(output, "solvent_on", False):
        model = getattr(output, "solvent_model", None)
        solvent = getattr(output, "solvent_id", None)
        if model not in (None, ""):
            level["solvent_model"] = str(model)
        if solvent not in (None, ""):
            level["solvent"] = str(solvent)
    return level


def _xtb_stationarity_gradient(output: Any) -> float | None:
    """max|g| at the geometry this xTB result's numbers belong to, Eh/Bohr.

    xTB is the second program whose reader can bind a gradient to one
    structure, and it is the one where the question bites hardest.
    ChemSmart's Hessian job is ``--hess``: the second derivatives are
    taken at the geometry the job was handed and that geometry is never
    relaxed, so a GFN2 Hessian at a structure optimised with another
    Hamiltonian -- or handed in from another program -- is the ordinary
    case rather than the exotic one.  One xTB invocation touches one
    geometry, so the gradient it writes and the spectrum it prints
    describe the same structure, which is exactly what ORCA's and
    Gaussian's per-step ``forces`` cannot promise.

    The number is read from a gradient **vector** -- ``.engrad`` first,
    then the Turbomole ``gradient`` sidecar -- and never from the
    ``GRADIENT NORM`` the summary prints.  A norm over 3N components is
    an upper bound on the largest one, so reporting it here would raise
    the anomaly at geometries whose largest component is under the
    criterion: a wrong number under a right name.  A result that wrote
    no gradient answers nothing, which is what "the host cannot say"
    means, and the receipt records ``unmeasured``.
    """

    import numpy as np

    vectors = getattr(output, "final_forces", None)
    if vectors is None:
        gradient_file = getattr(output, "gradient_file", None)
        forces = getattr(gradient_file, "forces", None)
        vectors = forces[-1] if forces else None
    if vectors is None:
        return None
    values = np.asarray(vectors, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        return None
    if not values.size or not bool(np.isfinite(values).all()):
        return None
    # The gradient belongs to this result's structure only if it counts
    # the same atoms.  A stale or foreign sidecar is an absence, never a
    # number read against the wrong molecule.
    molecule = getattr(output, "molecule", None)
    num_atoms = getattr(molecule, "num_atoms", None)
    if num_atoms is not None and int(num_atoms) != int(values.shape[0]):
        return None
    return float(np.max(np.abs(values)))


def _xtb_geometry_source_path(output: Any, selector: str) -> Path | None:
    """Name the xTB sidecar only for the reached-optimisation selector."""

    if selector != "reached_positions":
        return None
    geometry_file = getattr(output, "xtbopt_geometry_file", None)
    path = getattr(geometry_file, "filepath", None)
    return Path(str(path)) if path else None


def _xtb_native_evidence_paths(output: Any, selector: str) -> tuple[Path, ...]:
    """Return xTB sidecars whose bytes a selector directly consumes."""

    def _sidecar(source: Any) -> Path | None:
        """The path of one native table, unless it is the artifact itself."""

        if source is None or source is getattr(output, "main_out", None):
            return None
        path = getattr(source, "filepath", None) or getattr(
            source, "filename", None
        )
        return Path(str(path)) if path else None

    # The spectrum table carries the frequency column and the IR column as
    # one paired observation, and the degeneracy grouping is read from that
    # same column.  Bind their receipts to that exact table so a later
    # comparison can establish row-wise correspondence rather than merely
    # matching list lengths from two unrelated files.
    if selector in {
        "ir_intensities",
        "vibrational_frequencies",
        "vibrational_mode_degeneracy_group",
    }:
        spectrum = _sidecar(getattr(output, "ir_spectrum_source", None))
        return (spectrum,) if spectrum else ()

    # Per-atom participation reads both tables: the displacement vectors of
    # the mode table and the frequency list the row index is shared with.
    # Sealing one of the two would leave the delivered rows standing on
    # bytes no receipt names.  Seven real Hessians agreed row for row between
    # the two (at most 0.0050 cm^-1, print rounding), so the pairing is
    # sealed here and not guarded: no disagreement has been observed to
    # refuse.
    if selector == "vibrational_mode_atom_participation":
        paths = (
            _sidecar(getattr(output, "ir_spectrum_source", None)),
            _sidecar(getattr(output, "vibrational_mode_source", None)),
        )
        return tuple(path for path in paths if path)

    file_attribute = {
        "wiberg_bond_orders": "wbo_file",
        "xtb_scc_atomic_charges": "charges_file",
    }.get(selector)
    if file_attribute is None:
        return ()
    native_file = getattr(output, file_attribute, None)
    path = getattr(native_file, "filepath", None) or getattr(
        native_file, "filename", None
    )
    return (Path(str(path)),) if path else ()


def _gaussian_frontier(attribute: str) -> Callable[[Any], float]:
    """Closed-shell frontier value; open-shell refuses rather than collapses.

    The ORCA frontier repair defined HOMO/LUMO across both spin channels
    because the unrestricted extremum can straddle them.  The Gaussian
    parser's single-channel values have not been audited against an
    open-shell log here, so an unrestricted result refuses instead of
    silently reporting one channel as the frontier.
    """

    def _read(output: Any) -> float:
        if int(getattr(output, "multiplicity", 1) or 1) != 1:
            raise MissingQuantityError(
                "open-shell Gaussian frontier orbitals are not audited "
                "cross-spin; only closed-shell results resolve "
                f"{attribute!r} here"
            )
        return float(getattr(output, attribute))

    return _read


def _xtb_ir_intensities(output: Any) -> list[float]:
    """Ordered per-normal-mode IR absorption intensities in km/mol."""
    from chemsmart.analysis.result_quantities import QuantityContractError

    raw_frequencies = getattr(output, "vibrational_frequencies", None)
    if raw_frequencies is None:
        raise MissingQuantityError(
            "xTB IR intensities require corresponding vibrational frequencies, "
            "but vibrational_frequencies is absent"
        )
    try:
        frequencies = [float(item) for item in raw_frequencies]
    except (TypeError, ValueError) as err:
        raise QuantityContractError(
            f"malformed vibrational frequencies: {err}"
        ) from err

    if not frequencies:
        raise MissingQuantityError(
            "xTB IR intensities require corresponding vibrational frequencies, "
            "but vibrational_frequencies is empty"
        )

    raw_intensities = getattr(output, "ir_intensities", None)
    if raw_intensities is None:
        raise MissingQuantityError(
            "xTB result carries no IR intensity evidence"
        )

    try:
        intensities = [float(item) for item in raw_intensities]
    except (TypeError, ValueError) as err:
        raise QuantityContractError(
            f"malformed IR intensities: {err}"
        ) from err

    if not intensities:
        raise MissingQuantityError(
            "xTB result carries empty IR intensity evidence"
        )

    for val in intensities:
        if not math.isfinite(val):
            raise QuantityContractError(
                f"non-finite IR intensity value: {val}"
            )

    if len(intensities) != len(frequencies):
        raise QuantityContractError(
            f"xTB IR intensities count ({len(intensities)}) does not match "
            f"vibrational frequencies count ({len(frequencies)})"
        )

    return intensities


def _xtb_accessors() -> dict[str, Callable[[Any], Any]]:
    # The xTB parser resolves charge, multiplicity, and (after a Hessian)
    # the free energy; an expert review found them parsed but undeclared,
    # with this very comment claiming the opposite.
    return {
        "energy": _last_energy,
        "charge": _xtb_state_integer("charge"),
        "multiplicity": _xtb_state_integer("multiplicity"),
        "gibbs_free_energy": _xtb_gibbs,
        "reached_positions": _xtb_reached_positions,
        "vibrational_mode_atom_participation": (
            _vibrational_mode_atom_participation
        ),
        "vibrational_mode_degeneracy_group": (
            _vibrational_mode_degeneracy_group
        ),
        "vibrational_frequencies": lambda output: [
            float(item) for item in output.vibrational_frequencies
        ],
        "ir_intensities": _xtb_ir_intensities,
        "positions": _positions,
        "xtb_scc_atomic_charges": _xtb_scc_atomic_charges,
        "connectivity": lambda output: _connectivity_matrix(output.molecule),
        "symbols": _symbols,
        # xTB prints vector components in atomic units but the trailing total
        # in Debye.  Preserve both native tokens: QuantityValue then exposes
        # the real source unit and printed precision while its canonical value
        # remains available in Debye for dimensional arithmetic.
        "dipole_moment": lambda output: list(
            output.molecular_dipole_full_tokens
        ),
        "dipole_moment_magnitude": lambda output: (
            output.total_molecular_dipole_moment_token
        ),
        "homo": lambda output: float(output.homo_energy),
        "lumo": lambda output: float(output.lumo_energy),
        "gap": lambda output: float(output.fmo_gap),
        "dispersion_energy": _xtb_dispersion_energy,
        "wiberg_bond_orders": _xtb_wiberg_bond_orders,
        # What the solvent cost, and what it was.  The plan could already
        # switch ALPB or GBSA on -- ``solvent_model`` and ``solvent_id``
        # are advertised settings -- and no completed xTB result could say
        # which model ran, in what, or what it was worth.
        "solvation_model": _xtb_solvation_model,
        "solvent": _xtb_solvent,
        **{
            selector: _xtb_solvation_energy(selector)
            for selector, _attribute in _XTB_SOLVATION_TERMS
        },
    }


def _xyz_output(path: Path) -> Any:
    """Open a registered XYZ as molecular data, not as a program run."""

    from chemsmart.io.xyz.xyzfile import XYZFile

    return XYZFile(filename=str(path))


_XYZ_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
_XYZ_HARTREE_PATTERNS = (
    re.compile(
        rf"\bEnergy\s*\(\s*Hartree\s*\)\s*:\s*({_XYZ_NUMBER})\b",
        re.IGNORECASE,
    ),
    re.compile(
        rf"^Coordinates from ORCA-job\b.*\sE\s+({_XYZ_NUMBER})\s*$",
        re.IGNORECASE,
    ),
)


def _xyz_energy(output: Any) -> float:
    """Read only an explicitly Hartree-grounded XYZ comment energy."""

    comment = str(output.comments or "").strip()
    for pattern in _XYZ_HARTREE_PATTERNS:
        match = pattern.search(comment)
        if match is not None:
            return float(match.group(1).replace("D", "E").replace("d", "e"))
    raise MissingQuantityError(
        "XYZ energy requires an explicit Energy(Hartree) label or the "
        "ChemSmart/ORCA 'Coordinates from ORCA-job ... E <value>' form"
    )


def _xyz_accessors() -> dict[str, Callable[[Any], Any]]:
    """Expose the quantities ChemSmart's XYZ parser actually establishes."""

    return {
        "energy": _xyz_energy,
        "positions": _positions,
        "connectivity": lambda output: _connectivity_matrix(output.molecule),
        "symbols": _symbols,
        "trajectory_frame_count": lambda output: len(
            _xyz_trajectory_structures(output)
        ),
        "trajectory_start_positions": lambda output: [
            [float(value) for value in row]
            for row in _xyz_trajectory_structures(output)[0].positions
        ],
        "trajectory_end_positions": lambda output: [
            [float(value) for value in row]
            for row in _xyz_trajectory_structures(output)[-1].positions
        ],
        "trajectory_start_connectivity": lambda output: _connectivity_matrix(
            _xyz_trajectory_structures(output)[0]
        ),
        "trajectory_end_connectivity": lambda output: _connectivity_matrix(
            _xyz_trajectory_structures(output)[-1]
        ),
        "trajectory_connectivity_changed": (
            _xyz_trajectory_connectivity_changed
        ),
    }


#: Programs whose results can be read into typed quantities.  PySCF keeps its
#: dedicated structured-HDF5 path in ``result_quantities``; the entries here are
#: the log-parsing programs that had none.
def _pyscf_output(path: Path) -> Any:
    """Open a receipt-bound structured PySCF result, green or failed.

    The HDF5 path carries an admission guard no log format has: the file
    states its own contract version and a sibling run receipt binds these
    exact bytes and the whole ancestry of digests behind them.  That guard
    runs here at open time.  What used to run here as well -- normal
    termination, an empty findings list, a green receipt state -- made a
    failed PySCF result unopenable by any consumer of the shared plane,
    so the recovery route "restart from the geometry the run reached" was
    unreachable for the one program whose artifact records that geometry
    by construction.  Those requirements now live in
    ``admit_for_analysis``, which extraction and thermochemistry call, and
    the per-dataset unit check runs inside each accessor.
    """

    from chemsmart.analysis.result_quantities import (
        result_file_sha256,
        validate_pyscf_result_binding,
    )

    output, _receipt = validate_pyscf_result_binding(
        path, expected_sha256=result_file_sha256(path)
    )
    return output


def _pyscf_admit_for_analysis(path: Path) -> Any:
    """The green admission every PySCF quantity and free energy stands on."""

    from chemsmart.analysis.result_quantities import (
        result_file_sha256,
        validate_pyscf_analysis_artifact,
    )

    return validate_pyscf_analysis_artifact(
        path, expected_sha256=result_file_sha256(path)
    )


def _pyscf_require_units(selector: str, output: Any) -> None:
    """Refuse a dataset whose stored unit is not the one we read it as.

    Absence and disagreement are different answers.  A dataset this run never
    wrote has no unit at all, which is the ordinary absent quantity; a dataset
    that is present under a unit we did not expect means the writer's contract
    and this reader have diverged, and that is a defect to state rather than a
    value to convert.
    """

    from chemsmart.analysis import result_quantities as rq

    expected = rq._SELECTOR_RESULT_UNITS.get(selector, {})
    if not expected:
        # Job identity and the electronic state are read from the spec block
        # rather than a numeric dataset, so there is no stored unit to audit.
        return
    observed = output.result_units
    # A selector may have more than one home -- the spectrum of the
    # geometry a path stage was handed lives under an IRC's start and
    # under a saddle search's seed -- so it is absent only when every
    # home is, while a home that *is* present keeps its unit audited.
    if all(observed.get(path) is None for path in expected):
        raise MissingQuantityError(
            f"pyscf result contains no {selector!r} value "
            f"(datasets absent: {sorted(expected)})"
        )
    wrong = {
        path: {"expected": unit, "observed": observed.get(path)}
        for path, unit in expected.items()
        if observed.get(path) is not None and observed.get(path) != unit
    }
    if wrong:
        raise rq.QuantityExtractionError(
            f"PySCF result units are absent or incompatible: {wrong}"
        )


def _pyscf_spin(name: str, output: Any) -> float:
    """Return one member of the <S^2> diagnostic set.

    The target is derived from the bound multiplicity rather than read, so it
    exists for any result; the observed value and the effective multiplicity
    are written only for an unrestricted reference, and their absence is the
    ordinary meaning of a closed-shell run.
    """

    multiplicity = output.multiplicity
    if not isinstance(multiplicity, int) or multiplicity <= 0:
        raise MissingQuantityError(
            "pyscf result does not establish a positive multiplicity"
        )
    target = (float(multiplicity) ** 2 - 1.0) / 4.0
    if name == "spin_square_target":
        return target
    spin_value = output.results.get("spin_square")
    effective_value = output.results.get("spin_square_effective_multiplicity")
    if spin_value is None or effective_value is None:
        raise MissingQuantityError(
            "pyscf result has no complete <S^2> diagnostic"
        )
    spin_square = float(_first_scalar(spin_value))
    if name == "spin_square":
        return spin_square
    if name == "spin_square_deviation":
        return spin_square - target
    return float(_first_scalar(effective_value))


def _first_scalar(value: Any) -> float:
    """Return a stored scalar that HDF5 may have shaped as a 1-element array."""

    while isinstance(value, (list, tuple)) or hasattr(value, "shape"):
        try:
            value = value[0]
        except (IndexError, KeyError, TypeError):
            break
    return float(value)


def _pyscf_supplied_positions(output: Any) -> list[list[float]]:
    """The geometry the driver was handed, from ``spec/positions``.

    The unit is a spec field rather than a dataset attribute on older
    contracts, so it is checked here by name instead of through the
    dataset-unit guard: a supplied structure in any unit but Angstrom
    would be a writer change, not an absence.
    """

    from chemsmart.analysis import result_quantities as rq

    unit = getattr(output, "supplied_positions_unit", None)
    if unit != "Angstrom":
        raise rq.QuantityExtractionError(
            "PySCF supplied geometry is not in Angstrom: "
            f"spec/unit is {unit!r}"
        )
    values = output.supplied_positions
    if values is None:
        raise MissingQuantityError(
            "pyscf result records no supplied geometry (spec/positions)"
        )
    return [[float(value) for value in row] for row in values]


def _pyscf_optimization_converged(output: Any) -> int:
    """1 when the optimiser converged, 0 when it stopped short.

    None -- a fixed-geometry stage -- is refused rather than served as a
    0 that would read as an observed failure.
    """

    value = output.converged
    if value is None:
        raise MissingQuantityError(
            "this result ran no optimisation, so it has no convergence "
            "to report; 'converged' is declared for opt"
        )
    return int(bool(value))


def _pyscf_functional(output: Any) -> str:
    """The functional this run applied, in the ChemSmart vocabulary.

    ``spec/method`` is the name the project asked for and ``spec/xc`` the
    libxc name the writer resolved it to; where the name has a
    program-neutral meaning the answer is that literal (``b3lypg`` and
    ``b3lyp`` are one functional and answer ``b3lyp``, as ORCA's
    ``B3LYP/G`` and Gaussian's ``B3LYP`` do), and otherwise the requested
    name, as before.
    """

    from chemsmart.jobs.settings import canonical_functional_literal

    if not output.spec.get("xc"):
        raise MissingQuantityError(
            "this result ran a Hartree-Fock reference and names no "
            "functional; the method identity is read by 'ab_initio'"
        )
    value = output.method
    if not value:
        raise MissingQuantityError("pyscf result records no method name")
    return canonical_functional_literal(value) or str(value)


def _pyscf_surface_id(output: Any) -> str:
    """The electronic surface this result recorded, as one token."""

    surface = getattr(output, "surface", None)
    if not surface:
        raise MissingQuantityError(
            "this result records no surface identity; it was written "
            "under a contract older than v6, and the host states the "
            "surface as absent rather than rebuilding one from the "
            "settings it happens to recognise"
        )
    return surface_token(surface)


def _pyscf_ab_initio(output: Any) -> str:
    value = output.spec.get("ab_initio")
    if not value:
        raise MissingQuantityError(
            "this result ran a DFT functional and names no ab initio "
            "method; the functional identity is read by 'functional'"
        )
    return str(value)


def _pyscf_spin_populations(output: Any) -> list[float]:
    """Per-atom Mulliken spin populations of an open-shell result, in order.

    Written by the driver for open-shell references only; a closed shell
    carries the ``not_applicable`` property and is refused rather than
    served as zeros.  Checked against 2S = multiplicity - 1, the sum a
    complete vector has by construction.
    """

    values = output.mulliken_atomic_spin_populations
    if values is None:
        raise MissingQuantityError(
            "mulliken_atomic_spin_populations: this closed-shell result "
            "carries no spin to partition (property not_applicable)"
        )
    values = [float(item) for item in values]
    try:
        multiplicity = int(output.multiplicity)
    except (TypeError, ValueError):
        multiplicity = None
    if multiplicity is not None:
        expected = float(multiplicity - 1)
        total = sum(values)
        if abs(total - expected) > 0.05:
            raise MissingQuantityError(
                f"mulliken_atomic_spin_populations sums to {total:.3f} "
                f"where 2S for multiplicity {multiplicity} is "
                f"{expected:.1f}; the vector is not the complete molecule "
                "in order"
            )
    return values


def _pyscf_excited_records(output: Any) -> list[dict[str, Any]]:
    records = list(getattr(output, "excited_state_records", None) or ())
    if not records:
        raise MissingQuantityError(
            "this result ran no response stage; excited-state selectors "
            "are declared for td and for an opt carrying excited_state_root"
        )
    return records


def _pyscf_manifold_values(output: Any, multiplicity: int, key: str):
    """Values of one restricted manifold, refused for any other manifold.

    An artifact carries one manifold; the singlet selectors answer on a
    singlet manifold and refuse a triplet one, and an unrestricted
    reference has neither, which is an honest absence rather than a
    silent relabelling.
    """

    records = _pyscf_excited_records(output)
    values = [
        record[key]
        for record in records
        if record.get("multiplicity") == multiplicity
    ]
    if not values:
        manifold = getattr(output, "state_manifold", None)
        raise MissingQuantityError(
            f"this result carries the {manifold!r} manifold, not the "
            f"multiplicity-{multiplicity} one"
        )
    return values


def _pyscf_transition_dipoles(output: Any) -> list[list[float]]:
    values = output.transition_dipole_moments
    if values is None:
        raise MissingQuantityError(
            "this result records no transition dipole moments"
        )
    return [[float(item) for item in row] for row in values]


def _pyscf_excited_converged(output: Any) -> list[int]:
    values = output.excited_state_converged
    if values is None:
        raise MissingQuantityError(
            "this result records no per-root convergence"
        )
    return [int(bool(item)) for item in values]


def _pyscf_followed_root(output: Any) -> int:
    value = output.excited_state_followed_root
    if value is None:
        raise MissingQuantityError(
            "this result followed no excited root; the selector is declared "
            "for an opt carrying excited_state_root"
        )
    return int(value)


def _pyscf_correlated_scalar(name: str) -> Callable[[Any], float]:
    def accessor(output: Any) -> float:
        value = getattr(output, name)
        if value is None:
            method = getattr(output, "correlated_method", None)
            detail = (
                f"the {method} stage does not produce it"
                if method
                else "this result ran no correlated stage"
            )
            raise MissingQuantityError(f"{name}: {detail}")
        return float(value)

    return accessor


def _pyscf_channel_eigenvalues(
    output: Any, channel: str
) -> tuple[list[float], list[float]]:
    """``(occupied, virtual)`` orbital energies of one spin channel, in eV.

    The artifact's own ``results/mo_energy`` split by ``results/mo_occ``.
    A restricted reference carries the same array in both channels, and a
    restricted-open-shell one puts its singly occupied orbital in the
    alpha occupied list and the beta virtual list, which is what its one
    spatial orbital set means.
    """

    if str(channel) == "alpha":
        return (
            list(output.alpha_occ_eigenvalues or ()),
            list(output.alpha_virtual_eigenvalues or ()),
        )
    return (
        list(output.beta_occ_eigenvalues or ()),
        list(output.beta_virtual_eigenvalues or ()),
    )


def _pyscf_frontier_pair(output: Any) -> tuple[float, float]:
    """Return (HOMO, LUMO) in eV as extrema over every occupied channel.

    The same definition the ORCA reader states, for the same reason: for
    an unrestricted reference the frontier orbitals need not share a spin
    channel, and the extremum over both is what survives that case.  This
    program had grown its own -- ``homo`` and ``lumo`` refused outright
    for any open shell while ``gap`` was served as the lowest virtual of
    either channel minus the *highest SOMO*, which pairs a channel's
    occupied level with the other channel's virtual one.  On the archived
    hydroxyl radical that pairing is the alpha and beta halves of one
    singly occupied orbital, so 4.80 eV was reported under a name a
    session reads as a frontier separation when the beta channel's own
    separation is 4.07.  One selector name, one meaning, on both
    programs; a question about one channel is asked through the
    spin-resolved selectors beside these.
    """

    occupied: list[float] = []
    virtual: list[float] = []
    for channel in ("alpha", "beta"):
        channel_occupied, channel_virtual = _pyscf_channel_eigenvalues(
            output, channel
        )
        occupied.extend(channel_occupied)
        virtual.extend(channel_virtual)
    if not occupied or not virtual:
        raise MissingQuantityError(
            "pyscf result establishes no occupied and virtual orbital pair "
            "(results/mo_energy read through results/mo_occ)"
        )
    return max(occupied), min(virtual)


def _pyscf_channel_frontier(
    output: Any, channel: str, occupied: bool
) -> float:
    values = _pyscf_channel_eigenvalues(output, channel)[0 if occupied else 1]
    if not values:
        raise MissingQuantityError(
            f"pyscf result establishes no "
            f"{'occupied' if occupied else 'virtual'} {channel} orbital"
        )
    return max(values) if occupied else min(values)


def _pyscf_solvation_model(output: Any) -> str:
    """The continuum model the run applied, from what it recorded.

    Read from the applied spec rather than from the project's request:
    the driver resolves the PCM variant and the dielectric in the target
    environment, and the artifact records what was attached.
    """

    value = getattr(output, "solvent_model", None)
    if not value:
        raise MissingQuantityError(
            "this pyscf result applied no continuum solvent model"
        )
    return str(value)


def _pyscf_solvent(output: Any) -> str:
    """The solvent the continuum was parameterised for."""

    value = getattr(output, "solvent_id", None)
    if not value:
        raise MissingQuantityError(
            "this pyscf result names no solvent (no continuum model was "
            "applied)"
        )
    return str(value)


def _pyscf_static_dielectric(output: Any) -> float:
    """The permittivity the continuum polarised the density with.

    A solvent name is not a solvent: two legs of one cycle can both say
    "water" and run at different permittivities, and PySCF's PCM runs at
    water's whatever name it was given unless the driver sets the number.
    The artifact records what was set, so a claim can carry it.
    """

    value = getattr(output, "solvent_dielectric", None)
    if value is None:
        if not getattr(output, "solvent_model", None):
            raise MissingQuantityError(
                "this pyscf result applied no continuum solvent model, so "
                "there is no dielectric constant to report"
            )
        raise MissingQuantityError(
            "this pyscf result applied a continuum model and recorded no "
            "dielectric constant (spec/solvent_eps)"
        )
    return float(value)


def _pyscf_response_dielectric(output: Any) -> float:
    """The permittivity the excitation stage's fast term applied.

    Declared beside the excitation energies for the reason ``solvent``
    and ``solvation_model`` are declared beside the solvation terms: a
    number whose solvent a claim cannot carry is how two spectra come to
    be compared as a solvatochromic shift. Here the trap is sharper than
    a missing name, because the static permittivity a session *can* see
    is not the one the response ran on. PySCF 2.14 answers every
    non-equilibrium response with 1.78; a run named for toluene records
    2.3741 on its density and 1.78 on its spectrum, and only the second
    number bounds how far that spectrum can shift.
    """

    value = getattr(output, "response_dielectric", None)
    if value is None:
        if getattr(output, "td_stage", None) is None:
            raise MissingQuantityError(
                "this pyscf result ran no response stage, so no dielectric "
                "was applied to an excitation"
            )
        if not getattr(output, "solvent_model", None):
            raise MissingQuantityError(
                "this pyscf response stage ran in the gas phase, so no "
                "dielectric was applied to it"
            )
        raise MissingQuantityError(
            "this pyscf response stage recorded no applied dielectric "
            "(status/stages/td/solvent/response_eps_applied)"
        )
    return float(value)


def _pyscf_stability_eigenvalue_unit() -> str:
    """The unit the driver records stability eigenvalues in (its table)."""

    from chemsmart.jobs.pyscf.settings import PYSCF_STABILITY_EIGENVALUE_UNIT

    return PYSCF_STABILITY_EIGENVALUE_UNIT


def _pyscf_stability_record(output: Any) -> Mapping[str, Any] | None:
    """The stability record this result carries, or None."""

    record = getattr(output, "scf_stability", None)
    return record if isinstance(record, Mapping) else None


def _pyscf_printed_stability(output: Any, question: str) -> tuple[str, ...]:
    """Where this run's own PySCF log prints its answer to one question.

    A pointer, never a reading: the typed answer is the record's, and a
    record written before the driver listened to the analysis does not
    hold real -> complex or any eigenvalue although PySCF logged both.
    Saying "not determined" there without saying where PySCF said it would
    be a statement of absence the output contradicts.  The log is this
    run's only when the driver configuration it echoes carries this
    artifact's run nonce; the echoed driver script is skipped, and a
    Davidson's eigenvalue line is paired with the verdict PySCF notes
    after it, in PySCF's own words (``PYSCF_STABILITY_PRINTED_KINDS``).
    """

    from chemsmart.jobs.pyscf.settings import PYSCF_STABILITY_PRINTED_KINDS

    kinds = {
        kind
        for kind, name in PYSCF_STABILITY_PRINTED_KINDS.items()
        if name == question
    }
    nonce = str((getattr(output, "spec", None) or {}).get("run_nonce") or "")
    log = Path(str(getattr(output, "logfile", "") or ""))
    if not kinds or not nonce or not log.is_file() or log.is_symlink():
        return ()
    try:
        if log.stat().st_size > 64 * 1024 * 1024:
            return ()
        lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ()
    if not any(f'"run_nonce": "{nonce}"' in line for line in lines[:200]):
        return ()
    found: list[str] = []
    in_script = pending = False
    for number, line in enumerate(lines, 1):
        text = line.strip()
        if "#INFO: **** input file is" in text:
            in_script = True
            continue
        if in_script:
            in_script = "input file end" not in text
            continue
        head, sep, _rest = text.partition(": lowest eigs of H = ")
        if sep:
            pending = head.rsplit("_", 1)[-1] in kinds
            if pending:
                found.append(f"{log.name}:{number}: {text[:120]}")
            continue
        if pending and text.startswith("<class ") and "wavefunction" in text:
            found.append(f"{log.name}:{number}: {text[:120]}")
            pending = False
    return tuple(found[-2:])


def _printed_pointer(output: Any, question: str) -> str:
    printed = _pyscf_printed_stability(output, question)
    if not printed:
        return ""
    return "; this run's own PySCF log prints it: " + " | ".join(printed)


def _pyscf_stability_entry(output: Any, question: str) -> Mapping[str, Any]:
    """One question's entry from the recorded analysis, or why there is none.

    Four absences a session reads differently, so each names itself: the
    run was never asked, the artifact predates result contract v7, PySCF
    could not answer this question for this reference (an ROHF reference
    has no external answer at all), or it returned no answer. None of
    them is stability, which is why absence is spelled out rather than
    defaulted.  A fifth belongs to real -> complex alone: a record written
    before the driver listened to PySCF's analysis names it not
    determined, with the reason, and that reason is what is said.
    """

    record = _pyscf_stability_record(output)
    if not isinstance(record, Mapping):
        if getattr(output, "scf_stability_requested", None) is False:
            raise MissingQuantityError(
                "this run was not asked for a stability analysis (project "
                "key scf_stability); an absent analysis is not a stable "
                "reference"
            )
        raise MissingQuantityError(
            "this pyscf result records no stability analysis: it was "
            "either not asked or written before result contract v7, and "
            "neither says the reference is stable"
        )
    entry = (record.get("analyses") or {}).get(question)
    if not isinstance(entry, Mapping):
        undetermined = (record.get("not_determined") or {}).get(question)
        if isinstance(undetermined, Mapping):
            raise MissingQuantityError(
                f"this stability record names the {question!r} question "
                f"({undetermined.get('rotation_space') or question}) not "
                f"determined: {undetermined.get('reason') or 'no reason'}"
                + _printed_pointer(output, question)
            )
        raise MissingQuantityError(
            f"this stability analysis carries no {question!r} question"
        )
    if entry.get("unavailable"):
        raise MissingQuantityError(
            f"PySCF could not answer the {question!r} stability question "
            f"for this {record.get('reference_family') or 'reference'} "
            f"reference: {str(entry['unavailable']).strip()}"
        )
    return entry


def _pyscf_stability_verdict(question: str) -> Callable[[Any], str]:
    """``stable`` or ``unstable`` for one orbital-rotation question.

    Deliberately not served under Gaussian's ``wavefunction_stability_*``
    names. Gaussian prints one unnamed verdict in a vocabulary of its
    own; PySCF answers two questions separately and names the rotation
    space of each, and on triplet dioxygen at UKS the two answers differ.
    One name carrying both would make two programs look comparable where
    only one of them says which question it answered -- the same trap as
    two equal level strings that are not equal methods.
    """

    def accessor(output: Any) -> str:
        entry = _pyscf_stability_entry(output, question)
        answered = entry.get("stable")
        if answered is None:
            raise MissingQuantityError(
                f"PySCF returned no answer to the {question!r} stability "
                "question for this reference"
            )
        return "stable" if bool(answered) else "unstable"

    return accessor


def _pyscf_stability_rotation_space(output: Any) -> str:
    """Which rotation the external answer is an answer about.

    ``external`` is not one question: PySCF searches RHF/RKS -> UHF/UKS
    for a restricted reference and UHF/UKS -> GHF/GKS for an unrestricted
    one, so ``unstable`` means a different thing in each. The space
    travels as a value a claim can carry rather than as a display detail,
    because the last goal on this surface showed that a fact a session
    never has to ask for is a fact it may not read.
    """

    entry = _pyscf_stability_entry(output, "external")
    space = entry.get("rotation_space")
    if not space:
        raise MissingQuantityError(
            "this external stability answer names no rotation space"
        )
    return str(space)


def _pyscf_stability_lowest_eigenvalue(
    question: str,
) -> Callable[[Any], float]:
    """The lowest eigenvalue PySCF found for one rotation question, in Eh.

    The number the verdict is drawn from: PySCF calls the reference
    unstable when it lies below its threshold (-1e-5 Eh), so a value just
    above zero is a reference close to an instability and a value just
    below the threshold one barely across it -- what a bare word hides.
    It is a root of PySCF's orbital Hessian in PySCF's own normalisation,
    comparable across PySCF results and against that threshold, and not
    against another program's stability matrix, which is why the name is
    PySCF's own family and not Gaussian's.  The record states the unit
    it wrote; one that differs from the driver's table is a divergence,
    not an absence.
    """

    def accessor(output: Any) -> float:
        from chemsmart.analysis import result_quantities as rq

        entry = _pyscf_stability_entry(output, question)
        values = entry.get("lowest_eigenvalues")
        if not values:
            raise MissingQuantityError(
                f"this {question!r} stability answer records no "
                "eigenvalues: it was written before the driver kept the "
                "numbers PySCF's analysis logged"
                + _printed_pointer(output, question)
            )
        record = _pyscf_stability_record(output) or {}
        written = record.get("eigenvalue_unit")
        if written != _pyscf_stability_eigenvalue_unit():
            raise rq.QuantityExtractionError(
                f"this stability record states its eigenvalues in "
                f"{written!r}; the driver's table writes "
                f"{_pyscf_stability_eigenvalue_unit()!r}"
            )
        return float(values[0])

    return accessor


def _pyscf_absence_reason(
    read: Callable[[Any], Any],
    output: Any,
    fallback: MissingQuantityError,
) -> MissingQuantityError:
    """The accessor's own account of an absence, or the generic one.

    Reading a value that the unit audit has already found absent cannot
    return one: the caller raises whichever error comes back either way,
    so nothing this produces reaches a consumer as a quantity.
    """

    try:
        read(output)
    except MissingQuantityError as reason:
        return reason
    except Exception:  # noqa: BLE001 - the generic absence still stands
        return fallback
    return fallback


def _pyscf_decomposition_scalar(name: str) -> Callable[[Any], float]:
    """A term of the total, or the reason this result does not have it.

    Three absences that a single "not present" would blur, because each
    means something different to a session reading the number beside it:
    a gas-phase run has no continuum at all, a PCM-family run has
    electrostatics and no CDS term by construction, and an artifact
    written before contract v10 was never asked to record either.
    """

    def accessor(output: Any) -> float:
        value = getattr(output, name, None)
        if value is not None:
            return float(value)
        record = getattr(output, "solvation_property_status", None)
        if name == "dispersion_energy":
            detail = "this result applied no dispersion correction"
        elif not getattr(output, "solvent_model", None):
            detail = "this result applied no continuum solvent model"
        elif name == "solvation_nonelectrostatic_energy":
            detail = (
                "the %s model carries electrostatics only; the "
                "cavitation-dispersion-solvent-structure term exists "
                "for SMD alone" % str(output.solvent_model)
            )
        elif record is None:
            detail = (
                "this artifact was written before result contract v10, "
                "which is when the decomposition began to be recorded"
            )
        else:
            detail = (
                "the run applied a continuum model and recorded no "
                "electrostatic term"
            )
        raise MissingQuantityError(f"{name}: {detail}")

    return accessor


def _pyscf_scf_energy(output: Any) -> float:
    value = output.scf_energy
    if value is None:
        raise MissingQuantityError("this result records no SCF energy")
    return float(value)


def _pyscf_total_energy(output: Any) -> float:
    """The energy of the surface the job computed on (see ``total_energy``)."""

    value = output.total_energy
    if value is None:
        raise MissingQuantityError("this result records no total energy")
    return float(value)


def _pyscf_level(output: Any) -> dict[str, Any]:
    """The level this result computed at, from the artifact's own record.

    The method and basis the spec names and the solvent the SCF was
    attached to; the frozen-core count a correlated stage applied (PySCF
    correlates every electron unless told otherwise, so 0 is a level and
    never an absence); and the response an excited stage ran on -- its
    method, its manifold, the roots requested -- with the root an
    optimisation followed.  A matching functional string across programs
    is necessary and never sufficient, which is why the level is shown
    and never compared.
    """

    spec = output.spec if isinstance(output.spec, Mapping) else {}
    level: dict[str, Any] = {}
    ab_initio = spec.get("ab_initio")
    method = spec.get("method")
    if ab_initio not in (None, ""):
        level["ab_initio"] = ab_initio
    elif method not in (None, ""):
        # ``method`` is the functional the project asked for when no ab
        # initio method was named; the level states the literal it means,
        # as the ``functional`` selector does, so one functional is one
        # value on every program's level.
        from chemsmart.jobs.settings import canonical_functional_literal

        level["functional"] = canonical_functional_literal(method) or method
        level["dispersion"] = _normalized_dispersion(spec.get("dispersion"))
    if spec.get("basis") not in (None, ""):
        level["basis"] = spec["basis"]
    if getattr(output, "solvent_on", False):
        level["solvent_model"] = output.solvent_model
        level["solvent"] = output.solvent_id
        if output.solvent_dielectric is not None:
            level["solvent_dielectric"] = output.solvent_dielectric
    correlated = getattr(output, "correlated_method", None)
    if correlated:
        level["ab_initio"] = correlated
        frozen = output.frozen_core_applied
        if frozen is not None:
            level["frozen_core"] = int(frozen)
    stage = getattr(output, "td_stage", None)
    if isinstance(stage, Mapping):
        for key, name in (
            ("response_method_applied", "response_method"),
            ("state_manifold_applied", "state_manifold"),
            ("nstates_requested", "nstates"),
        ):
            if stage.get(key) is not None:
                level[name] = stage[key]
        # Which continuum the spectrum itself ran under. The equilibrium
        # word is what makes the second permittivity the operative one,
        # and it is a level fact rather than a quantity: it says how the
        # numbers beside it were computed, not what they are.
        if output.response_dielectric is not None:
            level["excitation_response_dielectric"] = (
                output.response_dielectric
            )
        if output.response_equilibrium_solvation is not None:
            level["excitation_response_solvation"] = (
                "equilibrium"
                if output.response_equilibrium_solvation
                else "non_equilibrium"
            )
    record = getattr(output, "excited_state_record", None)
    if isinstance(record, Mapping) and record.get("root") is not None:
        level["excited_state_root"] = int(record["root"])
        for name in ("response_method", "state_manifold", "nstates"):
            if record.get(name) is not None:
                level.setdefault(name, record[name])
    # What the basis was built of: the angular form the driver passes to
    # ``pyscf.M`` and the core potentials the molecule it built carried.
    from chemsmart.jobs.pyscf.writer import PYSCF_BASIS_FUNCTIONS

    if "basis" in level:
        level["basis_functions"] = PYSCF_BASIS_FUNCTIONS
    symbols = [str(symbol) for symbol in spec.get("symbols") or ()]
    try:
        cores = output.ecp_core_electrons
    except Exception:  # noqa: BLE001 - an unreadable record states nothing
        cores = None
    if cores is not None and symbols:
        level["ecp_core_electrons"] = dict(cores)
    requested = spec.get("frozen_core")
    rule = (
        "pyscf_default"
        if requested in (None, "")
        else "pyscf_auto" if str(requested).lower() == "auto" else None
    )
    _add_frozen_core_conventions(level, symbols, rule)
    # The family the driver built and the validator held the runtime class
    # to; the broken-symmetry singlet is the one uks/uhf at spin 0.
    family = str(spec.get("reference_family") or "").strip().lower()
    _add_reference_identity(
        level, family or None, spec.get("broken_symmetry") is True
    )
    return level


def _pyscf_irc_path_energies(output: Any) -> list[float]:
    """The energy at every accepted frame of an IRC branch, saddle first."""

    values = getattr(output, "irc_path_energies", None)
    if not values:
        raise MissingQuantityError(
            "pyscf result records no IRC path (results/irc/path_energies)"
        )
    return [float(item) for item in values]


def _pyscf_ts_seed_frequencies(output: Any) -> list[float]:
    """The spectrum at the geometry a saddle search was handed.

    Served under the same name an IRC's start spectrum is, because it is
    the same fact -- the harmonic spectrum of the structure this stage
    started from, on the surface it walked -- and a second name would buy
    a second vocabulary for one question.  What differs is the promise:
    ``START_POINT_PROMISES`` declares one imaginary mode for an ``irc``
    and declares nothing for a ``ts``, whose seed is a guess and is
    allowed to be anything.
    """

    values = getattr(output, "ts_seed_frequencies", None)
    if values is None:
        raise MissingQuantityError(
            "pyscf result records no transition-state seed spectrum "
            "(results/ts/seed_frequencies)"
        )
    return [float(item) for item in values]


def _pyscf_irc_start_frequencies(output: Any) -> list[float]:
    """The saddle's spectrum on the surface the IRC walked.

    Served under its own name, never as ``vibrational_frequencies``: the
    numbers belong to the supplied geometry, and every other property of
    an IRC artifact belongs to where the branch ended.
    """

    values = getattr(output, "irc_start_frequencies", None)
    if not values:
        raise MissingQuantityError(
            "pyscf result records no IRC start spectrum "
            "(results/irc/start_frequencies)"
        )
    return [float(item) for item in values]


def _pyscf_irc_direction(output: Any) -> str:
    """The branch the IRC was asked for, which the validator holds the
    measured first step to."""

    value = getattr(output, "irc_direction", None)
    if value not in {"forward", "backward"}:
        raise MissingQuantityError("result does not establish an IRC branch")
    return str(value)


def _pyscf_accessors() -> dict[str, Callable[[Any], Any]]:
    """Selector name to a callable reading it from a structured PySCF result.

    This was an ``if``/``elif`` chain on a second extraction path that built
    its own quantity values, carried its own unit table and had no job-type
    declaration gate -- so a single point could be asked for frequencies and
    only a runtime absence message stood in the way, while the capability
    query could not report that this program answers any selector at all.
    Reading through the shared plane also cross-checks each value's observed
    dimension against the declared one, which the separate path never did.
    """

    raw: dict[str, Callable[[Any], Any]] = {
        # The energy of the surface the job computed on: the SCF for an HF
        # or DFT job and for a spectrum, the followed root's total for an
        # excited-root optimisation, the correlated total otherwise --
        # ORCA's "FINAL SINGLE POINT ENERGY" semantics.  ``energies`` stays
        # the SCF trace in stage order and ``scf_energy`` names the
        # reference at the final geometry.
        "energy": _pyscf_total_energy,
        "scf_energy": _pyscf_scf_energy,
        "energies": lambda output: [float(item) for item in output.energies],
        "positions": lambda output: [
            [float(value) for value in row] for row in output.positions
        ],
        # The structure the optimiser stopped on, converged or not: the
        # same dataset as ``positions`` here, because every PySCF quantity
        # belongs to the final structure; declared for ``opt`` only, since
        # a fixed-geometry stage reaches nothing beyond what it was handed.
        "reached_positions": lambda output: [
            [float(value) for value in row] for row in output.positions
        ],
        "supplied_positions": _pyscf_supplied_positions,
        "converged": _pyscf_optimization_converged,
        "functional": _pyscf_functional,
        "ab_initio": _pyscf_ab_initio,
        "surface_id": _pyscf_surface_id,
        "mulliken_atomic_spin_populations": lambda output: (
            _pyscf_spin_populations(output)
        ),
        "symbols": lambda output: [
            str(symbol) for symbol in output.chemical_symbols
        ],
        "connectivity": lambda output: _connectivity_matrix(
            output.get_molecule()
        ),
        "charge": lambda output: int(output.charge),
        "multiplicity": lambda output: int(output.multiplicity),
        "method": lambda output: str(output.method),
        "basis": lambda output: str(output.basis),
        # The frontier pair and the four spin-resolved levels, all read
        # from this artifact's own orbital energies and occupations. An
        # open-shell result served no HOMO and no LUMO at all and served
        # a ``gap`` built from a different pairing; a radical in a redox
        # or hydrogen-transfer workflow needs the channel it is actually
        # asking about.
        "homo": lambda output: _pyscf_frontier_pair(output)[0],
        "lumo": lambda output: _pyscf_frontier_pair(output)[1],
        "gap": lambda output: (lambda pair: pair[1] - pair[0])(
            _pyscf_frontier_pair(output)
        ),
        "alpha_homo": lambda output: _pyscf_channel_frontier(
            output, "alpha", True
        ),
        "alpha_lumo": lambda output: _pyscf_channel_frontier(
            output, "alpha", False
        ),
        "beta_homo": lambda output: _pyscf_channel_frontier(
            output, "beta", True
        ),
        "beta_lumo": lambda output: _pyscf_channel_frontier(
            output, "beta", False
        ),
        "dipole_moment": lambda output: [
            float(value) for value in output.dipole_moment
        ],
        "dipole_moment_magnitude": lambda output: (
            sum(float(value) ** 2 for value in output.dipole_moment) ** 0.5
        ),
        "mulliken_atomic_charges": lambda output: _per_atom_vector(
            output.mulliken_atomic_charges,
            [str(symbol) for symbol in output.chemical_symbols],
            quantity="mulliken_atomic_charges",
        ),
        "vibrational_frequencies": lambda output: [
            float(item) for item in output.vibrational_frequencies
        ],
        "vibrational_mode_atom_participation": (
            _vibrational_mode_atom_participation
        ),
        "vibrational_mode_degeneracy_group": (
            _vibrational_mode_degeneracy_group
        ),
        # The response stage (contract v5): roots are ascending indices
        # within one manifold at this artifact's own geometry.  The
        # excitation energies are stored in hartree and read through the
        # declared source unit; the record shape is the log readers' own,
        # so a root is selected by index in an expression exactly as for
        # ORCA, and the manifold-specific selectors refuse the other
        # manifold and the unrestricted one rather than relabelling it.
        "excitation_energies": lambda output: [
            float(item) for item in output.excitation_energies
        ],
        "oscillator_strengths": lambda output: [
            float(item) for item in output.oscillator_strengths
        ],
        "excited_state_indices": lambda output: [
            int(item["state_index"]) for item in _pyscf_excited_records(output)
        ],
        "excited_state_manifold_roots": lambda output: [
            int(item["manifold_root"])
            for item in _pyscf_excited_records(output)
        ],
        "excited_state_multiplicities": lambda output: [
            int(item)
            for item in _required_record_values(
                _pyscf_excited_records(output),
                "multiplicity",
                "a spin multiplicity (an unrestricted manifold has none)",
            )
        ],
        "singlet_excitation_energies": lambda output: [
            float(item)
            for item in _pyscf_manifold_values(output, 1, "energy_eV")
        ],
        "triplet_excitation_energies": lambda output: [
            float(item)
            for item in _pyscf_manifold_values(output, 3, "energy_eV")
        ],
        "singlet_oscillator_strengths": lambda output: [
            float(item)
            for item in _pyscf_manifold_values(
                output, 1, "oscillator_strength"
            )
        ],
        "triplet_oscillator_strengths": lambda output: [
            float(item)
            for item in _pyscf_manifold_values(
                output, 3, "oscillator_strength"
            )
        ],
        "transition_dipole_moments": _pyscf_transition_dipoles,
        "excited_state_dominant_excitations": lambda output: (
            _dominant_excitation_values(
                _pyscf_excited_records(output), "labels"
            )
        ),
        "excited_state_dominant_weights": lambda output: (
            _dominant_excitation_values(
                _pyscf_excited_records(output), "weights"
            )
        ),
        "excited_state_converged": _pyscf_excited_converged,
        "excited_state_followed_root": _pyscf_followed_root,
        # The correlated stage: the program's own components at the final
        # geometry.  ``correlation_energy`` is the final method's whole
        # correlation, triples included, as the ORCA reader means it.
        # The decomposition of the total (contract v10): which continuum
        # was applied, in what solvent, and what the program itself put
        # into the energy it reports.  An energy whose solvation a
        # consumer cannot request is an energy two legs of a
        # thermodynamic cycle can disagree about silently.
        "solvation_model": _pyscf_solvation_model,
        "solvent": _pyscf_solvent,
        # The two permittivities the run applied. They are identities, as
        # the model and the solvent name are, and they are the identities
        # that decide what a difference between two solvated results is
        # allowed to be called.
        "solvent_dielectric": _pyscf_static_dielectric,
        "excitation_response_dielectric": _pyscf_response_dielectric,
        # What PySCF's own analysis said about the reference every number
        # above it stands on. The host sensor has read this record since
        # contract v7 and raised an anomaly on it; nothing could extract
        # it, so a session could not state, or claim, that the orbitals
        # its energy came from are not a minimum in rotation space.
        "scf_stability_internal": _pyscf_stability_verdict("internal"),
        "scf_stability_external": _pyscf_stability_verdict("external"),
        "scf_stability_external_rotation_space": (
            _pyscf_stability_rotation_space
        ),
        # The question PySCF solves inside its external analysis and does
        # not return, and the number behind each verdict. Recorded since
        # the driver listened to the analysis; before, the record said
        # real -> complex was not determined while the log beside it said.
        "scf_stability_real_to_complex": _pyscf_stability_verdict(
            "real_to_complex"
        ),
        "scf_stability_internal_lowest_eigenvalue": (
            _pyscf_stability_lowest_eigenvalue("internal")
        ),
        "scf_stability_external_lowest_eigenvalue": (
            _pyscf_stability_lowest_eigenvalue("external")
        ),
        "scf_stability_real_to_complex_lowest_eigenvalue": (
            _pyscf_stability_lowest_eigenvalue("real_to_complex")
        ),
        "solvation_electrostatic_energy": _pyscf_decomposition_scalar(
            "solvation_electrostatic_energy"
        ),
        "solvation_nonelectrostatic_energy": _pyscf_decomposition_scalar(
            "solvation_nonelectrostatic_energy"
        ),
        "dispersion_energy": _pyscf_decomposition_scalar("dispersion_energy"),
        "reference_energy": _pyscf_correlated_scalar("reference_energy"),
        "correlation_energy": _pyscf_correlated_scalar("correlation_energy"),
        "ccsd_correlation_energy": _pyscf_correlated_scalar(
            "ccsd_correlation_energy"
        ),
        "triples_correction": _pyscf_correlated_scalar("triples_correction"),
        # The IRC stage (contract v8): the accepted path, saddle first, read
        # through the same trajectory vocabulary the log readers answer for
        # a Gaussian IRC or an ORCA path sidecar -- here from the artifact
        # itself, so the start is the geometry the run was handed and the
        # end is the structure every other property belongs to.
        "trajectory_frame_count": lambda output: len(_irc_structures(output)),
        "trajectory_start_positions": lambda output: [
            [float(value) for value in row]
            for row in _irc_structures(output)[0].positions
        ],
        "trajectory_end_positions": lambda output: [
            [float(value) for value in row]
            for row in _irc_structures(output)[-1].positions
        ],
        "trajectory_start_connectivity": lambda output: (
            _connectivity_matrix(_irc_structures(output)[0])
        ),
        "trajectory_end_connectivity": lambda output: _connectivity_matrix(
            _irc_structures(output)[-1]
        ),
        "trajectory_connectivity_changed": _trajectory_connectivity_changed,
        "trajectory_energies": _pyscf_irc_path_energies,
        # One accessor for one question: the spectrum of the geometry a
        # path stage was handed, whichever stage ran.
        "trajectory_start_frequencies": lambda output: (
            _pyscf_ts_seed_frequencies(output)
            if getattr(output, "ts_stage", None) is not None
            else _pyscf_irc_start_frequencies(output)
        ),
        "irc_direction": _pyscf_irc_direction,
        "irc_converged": _irc_run_converged,
    }
    for name in (
        "spin_square",
        "spin_square_target",
        "spin_square_deviation",
        "effective_multiplicity",
    ):
        raw[name] = (lambda key: lambda output: _pyscf_spin(key, output))(name)

    def _guard(
        selector: str, read: Callable[[Any], Any]
    ) -> Callable[[Any], Any]:
        """Audit the stored unit, and let the accessor explain an absence.

        The unit check answers "is the dataset there, and is it stored as
        we read it".  It cannot answer *why* a selector is not there, and
        the absences differ: a gas-phase run applied no continuum at all,
        a PCM-family run has no cavitation term by construction, and an
        artifact written under an older contract was never asked to
        record either -- three facts a session reads differently and one
        "dataset absent" message flattens.  So when the dataset is gone
        the accessor gets to say which absence it is, and the generic
        message stands wherever the accessor has nothing to add.
        """

        def accessor(output: Any) -> Any:
            try:
                _pyscf_require_units(selector, output)
            except MissingQuantityError as unit_absence:
                raise _pyscf_absence_reason(
                    read, output, unit_absence
                ) from None
            return read(output)

        return accessor

    return {name: _guard(name, read) for name, read in raw.items()}


#: Selectors every executed PySCF stage writes.  The SCF block that follows
#: every stage stores energies, geometry, orbital energies and the population
#: and dipole properties, so a single point and an optimisation answer the
#: same set; a Hessian stage adds the vibrational quantities on top.
#: The environment set -- which continuum was applied, in what solvent, and
#: the terms the program put into its own total -- is declared here for the
#: same reason the spin populations are: every stage that converges an SCF
#: can answer it, and a run that applied no continuum refuses it as absent.
_PYSCF_SCF_SELECTORS = (
    "ab_initio",
    "alpha_homo",
    "alpha_lumo",
    "beta_homo",
    "beta_lumo",
    "surface_id",
    "basis",
    "charge",
    "connectivity",
    "dipole_moment",
    "dipole_moment_magnitude",
    "dispersion_energy",
    "effective_multiplicity",
    "energies",
    "energy",
    "functional",
    "gap",
    "homo",
    "lumo",
    "method",
    "mulliken_atomic_charges",
    "mulliken_atomic_spin_populations",
    "multiplicity",
    "positions",
    "scf_energy",
    "scf_stability_external",
    "scf_stability_external_lowest_eigenvalue",
    "scf_stability_external_rotation_space",
    "scf_stability_internal",
    "scf_stability_internal_lowest_eigenvalue",
    "scf_stability_real_to_complex",
    "scf_stability_real_to_complex_lowest_eigenvalue",
    "solvation_electrostatic_energy",
    "solvation_model",
    "solvation_nonelectrostatic_energy",
    "solvent",
    "solvent_dielectric",
    "spin_square",
    "spin_square_deviation",
    "spin_square_target",
    "supplied_positions",
    "symbols",
)

#: The response stage: declared for ``td`` and, because an excited-root
#: optimisation re-evaluates the spectrum at the reached geometry, for
#: ``opt``, where a ground-state optimisation refuses them as absent.
#: ``absorption_wavelengths`` is deliberately not here (``photon_wavelength``
#: is the one conversion authority), nor ``excited_state_labels`` (PySCF
#: prints none) nor ``excited_state_spin_square`` (PySCF prints no per-root
#: <S^2>).
_PYSCF_TD_SELECTORS = (
    "excitation_energies",
    "excited_state_converged",
    "excited_state_dominant_excitations",
    "excited_state_dominant_weights",
    "excited_state_indices",
    "excited_state_manifold_roots",
    "excited_state_multiplicities",
    "oscillator_strengths",
    "singlet_excitation_energies",
    "singlet_oscillator_strengths",
    "transition_dipole_moments",
    "triplet_excitation_energies",
    "triplet_oscillator_strengths",
)

#: The correlated stage: declared for ``sp`` and ``opt``, where an HF or
#: DFT result refuses them as absent.
_PYSCF_CORR_SELECTORS = (
    "ccsd_correlation_energy",
    "correlation_energy",
    "reference_energy",
    "triples_correction",
)

#: An optimisation additionally reports whether it converged and the
#: structure it reached, which for this program is the one every other
#: quantity belongs to; an excited-root optimisation also names its root.
_PYSCF_OPT_SELECTORS = tuple(
    sorted(
        _PYSCF_SCF_SELECTORS
        + _PYSCF_TD_SELECTORS
        + _PYSCF_CORR_SELECTORS
        + ("converged", "excited_state_followed_root", "reached_positions")
    )
)
_PYSCF_SP_SELECTORS = tuple(
    sorted(_PYSCF_SCF_SELECTORS + _PYSCF_CORR_SELECTORS)
)
#: An IRC branch (contract v8): the SCF set belongs to where the branch
#: ended, as for an optimisation, and the path adds the trajectory
#: vocabulary -- the start as supplied, the end as reached -- with the
#: energy at every frame and the saddle's own spectrum on the walked
#: surface. ``vibrational_frequencies`` is deliberately absent: the
#: artifact holds no Hessian at the endpoint, and serving the saddle's
#: under that name would describe the wrong structure.
_PYSCF_IRC_SELECTORS = tuple(
    sorted(
        _PYSCF_SCF_SELECTORS
        + (
            "irc_converged",
            "irc_direction",
            "reached_positions",
            "trajectory_connectivity_changed",
            "trajectory_end_connectivity",
            "trajectory_end_positions",
            "trajectory_energies",
            "trajectory_frame_count",
            "trajectory_start_connectivity",
            "trajectory_start_frequencies",
            "trajectory_start_positions",
        )
    )
)
#: ``excitation_response_dielectric`` is declared here and not in
#: ``_PYSCF_TD_SELECTORS``, which an ``opt`` inherits: an excited-root
#: optimisation is gas phase only -- PySCF has no solvated excited-state
#: gradient and the settings validator refuses one -- so on that job type
#: the selector could never be answered, and declaring a question a job
#: type cannot ever have is worse than not declaring it.
_PYSCF_TD_JOBTYPE_SELECTORS = tuple(
    sorted(
        _PYSCF_SCF_SELECTORS
        + _PYSCF_TD_SELECTORS
        + ("excitation_response_dielectric",)
    )
)
#: A saddle search (contract v9): the SCF set belongs to where the climb
#: ended, as for an optimisation, plus the spectrum of the seed it was
#: handed.  The trajectory vocabulary is deliberately absent: the frames
#: of a search are an optimiser's route to a structure and not a path on
#: the surface, and serving them under the name an IRC's path answers to
#: would invite reading a climb as a reaction coordinate.
#: ``vibrational_frequencies`` is absent for the reason it is on an
#: ``irc``: the artifact holds no Hessian where the search ended, and
#: serving the seed's under that name would describe the wrong structure.
_PYSCF_TS_SELECTORS = tuple(
    sorted(
        _PYSCF_SCF_SELECTORS
        + (
            "converged",
            "reached_positions",
            "trajectory_start_frequencies",
        )
    )
)

#: What each PySCF selector's value belongs to.  One structure per
#: result: the driver re-converges the SCF on the final geometry before
#: any property is read, so everything but the supplied structure and
#: the structure-free identities is ``as_reached`` -- and for ``sp`` and
#: ``hess`` the reached structure is the supplied one by construction,
#: which the validator enforces to 1e-8 A.  ``energies`` stays stateless
#: on purpose: on an opt it is [E(supplied), E(reached)].
_PYSCF_STRUCTURAL_STATES = tuple(
    sorted(
        [
            ("ccsd_correlation_energy", "as_reached"),
            ("connectivity", "as_reached"),
            ("converged", "as_reached"),
            ("correlation_energy", "as_reached"),
            ("alpha_homo", "as_reached"),
            ("alpha_lumo", "as_reached"),
            ("beta_homo", "as_reached"),
            ("beta_lumo", "as_reached"),
            ("dipole_moment", "as_reached"),
            ("dipole_moment_magnitude", "as_reached"),
            ("dispersion_energy", "as_reached"),
            ("effective_multiplicity", "as_reached"),
            ("energy", "as_reached"),
            ("excitation_energies", "as_reached"),
            ("excited_state_converged", "as_reached"),
            ("excited_state_dominant_excitations", "as_reached"),
            ("excited_state_dominant_weights", "as_reached"),
            ("excited_state_indices", "as_reached"),
            ("excited_state_manifold_roots", "as_reached"),
            ("excited_state_multiplicities", "as_reached"),
            ("gap", "as_reached"),
            ("homo", "as_reached"),
            ("irc_converged", "as_reached"),
            ("lumo", "as_reached"),
            ("mulliken_atomic_charges", "as_reached"),
            ("mulliken_atomic_spin_populations", "as_reached"),
            ("oscillator_strengths", "as_reached"),
            ("positions", "as_reached"),
            ("reached_positions", "as_reached"),
            ("reference_energy", "as_reached"),
            ("scf_energy", "as_reached"),
            # An answer about the density the final SCF converged, so it
            # belongs to the structure that SCF ran on, as every other
            # mean-field property here does.
            ("scf_stability_external", "as_reached"),
            ("scf_stability_external_lowest_eigenvalue", "as_reached"),
            ("scf_stability_external_rotation_space", "as_reached"),
            ("scf_stability_internal", "as_reached"),
            ("scf_stability_internal_lowest_eigenvalue", "as_reached"),
            ("scf_stability_real_to_complex", "as_reached"),
            (
                "scf_stability_real_to_complex_lowest_eigenvalue",
                "as_reached",
            ),
            ("singlet_excitation_energies", "as_reached"),
            ("singlet_oscillator_strengths", "as_reached"),
            ("solvation_electrostatic_energy", "as_reached"),
            ("solvation_nonelectrostatic_energy", "as_reached"),
            ("spin_square", "as_reached"),
            ("spin_square_deviation", "as_reached"),
            ("spin_square_target", "as_reached"),
            ("supplied_positions", "as_supplied"),
            ("surface_id", "stateless"),
            ("trajectory_end_connectivity", "as_reached"),
            ("trajectory_end_positions", "as_reached"),
            ("trajectory_start_connectivity", "as_supplied"),
            ("trajectory_start_frequencies", "as_supplied"),
            ("trajectory_start_positions", "as_supplied"),
            ("transition_dipole_moments", "as_reached"),
            ("triples_correction", "as_reached"),
            ("triplet_excitation_energies", "as_reached"),
            ("triplet_oscillator_strengths", "as_reached"),
            ("vibrational_frequencies", "as_reached"),
            ("vibrational_mode_atom_participation", "as_reached"),
            ("vibrational_mode_degeneracy_group", "as_reached"),
        ]
    )
)

#: Whose density or method each PySCF value belongs to.  The mean-field
#: properties the driver reads off ``mf`` after the final SCF -- dipole,
#: populations, orbital energies, the spin diagnostic -- are the
#: reference's, on every configuration: an excited-root optimisation and
#: a correlated result carry them beside a total that is not the
#: reference's.  ``energy`` is the surface the job computed on and is
#: resolved per artifact.
_PYSCF_ELECTRONIC_PROVENANCE = tuple(
    sorted(
        [
            ("alpha_homo", "reference"),
            ("alpha_lumo", "reference"),
            ("beta_homo", "reference"),
            ("beta_lumo", "reference"),
            ("ccsd_correlation_energy", "correlated"),
            ("correlation_energy", "correlated"),
            ("dipole_moment", "reference"),
            ("dipole_moment_magnitude", "reference"),
            # The decomposition is the mean field's: PySCF adds both
            # solvation terms and the dispersion correction into the
            # reference's own total, whatever surface the job went on to
            # compute, so a correlated or excited-root artifact carries
            # them beside a total that is not the reference's.
            ("dispersion_energy", "reference"),
            ("solvation_electrostatic_energy", "reference"),
            ("solvation_nonelectrostatic_energy", "reference"),
            ("effective_multiplicity", "reference"),
            ("energies", "reference"),
            ("energy", "computed_surface"),
            ("excitation_energies", "excited_root"),
            ("excited_state_converged", "excited_root"),
            ("excited_state_dominant_excitations", "excited_root"),
            ("excited_state_dominant_weights", "excited_root"),
            ("excited_state_followed_root", "excited_root"),
            ("excited_state_indices", "excited_root"),
            ("excited_state_manifold_roots", "excited_root"),
            ("excited_state_multiplicities", "excited_root"),
            ("gap", "reference"),
            ("homo", "reference"),
            ("lumo", "reference"),
            ("mulliken_atomic_charges", "reference"),
            ("mulliken_atomic_spin_populations", "reference"),
            ("oscillator_strengths", "excited_root"),
            ("reference_energy", "reference"),
            ("scf_energy", "reference"),
            # The analysis is about the reference itself, on every
            # configuration: a correlated or excited-root artifact
            # carries a total that is not the reference's beside a
            # verdict that is.
            ("scf_stability_external", "reference"),
            ("scf_stability_external_lowest_eigenvalue", "reference"),
            ("scf_stability_external_rotation_space", "reference"),
            ("scf_stability_internal", "reference"),
            ("scf_stability_internal_lowest_eigenvalue", "reference"),
            ("scf_stability_real_to_complex", "reference"),
            ("scf_stability_real_to_complex_lowest_eigenvalue", "reference"),
            ("surface_id", "stateless"),
            ("trajectory_energies", "reference"),
            ("trajectory_start_frequencies", "reference"),
            ("singlet_excitation_energies", "excited_root"),
            ("singlet_oscillator_strengths", "excited_root"),
            ("spin_square", "reference"),
            ("spin_square_deviation", "reference"),
            ("spin_square_target", "reference"),
            ("transition_dipole_moments", "excited_root"),
            ("triples_correction", "correlated"),
            ("triplet_excitation_energies", "excited_root"),
            ("triplet_oscillator_strengths", "excited_root"),
            ("vibrational_frequencies", "reference"),
            ("vibrational_mode_atom_participation", "reference"),
            ("vibrational_mode_degeneracy_group", "reference"),
        ]
    )
)

#: ORCA's word for the same axis, over the selectors its reader implements:
#: the excitation set belongs to a root, the mean-field properties to the
#: reference, and ``energy`` (``FINAL SINGLE POINT ENERGY``) to the surface
#: the job computed on, resolved from the ``ab_initio`` selector.
_ORCA_ELECTRONIC_PROVENANCE_DECLARED = (
    ("absorption_wavelengths", "excited_root"),
    ("alpha_homo", "reference"),
    ("alpha_lumo", "reference"),
    ("beta_homo", "reference"),
    ("beta_lumo", "reference"),
    ("correlation_energy", "correlated"),
    ("dipole_moment", "reference"),
    ("dipole_moment_magnitude", "reference"),
    ("dispersion_energy", "reference"),
    ("effective_multiplicity", "reference"),
    ("energies", "computed_surface"),
    ("energy", "computed_surface"),
    ("entropy_times_temperature", "computed_surface"),
    ("excitation_energies", "excited_root"),
    ("excited_state_dominant_excitations", "excited_root"),
    ("excited_state_dominant_weights", "excited_root"),
    ("excited_state_indices", "excited_root"),
    ("excited_state_labels", "excited_root"),
    ("excited_state_manifold_roots", "excited_root"),
    ("excited_state_multiplicities", "excited_root"),
    ("excited_state_spin_square", "excited_root"),
    ("gap", "reference"),
    ("gibbs_free_energy", "computed_surface"),
    ("hirshfeld_atomic_charges", "reference"),
    ("hirshfeld_atomic_spin_populations", "reference"),
    ("homo", "reference"),
    ("loewdin_atomic_charges", "reference"),
    ("loewdin_atomic_spin_populations", "reference"),
    ("lumo", "reference"),
    ("mayer_bond_orders", "reference"),
    ("mayer_free_valence", "reference"),
    ("mulliken_atomic_charges", "reference"),
    ("mulliken_atomic_spin_populations", "reference"),
    ("oscillator_strengths", "excited_root"),
    ("reference_energy", "reference"),
    ("scf_energy", "reference"),
    ("singlet_excitation_energies", "excited_root"),
    ("singlet_oscillator_strengths", "excited_root"),
    ("spin_square", "reference"),
    ("spin_square_after_annihilation", "reference"),
    ("spin_square_deviation", "reference"),
    ("spin_square_target", "reference"),
    # The singles amplitudes of the coupled-cluster calculation itself.
    ("t1_diagnostic", "correlated"),
    # Every point of an ORCA IRC path is a total on the surface the job
    # computed on, exactly as ``energy`` is, and is resolved the same way.
    ("trajectory_energies", "computed_surface"),
    ("triplet_excitation_energies", "excited_root"),
    ("triplet_oscillator_strengths", "excited_root"),
)


RESULT_READERS: dict[str, ResultReaderV1] = {
    "orca": ResultReaderV1(
        program="orca",
        artifact_kind="orca_output",
        parser_id="chemsmart.io.orca.output.ORCAOutput",
        open_output=_orca_output,
        accessors=_orca_accessors(),
        resolve_cartesian_hessian=_orca_cartesian_hessian,
        resolve_torsional_scan=_orca_torsional_scan,
        resolve_convergence_check=_orca_convergence_check,
        resolve_frozen_atoms=_orca_frozen_atoms,
        # An ORCA IRC's product structure is a sidecar beside the log, as
        # xTB's reached optimisation frame is: the geometry handoff seals
        # that file's digest, and extraction carries its bytes on the
        # receipt of every selector that reads it.
        geometry_source_path_for_selector=_orca_geometry_source_path,
        native_evidence_paths_for_selector=_orca_native_evidence_paths,
        #: A held coordinate is a distance for a bond and an angle for a
        #: torsion, so the kinds answer under their own names and each
        #: keeps its true unit. ``scan_coordinate_values`` has to declare
        #: none because a scan drives one coordinate whose kind the
        #: declaration cannot know; a constraint table states the kind of
        #: every row, so nothing here has to be dimensionless.
        #: An atom row and a count are pure numbers, and this vocabulary
        #: spells that ``1`` rather than leaving the unit empty: an
        #: introduced selector says what it is measured in, and "" reads
        #: as nobody having said.
        selector_declarations=(
            _CONSTRAINED_COORDINATE_DECLARATIONS
            + _EXCITED_CHARACTER_DECLARATIONS
            + (
                ("t1_diagnostic", "1", "DIMENSIONLESS"),
                ("mayer_bond_orders", "1", "DIMENSIONLESS"),
                ("mayer_free_valence", "1", "DIMENSIONLESS"),
            )
            + _HIRSHFELD_SPIN_DECLARATION
        ),
        #: An atom index this plane delivers indexes the vectors this
        #: plane delivers -- symbols, positions, every population -- so it
        #: is zero-based like all of them, and it says so where the model
        #: reads it rather than where the reader remembers it. ORCA's own
        #: label counts from one; that is a label, never the index, and
        #: the two are converted at the accessor. What the Agent *writes*
        #: -- the constrained coordinate on a modred node -- is one-based,
        #: which is the round trip this record exists to keep honest.
        atom_resolved_declarations=(
            _CONSTRAINED_COORDINATE_ATOM_DECLARATIONS
            + _HIRSHFELD_SPIN_ATOM_DECLARATION
            + (
                (
                    "mayer_bond_orders",
                    tuple(
                        {
                            "semantic_quantity": "bond_order",
                            "population_scheme": "Mayer",
                            "atom_order": "zero-based molecular atom order",
                            "data_shape": "rows of [atom_i, atom_j, order]",
                            "sparsity": (
                                "ORCA prints only orders above 0.1; an "
                                "omitted pair has no printed value and is "
                                "not zero"
                            ),
                        }.items()
                    ),
                ),
                (
                    "mayer_free_valence",
                    (
                        ("semantic_quantity", "free_valence"),
                        ("population_scheme", "Mayer"),
                        ("atom_order", "zero-based molecular atom order"),
                    ),
                ),
            )
        ),
        # Coverage is ``parser_supported_when_emitted``: it states what a job
        # of this type can be asked for, while method and settings still
        # decide whether the engine prints it.  The spin family and the
        # dispersion term are settings-dependent rather than job-dependent --
        # an unrestricted reference prints <S^2> at any job type -- so they
        # belong to every declaration here.  Frequencies and the
        # thermochemistry derived from them appear only where ChemSmart lets a
        # frequency step run: ``sp`` and ``td`` force ``freq`` off, while
        # ``opt`` and ``ts`` inherit the project's setting.  Job types without
        # a declaration stay unknown rather than guessed.  The frontier-orbital
        # family follows the converged reference, so it is declared wherever an
        # SCF converges -- ``sp``, ``opt``, ``ts`` -- and deliberately not on
        # ``td``, where the ground-state orbitals are not what the job was run
        # to answer.
        # What each selector's value belongs to. Read from the
        # accessors, not asserted: `_orca_positions` calls
        # `output.thermochemistry_molecule`, which is the geometry the
        # Hessian was computed at -- step 0 for an `OptTS Freq`, because
        # ChemSmart's own ORCA ts writer computes the Hessian first. So
        # `positions` is `thermochemistry_reference`, and the structure
        # the optimiser actually reached is `output.molecule`, which the
        # four `trajectory_*` accessors expose and which no jobtype
        # declared.
        selector_structural_states=(
            ("connectivity", "thermochemistry_reference"),
            # A held coordinate is delivered as the structure ORCA
            # returned actually has it, cross-checked against the value
            # ORCA declared it was holding. Held means the two agree, and
            # the reader checks rather than repeats the declaration.
            ("constrained_bond_angles", "as_reached"),
            ("constrained_bond_lengths", "as_reached"),
            ("constrained_dihedral_angles", "as_reached"),
            ("correlation_energy", "as_reached"),
            ("dispersion_energy", "as_reached"),
            ("energy", "as_reached"),
            ("entropy_times_temperature", "thermochemistry_reference"),
            ("gibbs_free_energy", "thermochemistry_reference"),
            ("positions", "thermochemistry_reference"),
            ("reached_positions", "as_reached"),
            ("reference_energy", "as_reached"),
            ("scan_coordinate_values", "scan_point"),
            ("scan_energies", "scan_point"),
            ("scan_point_indices", "scan_point"),
            ("scf_energy", "as_reached"),
            ("solvation_cavity_surface_area", "as_reached"),
            ("solvation_electrostatic_energy", "as_reached"),
            ("solvation_nonelectrostatic_energy", "as_reached"),
            ("t1_diagnostic", "as_reached"),
            # A comparison across the two ends of one branch belongs to
            # neither of them alone.
            ("trajectory_connectivity_changed", "trajectory_endpoint"),
            # The structure an IRC branch reached is what "reached" means
            # for a run that walks rather than optimises, and it is the
            # word PySCF's own irc endpoint already carries. It used to
            # read the log's single printed structure -- the saddle -- and
            # was declared for no jobtype, so nothing consumed the word
            # ``trajectory_endpoint`` it used to hold; now the endpoint is
            # ORCA's own sidecar and the recovery route can carry it.
            ("trajectory_end_connectivity", "as_reached"),
            ("trajectory_end_positions", "as_reached"),
            ("trajectory_start_connectivity", "as_supplied"),
            ("trajectory_start_positions", "as_supplied"),
            ("vibrational_frequencies", "thermochemistry_reference"),
            (
                "vibrational_mode_atom_participation",
                "thermochemistry_reference",
            ),
            ("vibrational_mode_degeneracy_group", "thermochemistry_reference"),
            ("vpt2_fundamental_frequencies", "thermochemistry_reference"),
            ("vpt2_harmonic_frequencies", "thermochemistry_reference"),
            (
                "vpt2_zero_point_rovibrational_energy",
                "thermochemistry_reference",
            ),
        ),
        selector_electronic_provenance=_electronic_provenance_table(
            _orca_accessors(), _ORCA_ELECTRONIC_PROVENANCE_DECLARED
        ),
        resolve_electronic_provenance=_resolve_computed_surface,
        resolve_level=_orca_level,
        resolve_surface=lambda output: surface_from_accessors(
            RESULT_READERS["orca"], output
        ),
        jobtype_selectors=(
            (
                "freq",
                (
                    # A frequency job at a fixed geometry: everything an
                    # optimisation's output means minus the optimisation
                    # claim. ORCA reads `! Freq` without `Opt` as this
                    # jobtype, and a budget-bound session that folded its
                    # frequencies into a single-point stage had every
                    # extraction refused for want of it (live, 2026-09-03).
                    "ab_initio",
                    "alpha_homo",
                    "alpha_lumo",
                    "basis",
                    "beta_homo",
                    "beta_lumo",
                    "charge",
                    "connectivity",
                    "correlation_energy",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "dispersion_energy",
                    "effective_multiplicity",
                    "energies",
                    "energy",
                    "functional",
                    "gap",
                    "hirshfeld_atomic_charges",
                    "hirshfeld_atomic_spin_populations",
                    "homo",
                    "loewdin_atomic_charges",
                    "loewdin_atomic_spin_populations",
                    "lumo",
                    "mayer_bond_orders",
                    "mayer_free_valence",
                    "mulliken_atomic_charges",
                    "mulliken_atomic_spin_populations",
                    "multiplicity",
                    "positions",
                    "reference_energy",
                    "scf_energy",
                    "solvation_cavity_surface_area",
                    "solvation_electrostatic_energy",
                    "solvation_model",
                    "solvation_nonelectrostatic_energy",
                    "solvent",
                    "spin_square",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "vibrational_frequencies",
                    "vibrational_mode_atom_participation",
                    "vibrational_mode_degeneracy_group",
                ),
            ),
            (
                "irc",
                # An ORCA IRC log prints a single structure -- the
                # transition state the branch was handed -- so every
                # state-dependent value read from the log body describes
                # the saddle, not the path: the first Agent-executed IRC
                # delivered the saddle's own distances as both endpoints,
                # and its ``energy`` differed from the true endpoint by the
                # entire barrier (33.40 kcal/mol, measured again on job
                # 2142379).  ``energy``, ``positions`` and the orbital,
                # dipole and spin families therefore stay undeclared.
                #
                # The path is not absent, though: ORCA prints it as the
                # IRC PATH SUMMARY table and writes the branch's endpoint
                # to its own ``_IRC_F.xyz`` / ``_IRC_B.xyz`` sidecar.  The
                # trajectory family reads exactly those, each end bound to
                # the other by the path table's own energies, so a
                # completed branch delivers its profile and its product
                # structure -- and the latter is in the ``as_reached``
                # state, which is what lets the recovery route carry it
                # into the optimisation that identifies the minimum.  A
                # ``direction both`` run leaves two branch endpoints and
                # the geometry selectors refuse it by name.
                (
                    "ab_initio",
                    "basis",
                    "charge",
                    "functional",
                    "irc_converged",
                    "irc_direction",
                    "multiplicity",
                    "solvation_model",
                    "solvent",
                    "symbols",
                    "trajectory_connectivity_changed",
                    "trajectory_end_connectivity",
                    "trajectory_end_positions",
                    "trajectory_energies",
                    "trajectory_frame_count",
                    "trajectory_start_connectivity",
                    "trajectory_start_positions",
                ),
            ),
            (
                # A constrained optimisation is an optimisation that was
                # not allowed to finish the job: it relaxes every degree
                # of freedom except the ones it holds, so what it returns
                # is a structure on the surface at a chosen value of a
                # chosen coordinate. That is the calculation whose result
                # seeds a saddle search, and until now a completed one
                # answered nothing -- the classifier already called it
                # ``modred`` and no reader had declared the name, so every
                # selector on it was refused, including its energy.
                #
                # The vibrational family is deliberately absent. A
                # constrained optimum is a stationary point only in the
                # subspace orthogonal to what is held; the gradient along
                # the constraint is whatever it is, ORCA projects nothing
                # out of the Hessian it prints, and the spectrum of that
                # Hessian is not the spectrum of a stationary point. So no
                # frequency, no thermochemistry and no stationary-point
                # claim is promised here, even when the project's modred
                # section runs ``! Opt Freq`` and ORCA prints one. That is
                # the same rule ``scan`` already follows, stated for the
                # job type that can actually emit the numbers.
                "modred",
                (
                    "ab_initio",
                    "alpha_homo",
                    "alpha_lumo",
                    "basis",
                    "beta_homo",
                    "beta_lumo",
                    "charge",
                    "connectivity",
                    "constrained_angle_atoms",
                    "constrained_bond_angles",
                    "constrained_bond_atoms",
                    "constrained_bond_lengths",
                    "constrained_coordinate_count",
                    "constrained_dihedral_angles",
                    "constrained_dihedral_atoms",
                    # ORCA converges a constrained optimisation on the
                    # free coordinates and says so in the same words, so
                    # the flag keeps its meaning: the relaxation finished.
                    "converged",
                    "correlation_energy",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "dispersion_energy",
                    "effective_multiplicity",
                    "energies",
                    "energy",
                    "functional",
                    "gap",
                    "hirshfeld_atomic_charges",
                    "hirshfeld_atomic_spin_populations",
                    "homo",
                    "loewdin_atomic_charges",
                    "loewdin_atomic_spin_populations",
                    "lumo",
                    "mayer_bond_orders",
                    "mayer_free_valence",
                    "mulliken_atomic_charges",
                    "mulliken_atomic_spin_populations",
                    "multiplicity",
                    "positions",
                    # The structure at the held value is the whole point
                    # of running one, and it is what a later saddle search
                    # or single point consumes.
                    "reached_positions",
                    "reference_energy",
                    "scf_energy",
                    "solvation_cavity_surface_area",
                    "solvation_electrostatic_energy",
                    "solvation_model",
                    "solvation_nonelectrostatic_energy",
                    "solvent",
                    "spin_square",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                ),
            ),
            (
                "opt",
                (
                    # Printed thermochemistry is deliberately not
                    # declared for ORCA: 6.x applies quasi-RRHO entropy
                    # by default with no keyword, while the typed
                    # derive_thermochemistry route computes the
                    # convention its receipt states and refuses
                    # imaginary modes the printed value silently drops.
                    # One free-energy route, gated and visible.
                    "ab_initio",
                    "alpha_homo",
                    "alpha_lumo",
                    "basis",
                    "beta_homo",
                    "beta_lumo",
                    "charge",
                    "connectivity",
                    "converged",
                    "correlation_energy",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "dispersion_energy",
                    "effective_multiplicity",
                    "energies",
                    "energy",
                    "functional",
                    "gap",
                    "hirshfeld_atomic_charges",
                    "hirshfeld_atomic_spin_populations",
                    "homo",
                    "loewdin_atomic_charges",
                    "loewdin_atomic_spin_populations",
                    "lumo",
                    "mayer_bond_orders",
                    "mayer_free_valence",
                    "mulliken_atomic_charges",
                    "mulliken_atomic_spin_populations",
                    "multiplicity",
                    "positions",
                    # An optimiser that ran and stopped prints its
                    # last structure; that is what "reached" means,
                    # and it is not what ``positions`` answers.  Not
                    # declared for ``irc`` (the log prints only the
                    # starting point), ``scan`` (the last printed
                    # structure is the last scan point, a different
                    # state), or the fixed-geometry jobtypes, where
                    # there is no second structure to distinguish.
                    "reached_positions",
                    "reference_energy",
                    "scf_energy",
                    "solvation_cavity_surface_area",
                    "solvation_electrostatic_energy",
                    "solvation_model",
                    "solvation_nonelectrostatic_energy",
                    "solvent",
                    "spin_square",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "vibrational_frequencies",
                    "vibrational_mode_atom_participation",
                    "vibrational_mode_degeneracy_group",
                ),
            ),
            (
                # A relaxed scan is a surface, and the surface is what the job
                # was run to establish, so the scan family is declared here
                # rather than left unknown -- the one job type whose whole
                # purpose is producing a profile could not promise it.  The
                # scalars below are the last converged point rather than a
                # stationary state; frequencies and thermochemistry do not
                # appear because ChemSmart runs no frequency step in a scan.
                "scan",
                (
                    "alpha_homo",
                    "alpha_lumo",
                    "beta_homo",
                    "beta_lumo",
                    "charge",
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "dispersion_energy",
                    "effective_multiplicity",
                    # ``energies`` is deliberately not declared: on a
                    # scan result it is the optimizer trace over every
                    # micro-iteration (max-min spread 143.4 kcal/mol on a
                    # real butane scan whose surface spans 5.15), sitting
                    # one natural name away from ``scan_energies``, which
                    # is the surface.
                    "energy",
                    "gap",
                    "homo",
                    "lumo",
                    "multiplicity",
                    "positions",
                    "reference_energy",
                    "scan_coordinate_values",
                    "scan_energies",
                    "scan_point_indices",
                    "scan_steps_planned",
                    "scan_steps_reached",
                    "scf_energy",
                    "solvation_model",
                    "solvent",
                    "spin_square",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                ),
            ),
            (
                "sp",
                (
                    "ab_initio",
                    "alpha_homo",
                    "alpha_lumo",
                    "basis",
                    "beta_homo",
                    "beta_lumo",
                    "charge",
                    "connectivity",
                    "correlation_energy",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "dispersion_energy",
                    "effective_multiplicity",
                    "energies",
                    "energy",
                    "functional",
                    "gap",
                    "hirshfeld_atomic_charges",
                    "hirshfeld_atomic_spin_populations",
                    "homo",
                    "loewdin_atomic_charges",
                    "loewdin_atomic_spin_populations",
                    "lumo",
                    "mayer_bond_orders",
                    "mayer_free_valence",
                    "mulliken_atomic_charges",
                    "mulliken_atomic_spin_populations",
                    "multiplicity",
                    "positions",
                    "reference_energy",
                    "scf_energy",
                    "solvation_cavity_surface_area",
                    "solvation_electrostatic_energy",
                    "solvation_model",
                    "solvation_nonelectrostatic_energy",
                    "solvent",
                    "spin_square",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    # The single-reference check of a coupled-cluster
                    # number, printed beside it.
                    "t1_diagnostic",
                ),
            ),
            (
                "td",
                (
                    "ab_initio",
                    "absorption_wavelengths",
                    "basis",
                    "charge",
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "dispersion_energy",
                    "effective_multiplicity",
                    "energies",
                    "energy",
                    "excitation_energies",
                    # What each root is made of: its largest single
                    # excitation and that excitation's weight.
                    "excited_state_dominant_excitations",
                    "excited_state_dominant_weights",
                    "excited_state_indices",
                    "excited_state_manifold_roots",
                    "excited_state_multiplicities",
                    "excited_state_spin_square",
                    "functional",
                    "multiplicity",
                    "oscillator_strengths",
                    "positions",
                    "reference_energy",
                    "scf_energy",
                    "singlet_excitation_energies",
                    "singlet_oscillator_strengths",
                    "solvation_model",
                    "solvent",
                    "spin_square",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    # ORCA solves the spin-adapted triplets beside the
                    # singlets; they were read and never declared, so a
                    # triplet request could not be answered by name.
                    "triplet_excitation_energies",
                    "triplet_oscillator_strengths",
                ),
            ),
            (
                "ts",
                (
                    # Printed thermochemistry is deliberately not
                    # declared for ORCA: 6.x applies quasi-RRHO entropy
                    # by default with no keyword, while the typed
                    # derive_thermochemistry route computes the
                    # convention its receipt states and refuses
                    # imaginary modes the printed value silently drops.
                    # One free-energy route, gated and visible.
                    "ab_initio",
                    "alpha_homo",
                    "alpha_lumo",
                    "basis",
                    "beta_homo",
                    "beta_lumo",
                    "charge",
                    "connectivity",
                    "converged",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "dispersion_energy",
                    "effective_multiplicity",
                    "energies",
                    "energy",
                    "functional",
                    "gap",
                    "hirshfeld_atomic_charges",
                    "hirshfeld_atomic_spin_populations",
                    "homo",
                    "loewdin_atomic_charges",
                    "loewdin_atomic_spin_populations",
                    "lumo",
                    "mayer_bond_orders",
                    "mayer_free_valence",
                    "mulliken_atomic_charges",
                    "mulliken_atomic_spin_populations",
                    "multiplicity",
                    "positions",
                    # An optimiser that ran and stopped prints its
                    # last structure; that is what "reached" means,
                    # and it is not what ``positions`` answers.  Not
                    # declared for ``irc`` (the log prints only the
                    # starting point), ``scan`` (the last printed
                    # structure is the last scan point, a different
                    # state), or the fixed-geometry jobtypes, where
                    # there is no second structure to distinguish.
                    "reached_positions",
                    "reference_energy",
                    "scf_energy",
                    "solvation_cavity_surface_area",
                    "solvation_electrostatic_energy",
                    "solvation_model",
                    "solvation_nonelectrostatic_energy",
                    "solvent",
                    "spin_square",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "vibrational_frequencies",
                    "vibrational_mode_atom_participation",
                    "vibrational_mode_degeneracy_group",
                ),
            ),
        ),
    ),
    "gaussian": ResultReaderV1(
        program="gaussian",
        artifact_kind="gaussian_output",
        parser_id="chemsmart.io.gaussian.output.Gaussian16Output",
        open_output=_gaussian_output,
        accessors=_gaussian_accessors(),
        resolve_cartesian_hessian=_gaussian_cartesian_hessian,
        resolve_torsional_scan=_gaussian_torsional_scan,
        resolve_convergence_check=_gaussian_convergence_check,
        resolve_frozen_atoms=_gaussian_frozen_atoms,
        # Gaussian prints the SMD-CDS term in kcal/mol to two decimals.
        source_units={"solvation_nonelectrostatic_energy": "kcal/mol"},
        # A held coordinate keeps its unit and its atoms as ORCA's do: one
        # declaration for both programs' constrained optimisations.
        selector_declarations=(
            _CONSTRAINED_COORDINATE_DECLARATIONS
            + _EXCITED_CHARACTER_DECLARATIONS
            # <R**2> of the SCF density about the centre of nuclear
            # charge, as Gaussian prints it (atomic units, bohr^2).
            + (("electronic_spatial_extent", "bohr^2", "AREA"),)
            # The volume keyword's per-molecule figure, bohr^3.
            + (("molecular_volume", "bohr^3", "VOLUME"),)
            # The number and the space beside Gaussian's stability word.
            + (
                ("wavefunction_stability_lowest_eigenvalue", "Eh", "ENERGY"),
                ("wavefunction_stability_rotation_space", "", "DIMENSIONLESS"),
            )
            + _HIRSHFELD_SPIN_DECLARATION
        ),
        atom_resolved_declarations=(
            _CONSTRAINED_COORDINATE_ATOM_DECLARATIONS
            + _HIRSHFELD_SPIN_ATOM_DECLARATION
        ),
        # Coverage is ``parser_supported_when_emitted``, as for ORCA: it
        # states what a job of this type can be asked for, while route and
        # settings still decide what Gaussian prints.  The spin family, the
        # dipole and the wavefunction-stability pair are settings-dependent
        # rather than job-dependent -- an unrestricted reference or a
        # stability check can accompany any job -- so they belong to every
        # declaration.  Frequencies and the thermochemistry derived from
        # them appear where a frequency step runs, and the excited-state
        # family where a TD route does.  Job types with no grounding
        # artifact stay undeclared rather than guessed.
        jobtype_selectors=(
            (
                # A one-direction IRC branch.  ChemSmart writes a Gaussian
                # IRC as two native inputs -- forward and reverse -- and the
                # route of each log says which, so the job type this reader
                # is asked about is never the bare ``irc``: the branch words
                # are what a completed Gaussian IRC log answers to.  Unlike
                # ORCA, whose IRC log prints only where the path started,
                # Gaussian prints every accepted point, so the trajectory
                # family is real here and ``reached_positions`` is where the
                # branch's walk ended.
                "ircf",
                _GAUSSIAN_IRC_BRANCH_SELECTORS,
            ),
            (
                "ircr",
                _GAUSSIAN_IRC_BRANCH_SELECTORS,
            ),
            (
                # A constrained optimisation: the coordinates the
                # ModRedundant section froze are held and everything else
                # relaxes, so the run ends on one converged structure, and
                # the scan family is absent because no coordinate was
                # driven.
                #
                # No Gibbs energy is declared here, for the reason ORCA's
                # modred gives: a constrained optimum is stationary only
                # orthogonal to what it held, so the free energy Gaussian
                # prints after one belongs to no stationary point. The
                # archived fe_ch_quintet_modred_link.log printed a
                # "Sum of electronic and thermal Free Energies" beside a
                # -1380 cm-1 mode it silently left out, and this
                # declaration served that number as gibbs_free_energy.
                #
                # Nor the vibrational family, for the same reason and as
                # ORCA's modred has none: Gaussian's frequency step after a
                # constrained optimisation diagonalises the whole Hessian at
                # a point that is not stationary along what was held and
                # projects nothing out, so those numbers are not the
                # harmonic frequencies of any stationary point -- and a
                # served ``vibrational_frequencies`` is what the host's own
                # thermochemistry and order checks read. The Hessian still
                # runs when the project asks for it.
                "modred",
                (
                    "ab_initio",
                    "basis",
                    "charge",
                    "connectivity",
                    # What it held and whether it finished relaxing the
                    # rest, as ORCA's constrained optimisation answers.
                    "constrained_angle_atoms",
                    "constrained_bond_angles",
                    "constrained_bond_atoms",
                    "constrained_bond_lengths",
                    "constrained_coordinate_count",
                    "constrained_dihedral_angles",
                    "constrained_dihedral_atoms",
                    "converged",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "effective_multiplicity",
                    "electronic_spatial_extent",
                    "energies",
                    "energy",
                    "functional",
                    "gap",
                    "hirshfeld_atomic_charges",
                    "hirshfeld_atomic_spin_populations",
                    "homo",
                    "lumo",
                    "mulliken_atomic_charges",
                    "mulliken_atomic_spin_populations",
                    "multiplicity",
                    "positions",
                    "reached_positions",
                    "scf_energy",
                    "spin_square",
                    "spin_square_after_annihilation",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "wavefunction_stability_history",
                    "wavefunction_stability_lowest_eigenvalue",
                    "wavefunction_stability_rotation_space",
                    "wavefunction_stability_verdict",
                ),
            ),
            (
                "opt",
                (
                    "ab_initio",
                    "basis",
                    "charge",
                    "connectivity",
                    "converged",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "effective_multiplicity",
                    "electronic_spatial_extent",
                    "energies",
                    "energy",
                    "functional",
                    "gap",
                    "gibbs_free_energy",
                    "hirshfeld_atomic_charges",
                    "hirshfeld_atomic_spin_populations",
                    "homo",
                    "ir_intensities",
                    "lumo",
                    "mulliken_atomic_charges",
                    "mulliken_atomic_spin_populations",
                    "multiplicity",
                    "positions",
                    # The structure the optimiser stopped on, as distinct
                    # from ``positions``, which answers with the last
                    # parsed orientation whatever produced it.  For a
                    # Gaussian ``opt freq`` the two agree, because the
                    # spectrum is taken at the converged geometry; the
                    # roles are still different questions, and only this
                    # one is admissible as a structure to carry forward.
                    "reached_positions",
                    "scf_energy",
                    "solvation_model",
                    "solvation_nonelectrostatic_energy",
                    "solvent",
                    "spin_square",
                    "spin_square_after_annihilation",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "vibrational_frequencies",
                    "wavefunction_stability_history",
                    "wavefunction_stability_lowest_eigenvalue",
                    "wavefunction_stability_rotation_space",
                    "wavefunction_stability_verdict",
                ),
            ),
            (
                # A relaxed scan is a surface, and the surface is the thing
                # the job was run to establish.  Gaussian prints no profile
                # table: it prints the optimiser's trace, 29 energies for a
                # 13-point surface on the first scan this reader was
                # validated against, so the profile is the parser's
                # assembly of the converged points and their own driven
                # coordinate.  ``energies`` is deliberately absent for the
                # same reason it is absent from ORCA's scan: it is that
                # optimiser trace, one natural name away from the surface.
                # No frequency step runs in a scan, so nothing vibrational
                # or thermochemical is declared.
                "scan",
                (
                    "ab_initio",
                    "basis",
                    "charge",
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "effective_multiplicity",
                    "energy",
                    "functional",
                    "gap",
                    "homo",
                    "lumo",
                    "multiplicity",
                    "positions",
                    "scan_coordinate_values",
                    "scan_energies",
                    "scan_point_indices",
                    "scan_steps_planned",
                    "scan_steps_reached",
                    "scf_energy",
                    "spin_square",
                    "spin_square_after_annihilation",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "wavefunction_stability_history",
                    "wavefunction_stability_lowest_eigenvalue",
                    "wavefunction_stability_rotation_space",
                    "wavefunction_stability_verdict",
                ),
            ),
            (
                "sp",
                (
                    "ab_initio",
                    "basis",
                    "charge",
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "effective_multiplicity",
                    "electronic_spatial_extent",
                    "energies",
                    "energy",
                    "functional",
                    "gap",
                    "hirshfeld_atomic_charges",
                    "hirshfeld_atomic_spin_populations",
                    "homo",
                    "lumo",
                    "molecular_volume",
                    "mulliken_atomic_charges",
                    "mulliken_atomic_spin_populations",
                    "multiplicity",
                    "positions",
                    "scf_energy",
                    "solvation_model",
                    "solvation_nonelectrostatic_energy",
                    "solvent",
                    "spin_square",
                    "spin_square_after_annihilation",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "wavefunction_stability_history",
                    "wavefunction_stability_lowest_eigenvalue",
                    "wavefunction_stability_rotation_space",
                    "wavefunction_stability_verdict",
                ),
            ),
            (
                # A fixed-geometry response calculation.  Gaussian activates
                # it with a route keyword and nothing else, so the route
                # line alone reads as a single point; the excited-state
                # family therefore lives here rather than on every ``sp``,
                # where it promised a transition list for a ground-state
                # energy.
                "td",
                (
                    "ab_initio",
                    "absorption_wavelengths",
                    "basis",
                    "charge",
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "effective_multiplicity",
                    "electronic_spatial_extent",
                    "energies",
                    "energy",
                    "excitation_energies",
                    # What each root is made of: its largest single
                    # excitation and that excitation's weight.
                    "excited_state_dominant_excitations",
                    "excited_state_dominant_weights",
                    "excited_state_indices",
                    "excited_state_labels",
                    # Which manifold each state belongs to and its rank
                    # there, and each spin block by name: the words ORCA's
                    # and PySCF's td already answered.  The accessors read
                    # Gaussian's own spin labels; undeclared, a 50-50 run
                    # could be read only by list position, where Gaussian
                    # interleaves singlets and triplets by energy.
                    "excited_state_manifold_roots",
                    "excited_state_multiplicities",
                    "excited_state_spin_square",
                    "functional",
                    "hirshfeld_atomic_charges",
                    "hirshfeld_atomic_spin_populations",
                    "mulliken_atomic_charges",
                    "mulliken_atomic_spin_populations",
                    "multiplicity",
                    "oscillator_strengths",
                    "positions",
                    "scf_energy",
                    "singlet_excitation_energies",
                    "singlet_oscillator_strengths",
                    "spin_square",
                    "spin_square_after_annihilation",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "triplet_excitation_energies",
                    "triplet_oscillator_strengths",
                    "wavefunction_stability_history",
                    "wavefunction_stability_lowest_eigenvalue",
                    "wavefunction_stability_rotation_space",
                    "wavefunction_stability_verdict",
                ),
            ),
            (
                "ts",
                (
                    "ab_initio",
                    "basis",
                    "charge",
                    "connectivity",
                    "converged",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "effective_multiplicity",
                    "electronic_spatial_extent",
                    "energies",
                    "energy",
                    "functional",
                    "gibbs_free_energy",
                    "hirshfeld_atomic_charges",
                    "hirshfeld_atomic_spin_populations",
                    "ir_intensities",
                    "mulliken_atomic_charges",
                    "mulliken_atomic_spin_populations",
                    "multiplicity",
                    "positions",
                    "reached_positions",
                    "scf_energy",
                    "spin_square",
                    "spin_square_after_annihilation",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "vibrational_frequencies",
                    "wavefunction_stability_history",
                    "wavefunction_stability_lowest_eigenvalue",
                    "wavefunction_stability_rotation_space",
                    "wavefunction_stability_verdict",
                ),
            ),
        ),
        # ChemSmart runs one Gaussian ``irc`` job as two native
        # one-direction inputs, so the word a plan names the stage with
        # appears in no log it writes.  Nothing joined the two, and the
        # consequence was not a missing selector but a stage that could
        # not be used at all: an extraction reading an IRC producer was
        # refused while the plan was still being built, and no geometry
        # edge could leave the node, although both branch declarations
        # carry ``reached_positions`` in the ``as_reached`` state.
        #
        # The bare word is deliberately not declared above instead.  A
        # route that says ``irc`` with no direction is Gaussian's own
        # both-direction job: one log holding two legs, whose last frame
        # is the end of the second leg and whose first and last frames
        # are not the two ends of one walk.  ``reached_positions`` and
        # the trajectory pair mean something else there, and that log has
        # not been audited.
        stage_result_jobtypes=(("irc", ("ircf", "ircr")),),
        # Which molecular state each value belongs to.  The reader declared
        # none at all until this round, which is not the same as every
        # value being stateless: it meant no consumer could ask this reader
        # for a *role*, so every geometry route refused a Gaussian result
        # outright and the structure a converged optimisation reached could
        # not be carried into the next calculation through the host.
        #
        # For Gaussian the thermochemistry reference and the reached
        # structure coincide -- an ``opt freq`` takes its spectrum at the
        # geometry it converged on, and the log's last orientation is that
        # geometry -- so the same state word is honest for both, and the
        # ORCA case that forced this distinction (a ``OptTS Freq`` whose
        # ``positions`` is step 0) does not arise here.  A scan's values
        # belong to a sampled point, and an IRC branch's to where its walk
        # ended.
        selector_structural_states=tuple(
            sorted(
                [
                    ("charge", "as_reached"),
                    ("connectivity", "as_reached"),
                    # A held coordinate is measured in the structure the
                    # run returned, as ORCA's reader delivers it.
                    ("constrained_bond_angles", "as_reached"),
                    ("constrained_bond_lengths", "as_reached"),
                    ("constrained_dihedral_angles", "as_reached"),
                    ("dipole_moment", "as_reached"),
                    ("dipole_moment_magnitude", "as_reached"),
                    ("electronic_spatial_extent", "as_reached"),
                    ("energy", "as_reached"),
                    ("gap", "as_reached"),
                    ("gibbs_free_energy", "as_reached"),
                    ("hirshfeld_atomic_charges", "as_reached"),
                    ("homo", "as_reached"),
                    ("ir_intensities", "as_reached"),
                    ("lumo", "as_reached"),
                    ("mulliken_atomic_charges", "as_reached"),
                    ("molecular_volume", "as_reached"),
                    ("mulliken_atomic_spin_populations", "as_reached"),
                    ("multiplicity", "as_reached"),
                    ("positions", "as_reached"),
                    ("reached_positions", "as_reached"),
                    ("scan_coordinate_values", "scan_point"),
                    ("scan_energies", "scan_point"),
                    ("scan_point_indices", "scan_point"),
                    ("scf_energy", "as_reached"),
                    ("solvation_nonelectrostatic_energy", "as_reached"),
                    ("symbols", "stateless"),
                    ("trajectory_end_connectivity", "trajectory_endpoint"),
                    ("trajectory_end_positions", "trajectory_endpoint"),
                    (
                        "trajectory_connectivity_changed",
                        "trajectory_endpoint",
                    ),
                    ("trajectory_start_connectivity", "as_supplied"),
                    ("trajectory_start_positions", "as_supplied"),
                    ("vibrational_frequencies", "as_reached"),
                ]
            )
        ),
        resolve_reference_diagnostics=_gaussian_reference_diagnostics,
        resolve_level=_gaussian_level,
        selector_electronic_provenance=_electronic_provenance_table(
            _gaussian_accessors(), _GAUSSIAN_ELECTRONIC_PROVENANCE_DECLARED
        ),
        resolve_electronic_provenance=_resolve_computed_surface,
    ),
    "xtb": ResultReaderV1(
        program="xtb",
        artifact_kind="xtb_output",
        parser_id="chemsmart.io.xtb.output.XTBOutput",
        open_output=_xtb_output,
        accessors=_xtb_accessors(),
        source_units={"dipole_moment": "e bohr"},
        geometry_source_path_for_selector=_xtb_geometry_source_path,
        native_evidence_paths_for_selector=_xtb_native_evidence_paths,
        #: The three ALPB/GBSA terms the shared vocabulary has no true name
        #: for.  ``solvation_free_energy`` is generic on purpose: it is the
        #: whole solvation free energy the applied model charged, which any
        #: program can mean, and the model it belongs to travels beside it.
        selector_declarations=(
            ("solvation_free_energy", "Eh", "ENERGY"),
            ("xtb_solvation_sasa_energy", "Eh", "ENERGY"),
            ("xtb_solvation_hydrogen_bond_energy", "Eh", "ENERGY"),
            ("xtb_solvation_shift_energy", "Eh", "ENERGY"),
        ),
        jobtype_selectors=(
            (
                "hess",
                (
                    "charge",
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "dispersion_energy",
                    "energy",
                    "gap",
                    "gibbs_free_energy",
                    "homo",
                    "ir_intensities",
                    "lumo",
                    "multiplicity",
                    "positions",
                    "solvation_electrostatic_energy",
                    "solvation_free_energy",
                    "solvation_model",
                    "solvent",
                    "symbols",
                    "vibrational_frequencies",
                    "vibrational_mode_atom_participation",
                    "vibrational_mode_degeneracy_group",
                    "wiberg_bond_orders",
                    "xtb_scc_atomic_charges",
                    "xtb_solvation_hydrogen_bond_energy",
                    "xtb_solvation_sasa_energy",
                    "xtb_solvation_shift_energy",
                ),
            ),
            (
                "opt",
                (
                    "charge",
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "dispersion_energy",
                    "energy",
                    "gap",
                    "gibbs_free_energy",
                    "homo",
                    "lumo",
                    "multiplicity",
                    "positions",
                    "reached_positions",
                    "solvation_electrostatic_energy",
                    "solvation_free_energy",
                    "solvation_model",
                    "solvent",
                    "symbols",
                    "wiberg_bond_orders",
                    "xtb_scc_atomic_charges",
                    "xtb_solvation_hydrogen_bond_energy",
                    "xtb_solvation_sasa_energy",
                    "xtb_solvation_shift_energy",
                ),
            ),
            (
                "sp",
                (
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "dispersion_energy",
                    "energy",
                    "gap",
                    "homo",
                    "lumo",
                    "positions",
                    "solvation_electrostatic_energy",
                    "solvation_free_energy",
                    "solvation_model",
                    "solvent",
                    "symbols",
                    "wiberg_bond_orders",
                    "xtb_scc_atomic_charges",
                    "xtb_solvation_hydrogen_bond_energy",
                    "xtb_solvation_sasa_energy",
                    "xtb_solvation_shift_energy",
                ),
            ),
        ),
        selector_structural_states=tuple(
            sorted(
                [
                    ("charge", "as_reached"),
                    ("connectivity", "as_reached"),
                    ("dipole_moment", "as_reached"),
                    ("dipole_moment_magnitude", "as_reached"),
                    ("dispersion_energy", "as_reached"),
                    ("energy", "as_reached"),
                    ("gap", "as_reached"),
                    ("gibbs_free_energy", "as_reached"),
                    ("homo", "as_reached"),
                    ("lumo", "as_reached"),
                    ("multiplicity", "as_reached"),
                    ("positions", "as_reached"),
                    ("reached_positions", "as_reached"),
                    # The solvation terms are charged at the structure the
                    # result reached; the model and solvent names belong to
                    # the run rather than to any one geometry, so they stay
                    # stateless, exactly as ORCA declares them.
                    ("solvation_electrostatic_energy", "as_reached"),
                    ("solvation_free_energy", "as_reached"),
                    ("symbols", "stateless"),
                    ("vibrational_frequencies", "as_reached"),
                    ("ir_intensities", "as_reached"),
                    ("vibrational_mode_atom_participation", "as_reached"),
                    ("vibrational_mode_degeneracy_group", "as_reached"),
                    ("wiberg_bond_orders", "as_reached"),
                    ("xtb_scc_atomic_charges", "as_reached"),
                    ("xtb_solvation_hydrogen_bond_energy", "as_reached"),
                    ("xtb_solvation_sasa_energy", "as_reached"),
                    ("xtb_solvation_shift_energy", "as_reached"),
                ]
            )
        ),
        # xTB's parsed record names its Hamiltonian but does not yet carry a
        # program-neutral electronic-surface identity.  Leave that axis
        # absent rather than labelling every quantity ``computed_surface``
        # without an identity a consumer can resolve.
        resolve_level=_xtb_level,
        resolve_stationarity_gradient=_xtb_stationarity_gradient,
        resolve_convergence_check=_xtb_convergence_check,
    ),
    "pyscf": ResultReaderV1(
        program="pyscf",
        artifact_kind="pyscf_hdf5",
        parser_id="chemsmart.io.pyscf.output.PySCFOutput",
        open_output=_pyscf_output,
        accessors=_pyscf_accessors(),
        # PySCF stores excitation energies in hartree where the log-parsing
        # programs print electronvolts.  One entry closes a cross-program
        # disagreement this project had recorded and never reconciled.
        source_units={"excitation_energies": "Eh"},
        #: Two relative permittivities, both dimensionless, both applied
        #: rather than requested.  The names are program-neutral because
        #: the facts are: every continuum program has a dielectric, and
        #: every non-equilibrium response has a second one.  ORCA prints
        #: both (``Epsilon``, ``Refrac``) and Gaussian both (``Eps``,
        #: ``EpsInf``); neither reader serves them yet, and the shared
        #: vocabulary is where they will meet when one does.
        #: The stability triple is program-local on purpose: Gaussian's
        #: ``wavefunction_stability_verdict`` is one word about one
        #: unnamed question, and PySCF answers two questions and names
        #: the rotation space of each.
        selector_declarations=(
            ("solvent_dielectric", "", "DIMENSIONLESS"),
            ("excitation_response_dielectric", "", "DIMENSIONLESS"),
            ("scf_stability_internal", "", "DIMENSIONLESS"),
            ("scf_stability_external", "", "DIMENSIONLESS"),
            ("scf_stability_external_rotation_space", "", "DIMENSIONLESS"),
            ("scf_stability_real_to_complex", "", "DIMENSIONLESS"),
            # Roots of PySCF's orbital Hessian, in hartree and in PySCF's
            # own normalisation: comparable within PySCF, not across
            # programs, so the names stay this family's.
            ("scf_stability_internal_lowest_eigenvalue", "Eh", "ENERGY"),
            ("scf_stability_external_lowest_eigenvalue", "Eh", "ENERGY"),
            (
                "scf_stability_real_to_complex_lowest_eigenvalue",
                "Eh",
                "ENERGY",
            ),
            *_EXCITED_CHARACTER_DECLARATIONS,
        ),
        jobtype_selectors=(
            (
                "hess",
                tuple(
                    sorted(
                        _PYSCF_SCF_SELECTORS
                        + (
                            "vibrational_frequencies",
                            "vibrational_mode_atom_participation",
                            "vibrational_mode_degeneracy_group",
                        )
                    )
                ),
            ),
            ("irc", _PYSCF_IRC_SELECTORS),
            ("opt", _PYSCF_OPT_SELECTORS),
            ("sp", _PYSCF_SP_SELECTORS),
            # The response stage is executable (contract v5): a td result
            # is one structure -- the supplied geometry, which the validator
            # holds it to -- carrying the reference's mean-field properties
            # beside its roots, so the SCF set is declared with the
            # excitation set and the provenance axis says whose each is.
            ("td", _PYSCF_TD_JOBTYPE_SELECTORS),
            ("ts", _PYSCF_TS_SELECTORS),
        ),
        selector_structural_states=_PYSCF_STRUCTURAL_STATES,
        selector_electronic_provenance=_PYSCF_ELECTRONIC_PROVENANCE,
        resolve_electronic_provenance=_resolve_computed_surface,
        resolve_level=_pyscf_level,
        resolve_surface=lambda output: getattr(output, "surface", None),
        resolve_reference_diagnostics=_pyscf_reference_diagnostics,
        resolve_stationarity_gradient=_pyscf_stationarity_gradient,
        resolve_cartesian_hessian=_pyscf_cartesian_hessian,
        admit_for_analysis=_pyscf_admit_for_analysis,
    ),
    "xyz": ResultReaderV1(
        program="xyz",
        artifact_kind="geometry_xyz",
        parser_id="chemsmart.io.xyz.xyzfile.XYZFile",
        open_output=_xyz_output,
        accessors=_xyz_accessors(),
        requires_normal_termination=False,
    ),
}


#: What the structure a result's spectrum belongs to is, in the host's own
#: word: minimum, first-order saddle, ... or not a stationary point.
STATIONARY_POINT_KIND_SELECTOR = "stationary_point_kind"


def _with_stationary_point_kind(reader: ResultReaderV1) -> ResultReaderV1:
    """Serve the host's stationary-point word wherever the modes are served.

    "Is it a minimum?" was the most asked categorical question in the
    archive and had no word to be answered with: sessions counted modes by
    expression or wrote a validation rule of their own and delivered its
    0/1 verdict, each choosing a convention (R10 Q23). The word is the
    judgement the host already makes (``terminal_states
    .stationary_point_kind``: the stationary-point rule's -20 cm^-1
    convention, and ``structure_stationarity`` -- the one function the
    characterisation and a free energy ask whether a structure is
    stationary), read from the same printed modes, beside
    ``vibrational_frequencies`` on every job type that declares it and for
    the structure those modes belong to.
    """

    frequencies = "vibrational_frequencies"
    if frequencies not in reader.accessors:
        return reader

    def kind(output: Any) -> str | None:
        from chemsmart.agent.terminal_states import stationary_point_kind
        from chemsmart.analysis.result_quantities import (
            structure_stationarity,
        )

        return stationary_point_kind(
            tuple(getattr(output, "vibrational_frequencies", None) or ()),
            structure_stationarity(reader.program, output).stationarity,
        )

    selector = STATIONARY_POINT_KIND_SELECTOR
    state = reader.structural_state(frequencies)
    provenance = reader.electronic_provenance(frequencies)
    return replace(
        reader,
        accessors={**reader.accessors, selector: kind},
        jobtype_selectors=tuple(
            (
                jobtype,
                (
                    tuple(sorted({*selectors, selector}))
                    if frequencies in selectors
                    else selectors
                ),
            )
            for jobtype, selectors in reader.jobtype_selectors
        ),
        selector_structural_states=(
            tuple(
                sorted((*reader.selector_structural_states, (selector, state)))
            )
            if state != "stateless"
            else reader.selector_structural_states
        ),
        selector_electronic_provenance=(
            tuple(
                sorted(
                    (
                        *reader.selector_electronic_provenance,
                        (selector, provenance),
                    )
                )
            )
            if provenance != "stateless"
            else reader.selector_electronic_provenance
        ),
        selector_declarations=(
            *reader.selector_declarations,
            (selector, "", "DIMENSIONLESS"),
        ),
    )


for _program, _reader in tuple(RESULT_READERS.items()):
    RESULT_READERS[_program] = _with_stationary_point_kind(_reader)
del _program, _reader


#: Physical dimension of each selector, in the shared quantity vocabulary.
_SELECTOR_DIMENSIONS = {
    "surface_id": "DIMENSIONLESS",
    "mulliken_atomic_spin_populations": "DIMENSIONLESS",
    "loewdin_atomic_spin_populations": "DIMENSIONLESS",
    "functional": "DIMENSIONLESS",
    "method": "DIMENSIONLESS",
    "ab_initio": "DIMENSIONLESS",
    "basis": "DIMENSIONLESS",
    "converged": "DIMENSIONLESS",
    "irc_converged": "DIMENSIONLESS",
    "absorption_wavelengths": "LENGTH",
    "energy": "ENERGY",
    "energies": "ENERGY",
    "scan_energies": "ENERGY",
    "scan_coordinate_values": "DIMENSIONLESS",
    "scan_point_indices": "DIMENSIONLESS",
    "scan_steps_reached": "DIMENSIONLESS",
    "scan_steps_planned": "DIMENSIONLESS",
    "entropy_times_temperature": "ENERGY",
    "excitation_energies": "ENERGY",
    "singlet_excitation_energies": "ENERGY",
    "triplet_excitation_energies": "ENERGY",
    "excited_state_indices": "DIMENSIONLESS",
    "excited_state_manifold_roots": "DIMENSIONLESS",
    "excited_state_multiplicities": "DIMENSIONLESS",
    "excited_state_labels": "DIMENSIONLESS",
    "excited_state_spin_square": "DIMENSIONLESS",
    "gibbs_free_energy": "ENERGY",
    "oscillator_strengths": "DIMENSIONLESS",
    "singlet_oscillator_strengths": "DIMENSIONLESS",
    "triplet_oscillator_strengths": "DIMENSIONLESS",
    "transition_dipole_moments": "DIPOLE_MOMENT",
    "excited_state_converged": "DIMENSIONLESS",
    "excited_state_followed_root": "DIMENSIONLESS",
    "ccsd_correlation_energy": "ENERGY",
    "triples_correction": "ENERGY",
    "vibrational_frequencies": "FREQUENCY",
    "ir_intensities": "IR_INTENSITY",
    "vibrational_mode_atom_participation": "DIMENSIONLESS",
    "vibrational_mode_degeneracy_group": "DIMENSIONLESS",
    "vpt2_harmonic_frequencies": "FREQUENCY",
    "vpt2_fundamental_frequencies": "FREQUENCY",
    "vpt2_zero_point_rovibrational_energy": "FREQUENCY",
    "positions": "LENGTH",
    "reached_positions": "LENGTH",
    "supplied_positions": "LENGTH",
    "connectivity": "DIMENSIONLESS",
    "symbols": "DIMENSIONLESS",
    "charge": "DIMENSIONLESS",
    "multiplicity": "DIMENSIONLESS",
    "scf_energy": "ENERGY",
    "reference_energy": "ENERGY",
    "correlation_energy": "ENERGY",
    "dispersion_energy": "ENERGY",
    "auxiliary_basis": "DIMENSIONLESS",
    "auxiliary_basis_role": "DIMENSIONLESS",
    "dipole_moment": "DIPOLE_MOMENT",
    "dipole_moment_magnitude": "DIPOLE_MOMENT",
    "homo": "ENERGY",
    "lumo": "ENERGY",
    "gap": "ENERGY",
    "alpha_homo": "ENERGY",
    "alpha_lumo": "ENERGY",
    "beta_homo": "ENERGY",
    "beta_lumo": "ENERGY",
    "spin_square": "DIMENSIONLESS",
    "spin_square_after_annihilation": "DIMENSIONLESS",
    "spin_square_target": "DIMENSIONLESS",
    "spin_square_deviation": "DIMENSIONLESS",
    "effective_multiplicity": "DIMENSIONLESS",
    "wavefunction_stability_verdict": "DIMENSIONLESS",
    "wavefunction_stability_history": "DIMENSIONLESS",
    "trajectory_frame_count": "DIMENSIONLESS",
    "trajectory_energies": "ENERGY",
    "trajectory_start_frequencies": "FREQUENCY",
    "trajectory_start_positions": "LENGTH",
    "trajectory_end_positions": "LENGTH",
    "trajectory_start_connectivity": "DIMENSIONLESS",
    "trajectory_end_connectivity": "DIMENSIONLESS",
    "trajectory_connectivity_changed": "DIMENSIONLESS",
    "irc_direction": "DIMENSIONLESS",
    "solvation_model": "DIMENSIONLESS",
    "solvation_electrostatic_energy": "ENERGY",
    "solvation_nonelectrostatic_energy": "ENERGY",
    "solvation_cavity_surface_area": "AREA",
    "mulliken_atomic_charges": "CHARGE",
    "hirshfeld_atomic_charges": "CHARGE",
    "loewdin_atomic_charges": "CHARGE",
    "xtb_scc_atomic_charges": "CHARGE",
    "solvent": "DIMENSIONLESS",
    "wiberg_bond_orders": "DIMENSIONLESS",
}


def merge_selector_declarations(
    readers: Mapping[str, ResultReaderV1],
    *,
    units: dict[str, str],
    dimensions: dict[str, str],
    atom_metadata: dict[str, Mapping[str, str]],
) -> frozenset[str]:
    """Fold what each reader introduces into the shared selector vocabulary.

    A new selector used to be four edits in three flat tables that every
    program shares -- unit, dimension, atom metadata and the request
    gate's set -- so two programs gaining a selector each at the same
    time edited the same lines.  A reader now declares what it introduces
    beside the accessor that reads it, and this is the one place the
    declarations meet the shared tables.

    It refuses rather than overwrites: a selector name carries one unit
    and one dimension on every program, which the flat tables enforced
    only by having one line per name.
    """

    declared: set[str] = set()
    for program, reader in readers.items():
        for selector, unit, dimension in reader.selector_declarations:
            for table, value, word in (
                (units, unit, "unit"),
                (dimensions, dimension, "dimension"),
            ):
                known = table.get(selector)
                if known is not None and known != value:
                    raise ValueError(
                        f"{program} declares {word} {value!r} for selector "
                        f"{selector!r}, which the shared vocabulary already "
                        f"gives as {known!r}; one selector name has one "
                        f"{word} on every program"
                    )
                table[selector] = value
            declared.add(selector)
        for selector, pairs in reader.atom_resolved_declarations:
            record = dict(pairs)
            known_record = atom_metadata.get(selector)
            if known_record is not None and dict(known_record) != record:
                raise ValueError(
                    f"{program} declares atom-resolved metadata for "
                    f"{selector!r} that disagrees with the shared record"
                )
            atom_metadata[selector] = record
    return frozenset(declared)


#: Selector names the readers introduced themselves.  The request gate
#: admits them beside the shared set (``result_quantities.supported_selectors``).
DECLARED_SELECTORS = merge_selector_declarations(
    RESULT_READERS,
    units=SELECTOR_UNITS,
    dimensions=_SELECTOR_DIMENSIONS,
    atom_metadata=_ATOM_RESOLVED_SELECTOR_METADATA,  # type: ignore[arg-type]
)

_TEXT_SELECTORS = frozenset(
    {
        "irc_direction",
        "functional",
        "method",
        "ab_initio",
        "basis",
        "auxiliary_basis",
        "auxiliary_basis_role",
        "solvation_model",
        "solvent",
        "wavefunction_stability_verdict",
        "surface_id",
    }
)
_TEXT_VECTOR_SELECTORS = frozenset(
    {
        "excited_state_dominant_excitations",
        "excited_state_labels",
        "wavefunction_stability_history",
    }
)
_INTEGER_SELECTORS = frozenset(
    {
        "charge",
        "converged",
        "excited_state_followed_root",
        "irc_converged",
        "multiplicity",
        "trajectory_frame_count",
        "trajectory_connectivity_changed",
    }
)


def _derived_adjacency(reader: Any, output: Any) -> Any:
    """Bond list bundled onto a delivered positions read.

    Reading relations off raw Cartesian coordinates is the operation
    language models measurably fail at, so the host bundles the
    adjacency it already computes -- the same covalent-radius perception
    the ``connectivity`` selector serves -- whenever positions are
    actually delivered from a result. The hard line, stated once:
    derived adjacency is a measurement and is allowed; any perceived
    label -- isomer, conformer, species, ring class, stereo
    descriptor -- is the scientist's judgement and must never be
    attached here. A reader that serves no connectivity simply bundles
    nothing.
    """

    try:
        symbols, _unit = reader.read(output, "symbols")
        connectivity, _unit = reader.read(output, "connectivity")
    except Exception:
        return ()
    bonds = tuple(
        (int(row_index), int(column_index))
        for row_index, row in enumerate(connectivity)
        for column_index, bonded in enumerate(row)
        if bonded and column_index > row_index
    )
    counts: dict[str, int] = {}
    for symbol in symbols:
        counts[str(symbol)] = counts.get(str(symbol), 0) + 1
    formula = "".join(
        f"{symbol}{count if count > 1 else ''}"
        for symbol, count in sorted(counts.items())
    )
    if not formula:
        return ()
    delivered = {
        "formula": formula,
        "bond_atom_pairs": bonds,
        # Which convention decided, and by how much. A bond list alone is
        # a boolean per pair, and a boolean from a threshold cannot be
        # told from a structural fact: a converged formaldehyde lost both
        # C-H bonds by 1.5 mA and the reader had no way to see that. The
        # policy id makes the convention nameable; the margins make the
        # near-threshold decisions legible. Neither asserts chemistry --
        # the host says what its rule did and how narrowly, and the
        # scientist draws the conclusion.
        "adjacency_policy_id": BOND_PERCEPTION_POLICY_ID,
    }
    try:
        positions, _unit = reader.read(output, "positions")
    except Exception:
        return delivered
    try:
        pairs = perceive_pairs(symbols, positions, include_rejected=True)
    except Exception:
        return delivered
    rows = []
    for pair in pairs:
        near_miss = not pair.adjacent and abs(pair.relative_margin) <= 0.10
        if not pair.adjacent and not near_miss:
            continue
        rows.append(
            {
                "atoms": (pair.first_index, pair.second_index),
                "distance_angstrom": round(pair.distance_angstrom, 6),
                "cutoff_angstrom": round(pair.cutoff_angstrom, 6),
                "margin_angstrom": round(pair.margin_angstrom, 6),
                "adjacent": bool(pair.adjacent),
            }
        )
    if rows:
        delivered["adjacency_margins"] = tuple(rows)
    return delivered


def _verified_native_evidence(
    *, reader: ResultReaderV1, output: Any, artifact: Path, selectors: Any
) -> tuple[tuple[tuple[str, str, str], ...], tuple[tuple[Path, str], ...]]:
    """Seal the sidecars a requested selector will read.

    The primary output is already a trusted artifact.  A parser sidecar is
    trustworthy only when it is a regular, non-symlink file under that
    result's directory and its digest survives the extraction.  Keep the
    selector association because a filename alone has no scientific meaning.
    """

    from chemsmart.analysis import result_quantities as rq

    result_root = artifact.parent.resolve()
    records: set[tuple[str, str, str]] = set()
    observed: dict[Path, str] = {}
    for requested in selectors:
        selector = requested.selector
        for supplied_path in reader.native_evidence_paths_for_output(
            output, selector
        ):
            candidate = Path(supplied_path).expanduser()
            if not candidate.is_absolute():
                candidate = artifact.parent / candidate
            if candidate.is_symlink():
                raise rq.QuantityExtractionError(
                    f"native evidence sidecar for {selector!r} is a symlink"
                )
            resolved = candidate.resolve()
            try:
                filename = resolved.relative_to(result_root).as_posix()
            except ValueError as exc:
                raise rq.QuantityExtractionError(
                    f"native evidence sidecar for {selector!r} lies outside "
                    "the verified result directory"
                ) from exc
            if not resolved.is_file():
                raise rq.QuantityExtractionError(
                    f"native evidence sidecar for {selector!r} is not a "
                    "regular file"
                )
            digest = rq.result_file_sha256(resolved)
            prior = observed.setdefault(resolved, digest)
            if prior != digest:
                raise rq.QuantityExtractionError(
                    "one native evidence file resolved to inconsistent "
                    "digests before extraction"
                )
            records.add((selector, filename, digest))
    return tuple(sorted(records)), tuple(
        sorted(observed.items(), key=lambda item: str(item[0]))
    )


def extract_logged_quantities(
    *,
    request: Any,
    artifact_path: str | Path,
) -> Any:
    """Extract selected quantities from a log-parsing program's result.

    Mirrors the structured-PySCF path: the artifact bytes are verified before
    and after parsing, so a file that changes mid-extraction is refused rather
    than silently mixing two runs.
    """

    from chemsmart.analysis import result_quantities as rq

    reader = reader_for(request.program)
    if reader is None:
        raise rq.QuantityContractError(
            f"no result reader is registered for {request.program!r}; "
            f"registered: {list(registered_reader_programs())}"
        )
    artifact = rq._verify_artifact(artifact_path, request.artifact_sha256)
    output = reader.open_output(artifact)
    if (
        reader.requires_normal_termination
        and getattr(output, "normal_termination", None) is not True
    ):
        raise rq.QuantityExtractionError(
            f"{request.program} scientific quantities require a normally "
            "terminated program result"
        )
    if reader.admit_for_analysis is not None:
        # A reader that opens failed results for inspection keeps its
        # green admission for quantities here.
        reader.admit_for_analysis(artifact)
    evidence_ref = f"artifact:{request.artifact_id}#{request.artifact_sha256}"
    # The declaration gate.  A selector declared for a jobtype is a semantic
    # claim; until this check the registry was consulted only by the pre-plan
    # capability query, so an undeclared selector still extracted -- which is
    # how an IRC log's starting structure was delivered labeled as both path
    # endpoints, and how a scan's optimizer trace could stand in for its
    # surface.  A reader that declares no jobtype coverage at all (the xyz
    # reader over registered geometry artifacts) stays ungated: its values
    # carry no jobtype semantics to misread.
    if reader.jobtype_selectors:
        jobtype = str(getattr(output, "jobtype", "") or "").casefold()
        declared = reader.selectors_for_jobtype(jobtype)
        if declared is None:
            raise rq.QuantityExtractionError(
                f"{request.program} declares no selector coverage for "
                f"jobtype {jobtype or 'unknown'!r}; declared jobtypes: "
                f"{sorted(j for j, _ in reader.jobtype_selectors)}. What an "
                "undeclared jobtype's printed values mean has not been "
                "audited, so nothing is extracted from it."
            )
        undeclared = sorted(
            {
                item.selector
                for item in request.selectors
                if item.selector not in declared
            }
        )
        if undeclared:
            raise rq.QuantityExtractionError(
                f"selector(s) {undeclared} are not declared for "
                f"{request.program} jobtype {jobtype!r}; a declaration is a "
                "semantic claim about what the value means for this job "
                f"type. Declared here: {sorted(declared)}."
                + rq.thermochemistry_route_hint(undeclared)
            )
    native_evidence, native_evidence_digests = _verified_native_evidence(
        reader=reader,
        output=output,
        artifact=artifact,
        selectors=request.selectors,
    )
    quantities = []
    absent: list[tuple[str, str, str]] = []
    # Whose density or method each delivered value belongs to, resolved
    # against this artifact and carried on the receipt beside the value,
    # so a ground-state dipole read off an excited-root result says so
    # where the number is cited.  Present on the body only when a reader
    # declares the axis, so every receipt minted before it verifies.
    provenance: list[tuple[str, str]] = []
    # A model-selected quantity id has no scientific semantics by itself.
    # Keep its host selector and structural role in the digest-bearing record,
    # including explicit absences, so a later cycle need not recover them from
    # an ephemeral tool event.
    bindings: list[tuple[str, str]] = []
    states: list[tuple[str, str]] = []
    positions_delivered = False
    for selector in request.selectors:
        bindings.append((selector.quantity_id, selector.selector))
        states.append(
            (selector.quantity_id, reader.structural_state(selector.selector))
        )
        try:
            source_value, source_unit = reader.read(output, selector.selector)
        except MissingQuantityError as exc:
            # A quantity this run never produced is a stated gap, so it is
            # recorded rather than raised: one absence used to end the whole
            # extraction, and every consumer of every sibling quantity
            # skipped with it -- including, in the run that prompted this,
            # twelve nodes that named no absent value at all.  The reason
            # travels with the absence, and the artifact's own inventory
            # travels with it too, so a reader learns the shape of the
            # result here instead of one refused selector at a time.
            #
            # Only this error type collects.  QuantityExtractionError below
            # means the reader and the writer disagree about what a stored
            # value is, which is a defect rather than a gap, and it still
            # ends the extraction.
            available = reader.available_selectors(output)
            absent.append(
                (
                    selector.quantity_id,
                    selector.selector,
                    f"{exc}. This result resolves: {', '.join(available)}",
                )
            )
            continue
        dimension = getattr(rq, _SELECTOR_DIMENSIONS[selector.selector])
        if selector.selector in {"symbols", *_TEXT_VECTOR_SELECTORS}:
            value = source_value
            unit = source_unit
            data_kind = "text_vector"
        elif selector.selector in _TEXT_SELECTORS or isinstance(
            source_value, str
        ):
            # A word the reader returned is a word, whichever list
            # remembered to say so. PySCF's three stability words were
            # declared dimensionless and never listed above, so every
            # extraction of them died converting 'stable' to a float
            # (r9 g2-stability, Slurm 2145043).
            value = source_value
            unit = source_unit
            data_kind = "text"
        elif (
            isinstance(source_value, (list, tuple))
            and source_value
            and all(isinstance(item, str) for item in source_value)
        ):
            value = source_value
            unit = source_unit
            data_kind = "text_vector"
        elif selector.selector in _INTEGER_SELECTORS:
            value = int(source_value)
            unit = "1"
            data_kind = "integer"
        elif selector.selector in {"wiberg_bond_orders", "mayer_bond_orders"}:
            # Sparse bond-order rows, [atom_i, atom_j, order], either
            # scheme: the pair is an index, the order a pure number.
            value = tuple(
                (int(r[0]), int(r[1]), float(r[2])) for r in source_value
            )
            unit = "1"
            data_kind = "matrix"
        else:
            from chemsmart.analysis.quantity_expressions import (
                normalize_numeric_value,
            )

            value, unit, observed_dimension = normalize_numeric_value(
                source_value, source_unit
            )
            if observed_dimension != dimension:
                raise rq.QuantityExtractionError(
                    f"selector {selector.selector!r} produced an "
                    "incompatible unit"
                )
            data_kind = None
        quantities.append(
            rq.make_quantity_value(
                quantity_id=selector.quantity_id,
                source_value=source_value,
                source_unit=source_unit,
                value=value,
                unit=unit,
                dimension=dimension,
                evidence_ref=evidence_ref,
                data_kind=data_kind,
            )
        )
        if selector.selector == "positions":
            positions_delivered = True
        word = reader.electronic_provenance_for_output(
            output, selector.selector
        )
        if word != "stateless":
            provenance.append((selector.quantity_id, word))
    if rq.result_file_sha256(artifact) != request.artifact_sha256:
        raise rq.QuantityExtractionError(
            "result artifact changed during extraction"
        )
    for sidecar, expected_sha256 in native_evidence_digests:
        if rq.result_file_sha256(sidecar) != expected_sha256:
            raise rq.QuantityExtractionError(
                "native evidence sidecar changed during extraction"
            )
    body = rq.canonical_extraction_receipt_body(
        schema_version="chemsmart.quantity-extraction-receipt.v1",
        artifact_id=request.artifact_id,
        artifact_sha256=request.artifact_sha256,
        program=request.program,
        parser_id=reader.parser_id,
        quantities=tuple(quantities),
        status="partial" if absent else "extracted",
        absent=tuple(absent),
        derived_adjacency=(
            _derived_adjacency(reader, output) if positions_delivered else ()
        ),
        electronic_provenance=tuple(provenance),
        selector_bindings=tuple(bindings),
        structural_states=tuple(states),
        level=reader.level_for_output(output),
        native_evidence=native_evidence,
    )
    return rq.QuantityExtractionReceiptV1(
        **body, receipt_sha256=rq.canonical_quantity_sha256(body)
    )


def selector_dimension(selector: str):
    """The fixed dimension a selector's extracted value carries, or None.

    Exposed for the plan-time gate: an extraction output's declared
    unit is checkable against this table before any engine runs. A
    selector absent from the table yields None and the caller skips.
    """

    from chemsmart.analysis import result_quantities as rq

    name = _SELECTOR_DIMENSIONS.get(str(selector).strip())
    return getattr(rq, name) if name else None


def reader_for(program: str) -> ResultReaderV1 | None:
    """Return the registered reader for ``program``, or ``None``."""

    return RESULT_READERS.get(str(program).strip().lower())


def registered_reader_programs() -> tuple[str, ...]:
    """Return every program with a registered log reader, sorted."""

    return tuple(sorted(RESULT_READERS))


def registered_reader_jobtype_selectors(
    program: str, jobtype: str
) -> tuple[str, ...] | None:
    """Return job-scoped parser support, conditional on engine emission."""

    reader = reader_for(program)
    if reader is None:
        return None
    return reader.selectors_for_jobtype(jobtype)


def registered_reader_selectors() -> dict[str, tuple[str, ...]]:
    """Return the model-discoverable selector inventory for each reader."""

    return {
        program: tuple(sorted(reader.selectors))
        for program, reader in sorted(RESULT_READERS.items())
    }
