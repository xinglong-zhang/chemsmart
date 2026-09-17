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

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping

from chemsmart.io.molecules.perception import (
    BOND_PERCEPTION_POLICY_ID,
    molecule_adjacency_matrix,
    perceive_pairs,
)

__all__ = [
    "RESULT_READERS",
    "MissingQuantityError",
    "ResultReaderV1",
    "SELECTOR_UNITS",
    "reader_for",
    "registered_reader_jobtype_selectors",
    "registered_reader_programs",
    "registered_reader_selectors",
    "atom_resolved_selector_metadata",
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
    is no spin to partition; checked against 2S = multiplicity - 1, the
    sum ORCA prints beneath the block, so a dropped or duplicated atom
    cannot pass as a population.
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
    if multiplicity is not None:
        expected = float(multiplicity - 1)
        total = sum(values)
        if abs(total - expected) > 0.05:
            raise MissingQuantityError(
                f"{quantity} sums to {total:.3f} where 2S for multiplicity "
                f"{multiplicity} is {expected:.1f}; the vector is not the "
                "complete molecule in order"
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


def _orca_excitation_energies(output: Any) -> list[float]:
    records = list(output.excited_state_records or ())
    if records:
        return [float(item["energy_eV"]) for item in records]
    # ORCA 6 spectrum-only fragments remain useful for the legacy aggregate
    # selector, but they cannot support multiplicity-specific selectors.
    return [float(item) for item in output.excitation_energies_eV]


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
    #: Program-native units that differ from the shared selector display unit.
    source_units: Mapping[str, str] = field(default_factory=dict)
    #: Selectors this parser can extract from a program job type when the
    #: chosen method/settings emit them. Missing coverage means unknown, never
    #: that the job produces no quantities.
    jobtype_selectors: tuple[tuple[str, tuple[str, ...]], ...] = ()
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
    #: level.  A level is shown and never compared: whether two programs
    #: mean one thing by a keyword is a fact about the programs.
    resolve_level: Callable[[Any], Mapping[str, Any]] | None = None
    #: The electronic surface one opened result is on, as a mapping over
    #: ``SURFACE_IDENTITY_FIELDS``.  Unlike the level, a surface is *for*
    #: comparing: two results are one surface when every field agrees and
    #: neither reader wrote ``unknown``.  None means this reader cannot
    #: say, and the organs that ask treat that as "not comparable" rather
    #: than as agreement.
    resolve_surface: Callable[[Any], Mapping[str, Any] | None] | None = None

    def surface_for_output(self, output: Any) -> Mapping[str, Any] | None:
        """The surface this result is on, or None when unknowable."""

        if self.resolve_surface is None:
            return None
        return self.resolve_surface(output)

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
        return values[-1]
    return _last_energy(output)


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


def _orca_accessors() -> dict[str, Callable[[Any], Any]]:
    accessors = _text_output_accessors(mode_composition=True)
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
            # Reached versus planned is the whole diagnosis when a scan
            # dies partway. The parser has always known both -- the step
            # announcements it counted, and the step total ORCA stated
            # before the first point -- and threw them away: a truncated
            # scan reached the typed layer only as shorter vectors, with
            # nothing to say what the target was. A live scan died at
            # step 2 of 12 and no tool could state either number.
            "scan_steps_reached": _scan_steps_reached,
            "scan_steps_planned": _scan_steps_planned,
            "absorption_wavelengths": lambda output: [
                float(item) for item in output.absorption_wavelengths
            ],
            "excitation_energies": _orca_excitation_energies,
            "oscillator_strengths": lambda output: [
                float(item) for item in output.oscillator_strengths
            ],
            "excited_state_indices": lambda output: [
                int(item["state_index"])
                for item in output.excited_state_records
            ],
            "excited_state_manifold_roots": lambda output: [
                int(item["manifold_root"])
                for item in output.excited_state_records
            ],
            "excited_state_multiplicities": lambda output: [
                int(item)
                for item in _required_record_values(
                    output.excited_state_records,
                    "multiplicity",
                    "a spin multiplicity",
                )
            ],
            "excited_state_spin_square": lambda output: [
                float(item["spin_square"])
                for item in output.excited_state_records
            ],
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
                for item in output.electronic_absorption_transition_records
                if item["multiplicity"] == 1
            ],
            "triplet_oscillator_strengths": lambda output: [
                float(item["oscillator_strength"])
                for item in output.electronic_absorption_transition_records
                if item["multiplicity"] == 3
            ],
            "energy": _orca_total_energy,
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
            "functional": _route_functional,
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


def _gaussian_accessors() -> dict[str, Callable[[Any], Any]]:
    accessors = _text_output_accessors()
    accessors.update(
        {
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
            "excited_state_spin_square": lambda output: [
                float(item)
                for item in _required_record_values(
                    output.excited_state_records,
                    "spin_square",
                    "a printed excited-state <S^2>",
                )
            ],
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
            "functional": _route_functional,
            "homo": _gaussian_frontier("homo_energy"),
            "lumo": _gaussian_frontier("lumo_energy"),
            "gap": _gaussian_frontier("fmo_gap"),
            "ab_initio": _route_ab_initio,
            "basis": _route_basis,
        }
    )
    return accessors


def _xtb_output(path: Path) -> Any:
    from chemsmart.io.xtb.output import XTBOutput

    # Unlike the Gaussian and ORCA readers, XTBOutput represents the complete
    # calculation directory because xTB can distribute one result across the
    # main log, geometry, Hessian, and vibrational-spectrum files.  The host
    # has already resolved and verified the exact main output artifact; its
    # parent is therefore the corresponding calculation directory.
    return XTBOutput(folder=str(path.parent))


def _xtb_state_integer(attribute: str) -> Callable[[Any], int]:
    def _read(output: Any) -> int:
        value = getattr(output, attribute, None)
        if value is None:
            raise MissingQuantityError(
                f"this xtb result records no {attribute}"
            )
        return int(value)

    return _read


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


def _xtb_geometry_source_path(output: Any, selector: str) -> Path | None:
    """Name the xTB sidecar only for the reached-optimisation selector."""

    if selector != "reached_positions":
        return None
    geometry_file = getattr(output, "xtbopt_geometry_file", None)
    path = getattr(geometry_file, "filepath", None)
    return Path(str(path)) if path else None


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
    absent = [path for path in expected if observed.get(path) is None]
    if absent:
        raise MissingQuantityError(
            f"pyscf result contains no {selector!r} value "
            f"(datasets absent: {sorted(absent)})"
        )
    wrong = {
        path: {"expected": unit, "observed": observed.get(path)}
        for path, unit in expected.items()
        if observed.get(path) != unit
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
    """The functional the project asked for, as ``spec/method`` names it.

    The literal libxc ran (``spec/xc``) is the same functional under a
    different spelling wherever the alias table rewrote it -- ``b3lyp``
    and ``b3lypg`` are one libxc code -- so the requested name is what a
    program-neutral identity means, as for ORCA and Gaussian; the applied
    materialisation rides the review as the functional-resolution receipt.
    """

    if not output.spec.get("xc"):
        raise MissingQuantityError(
            "this result ran a Hartree-Fock reference and names no "
            "functional; the method identity is read by 'ab_initio'"
        )
    value = output.method
    if not value:
        raise MissingQuantityError("pyscf result records no method name")
    return str(value)


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
        # initio method was named; the name is the project's, not libxc's.
        level["functional"] = method
    if spec.get("basis") not in (None, ""):
        level["basis"] = spec["basis"]
    if getattr(output, "solvent_on", False):
        level["solvent_model"] = output.solvent_model
        level["solvent"] = output.solvent_id
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
    record = getattr(output, "excited_state_record", None)
    if isinstance(record, Mapping) and record.get("root") is not None:
        level["excited_state_root"] = int(record["root"])
        for name in ("response_method", "state_manifold", "nstates"):
            if record.get(name) is not None:
                level.setdefault(name, record[name])
    return level


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
        "homo": lambda output: float(output.homo_energy),
        "lumo": lambda output: float(output.lumo_energy),
        "gap": lambda output: float(output.fmo_gap),
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
        "excited_state_converged": _pyscf_excited_converged,
        "excited_state_followed_root": _pyscf_followed_root,
        # The correlated stage: the program's own components at the final
        # geometry.  ``correlation_energy`` is the final method's whole
        # correlation, triples included, as the ORCA reader means it.
        "reference_energy": _pyscf_correlated_scalar("reference_energy"),
        "correlation_energy": _pyscf_correlated_scalar("correlation_energy"),
        "ccsd_correlation_energy": _pyscf_correlated_scalar(
            "ccsd_correlation_energy"
        ),
        "triples_correction": _pyscf_correlated_scalar("triples_correction"),
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
        def accessor(output: Any) -> Any:
            _pyscf_require_units(selector, output)
            return read(output)

        return accessor

    return {name: _guard(name, read) for name, read in raw.items()}


#: Selectors every executed PySCF stage writes.  The SCF block that follows
#: every stage stores energies, geometry, orbital energies and the population
#: and dipole properties, so a single point and an optimisation answer the
#: same set; a Hessian stage adds the vibrational quantities on top.
_PYSCF_SCF_SELECTORS = (
    "ab_initio",
    "surface_id",
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
    "lumo",
    "method",
    "mulliken_atomic_charges",
    "mulliken_atomic_spin_populations",
    "multiplicity",
    "positions",
    "scf_energy",
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
_PYSCF_TD_JOBTYPE_SELECTORS = tuple(
    sorted(_PYSCF_SCF_SELECTORS + _PYSCF_TD_SELECTORS)
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
            ("dipole_moment", "as_reached"),
            ("dipole_moment_magnitude", "as_reached"),
            ("effective_multiplicity", "as_reached"),
            ("energy", "as_reached"),
            ("excitation_energies", "as_reached"),
            ("excited_state_converged", "as_reached"),
            ("excited_state_indices", "as_reached"),
            ("excited_state_manifold_roots", "as_reached"),
            ("excited_state_multiplicities", "as_reached"),
            ("gap", "as_reached"),
            ("homo", "as_reached"),
            ("lumo", "as_reached"),
            ("mulliken_atomic_charges", "as_reached"),
            ("mulliken_atomic_spin_populations", "as_reached"),
            ("oscillator_strengths", "as_reached"),
            ("positions", "as_reached"),
            ("reached_positions", "as_reached"),
            ("reference_energy", "as_reached"),
            ("scf_energy", "as_reached"),
            ("singlet_excitation_energies", "as_reached"),
            ("singlet_oscillator_strengths", "as_reached"),
            ("spin_square", "as_reached"),
            ("spin_square_deviation", "as_reached"),
            ("spin_square_target", "as_reached"),
            ("supplied_positions", "as_supplied"),
            ("surface_id", "stateless"),
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
            ("ccsd_correlation_energy", "correlated"),
            ("correlation_energy", "correlated"),
            ("dipole_moment", "reference"),
            ("dipole_moment_magnitude", "reference"),
            ("effective_multiplicity", "reference"),
            ("energies", "reference"),
            ("energy", "computed_surface"),
            ("excitation_energies", "excited_root"),
            ("excited_state_converged", "excited_root"),
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
            ("surface_id", "stateless"),
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
    ("excited_state_indices", "excited_root"),
    ("excited_state_labels", "excited_root"),
    ("excited_state_manifold_roots", "excited_root"),
    ("excited_state_multiplicities", "excited_root"),
    ("excited_state_spin_square", "excited_root"),
    ("gap", "reference"),
    ("gibbs_free_energy", "computed_surface"),
    ("hirshfeld_atomic_charges", "reference"),
    ("homo", "reference"),
    ("loewdin_atomic_charges", "reference"),
    ("loewdin_atomic_spin_populations", "reference"),
    ("lumo", "reference"),
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
            ("trajectory_connectivity_changed", "trajectory_endpoint"),
            ("trajectory_end_connectivity", "trajectory_endpoint"),
            ("trajectory_end_positions", "trajectory_endpoint"),
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
                    "homo",
                    "loewdin_atomic_charges",
                    "loewdin_atomic_spin_populations",
                    "lumo",
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
                # ORCA writes the reaction path to XYZ sidecars rather than
                # into the log, so the ``trajectory_*`` family is deliberately
                # absent here: an IRC log parses to a single structure, and
                # the path is read from the registered ``_IRC_Full_trj.xyz``
                # artifact through the ``xyz`` reader.  What the log itself
                # establishes is the endpoint it converged to and the
                # direction the ``%irc`` block explicitly declared.
                "irc",
                # Only job-level facts are declared.  An ORCA IRC log prints
                # a single structure -- the starting point -- so every
                # state-dependent value (geometry, energies, orbitals,
                # dipoles, spin) read from it describes the transition
                # state, not the path: the first Agent-executed IRC
                # delivered the saddle's own distances as both endpoints,
                # and its ``energy`` differed from the true endpoint by the
                # entire barrier.  The path lives in the trajectory sidecar,
                # which enters the typed layer as a registered geometry
                # artifact and is read by the ungated xyz reader; a
                # log-native path route (ORCA's IRC PATH SUMMARY table) is
                # future parser work.
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
                    "homo",
                    "loewdin_atomic_charges",
                    "loewdin_atomic_spin_populations",
                    "lumo",
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
                    "homo",
                    "loewdin_atomic_charges",
                    "loewdin_atomic_spin_populations",
                    "lumo",
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
                    "homo",
                    "loewdin_atomic_charges",
                    "loewdin_atomic_spin_populations",
                    "lumo",
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
                "opt",
                (
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
                    "gibbs_free_energy",
                    "homo",
                    "lumo",
                    "multiplicity",
                    "positions",
                    "spin_square",
                    "spin_square_after_annihilation",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "vibrational_frequencies",
                    "wavefunction_stability_history",
                    "wavefunction_stability_verdict",
                ),
            ),
            (
                "sp",
                (
                    "ab_initio",
                    "absorption_wavelengths",
                    "basis",
                    "charge",
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "effective_multiplicity",
                    "energies",
                    "energy",
                    "excitation_energies",
                    "excited_state_indices",
                    "excited_state_labels",
                    "excited_state_spin_square",
                    "functional",
                    "multiplicity",
                    "oscillator_strengths",
                    "positions",
                    "spin_square",
                    "spin_square_after_annihilation",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "wavefunction_stability_history",
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
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "effective_multiplicity",
                    "energies",
                    "energy",
                    "functional",
                    "gibbs_free_energy",
                    "multiplicity",
                    "positions",
                    "spin_square",
                    "spin_square_after_annihilation",
                    "spin_square_deviation",
                    "spin_square_target",
                    "symbols",
                    "vibrational_frequencies",
                    "wavefunction_stability_history",
                    "wavefunction_stability_verdict",
                ),
            ),
        ),
    ),
    "xtb": ResultReaderV1(
        program="xtb",
        artifact_kind="xtb_output",
        parser_id="chemsmart.io.xtb.output.XTBOutput",
        open_output=_xtb_output,
        accessors=_xtb_accessors(),
        source_units={"dipole_moment": "e bohr"},
        geometry_source_path_for_selector=_xtb_geometry_source_path,
        jobtype_selectors=(
            (
                "hess",
                (
                    "charge",
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "energy",
                    "gap",
                    "gibbs_free_energy",
                    "homo",
                    "lumo",
                    "multiplicity",
                    "positions",
                    "symbols",
                    "vibrational_frequencies",
                    "vibrational_mode_atom_participation",
                    "vibrational_mode_degeneracy_group",
                    "xtb_scc_atomic_charges",
                ),
            ),
            (
                "opt",
                (
                    "charge",
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "energy",
                    "gap",
                    "gibbs_free_energy",
                    "homo",
                    "lumo",
                    "multiplicity",
                    "positions",
                    "reached_positions",
                    "symbols",
                    "xtb_scc_atomic_charges",
                ),
            ),
            (
                "sp",
                (
                    "connectivity",
                    "dipole_moment",
                    "dipole_moment_magnitude",
                    "energy",
                    "gap",
                    "homo",
                    "lumo",
                    "positions",
                    "symbols",
                    "xtb_scc_atomic_charges",
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
                    ("energy", "as_reached"),
                    ("gap", "as_reached"),
                    ("gibbs_free_energy", "as_reached"),
                    ("homo", "as_reached"),
                    ("lumo", "as_reached"),
                    ("multiplicity", "as_reached"),
                    ("positions", "as_reached"),
                    ("reached_positions", "as_reached"),
                    ("symbols", "stateless"),
                    ("vibrational_frequencies", "as_reached"),
                    ("vibrational_mode_atom_participation", "as_reached"),
                    ("vibrational_mode_degeneracy_group", "as_reached"),
                    ("xtb_scc_atomic_charges", "as_reached"),
                ]
            )
        ),
        # xTB's parsed record names its Hamiltonian but does not yet carry a
        # program-neutral electronic-surface identity.  Leave that axis
        # absent rather than labelling every quantity ``computed_surface``
        # without an identity a consumer can resolve.
        resolve_level=_xtb_level,
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
            ("opt", _PYSCF_OPT_SELECTORS),
            ("sp", _PYSCF_SP_SELECTORS),
            # The response stage is executable (contract v5): a td result
            # is one structure -- the supplied geometry, which the validator
            # holds it to -- carrying the reference's mean-field properties
            # beside its roots, so the SCF set is declared with the
            # excitation set and the provenance axis says whose each is.
            ("td", _PYSCF_TD_JOBTYPE_SELECTORS),
        ),
        selector_structural_states=_PYSCF_STRUCTURAL_STATES,
        selector_electronic_provenance=_PYSCF_ELECTRONIC_PROVENANCE,
        resolve_electronic_provenance=_resolve_computed_surface,
        resolve_level=_pyscf_level,
        resolve_surface=lambda output: getattr(output, "surface", None),
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
}

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
    {"excited_state_labels", "wavefunction_stability_history"}
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
        elif selector.selector in _TEXT_SELECTORS:
            value = source_value
            unit = source_unit
            data_kind = "text"
        elif selector.selector in _INTEGER_SELECTORS:
            value = int(source_value)
            unit = "1"
            data_kind = "integer"
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
