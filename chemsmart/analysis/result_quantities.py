"""Typed, hash-bound extraction of scientific quantities from result files.

This module is deliberately smaller than a general file-query language.  A
caller selects from a finite semantic vocabulary and supplies an artifact path
that has already been resolved by the host.  The path is never interpreted as
model-authored input, and the exact bytes are checked before and after parsing.

PySCF is the first registered reader because ChemSmart controls its structured
HDF5 schema.  The contracts are program-neutral so that Gaussian, ORCA, and xTB
readers can be registered later without changing the expression evaluator.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
from dataclasses import asdict, dataclass, is_dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from chemsmart.analysis.thermochemistry import Thermochemistry
from chemsmart.io.pyscf.output import PySCFOutput
from chemsmart.jobs.pyscf.environment import canonical_sha256
from chemsmart.jobs.pyscf.writer import (
    RESULT_UNITS,
    SUPPORTED_RESULT_CONTRACT_VERSIONS,
)
from chemsmart.utils.constants import energy_conversion

# Historical quantities use six bases in the order energy, length,
# temperature, angle, frequency, pressure.  New independent physical
# dimensions may append a component without rewriting old receipts.
Dimension = tuple[int, ...]
DIMENSIONLESS: Dimension = (0, 0, 0, 0, 0, 0)
ENERGY: Dimension = (1, 0, 0, 0, 0, 0)
LENGTH: Dimension = (0, 1, 0, 0, 0, 0)
TEMPERATURE: Dimension = (0, 0, 1, 0, 0, 0)
ANGLE: Dimension = (0, 0, 0, 1, 0, 0)
FREQUENCY: Dimension = (0, 0, 0, 0, 1, 0)
PRESSURE: Dimension = (0, 0, 0, 0, 0, 1)
ENTROPY: Dimension = (1, 0, -1, 0, 0, 0)
# Electric charge is not otherwise a base quantity in the current expression
# vocabulary.  Keep dipole moment independent from length rather than calling a
# Debye an Angstrom or a dimensionless number.
DIPOLE_MOMENT: Dimension = (0, 0, 0, 0, 0, 0, 1)
# Atomic mass is independent from the historical bases and from electric
# dipole moment.  Appending it preserves every existing six- and seven-entry
# dimension tuple while allowing coordinate-derived inertia to remain typed.
MASS: Dimension = (0, 0, 0, 0, 0, 0, 0, 1)
MOMENT_OF_INERTIA: Dimension = (0, 2, 0, 0, 0, 0, 0, 1)

#: Electric charge, the ninth base, in units of the elementary charge.  It
#: exists so that an electrode potential is a dimension rather than a number
#: with a hopeful label: potential is energy per charge, which makes
#: dG = -nFE dimensionally checkable instead of asserted, and dissolves the
#: Faraday constant into the unit system where it belongs -- it is a
#: definition, not a value taken from the literature.
#:
#: DIPOLE_MOMENT above stays its own base rather than becoming charge times
#: length.  Rewriting it would change the dimension recorded in every dipole
#: receipt already written, and a stored receipt's dimension is part of its
#: identity.
#: Area, i.e. length squared.  Named so the selector dimension table can
#: refer to it; it introduces no new base and composes exactly as a product
#: of lengths, which is what a molecular cavity surface is.
AREA: Dimension = (0, 2, 0, 0, 0, 0)
#: Volume, length cubed, named for the same reason: the unit system
#: already composes ``bohr^3`` and ``angstrom^3``; a selector declaring a
#: molecular volume needs the name.
VOLUME: Dimension = (0, 3, 0, 0, 0, 0)

CHARGE: Dimension = (0, 0, 0, 0, 0, 0, 0, 0, 1)
ELECTRIC_POTENTIAL: Dimension = (1, 0, 0, 0, 0, 0, 0, 0, -1)
#: Molar infrared absorption intensity, the tenth base, in units of km/mol.
#: It exists so that vibrational mode IR absorption strength is a typed
#: physical dimension rather than an unlabeled vector.
IR_INTENSITY: Dimension = (0, 0, 0, 0, 0, 0, 0, 0, 0, 1)

SUPPORTED_PYSCF_SELECTORS = frozenset(
    {
        "energy",
        "energies",
        "positions",
        # The artifact carries two structures -- what the driver was
        # handed (spec/positions) and where it ended (results/positions)
        # -- and a consumer asks for the role it needs.
        "supplied_positions",
        "reached_positions",
        "converged",
        "functional",
        "ab_initio",
        "mulliken_atomic_spin_populations",
        "connectivity",
        "surface_id",
        "symbols",
        "vibrational_frequencies",
        "vibrational_mode_atom_participation",
        "vibrational_mode_degeneracy_group",
        "homo",
        "lumo",
        "gap",
        "charge",
        "multiplicity",
        "method",
        "basis",
        "dipole_moment",
        "dipole_moment_magnitude",
        "excitation_energies",
        "oscillator_strengths",
        "spin_square",
        "spin_square_target",
        "spin_square_deviation",
        "effective_multiplicity",
        # The response stage (contract v5): per-root transition dipoles and
        # convergence, the root an excited-surface optimisation followed.
        "transition_dipole_moments",
        "excited_state_converged",
        "excited_state_followed_root",
        # The coupled-cluster components beside ``correlation_energy``.
        "ccsd_correlation_energy",
        "triples_correction",
    }
)

#: Selectors the log-parsing readers add on top of the structured PySCF set.
#: A selector being in this union does not mean every program answers it: each
#: reader declares what it provides, and a run that produced no such value is
#: refused as absent rather than guessed at.
SUPPORTED_SELECTORS = SUPPORTED_PYSCF_SELECTORS | frozenset(
    {
        "absorption_wavelengths",
        "excited_state_indices",
        "excited_state_labels",
        "excited_state_manifold_roots",
        "excited_state_multiplicities",
        "excited_state_spin_square",
        "excitation_energies",
        "entropy_times_temperature",
        "gibbs_free_energy",
        "oscillator_strengths",
        "singlet_excitation_energies",
        "triplet_excitation_energies",
        "singlet_oscillator_strengths",
        "triplet_oscillator_strengths",
        "spin_square",
        "spin_square_after_annihilation",
        "spin_square_target",
        "spin_square_deviation",
        "effective_multiplicity",
        "wavefunction_stability_verdict",
        "wavefunction_stability_history",
        "trajectory_frame_count",
        "trajectory_start_positions",
        "trajectory_end_positions",
        "trajectory_start_connectivity",
        "trajectory_end_connectivity",
        "trajectory_connectivity_changed",
        # A path's energy profile and the spectrum at its start: served by
        # the PySCF IRC, whose artifact carries the whole accepted path.
        "trajectory_energies",
        "trajectory_start_frequencies",
        "irc_direction",
        "solvation_model",
        "solvent",
        "scf_energy",
        "reference_energy",
        "correlation_energy",
        "dispersion_energy",
        "auxiliary_basis",
        "auxiliary_basis_role",
        "vpt2_harmonic_frequencies",
        "vpt2_fundamental_frequencies",
        "vpt2_zero_point_rovibrational_energy",
        "ir_intensities",
        # Method identity and convergence.  A composed workflow spanning many
        # nodes has to be able to assert that one functional and one method
        # ran across all of them -- ``all_equal_text`` is the predicate that
        # exists for exactly that -- and to ask whether an optimization
        # converged.  All four have had accessors, jobtype declarations, units
        # and dimensions for some time; only the request gate was missing, so
        # the capability was unreachable rather than absent.
        "ab_initio",
        "functional",
        "converged",
        "irc_converged",
        # Spin-resolved frontier orbitals, declared beside the restricted
        # pair.  An open-shell species has no single HOMO, so a radical in a
        # redox or hydrogen-transfer workflow needs these rather than ``homo``.
        "alpha_homo",
        "alpha_lumo",
        "beta_homo",
        "beta_lumo",
        # The relaxed-scan surface.  Read by the ORCA reader; deliberately not
        # in the PySCF set above, which has no scan jobtype and no extraction
        # branch for these names.
        "scan_coordinate_values",
        "scan_energies",
        "scan_point_indices",
        # Reached versus planned is the whole diagnosis when a scan dies
        # partway; both were parsed and unconsumed until a live scan died
        # at step 2 of 12 and no tool could state either number.
        "scan_steps_planned",
        "scan_steps_reached",
        # The structure an optimiser stopped on, which is not the
        # structure ``positions`` answers with: for an ORCA ``OptTS Freq``
        # that one is the geometry the Hessian was computed at, step 0.
        # Declared for the jobtypes that print a second structure, so a
        # session diagnosing an unconverged run can read the one it
        # reached instead of re-reading its own seed.
        "reached_positions",
        # The continuum solvation decomposition ORCA prints whenever a
        # continuum model is active.  The names are scheme-neutral but the
        # meanings are not interchangeable between programs, so only ORCA
        # declares them; see the reader for why Gaussian, PySCF and xTB do
        # not.
        "solvation_electrostatic_energy",
        "solvation_nonelectrostatic_energy",
        "solvation_cavity_surface_area",
        # Per-atom populations, named by the scheme that produced them
        # because Mulliken, Loewdin, Hirshfeld, CM5 and a tight-binding
        # population are different quantities rather than different
        # spellings of one, and no renormalisation removes the difference.
        "mulliken_atomic_charges",
        "loewdin_atomic_charges",
        "hirshfeld_atomic_charges",
        # xTB's ``charges`` sidecar is a self-consistent-charge population
        # from its tight-binding density.  It is deliberately *not* called
        # Mulliken: a common name would claim the two partitions agree.
        "xtb_scc_atomic_charges",
        # The second column of the same population block, read for years
        # and discarded one index from where a session needed it
        # (NOVEL-3 ino3, 2026-09-05); open shells only, sum 2S.
        "mulliken_atomic_spin_populations",
        "loewdin_atomic_spin_populations",
        # Symmetric atom-pair electronic bond-order matrix from population
        # analysis, in zero-based molecular atom order.
        "wiberg_bond_orders",
    }
)

#: Quantities a caller will reasonably ask this plane for that another tool
#: owns.  Refusing them by name alone sends the caller looking for a selector
#: that does not exist, when the quantity is available one tool away: the RRHO
#: engine computes all of these from a Hessian result and the extraction plane
#: reads structured fields.  Naming the producer is the same courtesy the host
#: registries already extend when an ID is unknown.
QUANTITIES_FROM_ANOTHER_TOOL: Mapping[str, str] = {
    "electronic_energy": "derive_thermochemistry",
    "enthalpy": "derive_thermochemistry",
    "entropy": "derive_thermochemistry",
    "internal_energy": "derive_thermochemistry",
    "quasi_harmonic_enthalpy": "derive_thermochemistry",
    "quasi_harmonic_entropy": "derive_thermochemistry",
    "quasi_harmonic_entropy_times_temperature": "derive_thermochemistry",
    "quasi_harmonic_gibbs_free_energy": "derive_thermochemistry",
    "quasi_harmonic_thermal_gibbs_correction": "derive_thermochemistry",
    "thermal_enthalpy_correction": "derive_thermochemistry",
    "enthalpy_increment_above_zero_point": "derive_thermochemistry",
    "thermal_gibbs_correction": "derive_thermochemistry",
    "thermal_internal_energy_correction": "derive_thermochemistry",
    "zero_point_energy": "derive_thermochemistry",
}

#: Quantity IDs ``derive_result_thermochemistry`` writes into every receipt.
#: A planned thermochemistry node names its outputs from this vocabulary, so
#: the list lives beside the builder where drift between the two is visible.
DERIVABLE_THERMOCHEMISTRY_QUANTITIES: tuple[str, ...] = (
    "electronic_energy",
    "enthalpy",
    "enthalpy_increment_above_zero_point",
    "entropy",
    "entropy_times_temperature",
    "gibbs_free_energy",
    "heat_capacity_cv",
    "internal_energy",
    "near_zero_mode_count",
    "pressure",
    "temperature",
    "thermal_enthalpy_correction",
    "thermal_gibbs_correction",
    "thermal_internal_energy_correction",
    "zero_point_energy",
)

#: Written only when a quasi-harmonic entropy method is requested.  A strict
#: RRHO request never produces them, so a node that asks for one under 'rrho'
#: is asking for a quantity its own receipt will not carry.
QUASI_HARMONIC_THERMOCHEMISTRY_QUANTITIES: tuple[str, ...] = (
    "quasi_harmonic_entropy",
    "quasi_harmonic_entropy_times_temperature",
    "quasi_harmonic_gibbs_free_energy",
    "quasi_harmonic_thermal_gibbs_correction",
)

#: A quasi-harmonic treatment never changes what a harmonic name means: it
#: adds a counterpart beside it.  A Grimme receipt's ``gibbs_free_energy``
#: is the RRHO value and its Grimme value is
#: ``quasi_harmonic_gibbs_free_energy``; the receipt carries both, so one
#: receipt can measure the spread between them.  Keyed by the harmonic
#: name: the counterpart, and which requested treatment writes it --
#: ``entropy`` (Grimme or Truhlar), ``enthalpy`` (a Head-Gordon enthalpy
#: cutoff) or ``either``.
#:
#: The writer (``derive_result_thermochemistry``), the plan-time contract
#: and the planning schema all read this one table.  They were three: the
#: schema listed only the harmonic kinds, the executor binds a planned
#: output by its kind, and two live campaigns planned Grimme nodes whose
#: outputs were the harmonic ``gibbs_free_energy`` -- one delivered an
#: entropy-model uncertainty of exactly 0.0 kcal/mol where the same
#: receipts give 0.3564 (po3-r19 cycle 5), the other fed ten harmonic Gibbs
#: energies into a pKa under a review that said Grimme.
QUASI_HARMONIC_COUNTERPARTS: Mapping[str, tuple[str, str]] = {
    "entropy": ("quasi_harmonic_entropy", "entropy"),
    "entropy_times_temperature": (
        "quasi_harmonic_entropy_times_temperature",
        "entropy",
    ),
    "enthalpy": ("quasi_harmonic_enthalpy", "enthalpy"),
    "gibbs_free_energy": ("quasi_harmonic_gibbs_free_energy", "either"),
    "thermal_gibbs_correction": (
        "quasi_harmonic_thermal_gibbs_correction",
        "either",
    ),
}


def quasi_harmonic_counterparts_for_treatment(
    entropy_method: str | None = "rrho",
    enthalpy_cutoff_cm1: float | None = None,
) -> dict[str, str]:
    """Harmonic name -> the quasi-harmonic counterpart this treatment writes.

    Empty for a strictly harmonic request, which writes no counterpart.
    """

    entropy = str(entropy_method or "rrho").strip().lower() != "rrho"
    enthalpy = enthalpy_cutoff_cm1 is not None
    written = {
        "entropy": entropy,
        "enthalpy": enthalpy,
        "either": entropy or enthalpy,
    }
    return {
        harmonic: counterpart
        for harmonic, (
            counterpart,
            needs,
        ) in QUASI_HARMONIC_COUNTERPARTS.items()
        if written[needs]
    }


def thermochemistry_quantities_for_treatment(
    entropy_method: str | None = "rrho",
    enthalpy_cutoff_cm1: float | None = None,
) -> tuple[str, ...]:
    """Every quantity id one receipt carries under this treatment."""

    names = set(DERIVABLE_THERMOCHEMISTRY_QUANTITIES)
    names.update(
        quasi_harmonic_counterparts_for_treatment(
            entropy_method, enthalpy_cutoff_cm1
        ).values()
    )
    return tuple(sorted(names))


#: Names a scientist may reasonably use for one of the canonical IDs.  The
#: receipt writes the canonical name, so the plan-time contract and the
#: completion matcher both resolve through this one table.
THERMOCHEMISTRY_QUANTITY_ALIASES: Mapping[str, str] = {
    # Gibbs free-energy correction and thermal free-energy correction are the
    # same G(T)-E_electronic quantity; the typed engine uses the former label.
    "thermal_free_energy_correction": "thermal_gibbs_correction",
}


def canonical_thermochemistry_quantity(name: str) -> str:
    """Return the receipt's own ID for a declared thermochemistry quantity."""

    key = str(name).strip().lower()
    return THERMOCHEMISTRY_QUANTITY_ALIASES.get(key, key)


def derivable_thermochemistry_quantities(
    entropy_method: str | None = "rrho",
    enthalpy_cutoff_cm1: float | None = None,
) -> tuple[str, ...]:
    """Return the quantity IDs a receipt will carry for this treatment.

    The enthalpy cutoff is part of the treatment: a Head-Gordon request
    writes ``quasi_harmonic_enthalpy`` and a quasi-harmonic Gibbs energy
    even under harmonic entropy, and a table keyed on the entropy method
    alone told a plan those did not exist.
    """

    return thermochemistry_quantities_for_treatment(
        entropy_method, enthalpy_cutoff_cm1
    )


def thermochemistry_route_hint(selectors) -> str:
    """One sentence naming the open door, when refused selectors have one.

    A refusal that only lists what is declared teaches a session what it
    cannot have; it does not say where the quantity actually lives. A
    live session asked a result for its free energy, was correctly told
    no jobtype declares that selector, reported the closed gate as a
    limitation -- and never reached the thermochemistry stage that was
    open, while a sibling session used it on identical results the same
    hour. The refusal is the moment the route must be named.

    Returns "" when none of the refused selectors is a thermochemistry
    quantity, so callers can append unconditionally.
    """

    derivable = set(derivable_thermochemistry_quantities(None))
    derivable.update(
        counterpart for counterpart, _ in QUASI_HARMONIC_COUNTERPARTS.values()
    )
    matched = sorted(
        {
            str(item)
            for item in selectors
            if canonical_thermochemistry_quantity(item) in derivable
        }
    )
    if not matched:
        return ""
    return (
        f" A thermochemistry stage derives {matched} from a "
        "frequency-bearing result at a stationary point -- a converged "
        "optimisation or saddle search, or a Hessian taken at one -- under "
        "an explicitly bound temperature, pressure, and standard state; "
        "plan one on such a calculation instead of selecting these from the "
        "log. A structure held or driven along a coordinate (modred, scan) "
        "is not a stationary point: for its free energy as one, relax it "
        "without the constraint first; for the free energy of the surface "
        "at the held value, derive on the held (modred) result with "
        "projected_coordinates naming the coordinates it held."
    )


#: The layers one structure's thermochemical state functions are sums of,
#: in the order they accumulate: E(el), + ZPE = E0, + thermal energy above
#: the zero-point level = U(T), + pV = H(T), - TS = G(T).  An energy a
#: receipt carries is a fixed block of these for one structure, which is
#: what lets the host read a combination of energies as a reaction.
ENERGY_LAYERS = ("electronic", "zero_point", "thermal", "pV", "minus_TS")

#: What a contiguous block of layers is called when one species enters an
#: output with it.
ENERGY_LAYER_BLOCKS: Mapping[tuple[int, ...], str] = {
    (0,): "E",
    (0, 1): "E0 (E + ZPE)",
    (0, 1, 2): "U",
    (0, 1, 2, 3): "H",
    (0, 1, 2, 3, 4): "G",
}


@dataclass(frozen=True)
class EnergyKindV1:
    """What an energy-valued number is, beyond its dimension.

    Dimensional analysis checks that two numbers share a unit; it cannot
    tell a total energy from an orbital eigenvalue or an orbital-rotation
    curvature, all of which the vocabulary serves in hartree.  ``kind`` is
    one of:

    - ``state_energy``: the energy of one structure in one state, as the
      block of ``ENERGY_LAYERS`` in ``layers`` (``E``, ``E0``, ``U``, ``H``
      or ``G``);
    - ``correction``: a block of layers without the electronic one (a
      zero-point energy, a thermal correction, T*S with ``sign`` -1);
    - ``electronic_component``: a part of one electronic energy (a
      correlation, dispersion or solvation term) -- a decomposition, with
      no layer of its own;
    - ``orbital_energy``: a one-electron eigenvalue or a gap between two;
    - ``excitation_energy``: a vertical energy between two states of one
      structure, positive by convention;
    - ``orbital_rotation_curvature``: an eigenvalue of an SCF stability
      matrix -- the curvature of the energy along a rotation of the
      orbitals, not an energy difference between states -- under the
      ``normalisation`` of the matrix it belongs to.

    ``indexes_structures`` marks a vector whose elements are the energies
    of different structures -- the points of a scan, an IRC branch or an
    optimisation -- so an element picked by index is its own structure,
    not the one the result reached.
    """

    kind: str
    layers: tuple[int, ...] = ()
    sign: int = 1
    normalisation: str = ""
    indexes_structures: bool = False


_E = EnergyKindV1("state_energy", (0,))
_E_ALONG_A_PATH = EnergyKindV1("state_energy", (0,), indexes_structures=True)
_COMPONENT = EnergyKindV1("electronic_component")
_ORBITAL = EnergyKindV1("orbital_energy")
_EXCITATION = EnergyKindV1("excitation_energy")

#: One kind per energy-valued selector or thermochemistry quantity id.  The
#: two name spaces agree wherever they share a name (``gibbs_free_energy``
#: printed by a program and derived by the host are both a G), so they are
#: one table; ``test_every_energy_names_its_kind`` holds every energy a
#: reader declares, and every one a thermochemistry receipt writes, to it.
ENERGY_KINDS: Mapping[str, EnergyKindV1] = {
    # One structure's electronic energy, however the program reached it.
    "energy": _E,
    "energies": _E_ALONG_A_PATH,
    "scf_energy": _E,
    "reference_energy": _E,
    "electronic_energy": _E,
    "scan_energies": _E_ALONG_A_PATH,
    "trajectory_energies": _E_ALONG_A_PATH,
    # Thermochemical state functions of one structure.
    "internal_energy": EnergyKindV1("state_energy", (0, 1, 2)),
    "enthalpy": EnergyKindV1("state_energy", (0, 1, 2, 3)),
    "quasi_harmonic_enthalpy": EnergyKindV1("state_energy", (0, 1, 2, 3)),
    "gibbs_free_energy": EnergyKindV1("state_energy", (0, 1, 2, 3, 4)),
    "quasi_harmonic_gibbs_free_energy": EnergyKindV1(
        "state_energy", (0, 1, 2, 3, 4)
    ),
    # Corrections: layers above the electronic energy of one structure.
    "zero_point_energy": EnergyKindV1("correction", (1,)),
    "thermal_internal_energy_correction": EnergyKindV1("correction", (1, 2)),
    "thermal_enthalpy_correction": EnergyKindV1("correction", (1, 2, 3)),
    "enthalpy_increment_above_zero_point": EnergyKindV1("correction", (2, 3)),
    "thermal_gibbs_correction": EnergyKindV1("correction", (1, 2, 3, 4)),
    "quasi_harmonic_thermal_gibbs_correction": EnergyKindV1(
        "correction", (1, 2, 3, 4)
    ),
    "entropy_times_temperature": EnergyKindV1("correction", (4,), sign=-1),
    "quasi_harmonic_entropy_times_temperature": EnergyKindV1(
        "correction", (4,), sign=-1
    ),
    # Parts of one electronic energy.
    "correlation_energy": _COMPONENT,
    "ccsd_correlation_energy": _COMPONENT,
    "triples_correction": _COMPONENT,
    "dispersion_energy": _COMPONENT,
    "solvation_electrostatic_energy": _COMPONENT,
    "solvation_nonelectrostatic_energy": _COMPONENT,
    "solvation_free_energy": _COMPONENT,
    "xtb_solvation_sasa_energy": _COMPONENT,
    "xtb_solvation_hydrogen_bond_energy": _COMPONENT,
    "xtb_solvation_shift_energy": _COMPONENT,
    # One-electron eigenvalues.
    "homo": _ORBITAL,
    "lumo": _ORBITAL,
    "gap": _ORBITAL,
    "alpha_homo": _ORBITAL,
    "alpha_lumo": _ORBITAL,
    "beta_homo": _ORBITAL,
    "beta_lumo": _ORBITAL,
    # Vertical state-to-state energies.
    "excitation_energies": _EXCITATION,
    "singlet_excitation_energies": _EXCITATION,
    "triplet_excitation_energies": _EXCITATION,
    # Stability-matrix eigenvalues: curvatures, each of its own matrix
    # (R10 Q13: PySCF's internal root is 4 x Gaussian's singlet (A+B) root,
    # its external root equals Gaussian's triplet root to 1e-6 Eh, and
    # real -> complex is (A-B)).
    "scf_stability_internal_lowest_eigenvalue": EnergyKindV1(
        "orbital_rotation_curvature",
        normalisation=(
            "PySCF internal: real rotations that keep the reference's spin "
            "form, eigenvalue of its orbital Hessian 4(A+B)"
        ),
    ),
    "scf_stability_external_lowest_eigenvalue": EnergyKindV1(
        "orbital_rotation_curvature",
        normalisation=(
            "PySCF external: rotations into the next space (RHF/RKS -> "
            "UHF/UKS triplet block, UHF/UKS -> GHF/GKS), eigenvalue of that "
            "block's (A+B)"
        ),
    ),
    "scf_stability_real_to_complex_lowest_eigenvalue": EnergyKindV1(
        "orbital_rotation_curvature",
        normalisation=(
            "PySCF real -> complex: imaginary rotations, eigenvalue of (A-B)"
        ),
    ),
    "wavefunction_stability_lowest_eigenvalue": EnergyKindV1(
        "orbital_rotation_curvature",
        normalisation=(
            "Gaussian: lowest (A+B) root over the blocks it tested (singlet "
            "and RHF -> UHF triplet for a restricted reference)"
        ),
    ),
}


def energy_kind(name: str) -> EnergyKindV1 | None:
    """The kind of an energy selector or thermochemistry quantity id."""

    return ENERGY_KINDS.get(str(name))


def exists_only_at_a_stationary_point(name: str) -> bool:
    """Whether a quantity is defined only at a stationary point.

    Everything that carries a vibrational or thermal layer -- a zero-point
    energy, a thermal correction, T*S, U, H, G -- is a harmonic expansion
    about a point where the gradient vanishes; an electronic energy, an
    orbital eigenvalue or a curvature is defined anywhere.
    """

    kind = energy_kind(name)
    return kind is not None and any(layer > 0 for layer in kind.layers)


_IDENTIFIER = re.compile(r"^[A-Za-z][A-Za-z0-9_.:-]{0,127}$")
#: Executed-evidence contracts the analysis plane admits, owned by the
#: writer: a previous supported contract is a subset of the current one.
_SUPPORTED_PYSCF_RESULT_CONTRACTS = SUPPORTED_RESULT_CONTRACT_VERSIONS


def supported_selectors() -> frozenset[str]:
    """Every selector a request may name: the shared set and what readers add.

    ``SUPPORTED_SELECTORS`` is the vocabulary several programs share.  A
    selector only one program's parser can answer is declared on that
    program's reader (``ResultReaderV1.selector_declarations``) and joins the
    gate here, so introducing it is not an edit to a set every program owns.
    Imported lazily, as the other reader lookups in this module are: the
    readers import this module's types while they are being built.
    """

    from chemsmart.analysis.result_readers import DECLARED_SELECTORS

    return SUPPORTED_SELECTORS | DECLARED_SELECTORS


class QuantityContractError(ValueError):
    """Raised when a quantity request or result violates its typed contract."""


class QuantityExtractionError(QuantityContractError):
    """Raised when trusted result evidence cannot be extracted safely."""


#: Which HDF5 dataset(s) each selector reads.  The *unit* each dataset
#: carries is the writer's declaration (``RESULT_UNITS``), never a second
#: table here: two hand-written tables once disagreed on ``normal_modes``
#: (writer amu^-1/2, reader "dimensionless") and the declared selector was
#: refused on every real PySCF Hessian while every synthetic test passed.
#: Deriving the expectation from the writer means a disagreement can only
#: be an artifact written under another contract, which is a fact to state.
_SELECTOR_RESULT_DATASETS: dict[str, tuple[str, ...]] = {
    # ``energy`` reads ``results/total_energy`` on a v5 artifact and the
    # SCF trace before it, so the audited dataset is the one every
    # contract writes; the v5 scalars audit their own datasets.
    "energy": ("results/energies",),
    "energies": ("results/energies",),
    "scf_energy": ("results/scf_energy",),
    "excitation_energies": ("results/excitation_energies",),
    "oscillator_strengths": ("results/oscillator_strengths",),
    "singlet_excitation_energies": ("results/excitation_energies",),
    "triplet_excitation_energies": ("results/excitation_energies",),
    "singlet_oscillator_strengths": ("results/oscillator_strengths",),
    "triplet_oscillator_strengths": ("results/oscillator_strengths",),
    "excited_state_indices": ("results/excitation_energies",),
    "excited_state_manifold_roots": ("results/excitation_energies",),
    "excited_state_multiplicities": ("results/excited_state_multiplicities",),
    "excited_state_converged": ("results/excited_state_converged",),
    "transition_dipole_moments": ("results/transition_dipole_moments",),
    "reference_energy": ("results/reference_energy",),
    "correlation_energy": ("results/correlation_energy",),
    "ccsd_correlation_energy": ("results/ccsd_correlation_energy",),
    "triples_correction": ("results/triples_correction",),
    # The decomposition of the total (contract v10).  The model and the
    # solvent are spec fields and have no stored unit to audit; the three
    # terms are datasets and do.
    "dispersion_energy": ("results/dispersion_energy",),
    "solvation_electrostatic_energy": (
        "results/solvation_electrostatic_energy",
    ),
    "solvation_nonelectrostatic_energy": (
        "results/solvation_nonelectrostatic_energy",
    ),
    "dipole_moment": ("results/dipole_moment",),
    "dipole_moment_magnitude": ("results/dipole_moment",),
    "mulliken_atomic_charges": ("results/mulliken_charges",),
    "mulliken_atomic_spin_populations": ("results/mulliken_spin_populations",),
    "positions": ("results/positions",),
    "reached_positions": ("results/positions",),
    "connectivity": ("results/positions",),
    "vibrational_frequencies": ("results/vibrational_frequencies",),
    "vibrational_mode_atom_participation": ("results/normal_modes",),
    "vibrational_mode_degeneracy_group": ("results/vibrational_frequencies",),
    "homo": ("results/mo_energy",),
    "lumo": ("results/mo_energy",),
    "gap": ("results/mo_energy",),
    # The same dataset, read through ``results/mo_occ`` per spin channel.
    "alpha_homo": ("results/mo_energy",),
    "alpha_lumo": ("results/mo_energy",),
    "beta_homo": ("results/mo_energy",),
    "beta_lumo": ("results/mo_energy",),
    "spin_square": ("results/spin_square",),
    "spin_square_deviation": ("results/spin_square",),
    "effective_multiplicity": ("results/spin_square_effective_multiplicity",),
    # The IRC stage (contract v8), under ``results/irc/``.
    "trajectory_frame_count": ("results/irc/path_positions",),
    "trajectory_start_positions": ("results/irc/path_positions",),
    "trajectory_end_positions": ("results/irc/path_positions",),
    "trajectory_start_connectivity": ("results/irc/path_positions",),
    "trajectory_end_connectivity": ("results/irc/path_positions",),
    "trajectory_connectivity_changed": ("results/irc/path_positions",),
    "trajectory_energies": ("results/irc/path_energies",),
    # One question with two homes: the spectrum of the geometry a path
    # stage was handed, which an ``irc`` files under its start and a
    # ``ts`` under its seed. A selector is absent when *no* home carries
    # it; a home that is present under the wrong unit is still a
    # divergence to state, which is why the two are asked separately.
    "trajectory_start_frequencies": (
        "results/irc/start_frequencies",
        "results/ts/seed_frequencies",
    ),
}


def _dataset_unit(path: str) -> str:
    leaf = path.rsplit("/", 1)[-1]
    try:
        return RESULT_UNITS[leaf]
    except (
        KeyError
    ) as exc:  # a selector naming a dataset the writer never writes
        raise QuantityContractError(
            f"selector dataset {path!r} has no unit in the writer's table"
        ) from exc


#: Derived, never hand-written: selector -> {dataset path: expected unit}.
_SELECTOR_RESULT_UNITS: dict[str, dict[str, str]] = {
    selector: {path: _dataset_unit(path) for path in paths}
    for selector, paths in _SELECTOR_RESULT_DATASETS.items()
}


def _freeze(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return _freeze(value.tolist())
    if isinstance(value, np.generic):
        return _freeze(value.item())
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, tuple):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, dict):
        return tuple(
            (str(key), _freeze(item))
            for key, item in sorted(
                value.items(), key=lambda pair: str(pair[0])
            )
        )
    if isinstance(value, float) and not math.isfinite(value):
        raise QuantityContractError("quantity values must be finite")
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise QuantityContractError(
        f"unsupported immutable quantity value: {type(value).__name__}"
    )


def _canonical_data(value: Any) -> Any:
    if is_dataclass(value):
        return _canonical_data(asdict(value))
    if isinstance(value, dict):
        return {
            str(key): _canonical_data(item)
            for key, item in sorted(
                value.items(), key=lambda pair: str(pair[0])
            )
        }
    if isinstance(value, (tuple, list)):
        return [_canonical_data(item) for item in value]
    if isinstance(value, float):
        if not math.isfinite(value):
            raise QuantityContractError("canonical records must be finite")
        return value
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise QuantityContractError(
        f"unsupported canonical value: {type(value).__name__}"
    )


def canonical_quantity_sha256(value: Any) -> str:
    encoded = json.dumps(
        _canonical_data(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def result_file_sha256(path: str | os.PathLike[str]) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_identifier(value: str, field: str) -> str:
    normalized = str(value).strip()
    if not _IDENTIFIER.fullmatch(normalized):
        raise QuantityContractError(f"{field} is not a stable identifier")
    return normalized


def _require_sha256(value: str) -> str:
    normalized = str(value).strip().lower()
    if len(normalized) != 64:
        raise QuantityContractError(
            "artifact_sha256 must contain 64 hex digits"
        )
    try:
        int(normalized, 16)
    except ValueError as exc:
        raise QuantityContractError(
            "artifact_sha256 must contain 64 hex digits"
        ) from exc
    return normalized


@dataclass(frozen=True)
class QuantitySelectorV1:
    """Select one semantic quantity without exposing a file-query language."""

    quantity_id: str
    selector: str

    def __post_init__(self) -> None:
        _require_identifier(self.quantity_id, "quantity_id")
        if self.selector not in supported_selectors():
            elsewhere = QUANTITIES_FROM_ANOTHER_TOOL.get(self.selector)
            detail = (
                f"; that quantity is produced by {elsewhere}, not by result "
                "extraction"
                if elsewhere
                else f"; supported selectors: {sorted(supported_selectors())}"
            )
            raise QuantityContractError(
                f"unsupported quantity selector: {self.selector!r}{detail}"
            )


@dataclass(frozen=True)
class ResultQuantityExtractionRequestV1:
    schema_version: str
    artifact_id: str
    artifact_sha256: str
    program: str
    selectors: tuple[QuantitySelectorV1, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "selectors", tuple(self.selectors))
        if self.schema_version != "chemsmart.quantity-extraction-request.v1":
            raise QuantityContractError(
                "unsupported extraction request schema"
            )
        _require_identifier(self.artifact_id, "artifact_id")
        _require_sha256(self.artifact_sha256)
        from chemsmart.analysis.result_readers import reader_for

        if reader_for(self.program) is None:
            raise QuantityContractError(
                f"no result reader is registered for {self.program!r}"
            )
        if not self.selectors:
            raise QuantityContractError(
                "at least one quantity selector is required"
            )
        quantity_ids = [selector.quantity_id for selector in self.selectors]
        if len(quantity_ids) != len(set(quantity_ids)):
            raise QuantityContractError("quantity_id values must be unique")


@dataclass(frozen=True)
class QuantityValueV1:
    """An immutable value with its parser unit and canonical arithmetic unit."""

    schema_version: str
    quantity_id: str
    data_kind: str
    source_value: Any
    source_unit: str
    value: Any
    unit: str
    dimension: Dimension
    evidence_ref: str
    value_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.quantity-value.v1":
            raise QuantityContractError("unsupported quantity value schema")
        _require_identifier(self.quantity_id, "quantity_id")
        if self.data_kind not in {
            "scalar",
            "vector",
            "matrix",
            "integer",
            "text",
            "text_vector",
        }:
            raise QuantityContractError("unsupported quantity data kind")
        if len(self.dimension) not in {6, 7, 8, 9, 10} or not all(
            isinstance(exponent, int) for exponent in self.dimension
        ):
            raise QuantityContractError(
                "dimension must contain six legacy, seven dipole-extended, "
                "eight mass-extended, nine charge-extended, or ten "
                "ir-intensity-extended integers"
            )
        object.__setattr__(self, "source_value", _freeze(self.source_value))
        object.__setattr__(self, "value", _freeze(self.value))
        body = {
            "schema_version": self.schema_version,
            "quantity_id": self.quantity_id,
            "data_kind": self.data_kind,
            "source_value": self.source_value,
            "source_unit": self.source_unit,
            "value": self.value,
            "unit": self.unit,
            "dimension": self.dimension,
            "evidence_ref": self.evidence_ref,
        }
        if self.value_sha256 != canonical_quantity_sha256(body):
            raise QuantityContractError("quantity value digest mismatch")


def canonical_extraction_receipt_body(
    *,
    schema_version: str,
    artifact_id: str,
    artifact_sha256: str,
    program: str,
    parser_id: str,
    quantities: Any,
    status: str,
    absent: Any = (),
    derived_adjacency: Any = (),
    electronic_provenance: Any = (),
    selector_bindings: Any = (),
    structural_states: Any = (),
    level: Any = (),
    native_evidence: Any = (),
) -> dict[str, Any]:
    """Return the one body an extraction receipt is digested over.

    Two digests used to be computed over this record -- the receipt's own,
    from a hand-built dict, and the durable event's, from the dataclass
    minus its digest field -- and they agreed only because the two happened
    to contain the same keys.  Adding a field that is present on the
    dataclass and absent from a complete receipt's body broke that
    coincidence and settled a node as failed with a digest mismatch, which
    is why the rule now lives in one place and is called from all three.

    ``absent`` enters the body only when it is non-empty, so every receipt
    written before absences existed stays verifiable.  Status and absences
    are locked together by the receipt, so a body without the key is
    unambiguously a complete extraction.  ``derived_adjacency`` follows the
    same rule.
    """

    body: dict[str, Any] = {
        "schema_version": schema_version,
        "artifact_id": artifact_id,
        "artifact_sha256": artifact_sha256,
        "program": program,
        "parser_id": parser_id,
        "quantities": quantities,
        "status": status,
    }
    if absent:
        body["absent"] = absent
    if derived_adjacency:
        # Present only when a positions read bundled the perceived bond
        # list, so every receipt minted before the field existed -- and
        # every read that carries no adjacency -- verifies under one
        # arithmetic.
        body["derived_adjacency"] = derived_adjacency
    if electronic_provenance:
        # Whose density or method each delivered value belongs to, as
        # ``((quantity_id, word), ...)``; present only when a reader
        # declares the axis, under the same rule as the adjacency.
        body["electronic_provenance"] = electronic_provenance
    if selector_bindings:
        # The quantity id is model-chosen, while the selector is the host's
        # semantic contract. Preserve their exact pairing durably rather than
        # only on the transient tool event.
        body["selector_bindings"] = selector_bindings
    if structural_states:
        # State answers which geometry a value belongs to, including
        # ``stateless`` values whose density provenance is absent.
        body["structural_states"] = structural_states
    if level:
        # This is the producer's own level record. It is displayed for a
        # scientist to compare, never used to infer equivalence.
        body["level"] = level
    if native_evidence:
        # Selectors such as xTB Wiberg bond orders are read from a native
        # sidecar rather than the result's primary log.  Bind the selector,
        # result-relative filename, and exact bytes into the same receipt as
        # the delivered quantity.  This is deliberately selector-scoped:
        # the quantity id is model-authored and has no parser semantics.
        body["native_evidence"] = native_evidence
    return body


@dataclass(frozen=True)
class QuantityExtractionReceiptV1:
    schema_version: str
    artifact_id: str
    artifact_sha256: str
    program: str
    parser_id: str
    quantities: tuple[QuantityValueV1, ...]
    status: str
    receipt_sha256: str
    #: Quantities this request asked for that the result does not carry, as
    #: ``(quantity_id, selector, reason)``.  A receipt used to record only
    #: what it held, which made a partial extraction indistinguishable from
    #: a complete one for a smaller request -- so the evidence chain could
    #: weaken without saying so.  An absence is a stated gap with the
    #: reader's own reason, never a value the host chose to substitute.
    absent: tuple[tuple[str, str, str], ...] = ()
    #: Bond list bundled onto a delivered positions read: the perceived
    #: covalent adjacency of the same structure, as
    #: ``{"formula": str, "bond_atom_pairs": ((i, j), ...)}``.  The hard
    #: line, stated once and forever: derived adjacency is a measurement
    #: and is allowed; any perceived label -- isomer, conformer, species,
    #: ring class, stereo descriptor -- is the scientist's judgement and
    #: must never be attached here.
    derived_adjacency: Any = ()
    #: Whose density or method each delivered value belongs to, as
    #: ``((quantity_id, word), ...)`` over the reader's declared
    #: electronic-provenance axis: ``reference`` for a mean-field property,
    #: ``excited_root`` for a root's quantity, ``correlated`` for a
    #: correlated component, resolved per artifact.  A geometry identity
    #: says nothing about whose density a dipole is, and this is where the
    #: receipt says it.  Empty for a reader that declares no axis.
    electronic_provenance: Any = ()
    #: Exact ``(quantity_id, selector)`` pairs requested from the reader.
    #: A selector carries scheme semantics, so it cannot remain event-only.
    selector_bindings: Any = ()
    #: Exact structural role per delivered or explicitly absent quantity.
    structural_states: Any = ()
    #: Producer level read from the same verified result artifact.
    level: Any = ()
    #: Exact native sidecars consumed by requested selectors, as
    #: ``(selector, result-relative filename, sha256)``.  The primary result
    #: remains ``artifact_id``/``artifact_sha256`` above; this field prevents
    #: a sidecar-derived value from being represented as evidence from that
    #: log alone.
    native_evidence: Any = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "quantities", tuple(self.quantities))
        object.__setattr__(
            self, "absent", tuple(tuple(item) for item in self.absent)
        )
        if self.derived_adjacency:
            adjacency = dict(self.derived_adjacency)
            pairs = tuple(
                tuple(int(index) for index in pair)
                for pair in adjacency.get("bond_atom_pairs", ())
            )
            # An allow-list, not an exact set. The line this enforces is
            # that no perceived *label* may ride here; it is not that a
            # measurement may carry nothing about itself. A bond list
            # alone is one boolean per pair, and a boolean produced by a
            # threshold cannot be distinguished from a structural fact --
            # a converged formaldehyde lost both C-H bonds by 1.5 mA and
            # the reader had no way to see it (sm1, 2026-09-11). The
            # convention's id and the per-pair margin are provenance and
            # arithmetic about the same measurement, so they are
            # admitted; anything not named here is still refused.
            permitted = {
                "formula",
                "bond_atom_pairs",
                "adjacency_policy_id",
                "adjacency_margins",
            }
            unknown = sorted(set(adjacency) - permitted)
            if (
                unknown
                or "formula" not in adjacency
                or "bond_atom_pairs" not in adjacency
                or not (
                    isinstance(adjacency.get("formula"), str)
                    and adjacency["formula"]
                    and all(len(pair) == 2 for pair in pairs)
                )
            ):
                raise QuantityContractError(
                    "derived adjacency carries a formula, bond atom pairs, "
                    "and optionally the perception policy id and per-pair "
                    "margins; a perceived label -- isomer, conformer, "
                    "species, ring class, stereo descriptor -- is the "
                    "scientist's judgement and is never attached here"
                    + (f" (offending field(s): {unknown})" if unknown else "")
                )
            adjacency["bond_atom_pairs"] = pairs
            margins = adjacency.get("adjacency_margins")
            if margins:
                adjacency["adjacency_margins"] = tuple(
                    {
                        "atoms": tuple(int(i) for i in row["atoms"]),
                        "distance_angstrom": float(row["distance_angstrom"]),
                        "cutoff_angstrom": float(row["cutoff_angstrom"]),
                        "margin_angstrom": float(row["margin_angstrom"]),
                        "adjacent": bool(row["adjacent"]),
                    }
                    for row in margins
                )
            object.__setattr__(self, "derived_adjacency", adjacency)
        if self.schema_version != "chemsmart.quantity-extraction-receipt.v1":
            raise QuantityContractError(
                "unsupported extraction receipt schema"
            )
        if self.status not in {"extracted", "partial"}:
            raise QuantityContractError("invalid extraction receipt status")
        if bool(self.absent) != (self.status == "partial"):
            # The status is what an older reader checks, so it must be the
            # thing that changes.  A partial receipt with no absence named,
            # or a complete one carrying absences, would let the two drift.
            raise QuantityContractError(
                "extraction receipt status must be 'partial' exactly when "
                "the request named quantities the result does not carry"
            )
        for item in self.absent:
            if len(item) != 3 or not all(
                isinstance(field, str) and field for field in item
            ):
                raise QuantityContractError(
                    "each absence records a quantity id, a selector and a "
                    "reason"
                )
        provenance = tuple(
            (str(quantity_id), str(word))
            for quantity_id, word in (self.electronic_provenance or ())
        )
        for _quantity_id, word in provenance:
            if not word or word == "stateless":
                raise QuantityContractError(
                    "electronic provenance names a state or method per "
                    "delivered quantity; a stateless value is not recorded"
                )
        object.__setattr__(self, "electronic_provenance", provenance)
        bindings = tuple(
            (str(quantity_id), str(selector))
            for quantity_id, selector in (self.selector_bindings or ())
        )
        if bindings and len({item[0] for item in bindings}) != len(bindings):
            raise QuantityContractError(
                "selector bindings must name each quantity id at most once"
            )
        object.__setattr__(self, "selector_bindings", bindings)
        states = tuple(
            (str(quantity_id), str(state))
            for quantity_id, state in (self.structural_states or ())
        )
        if states and len({item[0] for item in states}) != len(states):
            raise QuantityContractError(
                "structural states must name each quantity id at most once"
            )
        object.__setattr__(self, "structural_states", states)
        object.__setattr__(self, "level", dict(self.level or {}))
        native_evidence = tuple(
            (str(selector), str(filename), str(sha256))
            for selector, filename, sha256 in (self.native_evidence or ())
        )
        if native_evidence != tuple(sorted(set(native_evidence))):
            raise QuantityContractError(
                "native evidence must be sorted, unique selector/path/digest "
                "records"
            )
        for selector, filename, sha256 in native_evidence:
            if selector not in supported_selectors():
                raise QuantityContractError(
                    f"native evidence names unsupported selector {selector!r}"
                )
            if (
                not filename
                or filename.startswith("/")
                or filename in {".", ".."}
                or ".." in filename.split("/")
            ):
                raise QuantityContractError(
                    "native evidence filename must be a non-empty "
                    "result-relative path"
                )
            _require_sha256(sha256)
        object.__setattr__(self, "native_evidence", native_evidence)
        body = canonical_extraction_receipt_body(
            schema_version=self.schema_version,
            artifact_id=self.artifact_id,
            artifact_sha256=self.artifact_sha256,
            program=self.program,
            parser_id=self.parser_id,
            quantities=self.quantities,
            status=self.status,
            absent=self.absent,
            derived_adjacency=self.derived_adjacency,
            electronic_provenance=self.electronic_provenance,
            selector_bindings=self.selector_bindings,
            structural_states=self.structural_states,
            level=self.level,
            native_evidence=self.native_evidence,
        )
        if self.receipt_sha256 != canonical_quantity_sha256(body):
            raise QuantityContractError(
                "quantity extraction receipt digest mismatch"
            )


@dataclass(frozen=True)
class ThermochemistryRequestV1:
    schema_version: str
    artifact_id: str
    artifact_sha256: str
    program: str
    temperature_k: float
    pressure_atm: float
    concentration_mol_l: float | None = None
    entropy_method: str = "rrho"
    entropy_cutoff_cm1: float | None = None
    enthalpy_cutoff_cm1: float | None = None
    alpha: int = 4
    use_weighted_mass: bool = False
    frequency_scale_factor: float = 1.0
    #: Which printed mode the session names as the reaction coordinate,
    #: 1-based, when the structure is a characterised saddle rather than
    #: a minimum. Naming it is the scientific act; the kernel already
    #: knows how to remove it and how to treat the remaining low modes
    #: through ``entropy_method`` and the two cutoffs. Zero means the
    #: structure is claimed as a minimum and the ordinary check applies.
    reaction_coordinate_mode: int = 0
    #: Internal coordinates to remove from the Hessian before the partition
    #: functions are formed, each as the 1-based atoms ``modred`` takes (two
    #: for a bond, three for an angle, four for a dihedral).  Naming them is
    #: the request for the free energy of the surface on which they keep
    #: their values -- a point of a free-energy profile along them -- rather
    #: than of a stationary point; see ``derive_result_thermochemistry``.
    projected_coordinates: tuple[tuple[int, ...], ...] = ()
    #: Torsions to count as one-dimensional hindered rotors instead of
    #: harmonic modes (:class:`InternalRotorRequestV1`), each with the
    #: relaxed scan whose energies are its potential.
    internal_rotors: tuple[InternalRotorRequestV1, ...] = ()

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.thermochemistry-request.v1":
            raise QuantityContractError(
                "unsupported thermochemistry request schema"
            )
        _require_identifier(self.artifact_id, "artifact_id")
        _require_sha256(self.artifact_sha256)
        if int(self.reaction_coordinate_mode) < 0:
            raise QuantityContractError(
                "reaction_coordinate_mode is a 1-based mode index"
            )
        object.__setattr__(
            self,
            "projected_coordinates",
            normalized_projected_coordinates(self.projected_coordinates),
        )
        object.__setattr__(
            self,
            "internal_rotors",
            normalized_internal_rotors(self.internal_rotors),
        )
        if self.projected_coordinates and int(self.reaction_coordinate_mode):
            raise QuantityContractError(
                "name a reaction coordinate or project held coordinates, "
                "not both: projecting the coordinate a saddle moves along "
                "removes the mode reaction_coordinate_mode would name"
            )
        if self.internal_rotors and (
            self.projected_coordinates or int(self.reaction_coordinate_mode)
        ):
            raise QuantityContractError(
                "internal_rotors treats torsions of a minimum as hindered "
                "rotors; it is not served together with projected_coordinates "
                "(a held surface) or reaction_coordinate_mode (a saddle). "
                "Route: derive the rotor treatment on the minimum's own "
                "frequency result, and the held or saddle free energy in a "
                "separate request"
            )
        normalized_program = str(self.program).strip().lower()
        object.__setattr__(self, "program", normalized_program)
        from chemsmart.analysis.result_readers import reader_for

        if reader_for(normalized_program) is None:
            raise QuantityContractError(
                "no thermochemistry result reader is registered for "
                f"{normalized_program!r}"
            )
        if not math.isfinite(self.temperature_k) or self.temperature_k <= 0.0:
            raise QuantityContractError(
                "temperature_k must be finite and positive"
            )
        if not math.isfinite(self.pressure_atm) or self.pressure_atm <= 0.0:
            raise QuantityContractError(
                "pressure_atm must be finite and positive"
            )
        _validate_thermochemistry_controls(
            concentration_mol_l=self.concentration_mol_l,
            entropy_method=self.entropy_method,
            entropy_cutoff_cm1=self.entropy_cutoff_cm1,
            enthalpy_cutoff_cm1=self.enthalpy_cutoff_cm1,
            alpha=self.alpha,
            use_weighted_mass=self.use_weighted_mass,
            frequency_scale_factor=self.frequency_scale_factor,
        )
        object.__setattr__(
            self, "entropy_method", str(self.entropy_method).strip().lower()
        )
        for field in (
            "concentration_mol_l",
            "entropy_cutoff_cm1",
            "enthalpy_cutoff_cm1",
        ):
            value = getattr(self, field)
            if value is not None:
                object.__setattr__(self, field, float(value))
        object.__setattr__(
            self, "frequency_scale_factor", float(self.frequency_scale_factor)
        )


@dataclass(frozen=True)
class ThermochemistryReceiptV1:
    schema_version: str
    artifact_id: str
    artifact_sha256: str
    program: str
    engine_id: str
    temperature_k: float
    pressure_atm: float
    quantities: tuple[QuantityValueV1, ...]
    assumptions: tuple[str, ...]
    status: str
    receipt_sha256: str
    concentration_mol_l: float | None = None
    entropy_method: str = "rrho"
    #: Which printed mode the session named as the reaction coordinate,
    #: 1-based; 0 means the host removed the first genuine imaginary one.
    #: The receipt carried no record of the selection, so two different
    #: selections minted byte-identical receipts and a reader could not
    #: ask which mode a free energy had been computed without.
    entropy_cutoff_cm1: float | None = None
    enthalpy_cutoff_cm1: float | None = None
    alpha: int = 4
    use_weighted_mass: bool = False
    frequency_scale_factor: float = 1.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "quantities", tuple(self.quantities))
        object.__setattr__(self, "assumptions", tuple(self.assumptions))
        if self.schema_version != "chemsmart.thermochemistry-receipt.v1":
            raise QuantityContractError(
                "unsupported thermochemistry receipt schema"
            )
        if self.status != "derived":
            raise QuantityContractError(
                "invalid thermochemistry receipt status"
            )
        _validate_thermochemistry_controls(
            concentration_mol_l=self.concentration_mol_l,
            entropy_method=self.entropy_method,
            entropy_cutoff_cm1=self.entropy_cutoff_cm1,
            enthalpy_cutoff_cm1=self.enthalpy_cutoff_cm1,
            alpha=self.alpha,
            use_weighted_mass=self.use_weighted_mass,
            frequency_scale_factor=self.frequency_scale_factor,
        )
        legacy_body = {
            "schema_version": self.schema_version,
            "artifact_id": self.artifact_id,
            "artifact_sha256": self.artifact_sha256,
            "program": self.program,
            "engine_id": self.engine_id,
            "temperature_k": self.temperature_k,
            "pressure_atm": self.pressure_atm,
            "quantities": self.quantities,
            "assumptions": self.assumptions,
            "status": self.status,
        }
        extended_body = _extended_thermochemistry_body(legacy_body, self)
        extended_matches = self.receipt_sha256 == canonical_quantity_sha256(
            extended_body
        )
        legacy_matches = (
            self.program == "pyscf"
            and self.concentration_mol_l is None
            and self.entropy_method == "rrho"
            and self.entropy_cutoff_cm1 is None
            and self.enthalpy_cutoff_cm1 is None
            and self.alpha == 4
            and self.use_weighted_mass is False
            and self.frequency_scale_factor == 1.0
            and self.receipt_sha256 == canonical_quantity_sha256(legacy_body)
        )
        if not extended_matches and not legacy_matches:
            raise QuantityContractError(
                "thermochemistry receipt digest mismatch"
            )


def _extended_thermochemistry_body(
    legacy_body: Mapping[str, Any], controls: Any
) -> dict[str, Any]:
    """The digest body beyond the legacy contract, built in one place.

    Two organs answer this question -- the mint and the receipt's own
    revalidation -- so they call one function; they had already drifted
    once, and the drift is only visible when a stored receipt is read
    back.

    Which vibrational mode was taken as the reaction coordinate is
    deliberately *not* a key here, and not a field on the receipt. The
    ``record`` a receipt event carries **is** this body, and the event
    validator requires the two to hash alike, so a new key changes the
    digest of every receipt already written to disk and the run streams
    this laboratory has produced stop being evidence. The selection
    rides ``assumptions`` instead, beside every other control that
    changes what the numbers mean, and ``assumptions`` is inside this
    body already: no mode selected leaves the digest exactly as history
    recorded it, and two different selections still mint two different
    receipts.
    """

    body = {
        **legacy_body,
        "concentration_mol_l": controls.concentration_mol_l,
        "entropy_method": controls.entropy_method,
        "entropy_cutoff_cm1": controls.entropy_cutoff_cm1,
        "enthalpy_cutoff_cm1": controls.enthalpy_cutoff_cm1,
        "alpha": controls.alpha,
        "use_weighted_mass": controls.use_weighted_mass,
        "frequency_scale_factor": controls.frequency_scale_factor,
    }
    return body


def _validate_thermochemistry_controls(
    *,
    concentration_mol_l: float | None,
    entropy_method: str,
    entropy_cutoff_cm1: float | None,
    enthalpy_cutoff_cm1: float | None,
    alpha: int,
    use_weighted_mass: bool,
    frequency_scale_factor: float,
) -> None:
    method = str(entropy_method).strip().lower()
    if method not in {"rrho", "grimme", "truhlar"}:
        raise QuantityContractError(
            "entropy_method must be one of 'rrho', 'grimme', or 'truhlar'"
        )
    if concentration_mol_l is not None and (
        not math.isfinite(float(concentration_mol_l))
        or float(concentration_mol_l) <= 0.0
    ):
        raise QuantityContractError(
            "concentration_mol_l must be finite and positive"
        )
    for field, value in (
        ("entropy_cutoff_cm1", entropy_cutoff_cm1),
        ("enthalpy_cutoff_cm1", enthalpy_cutoff_cm1),
    ):
        if value is not None and (
            not math.isfinite(float(value)) or float(value) <= 0.0
        ):
            raise QuantityContractError(f"{field} must be finite and positive")
    if method == "rrho" and entropy_cutoff_cm1 is not None:
        raise QuantityContractError(
            "entropy_cutoff_cm1 requires entropy_method 'grimme' or 'truhlar'"
        )
    if method in {"grimme", "truhlar"} and entropy_cutoff_cm1 is None:
        raise QuantityContractError(
            f"entropy_method {method!r} requires entropy_cutoff_cm1"
        )
    if isinstance(alpha, bool) or not isinstance(alpha, int) or alpha <= 0:
        raise QuantityContractError("alpha must be a positive integer")
    if not isinstance(use_weighted_mass, bool):
        raise QuantityContractError("use_weighted_mass must be boolean")
    if (
        not math.isfinite(float(frequency_scale_factor))
        or float(frequency_scale_factor) <= 0.0
    ):
        raise QuantityContractError(
            "frequency_scale_factor must be finite and positive"
        )


def _numeric_kind(value: Any) -> str:
    array = np.asarray(value)
    if array.ndim == 0:
        return "scalar"
    if array.ndim == 1:
        return "vector"
    if array.ndim == 2:
        return "matrix"
    raise QuantityExtractionError(
        "quantities with rank greater than two are unsupported"
    )


def _make_quantity(
    *,
    quantity_id: str,
    source_value: Any,
    source_unit: str,
    value: Any,
    unit: str,
    dimension: Dimension,
    evidence_ref: str,
    data_kind: str | None = None,
) -> QuantityValueV1:
    frozen_source = _freeze(source_value)
    frozen_value = _freeze(value)
    kind = data_kind or _numeric_kind(frozen_value)
    body = {
        "schema_version": "chemsmart.quantity-value.v1",
        "quantity_id": quantity_id,
        "data_kind": kind,
        "source_value": frozen_source,
        "source_unit": source_unit,
        "value": frozen_value,
        "unit": unit,
        "dimension": dimension,
        "evidence_ref": evidence_ref,
    }
    return QuantityValueV1(
        **body, value_sha256=canonical_quantity_sha256(body)
    )


def make_quantity_value(
    *,
    quantity_id: str,
    source_value: Any,
    source_unit: str,
    value: Any,
    unit: str,
    dimension: Dimension,
    evidence_ref: str,
    data_kind: str | None = None,
) -> QuantityValueV1:
    """Build a validated immutable quantity for a deterministic derivation."""

    return _make_quantity(
        quantity_id=quantity_id,
        source_value=source_value,
        source_unit=source_unit,
        value=value,
        unit=unit,
        dimension=dimension,
        evidence_ref=evidence_ref,
        data_kind=data_kind,
    )


def _verify_artifact(
    path: str | os.PathLike[str], expected_sha256: str
) -> Path:
    artifact = Path(path).expanduser().resolve()
    if not artifact.is_file():
        raise QuantityExtractionError("trusted result artifact does not exist")
    expected = _require_sha256(expected_sha256)
    observed = result_file_sha256(artifact)
    if observed != expected:
        raise QuantityExtractionError(
            "trusted result artifact digest differs from the requested digest"
        )
    return artifact


def _require_receipt_bound_pyscf_result(
    *,
    artifact: Path,
    expected_sha256: str,
    output: PySCFOutput,
) -> dict[str, Any]:
    """Require a real, receipt-bound PySCF artifact -- green or not.

    This is the admission to *open* a result: the file states a supported
    contract, it is not a preview, and the sibling run receipt is
    digest-valid, binds these exact bytes and carries the ancestry of
    digests behind them.  It says nothing about whether the run succeeded.
    A failed optimisation opened through it is inspectable -- its terminal
    record, its last geometry, its printed numbers -- exactly as a failed
    ORCA log is, while ``_require_analysis_ready_pyscf_result`` keeps the
    green requirements every quantity extraction and free energy demands.
    """

    if (
        output.spec.get("preview_only") is not False
        or output.spec.get("result_contract_version")
        not in _SUPPORTED_PYSCF_RESULT_CONTRACTS
    ):
        raise QuantityExtractionError(
            "scientific quantities require a current executed PySCF result"
        )
    receipt_path = artifact.with_suffix(".receipt.json")
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise QuantityExtractionError(
            "scientific quantities require a readable sibling run receipt"
        ) from exc
    if not isinstance(receipt, dict):
        raise QuantityExtractionError(
            "PySCF run receipt must be a JSON object"
        )
    embedded_receipt_sha256 = receipt.get("receipt_sha256")
    receipt_body = dict(receipt)
    receipt_body.pop("receipt_sha256", None)
    if not isinstance(
        embedded_receipt_sha256, str
    ) or embedded_receipt_sha256 != canonical_sha256(receipt_body):
        raise QuantityExtractionError(
            "PySCF run receipt digest is absent or invalid"
        )
    if (
        receipt.get("fake") is not False
        or receipt.get("result_sha256") != expected_sha256
    ):
        raise QuantityExtractionError(
            "PySCF run receipt does not bind this exact result"
        )
    spec = output.spec
    provenance = output.provenance
    digest_bindings = {
        "script_sha256": (receipt, provenance, "script_sha256"),
        "input_receipt_sha256": (
            receipt,
            provenance,
            "input_receipt_sha256",
        ),
        "environment_receipt_sha256": (
            receipt,
            provenance,
            "environment_receipt_sha256",
        ),
        "input_geometry_sha256": (
            receipt,
            spec,
            "input_geometry_sha256",
        ),
        "requested_settings_sha256": (
            receipt,
            spec,
            "requested_settings_sha256",
        ),
        "project_yaml_sha256": (
            receipt,
            provenance,
            "project_yaml_digest",
        ),
        "input_artifact_sha256": (
            receipt,
            spec,
            "input_artifact_sha256",
        ),
        "applied_settings_sha256": (
            receipt,
            spec,
            "applied_settings_sha256",
        ),
    }
    for field, (
        receipt_source,
        hdf5_source,
        hdf5_field,
    ) in digest_bindings.items():
        expected = receipt_source.get(field)
        observed = hdf5_source.get(hdf5_field)
        # PySCF's CLI can receive an exact geometry value without a source
        # artifact binding.  In that normal path both records deliberately
        # carry a null ``input_artifact_sha256`` while the canonical
        # ``input_geometry_sha256`` above still binds atom order, coordinates,
        # units, charge and multiplicity.  Do not reject a matching absence as
        # if it were a broken digest.  The applied-settings digest is written
        # by the driver after the mean-field object exists, so a run that
        # died earlier carries a matching absence there too.
        if field in {"input_artifact_sha256", "applied_settings_sha256"} and (
            expected in (None, "")
        ):
            if observed in (None, ""):
                continue
        try:
            valid_digest = (
                isinstance(expected, str)
                and len(expected) == 64
                and int(expected, 16) >= 0
            )
        except ValueError:
            valid_digest = False
        if not valid_digest or observed != expected:
            raise QuantityExtractionError(
                "PySCF run receipt ancestry differs from structured result: "
                + field
            )
    for field in ("run_id", "run_nonce"):
        expected = receipt.get(field)
        if (
            not isinstance(expected, str)
            or not expected
            or spec.get(field) != expected
            or provenance.get(field) != expected
        ):
            raise QuantityExtractionError(
                "PySCF run identity differs from structured result: " + field
            )
    input_artifact_kind = receipt.get("input_artifact_kind")
    no_source_artifact = input_artifact_kind in (None, "")
    matching_absence = no_source_artifact and (
        spec.get("input_artifact_kind") in (None, "")
        and provenance.get("input_artifact_kind") in (None, "")
    )
    if not matching_absence and (
        not isinstance(input_artifact_kind, str)
        or not input_artifact_kind
        or spec.get("input_artifact_kind") != input_artifact_kind
        or provenance.get("input_artifact_kind") != input_artifact_kind
    ):
        raise QuantityExtractionError(
            "PySCF input artifact kind is absent or inconsistent"
        )
    if provenance.get("applied_settings_sha256") != receipt.get(
        "applied_settings_sha256"
    ):
        raise QuantityExtractionError(
            "PySCF applied settings digest differs across result provenance"
        )
    return receipt


def _require_analysis_ready_pyscf_result(
    *,
    artifact: Path,
    expected_sha256: str,
    output: PySCFOutput,
    required_units: dict[str, str],
) -> dict[str, Any]:
    """Require current, executed, receipt-bound, *green* result evidence.

    Fake-preview HDF5 files intentionally resemble real results so downstream
    readers can be exercised.  They are not numerical evidence.  Admission
    therefore requires both the machine contract inside HDF5 and the sibling
    ChemSmart run receipt that binds deterministic checks to these bytes,
    and every invariant green: this is what a quantity extraction and a free
    energy stand on.  A historical receipt whose Hessian the runner left
    ``unclassified`` stays admissible; the host applies the stationary-point
    promise itself.
    """

    receipt = _require_receipt_bound_pyscf_result(
        artifact=artifact, expected_sha256=expected_sha256, output=output
    )
    if (
        not output.normal_termination
        or not output.engine_complete
        or output.failure is not None
    ):
        raise QuantityExtractionError(
            "scientific quantities require a current executed PySCF result"
        )
    receipt_state = receipt.get("state")
    scientific_state = receipt.get("scientific_validation_state")
    if (
        receipt.get("engine_complete") is not True
        or receipt.get("child_returncode") != 0
        or receipt.get("findings") not in ([], ())
        or receipt_state not in {"validated", "engine_complete"}
        or scientific_state not in {"validated", "unclassified"}
    ):
        raise QuantityExtractionError(
            "PySCF run receipt does not admit this exact result for analysis"
        )
    observed_units = output.result_units
    mismatches = {
        path: {"expected": unit, "observed": observed_units.get(path)}
        for path, unit in required_units.items()
        if observed_units.get(path) != unit
    }
    if mismatches:
        raise QuantityExtractionError(
            f"PySCF result units are absent or incompatible: {mismatches}"
        )
    return receipt


def validate_pyscf_result_binding(
    artifact_path: str | os.PathLike[str],
    *,
    expected_sha256: str,
) -> tuple[PySCFOutput, dict[str, Any]]:
    """Open a receipt-bound structured result, green or failed.

    The open-time admission every reader consumer shares: the bytes are the
    receipt's bytes and the ancestry holds.  Whether the run succeeded is
    the result's own word (``normal_termination``), which extraction and
    thermochemistry check through ``validate_pyscf_analysis_artifact``.
    """

    artifact = _verify_artifact(artifact_path, expected_sha256)
    output = PySCFOutput(artifact)
    receipt = _require_receipt_bound_pyscf_result(
        artifact=artifact, expected_sha256=expected_sha256, output=output
    )
    return output, receipt


def validate_pyscf_analysis_artifact(
    artifact_path: str | os.PathLike[str],
    *,
    expected_sha256: str,
    required_units: dict[str, str] | None = None,
) -> tuple[PySCFOutput, dict[str, Any]]:
    """Validate a current structured result and its immutable run ancestry."""

    artifact = _verify_artifact(artifact_path, expected_sha256)
    output = PySCFOutput(artifact)
    receipt = _require_analysis_ready_pyscf_result(
        artifact=artifact,
        expected_sha256=expected_sha256,
        output=output,
        required_units=dict(required_units or {}),
    )
    return output, receipt


def _require_finite_numeric(value: Any, selector: str) -> Any:
    array = np.asarray(value, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise QuantityExtractionError(
            f"selector {selector!r} did not produce finite numeric data"
        )
    return value


def extract_pyscf_quantities(
    *,
    request: ResultQuantityExtractionRequestV1,
    artifact_path: str | os.PathLike[str],
) -> QuantityExtractionReceiptV1:
    """Extract selected values from a host-resolved structured PySCF result.

    PySCF now answers the shared selector vocabulary through the same
    registered reader, job-type declaration gate and unit-and-dimension
    cross-check as every log-parsing program, so this name survives only for
    callers that reach for it directly.
    """

    from chemsmart.analysis.result_readers import extract_logged_quantities

    if request.program != "pyscf":
        raise QuantityContractError(
            "extract_pyscf_quantities requires program 'pyscf'"
        )
    return extract_logged_quantities(
        request=request, artifact_path=artifact_path
    )


def _thermo_quantity(
    *,
    quantity_id: str,
    value: float,
    source_unit: str,
    normalized_unit: str,
    dimension: Dimension,
    evidence_ref: str,
) -> QuantityValueV1:
    finite_value = float(_require_finite_numeric(value, quantity_id))
    if dimension == ENERGY:
        normalized = energy_conversion("J/mol", "hartree", finite_value)
    elif dimension == ENTROPY:
        normalized = energy_conversion("J/mol", "hartree", finite_value)
    else:
        normalized = finite_value
    return _make_quantity(
        quantity_id=quantity_id,
        source_value=finite_value,
        source_unit=source_unit,
        value=normalized,
        unit=normalized_unit,
        dimension=dimension,
        evidence_ref=evidence_ref,
    )


def _thermochemistry_assumptions(
    request: ThermochemistryRequestV1,
    engine_statements: Sequence[str] = (),
) -> tuple[str, ...]:
    # PySCF receipts used to keep a shorter "legacy" assumption list at
    # default settings that never named the standard state; every program's
    # receipt now says the same things, and a PySCF receipt additionally
    # says which mass table produced its frequencies (below).
    #
    # The symmetry number used to be announced here as "derived by the
    # shared ChemSmart engine" while the engine read each program's own
    # printed value -- 1 from ORCA 6.0.1 for D-infinity-h CO2, 2 from
    # Gaussian and xTB for the same molecule. The engine now counts it and
    # says which number it used, and what the program had said, itself.
    assumptions = [
        (
            "rigid-rotor harmonic-oscillator thermochemistry for harmonic "
            "quantities, except the torsion(s) counted as hindered rotors "
            "below"
            if getattr(request, "internal_rotors", ())
            else "rigid-rotor harmonic-oscillator thermochemistry for "
            "harmonic quantities"
        ),
        "ground-state electronic degeneracy equals spin multiplicity",
        (
            "natural-abundance weighted isotopic masses"
            if request.use_weighted_mass
            else "most-abundant isotopic masses"
        ),
        *engine_statements,
        (
            "frequency scale factor 1.0; no frequency scaling"
            if request.frequency_scale_factor == 1.0
            else "vibrational frequencies multiplied by "
            f"{request.frequency_scale_factor:g} before thermochemistry"
        ),
    ]
    if request.concentration_mol_l is None:
        assumptions.append(
            f"ideal-gas translational standard state at {request.pressure_atm:g} atm"
        )
    else:
        assumptions.append(
            "solution translational standard state at "
            f"{request.concentration_mol_l:g} mol L^-1; pressure is recorded "
            "but not used in the translational partition function"
        )
    if request.entropy_method == "grimme":
        assumptions.append(
            "Grimme quasi-RRHO vibrational entropy with "
            f"{request.entropy_cutoff_cm1:g} cm^-1 cutoff and alpha "
            f"{request.alpha}"
        )
    elif request.entropy_method == "truhlar":
        assumptions.append(
            "Truhlar quasi-harmonic vibrational entropy with frequencies "
            f"below {request.entropy_cutoff_cm1:g} cm^-1 raised to the cutoff"
        )
    monoatomic = any(
        str(item).startswith("monoatomic:") for item in engine_statements
    )
    if request.program == "pyscf" and not monoatomic:
        # The frequencies came from PySCF's harmonic analysis under its
        # isotope-averaged masses while the rotational and translational
        # terms above use the table this engine names; both conventions
        # are on the receipt because the artifact records the first.
        # An atom has no frequencies, so there is no convention to state.
        assumptions.append(
            "vibrational frequencies from PySCF harmonic analysis under "
            "isotope-averaged atomic masses (artifact mass_convention)"
        )
    mode = int(getattr(request, "reaction_coordinate_mode", 0) or 0)
    if mode:
        # Which mode is the reaction coordinate is the scientist's
        # judgement; that a judgement was made, and which way, is
        # provenance the receipt owes its reader. It rides here rather
        # than in a field of its own because ``assumptions`` is already
        # inside the digest and inside the recorded record, so a
        # selection changes the receipt without invalidating every
        # receipt minted before the choice existed.
        assumptions.append(
            f"mode {mode} taken as the reaction coordinate and excluded "
            "from the vibrational partition function"
        )
    if request.enthalpy_cutoff_cm1 is not None:
        assumptions.append(
            "Head-Gordon quasi-RRHO vibrational enthalpy with "
            f"{request.enthalpy_cutoff_cm1:g} cm^-1 cutoff and alpha "
            f"{request.alpha}"
        )
    elif request.entropy_method != "rrho":
        assumptions.append(
            "quasi-harmonic correction applies to entropy only; enthalpy "
            "remains harmonic"
        )
    return tuple(assumptions)


#: Below this, a harmonic oscillator's entropy is dominated by a mode
#: the harmonic model describes worst; the observation names them.
LOW_FREQUENCY_MODE_THRESHOLD_CM1 = 50.0


def low_frequency_mode_entropy(
    frequencies_cm1: Sequence[float],
    *,
    temperature_k: float,
    threshold_cm1: float = LOW_FREQUENCY_MODE_THRESHOLD_CM1,
) -> dict[str, Any]:
    """The RRHO entropy carried by the real modes below a threshold.

    Two enantiomeric gauche rotamers, run in two sealed windows, differed
    by 0.06 kJ/mol in electronic energy and 0.40 kJ/mol in Gibbs energy,
    the difference coming entirely from harmonic-oscillator entropy of
    modes between 8 and 16 cm-1 -- on a question whose design threshold
    was 2 kJ/mol (NOVEL-1/2 po2, 2026-09-04). Per mode, with
    x = h c nu / k T:  S = R [ x / (e^x - 1) - ln(1 - e^-x) ]. An
    observation and never a verdict: the quasi-harmonic treatments the
    derivation already exposes are the scientist's to choose.
    """

    import math

    gas_constant = 8.314462618  # J mol^-1 K^-1
    second_radiation = 1.438776877  # h c / k in cm K
    temperature = float(temperature_k)
    low = sorted(
        float(value)
        for value in frequencies_cm1
        if 0.0 < float(value) < float(threshold_cm1)
    )
    entropy = 0.0
    for nu in low:
        x = second_radiation * nu / temperature
        entropy += gas_constant * (
            x / math.expm1(x) - math.log1p(-math.exp(-x))
        )
    return {
        "threshold_cm1": float(threshold_cm1),
        "temperature_k": temperature,
        "low_modes_cm1": tuple(round(value, 2) for value in low),
        "entropy_j_per_mol_k": round(entropy, 4),
        "entropy_term_kj_per_mol": round(temperature * entropy / 1000.0, 4),
    }


#: The words a stationarity reading can say, and what establishes each.
STATIONARITY_WORDS = ("stationary", "not_stationary", "unmeasured")


@dataclass(frozen=True)
class StructureStationarityV1:
    """Whether the structure a result's modes belong to is stationary.

    A harmonic free energy, like the order of a stationary point, is a
    property of a point where the gradient vanishes: the partition
    function expands the energy about it and counts every remaining
    motion as a vibration about it.  The answer is read from the result
    through its reader, in this order: one atom has no internal
    coordinate; a coordinate the result held, an atom it froze in space,
    or a coordinate it drove is not relaxed; a geometry search that printed
    its own non-convergence did not reach a stationary point; the
    program's own convergence check at the structure the modes belong to
    judges it by the criterion it was converged with, in the program's own
    coordinates (``check``); a gradient the reader binds to that structure,
    where no program judged it, is compared with the host's criterion; a
    search's own convergence marker; and a result handed its geometry,
    with none of these, is ``unmeasured`` -- neither shown nor refuted.

    A measurement outranks a label, and a criterion belongs to the
    coordinates it was applied in: Gaussian and ORCA converge forces on
    redundant internal coordinates, xTB the norm of the Cartesian
    gradient, geomeTRIC each Cartesian component, and a structure one of
    them converged can sit above another's threshold (R10 Q33: an xTB
    ``--opt loose`` result, gradient norm 1.2e-3 under its 4e-3, was called
    not stationary by its largest component against geomeTRIC's 4.5e-4).
    """

    stationarity: str
    basis: str
    program: str
    jobtype: str = ""
    max_abs_gradient_eh_per_bohr: float | None = None
    criterion_eh_per_bohr: float | None = None
    held_coordinates: int = 0
    driven_points: int = 0
    frozen_atoms: int = 0
    check: Any = None

    def __post_init__(self) -> None:
        if self.stationarity not in STATIONARITY_WORDS:
            raise QuantityContractError(
                f"stationarity is one of {list(STATIONARITY_WORDS)}"
            )

    def sentence(self) -> str:
        """What the structure is and what says so, in one clause."""

        from chemsmart.agent.terminal_states import GEOMETRY_SEARCH_JOBTYPES

        job = f"{self.program} {self.jobtype}".strip()
        noun = (
            "search" if self.jobtype in GEOMETRY_SEARCH_JOBTYPES else "result"
        )
        if self.basis == "atom":
            return "stationary point: one atom has no internal coordinate"
        if self.basis == "program_check":
            check = self.check
            stationary = self.stationarity == "stationary"
            others = check.words(force_rows=False)
            measured = (
                "; the largest Cartesian gradient component the host reads "
                f"there is {self.max_abs_gradient_eh_per_bohr:.3g} Eh/Bohr, "
                "a number on other coordinates than this criterion's"
                if self.max_abs_gradient_eh_per_bohr is not None
                else ""
            )
            return (
                (
                    "stationary point"
                    if stationary
                    else "not a stationary point"
                )
                + f": the {job} {noun}'s own convergence check at this "
                f"structure ({check.source}; {check.criterion}) gives "
                f"{check.words()}"
                + (f"; {check.verdict()}" if check.verdict() else "")
                + (f"; its other criteria: {others}" if others else "")
                + (
                    "; its predicted energy change "
                    f"{check.predicted_energy_change_eh:.2g} Eh"
                    if check.predicted_energy_change_eh is not None
                    else ""
                )
                + measured
            )
        if self.basis == "measured_gradient":
            relation = (
                "at or below" if self.stationarity == "stationary" else "above"
            )
            prefix = (
                "stationary point"
                if self.stationarity == "stationary"
                else "not a stationary point"
            )
            return (
                f"{prefix}: the largest gradient component at this "
                "result's structure is "
                f"{self.max_abs_gradient_eh_per_bohr:.3g} Eh/Bohr, "
                f"{relation} the optimiser's criterion of "
                f"{self.criterion_eh_per_bohr:g} (geomeTRIC convergence_gmax"
                + (
                    ", which PySCF's optimiser applies"
                    if self.program == "pyscf"
                    else ", the host's criterion where no optimiser judged "
                    "this structure"
                )
                + ")"
            )
        if self.basis == "held_coordinate":
            held = []
            if self.held_coordinates:
                held.append(
                    f"{self.held_coordinates} internal coordinate(s) fixed"
                )
            if self.frozen_atoms:
                held.append(f"{self.frozen_atoms} atom(s) fixed in space")
            return (
                f"not a stationary point: this {job} result held "
                + " and ".join(held)
                + " while the rest relaxed, so the energy still slopes along "
                "the held motion and its modes count that motion as a "
                "vibration"
                + (" or leave the held atoms out" if self.frozen_atoms else "")
            )
        if self.basis == "driven_coordinate":
            return (
                f"not a stationary point: this {job} result drove a "
                f"coordinate over {self.driven_points} point(s); no point "
                "of a scan is claimed stationary"
            )
        if self.basis == "search_not_converged":
            return (
                f"not a stationary point: this {job} search printed the "
                "program's own marker that it did not converge, so the "
                "structure its modes belong to is not one the program "
                "found stationary"
            )
        if self.basis == "search_converged":
            return (
                f"stationary point: the {job} search printed the "
                "program's own convergence marker, and its modes belong to "
                "the structure it reached"
            )
        if self.basis == "search_unjudged":
            return (
                f"stationarity unmeasured: this {job} search ended before "
                "its program judged any structure (it printed no "
                "convergence check and no verdict), so nothing shows or "
                "refutes that the geometry it stopped at is a stationary "
                "point of this surface; a free energy describes a state "
                "only if it is one"
            )
        return (
            f"stationarity unmeasured: this {job} result was handed its "
            "geometry and binds no gradient to it, so it neither shows "
            "nor refutes that the geometry is a stationary point of this "
            "surface; a free energy describes a state only if it is one"
        )


def _reader_answer(reader: Any, output: Any, selector: str) -> Any:
    """The reader's value for a selector, or None when it says nothing."""

    try:
        value, _unit = reader.read(output, selector)
    except Exception:  # noqa: BLE001 - an absence or an undeclared selector
        return None
    return value


def structure_stationarity(
    program: str, output: Any
) -> StructureStationarityV1:
    """Read whether the structure this result's modes belong to is stationary.

    One function, because every organ that says what a structure *is*
    asks it: a free energy (``derive_result_thermochemistry``) and the
    order of a stationary point both stand on a stationary point.  The
    characterisation asked only the gradient and the free energy asked
    nothing, so on the base of R10 Q21 a free energy was derived from a
    Gaussian ``modred`` held at HOOH = 90 deg, from an ORCA OptTS that
    printed its own non-convergence (po3-r19: a delivered 23.19 kcal/mol
    free energy of activation) and from the PySCF ``water_stretched_hess``
    Hessian at 41 times the gradient criterion -- the very artifact the
    characterisation refuses an order on.

    Every fact is read through the program's reader: the structure's
    atoms, a held, frozen or driven coordinate
    (``constrained_coordinate_count``, ``frozen_atoms_for_output``,
    ``scan_steps_planned``), a geometry search's own convergence marker
    (``converged``), the program's own convergence check at the structure
    (``convergence_check_for_output``) and the gradient the reader binds
    to that structure (``stationarity_gradient_for_output``).
    """

    from chemsmart.agent.terminal_states import (
        GEOMETRY_SEARCH_JOBTYPES,
        HESS_STATIONARITY_GRADIENT_EH_PER_BOHR,
    )
    from chemsmart.analysis.result_readers import reader_for

    normalized = str(program).strip().lower()
    reader = reader_for(normalized)
    if reader is None:
        raise QuantityContractError(
            f"no result reader is registered for {normalized!r}"
        )
    jobtype = str(getattr(output, "jobtype", "") or "").strip().lower()
    common = {"program": normalized, "jobtype": jobtype}
    symbols = _reader_answer(reader, output, "symbols")
    if symbols is not None and len(tuple(symbols)) == 1:
        return StructureStationarityV1(
            stationarity="stationary", basis="atom", **common
        )
    gradient = reader.stationarity_gradient_for_output(output)
    measured = float(f"{gradient:.6g}") if gradient is not None else None
    held = _reader_answer(reader, output, "constrained_coordinate_count")
    frozen = reader.frozen_atoms_for_output(output) or ()
    if held or frozen:
        # What a result held it was not relaxed along, whatever its
        # program's check on the rest says: that check is taken on the
        # surface the held coordinates define (Gaussian leaves a frozen
        # coordinate out of its forces and prints a zero force on a frozen
        # atom), so it cannot show the structure stationary on the full one.
        return StructureStationarityV1(
            stationarity="not_stationary",
            basis="held_coordinate",
            held_coordinates=int(held or 0),
            frozen_atoms=len(frozen),
            **common,
        )
    driven = _reader_answer(reader, output, "scan_steps_planned")
    if driven:
        return StructureStationarityV1(
            stationarity="not_stationary",
            basis="driven_coordinate",
            driven_points=int(driven),
            **common,
        )
    converged = None
    if jobtype in GEOMETRY_SEARCH_JOBTYPES:
        converged = _reader_answer(reader, output, "converged")
        if converged is None and "converged" not in reader.accessors:
            # A reader that declares no convergence selector (xTB) still
            # names the parser's own marker for the host's sensors.
            converged = getattr(output, "converged", None)
        if converged is not None and not bool(converged):
            # The search's own verdict outranks its last check: an ORCA
            # OptTS that ran out of cycles printed gradient rows within
            # tolerance at its last step and its modes where its Hessian
            # was taken, 1.23 A away (po3-r19).
            return StructureStationarityV1(
                stationarity="not_stationary",
                basis="search_not_converged",
                **common,
            )
    check = reader.convergence_check_for_output(output)
    if check is not None and (
        jobtype in GEOMETRY_SEARCH_JOBTYPES or check.at_modes
    ):
        return StructureStationarityV1(
            stationarity=(
                "stationary" if check.stationary else "not_stationary"
            ),
            basis="program_check",
            check=check,
            max_abs_gradient_eh_per_bohr=measured,
            **common,
        )
    if measured is not None:
        criterion = float(HESS_STATIONARITY_GRADIENT_EH_PER_BOHR)
        return StructureStationarityV1(
            stationarity=(
                "stationary" if measured <= criterion else "not_stationary"
            ),
            basis="measured_gradient",
            max_abs_gradient_eh_per_bohr=measured,
            criterion_eh_per_bohr=criterion,
            **common,
        )
    if converged is not None:
        return StructureStationarityV1(
            stationarity="stationary",
            basis="search_converged",
            **common,
        )
    if jobtype in GEOMETRY_SEARCH_JOBTYPES:
        return StructureStationarityV1(
            stationarity="unmeasured", basis="search_unjudged", **common
        )
    return StructureStationarityV1(
        stationarity="unmeasured", basis="fixed_geometry", **common
    )


def result_species(program: str, output: Any) -> tuple[str, Any, Any] | None:
    """``(formula, charge, multiplicity)`` of the structure a result is of.

    Read through the program's reader, as every other fact about a
    result is; the formula is Hill-ordered, and a charge or multiplicity
    the reader does not serve is None rather than a guess. None when the
    reader serves no atoms.
    """

    from chemsmart.analysis.quantity_expressions import hill_formula
    from chemsmart.analysis.result_readers import reader_for

    reader = reader_for(str(program).strip().lower())
    if reader is None:
        return None
    symbols = _reader_answer(reader, output, "symbols")
    if not symbols:
        return None

    def _integer(value: Any) -> Any:
        try:
            return int(round(float(value)))
        except (TypeError, ValueError):
            return None

    return (
        hill_formula(symbols),
        _integer(_reader_answer(reader, output, "charge")),
        _integer(_reader_answer(reader, output, "multiplicity")),
    )


#: The selectors that serve one structure's positions. Each reader declares
#: which of its molecular states (``result_readers.STRUCTURAL_STATES``) each
#: belongs to, so the geometry a number describes is read by the state its
#: own selector declares -- never assumed to be "the" structure of a result
#: (an unconverged ORCA ``OptTS Freq`` reports its energy where the search
#: stopped and its modes where its Hessian was computed, 1.23 A apart on a
#: live po3 search).
POSITION_SELECTORS = (
    "positions",
    "reached_positions",
    "supplied_positions",
    "trajectory_end_positions",
    "trajectory_start_positions",
)


def result_geometries(program: str, output: Any) -> dict[str, tuple]:
    """Each structure one result carries, by the state its reader declares.

    ``{structural state: sorted interatomic distances in Angstrom}`` for
    every position selector the reader serves on this output; a sampled
    point (``scan_point``) and a value no structure changes (``stateless``)
    name no one geometry and are left out.
    """

    from chemsmart.analysis.quantity_expressions import interatomic_distances
    from chemsmart.analysis.result_readers import reader_for

    reader = reader_for(str(program).strip().lower())
    if reader is None:
        return {}
    geometries: dict[str, tuple] = {}
    for selector in POSITION_SELECTORS:
        if selector not in reader.selectors:
            continue
        state = reader.structural_state(selector)
        if state in geometries or state in {"scan_point", "stateless"}:
            continue
        distances = interatomic_distances(
            _reader_answer(reader, output, selector)
        )
        if distances is not None:
            geometries[state] = distances
    return geometries


def geometry_of_selector(
    program: str, geometries: Mapping[str, tuple], selector: str
) -> tuple | None:
    """The geometry a selector's value describes, from ``result_geometries``.

    None when the reader declares no state for the selector, or serves no
    positions in that state: the host then does not know which structure
    the number belongs to and says so rather than guessing.
    """

    from chemsmart.analysis.result_readers import reader_for

    reader = reader_for(str(program).strip().lower())
    if reader is None:
        return None
    return geometries.get(reader.structural_state(str(selector)))


#: How far, in cm^-1, the host's own harmonic analysis of the Hessian it
#: read may sit from the frequencies the program printed before the host
#: says the matrix is not the one those modes came from.  Measured on the
#: archived H2O2 results: ORCA's ``.hess`` 0.005 (the log prints two
#: decimals), Gaussian's archive entry 0.009 (its masses are printed to
#: five decimals), PySCF's ``results/hessian`` 0.0000.
HESSIAN_REPRODUCTION_TOLERANCE_CM1 = 0.5

#: The one line of a receipt's ``assumptions`` that names the coordinates a
#: derivation projected, in the canonical form ``projected_coordinates``
#: normalises to.  Written by ``_held_coordinate_projection`` and read back
#: by ``projected_coordinates_of``, the only two places that know it.
PROJECTED_COORDINATES_STATEMENT = "projected coordinates (one-based atoms): "


def projected_coordinates_of(
    assumptions: Sequence[str],
) -> tuple[tuple[int, ...], ...]:
    """The coordinates a thermochemistry receipt projected; () if none."""

    for line in assumptions or ():
        text = str(line)
        if text.startswith(PROJECTED_COORDINATES_STATEMENT):
            return normalized_projected_coordinates(
                json.loads(text[len(PROJECTED_COORDINATES_STATEMENT) :])
            )
    return ()


#: What a held coordinate is called in a sentence, by atom count.
_COORDINATE_UNITS = {
    2: ("bond", "A"),
    3: ("angle", "deg"),
    4: ("dihedral", "deg"),
}


def normalized_projected_coordinates(
    value: Any,
) -> tuple[tuple[int, ...], ...]:
    """Coordinates to project, as sorted canonical tuples of 1-based atoms.

    Each is two atoms (a bond), three (an angle, vertex in the middle) or
    four (a dihedral about the middle pair), the rows ``modred`` takes.  A
    coordinate and its reverse are one coordinate; naming one twice is
    refused, because a direction cannot be removed twice.
    """

    if value is None:
        return ()
    if isinstance(value, (str, bytes)) or not isinstance(value, (list, tuple)):
        raise QuantityContractError(
            "projected_coordinates is a list of coordinates, each a list "
            "of 2 to 4 one-based atom indices"
        )
    coordinates = []
    for item in value:
        if isinstance(item, (str, bytes)) or not isinstance(
            item, (list, tuple)
        ):
            raise QuantityContractError(
                "each projected coordinate is a list of 2 to 4 one-based "
                f"atom indices, not {item!r}"
            )
        atoms = []
        for index in item:
            if isinstance(index, bool) or float(index) != int(float(index)):
                raise QuantityContractError(
                    f"atom index {index!r} is not an integer"
                )
            atoms.append(int(float(index)))
        if len(atoms) not in _COORDINATE_UNITS:
            raise QuantityContractError(
                f"a projected coordinate names 2, 3 or 4 atoms, not {atoms}"
            )
        if min(atoms) < 1 or len(set(atoms)) != len(atoms):
            raise QuantityContractError(
                f"projected coordinate {atoms} must name distinct one-based "
                "atoms"
            )
        coordinates.append(_canonical_coordinate(atoms))
    if len(set(coordinates)) != len(coordinates):
        raise QuantityContractError(
            "projected_coordinates names one coordinate twice"
        )
    return tuple(sorted(coordinates))


def _canonical_coordinate(atoms: Sequence[int]) -> tuple[int, ...]:
    atoms = tuple(int(index) for index in atoms)
    if len(atoms) == 2:
        return tuple(sorted(atoms))
    if len(atoms) == 3:
        return (min(atoms[0], atoms[2]), atoms[1], max(atoms[0], atoms[2]))
    return min(atoms, tuple(reversed(atoms)))


@dataclass(frozen=True)
class InternalRotorRequestV1:
    """One torsion counted as a hindered rotor, and where its potential is.

    ``torsion`` is the dihedral a-b-c-d as the one-based atoms ``modred``
    and ``scan`` take; the rotor turns about the b-c bond.  The potential
    is the relaxed scan ``scan_artifact_id`` (bytes ``scan_artifact_sha256``,
    read by the ``scan_program`` reader), which must drive a dihedral about
    the same bond over one full period of the rotor.
    """

    torsion: tuple[int, ...]
    scan_artifact_id: str
    scan_artifact_sha256: str
    scan_program: str

    def __post_init__(self) -> None:
        atoms = []
        for index in tuple(self.torsion or ()):
            if isinstance(index, bool) or float(index) != int(float(index)):
                raise QuantityContractError(
                    f"atom index {index!r} is not an integer"
                )
            atoms.append(int(float(index)))
        if len(atoms) != 4 or min(atoms) < 1 or len(set(atoms)) != 4:
            raise QuantityContractError(
                "an internal rotor names its torsion as four distinct "
                f"one-based atoms a-b-c-d (it turns about b-c), not {atoms}"
            )
        object.__setattr__(self, "torsion", tuple(atoms))
        _require_identifier(self.scan_artifact_id, "scan_artifact_id")
        _require_sha256(self.scan_artifact_sha256)
        program = str(self.scan_program).strip().lower()
        object.__setattr__(self, "scan_program", program)
        from chemsmart.analysis.result_readers import reader_for

        reader = reader_for(program)
        if reader is None or reader.resolve_torsional_scan is None:
            from chemsmart.analysis.result_readers import RESULT_READERS

            serving = sorted(
                name
                for name, item in RESULT_READERS.items()
                if item.resolve_torsional_scan is not None
            )
            raise QuantityContractError(
                f"no {program!r} reader serves a relaxed dihedral scan; "
                f"the programs whose scans are read are {serving}"
            )

    @property
    def axis(self) -> tuple[int, int]:
        """The bond turned about, one-based and sorted."""

        return tuple(sorted(self.torsion[1:3]))

    def record(self) -> dict[str, Any]:
        return {
            "torsion": list(self.torsion),
            "scan_artifact_id": self.scan_artifact_id,
            "scan_artifact_sha256": self.scan_artifact_sha256,
            "scan_program": self.scan_program,
        }


def normalized_internal_rotors(
    value: Any,
) -> tuple[InternalRotorRequestV1, ...]:
    """Internal rotors as typed requests, one per bond, in bond order."""

    if value is None:
        return ()
    if isinstance(value, (str, bytes)) or not isinstance(value, (list, tuple)):
        raise QuantityContractError(
            "internal_rotors is a list of rotors, each {torsion: [a, b, c, "
            "d], scan_artifact_id, ...}"
        )
    rotors = []
    for item in value:
        if isinstance(item, InternalRotorRequestV1):
            rotors.append(item)
            continue
        if not isinstance(item, Mapping):
            raise QuantityContractError(
                f"each internal rotor is a mapping, not {item!r}"
            )
        unknown = sorted(
            set(item)
            - {
                "torsion",
                "scan_artifact_id",
                "scan_artifact_sha256",
                "scan_program",
            }
        )
        if unknown:
            raise QuantityContractError(
                f"internal rotor fields {unknown} are not part of the request"
            )
        rotors.append(
            InternalRotorRequestV1(
                torsion=tuple(item.get("torsion") or ()),
                scan_artifact_id=str(item.get("scan_artifact_id", "")),
                scan_artifact_sha256=str(item.get("scan_artifact_sha256", "")),
                scan_program=str(item.get("scan_program", "")),
            )
        )
    axes = [rotor.axis for rotor in rotors]
    if len(set(axes)) != len(axes):
        raise QuantityContractError(
            "internal_rotors names one bond twice; a bond has one torsion"
        )
    return tuple(sorted(rotors, key=lambda rotor: rotor.axis))


#: The one line of a receipt's ``assumptions`` that names the rotors a
#: derivation treated, canonical JSON of their requests.  Written by
#: ``_internal_rotor_treatment`` and read back by ``internal_rotors_of``.
INTERNAL_ROTORS_STATEMENT = (
    "internal rotors (one-based torsion atoms, scan artifact): "
)


def internal_rotors_of(assumptions: Sequence[str]) -> tuple[dict, ...]:
    """The rotors a thermochemistry receipt treated; () if none."""

    for line in assumptions or ():
        text = str(line)
        if text.startswith(INTERNAL_ROTORS_STATEMENT):
            return tuple(json.loads(text[len(INTERNAL_ROTORS_STATEMENT) :]))
    return ()


def _coordinate_words(
    atoms: Sequence[int], symbols: Sequence[str], positions_bohr: Any
) -> str:
    """``dihedral H3-O1-O2-H4 at 90.00 deg``: one-based, measured."""

    from chemsmart.analysis.thermochemistry import internal_coordinate_value

    kind, unit = _COORDINATE_UNITS[len(atoms)]
    value = internal_coordinate_value(
        positions_bohr, [index - 1 for index in atoms]
    )
    if unit == "A":
        value *= 0.529177210903
    else:
        value = math.degrees(value)
    label = "-".join(f"{symbols[index - 1]}{index}" for index in atoms)
    return f"{kind} {label} at {value:.2f} {unit}"


#: The selectors through which a result says which coordinates it held. The
#: held-surface free energy is derived from exactly these (``_held_by_result``),
#: so a stage whose reader declares one of them is a stage the derivation
#: serves a free energy of -- the capability cell reads this table too
#: (``capabilities.coverage_for``).
HELD_COORDINATE_SELECTORS = (
    "constrained_bond_atoms",
    "constrained_angle_atoms",
    "constrained_dihedral_atoms",
)


def _held_by_result(reader: Any, output: Any) -> tuple[tuple[int, ...], ...]:
    """The coordinates this result itself held, canonical and one-based."""

    held = set()
    for selector in HELD_COORDINATE_SELECTORS:
        for row in _reader_answer(reader, output, selector) or ():
            held.add(
                _canonical_coordinate(
                    [int(round(float(index))) + 1 for index in row]
                )
            )
    return tuple(sorted(held))


@dataclass(frozen=True)
class _HeldCoordinateProjection:
    frequencies_cm1: tuple[float, ...]
    statements: tuple[str, ...]


def _torsion_rigid_turn(
    record: Any, atoms: Sequence[int], normal: Any
) -> dict[str, Any] | None:
    """The rigid turn a named dihedral stands for, where its normal is not.

    ``atoms`` is the one-based dihedral a-b-c-d and ``normal`` its
    Cartesian gradient.  None unless the dihedral turns a group about a
    bond outside any ring (``internal_rotor_tops``) and an end of that
    bond carries more than one atom off its axis: only then is the
    normal, which moves a and d alone, not the turn of the group.  The
    turn is the one ``_internal_rotor_treatment`` removes (the top's
    displacement per radian times the masses, so that its mass-weighted
    form is the top's mass-weighted turn).  Where the result records its
    gradient, the turn's surface is measured too: the gradient left once
    the turn and the rigid motions are removed, and the energy the host's
    own Hessian predicts the structure would still lose relaxing on that
    surface -- the strain one held dihedral leaves in the group it turns.
    """

    import numpy as np

    from chemsmart.analysis.thermochemistry import (
        ROTOR_AXIS_OFFSET_ANGSTROM,
        _rigid_motion_basis,
        internal_rotation_displacement,
        internal_rotor_tops,
    )

    if len(atoms) != 4:
        return None
    b, c = int(atoms[1]) - 1, int(atoms[2]) - 1
    x = np.asarray(record.positions_bohr, dtype=float)
    masses = np.asarray(record.masses_amu, dtype=float)
    symbols = [str(item) for item in record.symbols]
    angstrom = x * _BOHR_ANGSTROM
    try:
        tops = internal_rotor_tops(symbols, angstrom, (b, c))
    except ValueError:
        return None
    axis = angstrom[c] - angstrom[b]
    axis = axis / np.linalg.norm(axis)

    def off_axis(group: Sequence[int], pivot: int) -> int:
        count = 0
        for index in group:
            if index == pivot:
                continue
            relative = angstrom[index] - angstrom[pivot]
            radial = relative - (relative @ axis) * axis
            if float(np.linalg.norm(radial)) > ROTOR_AXIS_OFFSET_ANGSTROM:
                count += 1
        return count

    if max(off_axis(tops.top, c), off_axis(tops.frame, b)) < 2:
        return None
    turn = internal_rotation_displacement(x, (b, c), tops.top)
    rigid = _rigid_motion_basis(x, masses)
    root = np.sqrt(np.repeat(masses, 3))

    def unit(vector: Any) -> Any:
        vector = np.asarray(vector, dtype=float).ravel()
        vector = vector - rigid @ (rigid.T @ vector)
        return vector / np.linalg.norm(vector)

    along = unit(np.asarray(normal, dtype=float).ravel() / root)
    rotation = unit(turn.ravel() * root)
    found = {
        "direction": turn * masses[:, None],
        "overlap": float((along @ rotation) ** 2),
        "top": tuple(int(index) + 1 for index in tops.top),
        "bond": (
            f"{symbols[min(b, c)]}{min(b, c) + 1}-"
            f"{symbols[max(b, c)]}{max(b, c) + 1}"
        ),
        "off_axis": max(off_axis(tops.top, c), off_axis(tops.frame, b)),
        "residual_eh_per_bohr": None,
        "relaxation_eh": None,
    }
    if record.gradient is not None:
        size = 3 * x.shape[0]
        basis, _ = np.linalg.qr(np.column_stack([rigid, rotation]))
        projector = np.eye(size) - basis @ basis.T
        g_mw = np.asarray(record.gradient, dtype=float).ravel() / root
        left = projector @ g_mw
        hessian = np.asarray(record.hessian, dtype=float).reshape(size, size)
        weighted = 0.5 * (hessian + hessian.T) / np.outer(root, root)
        values, vectors = np.linalg.eigh(projector @ weighted @ projector)
        kept = np.abs(values) > 1e-6
        coefficients = vectors[:, kept].T @ left
        found["residual_eh_per_bohr"] = float(np.max(np.abs(left * root)))
        found["relaxation_eh"] = float(
            -0.5 * np.sum(coefficients**2 / np.abs(values[kept]))
        )
    return found


def _rigid_turn_statement(
    atoms: Sequence[int], turn: Mapping[str, Any], record: Any
) -> str:
    """One receipt clause for a dihedral removed as its group's rigid turn."""

    label = "-".join(f"{record.symbols[index - 1]}{index}" for index in atoms)
    strain = (
        "; on the surface where that turn is held the structure's "
        "gradient left is at most "
        f"{turn['residual_eh_per_bohr']:.2g} Eh/Bohr and the host's Hessian "
        f"predicts a relaxation of {turn['relaxation_eh'] * 627.509474:.2g} "
        "kcal/mol there: the strain one held dihedral leaves in the group "
        "it turns"
        if turn["residual_eh_per_bohr"] is not None
        else "; how far the structure is from stationary on the surface "
        "where that turn is held is unmeasured, since this result records "
        "no gradient at the structure"
    )
    return (
        f"dihedral {label} turns a group with {turn['off_axis']} atoms off "
        f"the {turn['bond']} bond (one-based atoms {list(turn['top'])}), so "
        f"its normal is only {turn['overlap']:.0%} that group's rigid turn "
        "and moves one atom against the others; the group's rigid turn "
        f"about {turn['bond']} -- the direction internal_rotors removes for "
        "this torsion -- is removed in its place" + strain
    )


def _projection_refusal(artifact_id: str, diagnosis: str) -> Exception:
    return QuantityExtractionError(
        "[thermochemistry.free_energy_needs_a_stationary_point] A free "
        "energy along held coordinates is the free energy of the surface "
        "on which they keep their values, and it exists where the "
        "structure is a minimum of that surface. Diagnosis: result "
        f"{artifact_id!r}: {diagnosis}"
    )


def _held_coordinate_projection(
    *,
    program: str,
    artifact_id: str,
    named: tuple[tuple[int, ...], ...],
    reader: Any,
    output: Any,
    stationarity: StructureStationarityV1,
    entropy_method: str = "rrho",
    entropy_cutoff_cm1: float | None = None,
) -> _HeldCoordinateProjection:
    """The spectrum of the surface the named coordinates are held on.

    Every fact is the result's own, read through its reader: the
    coordinates it held and whether its constrained search converged, the
    Cartesian Hessian at the structure its modes belong to (which must
    reproduce the program's printed spectrum before anything is removed
    from it), and the gradient where the result records one at that
    structure.  Refused, naming the numbers, when the named coordinates
    are not the ones the result held, when the structure is not
    stationary on their surface, when the Hessian cannot be read or is not
    the printed spectrum's, and when the structure is a saddle of the
    surface.
    """

    import numpy as np

    from chemsmart.agent.terminal_states import (
        HESS_STATIONARITY_GRADIENT_EH_PER_BOHR,
    )
    from chemsmart.analysis.quantity_expressions import (
        ONE_GEOMETRY_ANGSTROM,
        geometry_difference,
        interatomic_distances,
    )
    from chemsmart.analysis.thermochemistry import (
        NEAR_ZERO_FREQUENCY_TOLERANCE_CM,
        internal_coordinate_gradient,
        projected_harmonic_frequencies,
    )

    job = f"{program} {stationarity.jobtype}".strip()
    held = _held_by_result(reader, output)
    if held and set(named) != set(held):
        raise _projection_refusal(
            artifact_id,
            f"it held {[list(item) for item in held]} (one-based atoms) and "
            f"the request names {[list(item) for item in named]}. The "
            "structure is stationary only on the surface of what it held, "
            "and its energy still slopes along the held coordinates, so a "
            "free energy along any other set describes a surface it is not "
            "stationary on. Route: name exactly the coordinates it held.",
        )
    record = reader.cartesian_hessian_for_output(output)
    if record is None:
        raise _projection_refusal(
            artifact_id,
            f"the {program} reader serves no Cartesian Hessian for "
            "this result (an ORCA Freq run keeps it as the .hess sidecar "
            "beside its output, a Gaussian Freq job in its archive entry, a "
            "PySCF hess stage in results/hessian), and a coordinate can "
            "only be removed from a Hessian the host holds. Route: a "
            "frequency calculation on this structure with a program whose "
            "Hessian the host reads.",
        )
    atoms = len(record.symbols)
    if any(max(item) > atoms for item in named):
        raise _projection_refusal(
            artifact_id,
            f"{[list(item) for item in named]} names an atom beyond this "
            f"structure's {atoms}",
        )
    printed = sorted(
        float(value)
        for value in (
            _reader_answer(reader, output, "vibrational_frequencies") or ()
        )
        if float(value) != 0.0
    )
    rigid = projected_harmonic_frequencies(
        record.hessian, record.positions_bohr, record.masses_amu
    )
    if rigid.external != 6:
        raise _projection_refusal(
            artifact_id,
            "the structure is a linear rotor; a free energy along a held "
            "coordinate of a linear rotor is not served",
        )
    reproduced = sorted(rigid.frequencies_cm1)
    deviation = (
        max(abs(first - second) for first, second in zip(reproduced, printed))
        if len(printed) == len(reproduced) and printed
        else None
    )
    if deviation is None or deviation > HESSIAN_REPRODUCTION_TOLERANCE_CM1:
        raise _projection_refusal(
            artifact_id,
            f"the Hessian read from {record.source} does not reproduce the "
            f"{program} printed spectrum ("
            + (
                f"{len(reproduced)} modes against {len(printed)}"
                if deviation is None
                else f"largest difference {deviation:.3g} cm^-1"
            )
            + "), so it is not the matrix those modes came from",
        )
    read_structure = interatomic_distances(
        record.positions_bohr * 0.529177210903
    )
    printed_structure = interatomic_distances(
        _reader_answer(reader, output, "positions")
    )
    apart = geometry_difference(read_structure, printed_structure)
    if apart is None or apart > ONE_GEOMETRY_ANGSTROM:
        raise _projection_refusal(
            artifact_id,
            f"the Hessian read from {record.source} belongs to a structure "
            + (
                "the host could not compare with"
                if apart is None
                else f"{apart:.3g} A from"
            )
            + " the one this result's modes belong to",
        )
    directions = [
        internal_coordinate_gradient(
            record.positions_bohr, [index - 1 for index in item]
        )
        for item in named
    ]
    words = [
        _coordinate_words(item, record.symbols, record.positions_bohr)
        for item in named
    ]
    criterion = float(HESS_STATIONARITY_GRADIENT_EH_PER_BOHR)
    measured = None
    if record.gradient is not None:
        basis = np.array([direction.ravel() for direction in directions]).T
        gradient = np.asarray(record.gradient, dtype=float).ravel()
        slopes, *_ = np.linalg.lstsq(basis, gradient, rcond=None)
        residual = float(np.max(np.abs(gradient - basis @ slopes)))
        slope_words = "; ".join(
            f"dE/d({word.split(' at ')[0]}) = {slope:.3g} Eh per "
            + ("Bohr" if len(item) == 2 else "radian")
            for word, item, slope in zip(words, named, slopes)
        )
        measured = (
            "the gradient left after removing the held coordinates is at "
            f"most {residual:.2g} Eh/Bohr in Cartesian components; the energy "
            f"slopes along them ({slope_words}), which the projection removes"
        )
    # The program that held the coordinates judged the rest of the
    # structure by its own criterion, on the surface it held: Gaussian
    # leaves a frozen coordinate out of its forces, ORCA projects its
    # constraints out of its gradient.  That check, not a Cartesian
    # threshold from another optimiser, says whether the structure is
    # stationary there: a Gaussian-converged methanol held at 115 deg (max
    # internal force 4.3e-4 of 4.5e-4) was refused at a Cartesian residual
    # of 7.5e-4 against geomeTRIC's 4.5e-4, while the host's own Hessian
    # predicts it lies 1e-6 Eh from the held surface's stationary point
    # (R10 Q30 oracle O1; R10 Q33).
    check = reader.convergence_check_for_output(output) if held else None
    converged = _reader_answer(reader, output, "converged") if held else None
    if held and converged is not None and not bool(converged):
        raise _projection_refusal(
            artifact_id,
            f"this {job} search printed the program's own marker that "
            "it did not converge, so the structure is not a stationary "
            "point even of the surface it held. Route: continue the "
            "constrained search from the structure it reached.",
        )
    if check is not None:
        if not check.stationary:
            raise _projection_refusal(
                artifact_id,
                f"the {job} search that held these coordinates judged the "
                f"rest of the structure with its own check ({check.source}; "
                f"{check.criterion}) and it gives {check.words()}"
                + (f"; {check.verdict()}" if check.verdict() else "")
                + ", so the structure is not a stationary point even of the "
                "surface it held. Route: continue the constrained search "
                "from the structure it reached.",
            )
        surface = (
            "a stationary point of the held surface by the check of the "
            f"program that held it: the {job} search's own convergence "
            f"check ({check.source}; {check.criterion}) gives "
            f"{check.words()}"
            + (f"; {check.verdict()}" if check.verdict() else "")
            + (f"; {measured}" if measured else "")
            + "; not a stationary point of the full surface, whose slope "
            "along the held coordinates the projection removes"
        )
    elif held and measured is None:
        surface = (
            "a stationary point of the held surface: the "
            f"{job} search that held these coordinates "
            + (
                "printed the program's own convergence marker"
                if converged
                else "states no convergence marker, so this is unmeasured"
            )
            + "; not a stationary point of the full surface, whose slope "
            "along the held coordinates the projection removes"
        )
    elif not held and stationarity.stationarity == "stationary":
        # Stationary on the full surface by whatever judged it -- the
        # program's own check, or a gradient the host measured -- is
        # stationary on any surface through the structure: removing a
        # direction only removes gradient.  A Cartesian threshold from
        # another optimiser is not asked again here (R10 Q33: Gaussian's
        # converged ethane, largest Cartesian component 8.1e-4, was refused
        # the torsion's removal at its own minimum).
        surface = (
            f"{stationarity.sentence()}; the named coordinates are removed "
            "at that stationary point (a minimum's value is a profile's "
            "reference, a saddle's is its transition-state free energy when "
            "its imaginary mode is the removed coordinate)"
            + (f"; {measured}" if measured else "")
        )
    elif measured is not None:
        if residual > criterion:
            raise _projection_refusal(
                artifact_id,
                "the gradient left at this structure after the named "
                f"coordinates are removed is {residual:.3g} Eh/Bohr, above "
                f"the optimiser's criterion {criterion:g}, so the structure "
                "is not a stationary point of their surface. Route: hold "
                "those coordinates in a constrained optimisation (modred) "
                "and derive on its converged result.",
            )
        surface = (
            f"a stationary point of the held surface: {measured}, at or "
            f"below the optimiser's criterion {criterion:g} (geomeTRIC "
            "convergence_gmax)"
        )
    elif stationarity.stationarity == "not_stationary":
        raise _projection_refusal(
            artifact_id,
            f"it is {stationarity.sentence()}, and it held none of the "
            "named coordinates, so nothing says it is stationary on their "
            "surface either. Route: hold them in a constrained optimisation "
            "(modred) and derive on its converged result.",
        )
    else:
        surface = (
            f"{stationarity.sentence()}; the named coordinates are removed "
            "at that stationary point (a minimum's value is a profile's "
            "reference, a saddle's is its transition-state free energy when "
            "its imaginary mode is the removed coordinate)"
        )
    # A dihedral stands for the turn of the group it rotates.  Where that
    # group has one atom off the bond (H2O2's hydrogens) the dihedral's
    # normal is the group's rigid turn; where it has more, the normal moves
    # one of them against the rest, and removing it removes part of a rock
    # and keeps part of the turn (R10 Q30: ethane's CH3 rock 999.6 -> 723.5
    # cm^-1).  There the group's rigid turn -- the direction
    # ``internal_rotors`` removes for the same torsion -- is removed instead.
    removed, turns = [], []
    for item, direction in zip(named, directions):
        turn = _torsion_rigid_turn(record, item, direction)
        removed.append(direction if turn is None else turn["direction"])
        if turn is not None:
            turns.append((item, turn))
    spectrum = projected_harmonic_frequencies(
        record.hessian, record.positions_bohr, record.masses_amu, removed
    )
    imaginary = [
        value
        for value in spectrum.frequencies_cm1
        if value < -abs(NEAR_ZERO_FREQUENCY_TOLERANCE_CM)
    ]
    if imaginary:
        raise _projection_refusal(
            artifact_id,
            "with the held coordinates removed the structure still has "
            f"{len(imaginary)} imaginary mode(s) "
            f"({', '.join(f'{value:.1f}' for value in imaginary)} cm^-1): "
            "it is a saddle of the held surface, not a minimum of it",
        )
    total = 3 * atoms - spectrum.external
    projection = (
        "held-coordinate projection: "
        + "; ".join(words)
        + " (one-based atoms), removed from the mass-weighted Cartesian "
        f"Hessian together with the {spectrum.external} translations and "
        f"rotations; {spectrum.kept} of {total} vibrational modes kept "
        f"(3N-{spectrum.external + spectrum.internal}). The held "
        "coordinate carries no partition function: this is the free energy "
        "of the dividing surface at its held value, the generalized free "
        "energy of variational transition-state theory, which at a saddle "
        "whose imaginary mode is that coordinate is the transition-state "
        "free energy with that mode removed"
    )
    kind = (
        "rectilinear projection: the held coordinate's mass-weighted "
        "normal is the direction removed (Baboul & Schlegel, J. Chem. "
        "Phys. 107, 9413 (1997), Eq. 4, with that normal where they take "
        "the path tangent), so the curvature of the held surface is not "
        "included"
    )
    if turns:
        kind += "; " + "; ".join(
            _rigid_turn_statement(item, turn, record) for item, turn in turns
        )
    source = (
        f"Hessian read from {record.source}"
        + (f" (sha256 {record.source_sha256})" if record.source_sha256 else "")
        + ": its translation-rotation analysis reproduces the "
        f"{program} printed frequencies to {deviation:.2g} cm^-1 "
        f"with {record.mass_convention}, which the kept modes use too"
    )
    if entropy_method == "grimme":
        kept_modes = (
            "kept modes near and below "
            f"{entropy_cutoff_cm1:g} cm^-1 have their entropy "
            "interpolated toward a free rotor's (Grimme)"
        )
    elif entropy_method == "truhlar":
        kept_modes = (
            f"kept modes below {entropy_cutoff_cm1:g} cm^-1 are "
            "raised to it for the entropy (Truhlar)"
        )
    else:
        kept_modes = "every kept mode is a harmonic oscillator"
    rotor = (
        "rotor treatment: the held coordinate is removed, not treated as a "
        "rotor; the molecule rotates as a rigid rotor at this structure; "
        f"{kept_modes}; no kept mode is treated as a hindered internal rotor"
    )
    return _HeldCoordinateProjection(
        frequencies_cm1=tuple(spectrum.frequencies_cm1),
        statements=(
            PROJECTED_COORDINATES_STATEMENT
            + json.dumps([list(item) for item in named]),
            projection,
            kind,
            source,
            surface,
            rotor,
        ),
    )


#: The largest distance (kcal/mol) a scan point may sit from the Fourier
#: potential fitted through it: beyond it the samples are not one smooth
#: periodic function (a relaxed scan that jumped between conformers, or a
#: point whose constrained optimisation did not finish).
TORSIONAL_FIT_MAX_RESIDUAL_KCAL = 0.5

#: How far (kcal/mol) above the potential's lowest point the frequency
#: result's own torsion may sit and still be the well the rotor's levels
#: are counted from.
ROTOR_WELL_TOLERANCE_KCAL = 0.1

#: The widest gap between scan points, as a fraction of the rotor's period.
ROTOR_SCAN_MAX_GAP_FRACTION = 1.0 / 6.0

_EH_TO_CM1 = 219474.6313632
_KCAL_TO_CM1 = 349.7550882
_BOHR_ANGSTROM = 0.529177210903


@dataclass(frozen=True)
class _InternalRotorTreatment:
    frequencies_cm1: tuple[float, ...]
    rotors: tuple[Any, ...]
    statements: tuple[str, ...]
    rotor_bonds: tuple[str, ...] = ()


def _rotor_refusal(artifact_id: str, diagnosis: str) -> Exception:
    return QuantityExtractionError(
        "[thermochemistry.hindered_rotor_stands_on_its_scan] A hindered rotor "
        "is built from the "
        "frequency result's own Hessian and structure and from a relaxed "
        "scan of a dihedral about the same bond over one full period of the "
        f"rotor. Diagnosis: result {artifact_id!r}: {diagnosis}"
    )


def _normal_mode_overlap(record: Any, direction: Any) -> tuple[float, float]:
    """(overlap^2, wavenumber) of the normal mode nearest ``direction``.

    ``direction`` is mass-weighted; the translations and rotations are
    removed from it and from the Hessian first.
    """

    import numpy as np

    from chemsmart.analysis.thermochemistry import (
        HESSIAN_EIGENVALUE_TO_CM1,
        _rigid_motion_basis,
    )

    x = np.asarray(record.positions_bohr, dtype=float)
    masses = np.asarray(record.masses_amu, dtype=float)
    size = 3 * x.shape[0]
    h = np.asarray(record.hessian, dtype=float).reshape(size, size)
    inverse_root = 1.0 / np.sqrt(np.repeat(masses, 3))
    weighted = 0.5 * (h + h.T) * np.outer(inverse_root, inverse_root)
    rigid = _rigid_motion_basis(x, masses)
    projector = np.eye(size) - rigid @ rigid.T
    values, vectors = np.linalg.eigh(projector @ weighted @ projector)
    v = projector @ np.asarray(direction, dtype=float).ravel()
    v = v / np.linalg.norm(v)
    overlaps = (vectors.T @ v) ** 2
    best = int(np.argmax(overlaps))
    value = float(values[best])
    return float(overlaps[best]), float(
        np.sign(value) * math.sqrt(abs(value)) * HESSIAN_EIGENVALUE_TO_CM1
    )


#: How many torsions a harmonic receipt names, lowest mode first.
HARMONIC_TORSIONS_NAMED = 6


def harmonic_torsions(reader: Any, output: Any) -> tuple[dict[str, Any], ...]:
    """The torsions a harmonic derivation counts as oscillators.

    Every bond that is not in a ring and has atoms off its axis at both
    ends turns one part of the molecule against the rest; its rigid turn
    (mass-weighted, translations and rotations removed) is compared with
    the normal modes of the result's own Hessian, and the nearest mode is
    the harmonic oscillator that turn is counted as.  Empty when the reader
    serves no Hessian or the structure has no such bond.
    """

    import networkx as nx
    import numpy as np

    from chemsmart.analysis.thermochemistry import (
        internal_rotation_displacement,
        internal_rotor_tops,
    )
    from chemsmart.io.molecules.perception import adjacency_graph

    record = reader.cartesian_hessian_for_output(output)
    if record is None:
        return ()
    symbols = [str(item) for item in record.symbols]
    positions = np.asarray(record.positions_bohr, dtype=float) * _BOHR_ANGSTROM
    found = []
    for b, c in sorted(
        tuple(sorted(edge))
        for edge in nx.bridges(adjacency_graph(symbols, positions))
    ):
        try:
            tops = internal_rotor_tops(symbols, positions, (b, c))
        except ValueError:
            continue
        if len(tops.top) < 2 or len(tops.frame) < 2:
            continue
        displacement = internal_rotation_displacement(
            record.positions_bohr, (b, c), tops.top
        )
        try:
            overlap, mode = _normal_mode_overlap(
                record,
                displacement
                * np.sqrt(np.asarray(record.masses_amu, dtype=float))[:, None],
            )
        except (ValueError, np.linalg.LinAlgError):
            continue
        found.append(
            {
                "bond": f"{symbols[b]}{b + 1}-{symbols[c]}{c + 1}",
                "mode_cm1": mode,
                "overlap": overlap,
                "sigma_int": tops.symmetry_number,
            }
        )
    return tuple(sorted(found, key=lambda item: item["mode_cm1"]))


def _harmonic_torsion_statement(
    reader: Any, output: Any, rotor_bonds: Sequence[str] = ()
) -> str:
    """One receipt line naming the torsions counted as oscillators, or ''.

    ``rotor_bonds`` are bonds already treated as hindered rotors, which the
    line leaves out.
    """

    try:
        torsions = tuple(
            item
            for item in harmonic_torsions(reader, output)
            if item["bond"] not in set(rotor_bonds)
        )
    except Exception:  # noqa: BLE001 - a statement never fails a derivation
        return ""
    if not torsions:
        return ""
    named = [
        f"about {item['bond']} is {item['overlap']:.0%} the "
        f"{item['mode_cm1']:.1f} cm^-1 mode"
        for item in torsions[:HARMONIC_TORSIONS_NAMED]
    ]
    more = len(torsions) - len(named)
    return (
        "torsions counted as harmonic oscillators: the rigid turn "
        + "; ".join(named)
        + (f"; and {more} more" if more > 0 else "")
        + ". A torsion counted as a one-dimensional hindered rotor instead "
        "is internal_rotors, with a relaxed scan of a dihedral about its "
        "bond over one full period of the rotor"
    )


def _internal_rotor_treatment(
    *,
    program: str,
    artifact_id: str,
    rotors: tuple[InternalRotorRequestV1, ...],
    rotor_artifact_paths: Mapping[str, Any] | None,
    reader: Any,
    output: Any,
    temperature_k: float,
    use_weighted_mass: bool,
) -> _InternalRotorTreatment:
    """Kept modes, hindered rotors and receipt statements for ``rotors``.

    Every fact is read through a reader: the frequency result's Cartesian
    Hessian and structure (which must reproduce its printed spectrum), and
    each rotor's relaxed scan (which must be the same molecule in the same
    atom order, drive a dihedral about the rotor's bond, and cover one full
    period of the rotor).  Each rotor's rigid turn is projected from the
    Hessian, so the harmonic mode it replaces is not counted beside it.
    """

    import numpy as np

    from chemsmart.analysis.result_readers import reader_for, surfaces_agree
    from chemsmart.analysis.thermochemistry import (
        NEAR_ZERO_FREQUENCY_TOLERANCE_CM,
        fit_torsional_potential,
        hindered_rotor,
        internal_coordinate_value,
        internal_rotation_displacement,
        internal_rotation_moment,
        internal_rotor_tops,
        projected_harmonic_frequencies,
    )
    from chemsmart.utils.periodictable import PeriodicTable

    record = reader.cartesian_hessian_for_output(output)
    if record is None:
        raise _rotor_refusal(
            artifact_id,
            f"the {program} reader serves no Cartesian Hessian for this "
            "result (an ORCA Freq run keeps it as the .hess sidecar, a "
            "Gaussian Freq job in its archive entry, a PySCF hess stage in "
            "results/hessian), and a rotor's harmonic mode can only be "
            "removed from a Hessian the host holds. Route: a frequency "
            "calculation on this minimum with a program whose Hessian the "
            "host reads.",
        )
    printed = sorted(
        float(value)
        for value in (
            _reader_answer(reader, output, "vibrational_frequencies") or ()
        )
        if float(value) != 0.0
    )
    rigid = projected_harmonic_frequencies(
        record.hessian, record.positions_bohr, record.masses_amu
    )
    if rigid.external != 6:
        raise _rotor_refusal(
            artifact_id,
            "the structure is a linear rotor, which has no torsion",
        )
    reproduced = sorted(rigid.frequencies_cm1)
    deviation = (
        max(abs(first - second) for first, second in zip(reproduced, printed))
        if len(printed) == len(reproduced) and printed
        else None
    )
    if deviation is None or deviation > HESSIAN_REPRODUCTION_TOLERANCE_CM1:
        raise _rotor_refusal(
            artifact_id,
            f"the Hessian read from {record.source} does not reproduce the "
            f"{program} printed spectrum, so it is not the matrix those "
            "modes came from",
        )
    symbols = [str(item) for item in record.symbols]
    positions = np.asarray(record.positions_bohr, dtype=float) * _BOHR_ANGSTROM
    table = PeriodicTable()
    masses = np.array(
        [
            (
                table.to_weighted_atomic_mass_by_abundance(symbol)
                if use_weighted_mass
                else table.to_most_abundant_atomic_mass(symbol)
            )
            for symbol in symbols
        ],
        dtype=float,
    )
    state = (
        _reader_answer(reader, output, "charge"),
        _reader_answer(reader, output, "multiplicity"),
    )
    energy_eh = _reader_answer(reader, output, "energy")
    surface = reader.surface_for_output(output)
    built, directions, statements, records, bonds = [], [], [], [], []
    for number, rotor in enumerate(rotors, start=1):
        if max(rotor.torsion) > len(symbols):
            raise _rotor_refusal(
                artifact_id,
                f"torsion {list(rotor.torsion)} names an atom beyond this "
                f"structure's {len(symbols)}",
            )
        a, b, c, d = (index - 1 for index in rotor.torsion)
        try:
            tops = internal_rotor_tops(symbols, positions, (b, c))
        except ValueError as exc:
            raise _rotor_refusal(artifact_id, str(exc)) from exc
        sigma = tops.symmetry_number
        moment = internal_rotation_moment(positions, masses, (b, c), tops.top)
        moment_frame = internal_rotation_moment(
            positions, masses, (b, c), tops.frame
        )
        path = (rotor_artifact_paths or {}).get(rotor.scan_artifact_id)
        if path is None:
            raise _rotor_refusal(
                artifact_id,
                f"the scan {rotor.scan_artifact_id!r} named for rotor "
                f"{number} is not a result bound to this request",
            )
        scan_path = _verify_artifact(path, rotor.scan_artifact_sha256)
        scan_reader = reader_for(rotor.scan_program)
        scan_output = scan_reader.open_output(scan_path)
        scan = scan_reader.torsional_scan_for_output(scan_output)
        if scan is None:
            raise _rotor_refusal(
                artifact_id,
                f"the {rotor.scan_program} reader reads no relaxed scan of "
                f"one dihedral from {rotor.scan_artifact_id!r}. Route: a "
                "relaxed scan (scan jobtype) of a dihedral about the rotor's "
                "bond over one full period of the rotor.",
            )
        if list(scan.symbols) != symbols:
            raise _rotor_refusal(
                artifact_id,
                f"the scan {rotor.scan_artifact_id!r} is not this molecule in "
                "this atom order, so its energies are not this torsion's",
            )
        if set(scan.atoms[1:3]) != {b + 1, c + 1}:
            raise _rotor_refusal(
                artifact_id,
                f"the scan drove dihedral {list(scan.atoms)} about the "
                f"{scan.atoms[1]}-{scan.atoms[2]} bond, not about the rotor's "
                f"{b + 1}-{c + 1} bond",
            )
        scan_state = (
            _reader_answer(scan_reader, scan_output, "charge"),
            _reader_answer(scan_reader, scan_output, "multiplicity"),
        )
        if None not in state and None not in scan_state:
            if tuple(round(float(v)) for v in state) != tuple(
                round(float(v)) for v in scan_state
            ):
                raise _rotor_refusal(
                    artifact_id,
                    f"the scan's charge and multiplicity {scan_state} are "
                    f"not the frequency result's {state}",
                )
        period = 2.0 * math.pi / sigma
        phi = np.radians(np.asarray(scan.values_deg, dtype=float))
        reduced = np.sort(np.mod(phi, period))
        gaps = np.diff(np.concatenate([reduced, [reduced[0] + period]]))
        if float(gaps.max()) > ROTOR_SCAN_MAX_GAP_FRACTION * period + 1e-9:
            raise _rotor_refusal(
                artifact_id,
                f"the scan's {len(phi)} point(s) leave a gap of "
                f"{math.degrees(float(gaps.max())):.1f} deg in the "
                f"{math.degrees(period):g}-deg period of this rotor "
                f"(sigma_int {sigma}); the potential is fitted only where "
                "every gap is at most "
                f"{math.degrees(ROTOR_SCAN_MAX_GAP_FRACTION * period):.1f} "
                "deg. Route: a relaxed scan over one full period of the "
                "rotor, starting at the structure's own value and offset "
                "from a planar 0 or 180 deg point.",
            )
        energies = np.asarray(scan.energies_eh, dtype=float)
        zero = float(energies.min())
        samples = (energies - zero) * _EH_TO_CM1
        potential = fit_torsional_potential(phi, samples, period)
        residual = samples - potential.value(phi)
        largest = float(np.abs(residual).max()) / _KCAL_TO_CM1
        if largest > TORSIONAL_FIT_MAX_RESIDUAL_KCAL:
            raise _rotor_refusal(
                artifact_id,
                f"a scan point sits {largest:.2f} kcal/mol from the Fourier "
                f"potential fitted through the scan (at most "
                f"{TORSIONAL_FIT_MAX_RESIDUAL_KCAL:g}), so its points are not "
                "one smooth periodic surface: a relaxed scan that jumped "
                "between conformers, or a point that did not finish its "
                "constrained optimisation",
            )
        phi_min, v_min = potential.minimum
        phi_eq = internal_coordinate_value(
            positions, [index - 1 for index in scan.atoms]
        )
        above = (float(potential.value(phi_eq)) - v_min) / _KCAL_TO_CM1
        if above > ROTOR_WELL_TOLERANCE_KCAL:
            raise _rotor_refusal(
                artifact_id,
                "the frequency result's own dihedral "
                f"{math.degrees(phi_eq):.1f} deg sits {above:.2f} kcal/mol "
                "above the scan's lowest point (at "
                f"{math.degrees(phi_min):.1f} deg), so the rotor's levels "
                "would not be counted from the well this result's modes "
                "belong to. Route: derive on the frequency result at the "
                "lowest well of this torsion.",
            )
        rotor_levels = hindered_rotor(potential, moment, sigma)
        displacement = internal_rotation_displacement(
            record.positions_bohr, (b, c), tops.top
        )
        directions.append(
            displacement * np.asarray(record.masses_amu, dtype=float)[:, None]
        )
        overlap, mode = _normal_mode_overlap(
            record,
            displacement
            * np.sqrt(np.asarray(record.masses_amu, dtype=float))[:, None],
        )
        built.append(rotor_levels)
        records.append(rotor.record())
        low, high = sorted((b, c))
        bonds.append(f"{symbols[low]}{low + 1}-{symbols[high]}{high + 1}")
        minima, maxima = potential.stationary_points()
        barriers = ", ".join(
            f"{(value - v_min) / _KCAL_TO_CM1:.3f} at "
            f"{math.degrees(where):.1f}"
            for where, value in maxima
        )
        label = "-".join(f"{symbols[i]}{i + 1}" for i in (a, b, c, d))
        scan_label = "-".join(
            f"{symbols[i - 1]}{i}" for i in (int(v) for v in scan.atoms)
        )
        agreement = surfaces_agree(
            surface, scan_reader.surface_for_output(scan_output)
        )
        surface_words = {
            True: "on the frequency result's surface",
            False: "on a different surface from the frequency result's "
            "(the readers' surface identities differ)",
            None: "surface identity not comparable",
        }[agreement]
        if energy_eh is not None:
            # The one cross-check both readers can always make: the
            # potential's energy at this result's own dihedral against the
            # result's energy -- about zero on one surface, a basis or a
            # functional apart otherwise.
            offset = (
                zero
                + float(potential.value(phi_eq)) / _EH_TO_CM1
                - float(energy_eh)
            )
            surface_words += (
                "; the scan's energy at this result's dihedral is "
                f"{offset * 627.509474:+.4f} kcal/mol from the result's own"
            )
        s_rotor = rotor_levels.thermodynamics(temperature_k)
        s_harm = rotor_levels.harmonic_thermodynamics(temperature_k)
        statements.extend(
            [
                (
                    f"hindered rotor {number}: dihedral {label} (one-based "
                    f"atoms) turned about {symbols[b]}{b + 1}-{symbols[c]}"
                    f"{c + 1}, at {math.degrees(internal_coordinate_value(positions, [a, b, c, d])):.2f} "
                    f"deg in this result's structure; top "
                    f"{[i + 1 for i in tops.top]} against "
                    f"{[i + 1 for i in tops.frame]}"
                ),
                (
                    f"rotor {number} potential: {len(phi)} points of the "
                    f"relaxed scan {rotor.scan_artifact_id!r} "
                    f"({rotor.scan_program}; {scan.source}; dihedral "
                    f"{scan_label} "
                    + (
                        "measured in each point's structure"
                        if scan.measured
                        else "as the targets the program held"
                    )
                    + f"), {surface_words}; Fourier series of order "
                    f"{potential.order} over the "
                    f"{math.degrees(period):g}-deg period (rms residual "
                    f"{potential.rms_residual_cm1 / _KCAL_TO_CM1:.3f}, "
                    f"largest {largest:.3f} kcal/mol); {len(minima)} well(s) "
                    f"per period, the lowest at {math.degrees(phi_min):.1f} "
                    f"deg; barrier(s) in kcal/mol at deg: {barriers}"
                ),
                (
                    f"rotor {number} moment: I(3,4) = {moment:.4f} amu A^2 "
                    f"(B = {rotor_levels.rotational_constant_cm1:.3f} "
                    "cm^-1), the rigid top turned about the bond at zero "
                    "overall angular momentum (East & Radom, J. Chem. Phys. "
                    "106, 6655 (1997)); turning the other end gives "
                    f"{moment_frame:.4f}; "
                    + (
                        "natural-abundance"
                        if use_weighted_mass
                        else "most-abundant"
                    )
                    + " isotopic masses, at this result's structure"
                ),
                (
                    f"rotor {number} symmetry: sigma_int {sigma}, the lcm of "
                    "the two ends' rotational orders about the bond "
                    f"({tops.top_order} and {tops.frame_order}; largest "
                    "image miss "
                    f"{max(tops.top_deviation, tops.frame_deviation):.3f} A); "
                    "every level over one full turn is summed and divided by "
                    "sigma_int, beside the external symmetry number"
                ),
                (
                    f"rotor {number} replaces a harmonic mode: its rigid turn "
                    "is projected from the mass-weighted Hessian with the "
                    "translations and rotations, so the harmonic mode is not "
                    f"counted beside the rotor; the removed direction is "
                    f"{overlap:.0%} the {mode:.1f} cm^-1 normal mode"
                ),
                (
                    f"rotor {number} levels: -B d2/dphi2 + V(phi) in a "
                    f"free-rotor basis |m| <= {rotor_levels.basis}; lowest "
                    f"level {rotor_levels.zero_point_cm1:.1f} cm^-1 above the "
                    f"potential minimum; at {float(temperature_k):g} K "
                    f"S = {s_rotor[1]:.3f} and Cv = {s_rotor[3]:.3f} "
                    "J/(K mol), where the harmonic oscillator of the same "
                    "moment and curvature "
                    f"({rotor_levels.harmonic_frequency_cm1:.1f} cm^-1) "
                    f"gives S = {s_harm[0]:.3f} and Cv = {s_harm[2]:.3f}"
                ),
            ]
        )
    spectrum = projected_harmonic_frequencies(
        record.hessian, record.positions_bohr, record.masses_amu, directions
    )
    imaginary = [
        value
        for value in spectrum.frequencies_cm1
        if value < -abs(NEAR_ZERO_FREQUENCY_TOLERANCE_CM)
    ]
    if imaginary:
        raise _rotor_refusal(
            artifact_id,
            f"with the rotors' turns removed the Hessian still has "
            f"{len(imaginary)} imaginary mode(s) "
            f"({', '.join(f'{value:.1f}' for value in imaginary)} cm^-1): "
            "the structure is not a minimum",
        )
    total = 3 * len(symbols) - spectrum.external
    statements.append(
        f"{spectrum.kept} of {total} vibrational modes kept beside "
        f"{len(built)} hindered rotor(s) (3N-{spectrum.external}-"
        f"{len(built)}); Hessian read from {record.source}"
        + (f" (sha256 {record.source_sha256})" if record.source_sha256 else "")
        + f", reproducing the {program} printed frequencies to "
        f"{deviation:.2g} cm^-1 with {record.mass_convention}; rotors are "
        "independent one-dimensional rotors, with no rotor-rotor coupling"
    )
    return _InternalRotorTreatment(
        frequencies_cm1=tuple(spectrum.frequencies_cm1),
        rotors=tuple(built),
        rotor_bonds=tuple(bonds),
        statements=(
            INTERNAL_ROTORS_STATEMENT
            + json.dumps(records, sort_keys=True, separators=(",", ":")),
            *statements,
        ),
    )


#: What a result can have a free energy of.
FREE_ENERGY_SURFACES = ("stationary_point", "held_surface", "none")


@dataclass(frozen=True)
class FreeEnergySurfaceV1:
    """Whether a result has a free energy, and of which surface.

    ``stationary_point``: the structure is one, or nothing shows it is not
    (``stationarity`` says which); the ordinary derivation applies.
    ``held_surface``: the structure is not a stationary point because it
    held ``held_coordinates``, and the host derives the free energy of the
    surface they are held on (``projected_coordinates``) from this result.
    ``none``: shown not to be a stationary point of any surface the host
    can derive a free energy on; ``reason`` says what was read.
    """

    surface: str
    stationarity: StructureStationarityV1
    held_coordinates: tuple[tuple[int, ...], ...] = ()
    reason: str = ""

    def __post_init__(self) -> None:
        if self.surface not in FREE_ENERGY_SURFACES:
            raise QuantityContractError(
                f"surface is one of {list(FREE_ENERGY_SURFACES)}"
            )

    def route(self) -> str:
        """The request that derives it, for a held surface; else ''."""

        if self.surface != "held_surface":
            return ""
        return (
            "derive_thermochemistry with projected_coordinates "
            f"{[list(item) for item in self.held_coordinates]} (the "
            "one-based atoms it held)"
        )


def free_energy_surface(
    program: str, output: Any, *, artifact_id: str = "result"
) -> FreeEnergySurfaceV1:
    """Is there a free energy at this result, and of which surface.

    One function, because three organs ask it and they had drifted apart:
    the derivation, which with ``projected_coordinates`` derives the free
    energy of a held surface; the check that verifies a session's refusal
    of a free energy, which went on signing "absent" over such a result
    (R10 Q27, the held 90-deg H2O2 of goal g1, CUHK 2153714: verified while
    the same host derived 0.346 kcal/mol from it); and the refusal's own
    route.  A held result has a held-surface free energy exactly when the
    derivation would give one: the same checks run here with the
    coordinates the result itself held.
    """

    normalized = str(program).strip().lower()
    stationarity = structure_stationarity(normalized, output)
    if stationarity.stationarity != "not_stationary":
        return FreeEnergySurfaceV1(
            surface="stationary_point", stationarity=stationarity
        )
    from chemsmart.analysis.result_readers import reader_for

    reader = reader_for(normalized)
    held = _held_by_result(reader, output)
    if not held:
        return FreeEnergySurfaceV1(
            surface="none",
            stationarity=stationarity,
            reason=stationarity.sentence(),
        )
    try:
        _held_coordinate_projection(
            program=normalized,
            artifact_id=artifact_id,
            named=held,
            reader=reader,
            output=output,
            stationarity=stationarity,
        )
    except QuantityContractError as exc:
        return FreeEnergySurfaceV1(
            surface="none",
            stationarity=stationarity,
            held_coordinates=held,
            reason=f"{stationarity.sentence()}; and no free energy of the "
            f"surface it held is derivable: {exc}",
        )
    return FreeEnergySurfaceV1(
        surface="held_surface",
        stationarity=stationarity,
        held_coordinates=held,
        reason=stationarity.sentence(),
    )


def _modes_the_program_removed(engine: Any) -> str:
    """A sentence when the printed spectrum is short of the structure's modes.

    A harmonic partition function counts every vibration a structure has --
    3N-6, or 3N-5 for a linear rotor -- unless the host removed one itself
    and says which (``projected_coordinates``).  A program can remove one
    before it prints: Gaussian's ``freq=projected`` drops the direction of
    the gradient, which at H2O2 held at 0 or 180 deg (a symmetric saddle,
    gradient orthogonal to the torsion) is a stretch-bend mixture, not the
    held torsion, and Gaussian's own free energy there came out 5.4 and 1.7
    kcal/mol below the saddles' (R10 Q27 oracle O1b, CUHK 2153717); a
    Gaussian optimisation with frozen atoms prints the modes of the others
    only (12 of a 14-atom structure's 36 in an archived fixture).  Such a
    spectrum is used as printed, and the receipt says what it lacks.
    """

    frequencies = engine.vibrational_frequencies
    if frequencies is None or engine.molecule.is_monoatomic:
        return ""
    atoms = int(engine.molecule.num_atoms)
    expected = 3 * atoms - (5 if engine.is_linear_rotor else 6)
    printed = len(frequencies)
    if engine.quasi_linear_padded_mode_cm1 is not None:
        printed += 1
    if printed >= expected:
        return ""
    return (
        f"the program printed {printed} vibrational mode(s) where this "
        f"{atoms}-atom structure has {expected}: it removed "
        f"{expected - printed} direction(s) itself before its own analysis "
        "(Gaussian's freq=projected removes the direction of the gradient; "
        "frozen atoms take their own motions out of the analysis), so this "
        "free energy counts only the printed modes and lacks "
        "motion(s) the host neither chose nor can name; removing a named "
        "held coordinate is projected_coordinates"
    )


def _stationarity_refusal(
    surface: FreeEnergySurfaceV1, artifact_id: str
) -> QuantityExtractionError:
    """The routed refusal of a free energy at a non-stationary structure."""

    stationarity = surface.stationarity
    held_route = ""
    if surface.surface == "held_surface":
        held_route = (
            " For the free energy of the structure at the value it held -- a "
            "point of a free-energy profile along the held coordinate -- "
            f"call {surface.route()}: the host removes them from the "
            "Hessian and derives the free energy of the 3N-7 modes of the "
            "surface they are held on, the generalized free energy of "
            "variational transition-state theory."
        )
    return QuantityExtractionError(
        "[thermochemistry.free_energy_needs_a_stationary_point] A free "
        "energy, an enthalpy or a zero-point energy from harmonic modes is "
        "a property of a stationary point: the partition function expands "
        "the energy about a point where the gradient vanishes and counts "
        "every other motion as a vibration about it. Diagnosis: result "
        f"{artifact_id!r} is {stationarity.sentence()}. Route: its "
        "electronic energy, its frequencies and every other number stay "
        "readable and deliverable as what they are -- the energy and the "
        "curvature at a structure that is not stationary (a held or driven "
        "coordinate's energies are points on a constrained surface)."
        + held_route
        + " For the free energy of a stationary point, reach a stationary "
        "point of the same surface -- relax without the constraint (opt for a minimum, "
        "ts for a saddle; a structure held at a symmetric point often "
        "converges in a few steps) or continue the search from the "
        "structure this one reached -- and derive thermochemistry on that "
        "result. Cost: one engine call to reach a stationary point."
    )


def derive_result_thermochemistry(
    *,
    request: ThermochemistryRequestV1,
    artifact_path: str | os.PathLike[str],
    rotor_artifact_paths: Mapping[str, str | os.PathLike[str]] | None = None,
) -> ThermochemistryReceiptV1:
    """Derive RRHO or quasi-harmonic thermochemistry from a trusted result.

    The formulas and molecular conventions remain owned by ChemSmart's common
    :class:`Thermochemistry` engine.  This function only binds conditions and
    serializes the resulting values with explicit units and provenance.
    ``rotor_artifact_paths`` maps each ``internal_rotors`` scan artifact id
    to its host-resolved path.
    """

    artifact = _verify_artifact(artifact_path, request.artifact_sha256)
    if request.program == "pyscf":
        output = PySCFOutput(artifact)
        # An atom has no vibration to take a Hessian of: its partition
        # function is translational and electronic, from the energy, the
        # position and the multiplicity a single point already carries
        # (the engine's own rule, analysis/thermochemistry.py).
        monoatomic = len(tuple(output.chemical_symbols or ())) == 1
        required_units = {
            "results/energies": "Eh",
            "results/positions": "Angstrom",
        }
        if not monoatomic:
            required_units.update(
                {
                    "results/hessian": "Eh/Bohr^2",
                    "results/vibrational_frequencies": "cm^-1",
                }
            )
        _require_analysis_ready_pyscf_result(
            artifact=artifact,
            expected_sha256=request.artifact_sha256,
            output=output,
            required_units=required_units,
        )
        if not output.freq and not monoatomic:
            raise QuantityExtractionError(
                "thermochemistry requires a validated PySCF Hessian result"
            )
        if output.result_sha256 != request.artifact_sha256:
            raise QuantityExtractionError(
                "PySCF parser observed substituted bytes"
            )
    # Whether the structure these modes belong to is a stationary point,
    # asked before anything is derived: a free energy is a property of a
    # stationary point, so a structure that is shown not to be one has none
    # to give, whatever its modes are -- and a structure nobody measured
    # says so in the receipt rather than passing as one.
    from chemsmart.analysis.result_readers import reader_for

    reader = reader_for(request.program)
    reached = reader.open_output(artifact)
    stationarity = structure_stationarity(request.program, reached)
    # A free energy along held coordinates is asked for by naming them:
    # the host removes them from the Hessian and the free energy is that
    # of the surface they are held on, which the structure must be a
    # minimum of (R10 Q27).  Without the request the stationary-point rule
    # below is untouched.
    projection = (
        _held_coordinate_projection(
            program=request.program,
            artifact_id=request.artifact_id,
            named=request.projected_coordinates,
            reader=reader,
            output=reached,
            stationarity=stationarity,
            entropy_method=request.entropy_method,
            entropy_cutoff_cm1=request.entropy_cutoff_cm1,
        )
        if request.projected_coordinates
        else None
    )
    # A result that printed no modes at all has nothing to derive from,
    # and says so below in its own words; the stationarity refusal is for
    # a structure whose modes exist and describe no stationary point. Its
    # route is the one ``free_energy_surface`` -- the function the refusal
    # verification also reads -- says this result has.
    if (
        projection is None
        and stationarity.stationarity == "not_stationary"
        and _reader_answer(reader, reached, "vibrational_frequencies")
    ):
        raise _stationarity_refusal(
            free_energy_surface(
                request.program, reached, artifact_id=request.artifact_id
            ),
            request.artifact_id,
        )
    # Torsions counted as hindered rotors instead of harmonic modes: each
    # rotor's rigid turn leaves the Hessian and its levels on the scan's
    # potential take its place.  The stationary-point rule above has
    # already been asked; a rotor is a treatment of a minimum's torsion.
    rotor_treatment = (
        _internal_rotor_treatment(
            program=request.program,
            artifact_id=request.artifact_id,
            rotors=request.internal_rotors,
            rotor_artifact_paths=rotor_artifact_paths,
            reader=reader,
            output=reached,
            temperature_k=request.temperature_k,
            use_weighted_mass=request.use_weighted_mass,
        )
        if request.internal_rotors
        else None
    )
    kept_frequencies = (
        projection.frequencies_cm1
        if projection is not None
        else (
            rotor_treatment.frequencies_cm1
            if rotor_treatment is not None
            else None
        )
    )
    engine = Thermochemistry(
        projected_frequencies=kept_frequencies,
        internal_rotors=(
            rotor_treatment.rotors if rotor_treatment is not None else ()
        ),
        filename=str(artifact),
        temperature=request.temperature_k,
        concentration=request.concentration_mol_l,
        pressure=request.pressure_atm,
        use_weighted_mass=request.use_weighted_mass,
        alpha=request.alpha,
        s_freq_cutoff=request.entropy_cutoff_cm1,
        entropy_method=(
            None
            if request.entropy_method == "rrho"
            else request.entropy_method
        ),
        h_freq_cutoff=request.enthalpy_cutoff_cm1,
        frequency_scale_factor=request.frequency_scale_factor,
        # A session that names the reaction coordinate has said what the
        # structure is; the kernel then removes that mode and treats the
        # rest through `entropy_method` and the two cutoffs, which is
        # machinery it has always had and nothing could reach.
        #
        # The refusal this replaces rejected every free energy of
        # activation -- the commonest thermochemistry in mechanism work
        # -- and prevented nothing: a live session computed
        # delta-E + delta-ZPE by hand instead, changing the rigour of a
        # delivered number without changing whether it was delivered.
        # The charter says outright that a wrong-stationary-point
        # result's "energy, its modes and its thermochemistry are
        # exactly where a finding lives". It is the same class of
        # refusal f9c853d3 removed elsewhere (SUFFICIENCY-5,
        # 2026-09-10). Naming the mode is a scientific act, not a
        # permission bit: the default still refuses, because a saddle
        # nobody has characterised is a failed optimisation.
        check_imaginary_frequencies=not int(
            getattr(request, "reaction_coordinate_mode", 0) or 0
        ),
        # The selection itself, not merely whether one was made. The
        # index stopped here and the engine removed the first imaginary
        # mode whatever the session named.
        reaction_coordinate_mode=int(
            getattr(request, "reaction_coordinate_mode", 0) or 0
        ),
    )
    if engine.program != request.program:
        raise QuantityExtractionError(
            "trusted result program differs from the requested thermochemistry "
            f"program: expected {request.program!r}, observed {engine.program!r}"
        )
    missing_modes = (
        _modes_the_program_removed(engine)
        if projection is None and rotor_treatment is None
        else ""
    )
    # `check_frequencies` keys on the program job label: one imaginary
    # mode is a correct transition state for a `ts` job and a refusal
    # for the identical structure labelled `opt`. The label is the
    # arbitrary part -- a saddle is a saddle whatever the input asked
    # for -- and this refusal rejected every free energy of activation
    # while preventing nothing, since a session can compute
    # delta-E + delta-ZPE by hand and one did.
    #
    # A named reaction coordinate replaces the label with something
    # stronger: the agent layer admits it only over a stationary-point
    # characterisation this host minted, whose order it checked against
    # the program's own printed frequencies. The human CLI's validator
    # is untouched (SUFFICIENCY-5, 2026-09-10).
    if not int(getattr(request, "reaction_coordinate_mode", 0) or 0):
        engine.check_frequencies()
    # A geometry optimisation that hit its iteration cap never reached
    # its frequency step, so the result carries no Hessian and no
    # thermochemistry -- and the arithmetic below was reaching straight
    # into those absent values. `float + None` is a TypeError, which no
    # caller was expecting from a kernel, so it escaped the analysis
    # phase, killed the executor, and left a goal with three and a half
    # hours of validated engine work unsettled (NOVEL-2 ino1,
    # 2026-09-04: five converged spin states lost with the sixth).
    # A missing quantity is a refusal, and the node that produced it is
    # the finding.
    absent = tuple(
        name
        for name in (
            "electronic_energy",
            "zero_point_energy",
            "total_internal_energy",
            "enthalpy",
            "entropy_times_temperature",
            "gibbs_free_energy",
        )
        if getattr(engine, name, None) is None
    )
    if absent:
        # A ScanTS with Freq converged and ran its frequency step; ORCA
        # printed the table for its first Hessian only. The refusal names
        # what the output itself says, not a cause it did not have
        # (R10 Q31, CUHK 2154022).
        unprinted = getattr(
            getattr(engine, "file_object", None),
            "unprinted_frequency_table_reason",
            None,
        )
        cause = (
            f" {unprinted[0].upper()}{unprinted[1:]}."
            if isinstance(unprinted, str) and unprinted
            else " A run whose optimisation did not converge never reached "
            "its frequency step, so there is no Hessian to derive it from."
        )
        raise QuantityExtractionError(
            f"{request.program} result {request.artifact_id!r} carries no "
            "thermochemistry: " + ", ".join(absent) + "." + cause
        )
    evidence_ref = f"artifact:{request.artifact_id}#{request.artifact_sha256}"
    energy_values = {
        "electronic_energy": engine.electronic_energy,
        "zero_point_energy": engine.zero_point_energy,
        "internal_energy": engine.electronic_energy
        + engine.total_internal_energy,
        "enthalpy": engine.enthalpy,
        "entropy_times_temperature": engine.entropy_times_temperature,
        "gibbs_free_energy": engine.gibbs_free_energy,
        "thermal_internal_energy_correction": engine.total_internal_energy,
        "thermal_enthalpy_correction": engine.enthalpy
        - engine.electronic_energy,
        "enthalpy_increment_above_zero_point": (
            engine.enthalpy
            - engine.electronic_energy
            - engine.zero_point_energy
        ),
        "thermal_gibbs_correction": (
            engine.gibbs_free_energy - engine.electronic_energy
        ),
    }
    quantities = [
        _thermo_quantity(
            quantity_id=name,
            value=value,
            source_unit="J mol^-1",
            normalized_unit="hartree",
            dimension=ENERGY,
            evidence_ref=evidence_ref,
        )
        for name, value in energy_values.items()
    ]
    quantities.extend(
        [
            _thermo_quantity(
                quantity_id="entropy",
                value=engine.total_entropy,
                source_unit="J mol^-1 K^-1",
                normalized_unit="hartree K^-1",
                dimension=ENTROPY,
                evidence_ref=evidence_ref,
            ),
            _make_quantity(
                quantity_id="temperature",
                source_value=request.temperature_k,
                source_unit="K",
                value=request.temperature_k,
                unit="K",
                dimension=TEMPERATURE,
                evidence_ref=evidence_ref,
            ),
            _make_quantity(
                quantity_id="pressure",
                source_value=request.pressure_atm,
                source_unit="atm",
                value=request.pressure_atm,
                unit="atm",
                dimension=PRESSURE,
                evidence_ref=evidence_ref,
            ),
            _make_quantity(
                quantity_id="near_zero_mode_count",
                source_value=engine.near_zero_mode_count,
                source_unit="1",
                value=engine.near_zero_mode_count,
                unit="1",
                dimension=DIMENSIONLESS,
                evidence_ref=evidence_ref,
            ),
            _thermo_quantity(
                quantity_id="heat_capacity_cv",
                value=engine.total_heat_capacity,
                source_unit="J mol^-1 K^-1",
                normalized_unit="hartree K^-1",
                dimension=ENTROPY,
                evidence_ref=evidence_ref,
            ),
        ]
    )
    if request.entropy_method != "rrho":
        quantities.extend(
            [
                _thermo_quantity(
                    quantity_id="quasi_harmonic_entropy",
                    value=engine.qrrho_total_entropy,
                    source_unit="J mol^-1 K^-1",
                    normalized_unit="hartree K^-1",
                    dimension=ENTROPY,
                    evidence_ref=evidence_ref,
                ),
                _thermo_quantity(
                    quantity_id="quasi_harmonic_entropy_times_temperature",
                    value=engine.qrrho_entropy_times_temperature,
                    source_unit="J mol^-1",
                    normalized_unit="hartree",
                    dimension=ENERGY,
                    evidence_ref=evidence_ref,
                ),
            ]
        )
    if request.enthalpy_cutoff_cm1 is not None:
        quantities.append(
            _thermo_quantity(
                quantity_id="quasi_harmonic_enthalpy",
                value=engine.qrrho_enthalpy,
                source_unit="J mol^-1",
                normalized_unit="hartree",
                dimension=ENERGY,
                evidence_ref=evidence_ref,
            )
        )
    if (
        request.entropy_method != "rrho"
        or request.enthalpy_cutoff_cm1 is not None
    ):
        if (
            request.entropy_method != "rrho"
            and request.enthalpy_cutoff_cm1 is not None
        ):
            quasi_harmonic_gibbs = engine.qrrho_gibbs_free_energy
        elif request.entropy_method != "rrho":
            quasi_harmonic_gibbs = engine.qrrho_gibbs_free_energy_qs
        else:
            quasi_harmonic_gibbs = engine.qrrho_gibbs_free_energy_qh
        quantities.extend(
            [
                _thermo_quantity(
                    quantity_id="quasi_harmonic_gibbs_free_energy",
                    value=quasi_harmonic_gibbs,
                    source_unit="J mol^-1",
                    normalized_unit="hartree",
                    dimension=ENERGY,
                    evidence_ref=evidence_ref,
                ),
                _thermo_quantity(
                    quantity_id=("quasi_harmonic_thermal_gibbs_correction"),
                    value=quasi_harmonic_gibbs - engine.electronic_energy,
                    source_unit="J mol^-1",
                    normalized_unit="hartree",
                    dimension=ENERGY,
                    evidence_ref=evidence_ref,
                ),
            ]
        )
    # The table a plan is checked against is the table this writer keeps:
    # a quantity written here and missing there is one a plan cannot name,
    # and one listed there and not written here is one a plan is promised.
    written = tuple(sorted(item.quantity_id for item in quantities))
    promised = thermochemistry_quantities_for_treatment(
        request.entropy_method, request.enthalpy_cutoff_cm1
    )
    if written != promised:
        raise QuantityContractError(
            "thermochemistry writer and its quantity table disagree: wrote "
            f"{list(written)}, the table promises {list(promised)}"
        )
    if result_file_sha256(artifact) != request.artifact_sha256:
        raise QuantityExtractionError(
            "result artifact changed during thermochemistry derivation"
        )
    # Which torsions this number counts as harmonic oscillators, and the
    # request that would count one as a hindered rotor instead: a low
    # torsion is not a harmonic oscillator, and a receipt that says only
    # "rigid-rotor harmonic-oscillator" does not say which modes that
    # assumption is worst for (R10 Q30).  Not said of a held surface,
    # whose statements already name what was removed.
    torsion_words = (
        _harmonic_torsion_statement(
            reader,
            reached,
            rotor_bonds=tuple(
                rotor_treatment.rotor_bonds
                if rotor_treatment is not None
                else ()
            ),
        )
        if projection is None
        else ""
    )
    for rotor in request.internal_rotors:
        if (
            result_file_sha256(rotor_artifact_paths[rotor.scan_artifact_id])
            != rotor.scan_artifact_sha256
        ):
            raise QuantityExtractionError(
                "a rotor's scan artifact changed during thermochemistry "
                "derivation"
            )
    assumptions = (
        _thermochemistry_assumptions(request, engine.convention_statements)
        + (
            # What the free energy stands on, with the number or the marker
            # that says so. Inside ``assumptions`` (already inside the digest
            # and the recorded record) for the reason the reaction-coordinate
            # selection is: a receipt minted before this line keeps verifying.
            # A projected free energy stands on the surface it was projected
            # onto, and says which coordinates, how many modes and which rotor
            # treatment, in the same place.
            (stationarity.sentence(),)
            if projection is None
            else projection.statements
        )
        # Which torsions are hindered rotors, where each potential came
        # from, which moment and symmetry numbers, and what replaced the
        # harmonic mode: said in the same place, for the same reason.
        + (rotor_treatment.statements if rotor_treatment is not None else ())
        + ((missing_modes,) if missing_modes else ())
        + ((torsion_words,) if torsion_words else ())
    )
    body = {
        "schema_version": "chemsmart.thermochemistry-receipt.v1",
        "artifact_id": request.artifact_id,
        "artifact_sha256": request.artifact_sha256,
        "program": request.program,
        "engine_id": "chemsmart.analysis.thermochemistry.Thermochemistry",
        "temperature_k": request.temperature_k,
        "pressure_atm": request.pressure_atm,
        "quantities": tuple(quantities),
        "assumptions": assumptions,
        "status": "derived",
    }
    body = _extended_thermochemistry_body(body, request)
    return ThermochemistryReceiptV1(
        **body, receipt_sha256=canonical_quantity_sha256(body)
    )


def derive_pyscf_thermochemistry(
    *,
    request: ThermochemistryRequestV1,
    artifact_path: str | os.PathLike[str],
) -> ThermochemistryReceiptV1:
    """Backward-compatible PySCF entry point for the shared implementation."""

    if request.program != "pyscf":
        raise QuantityContractError(
            "derive_pyscf_thermochemistry requires program 'pyscf'"
        )
    return derive_result_thermochemistry(
        request=request,
        artifact_path=artifact_path,
    )


def quantity_map(
    receipts: Iterable[QuantityExtractionReceiptV1 | ThermochemistryReceiptV1],
) -> dict[str, QuantityValueV1]:
    """Return a unique ID-to-value mapping for expression evaluation."""

    values: dict[str, QuantityValueV1] = {}
    for receipt in receipts:
        for quantity in receipt.quantities:
            if quantity.quantity_id in values:
                raise QuantityContractError(
                    f"duplicate quantity_id across receipts: {quantity.quantity_id}"
                )
            values[quantity.quantity_id] = quantity
    return values


def quantity_value_from_record(record: Mapping[str, Any]) -> QuantityValueV1:
    """Reconstruct and revalidate one canonical quantity event record."""

    values = dict(record)
    values["dimension"] = tuple(values.get("dimension") or ())
    return QuantityValueV1(**values)


def quantity_extraction_receipt_from_record(
    record: Mapping[str, Any], *, receipt_sha256: str
) -> QuantityExtractionReceiptV1:
    """Rehydrate an extraction receipt persisted by Runtime V2."""

    values = dict(record)
    values["quantities"] = tuple(
        quantity_value_from_record(item)
        for item in values.get("quantities") or ()
    )
    return QuantityExtractionReceiptV1(**values, receipt_sha256=receipt_sha256)


def thermochemistry_receipt_from_record(
    record: Mapping[str, Any], *, receipt_sha256: str
) -> ThermochemistryReceiptV1:
    """Rehydrate a thermochemistry receipt persisted by Runtime V2."""

    values = dict(record)
    values["quantities"] = tuple(
        quantity_value_from_record(item)
        for item in values.get("quantities") or ()
    )
    values["assumptions"] = tuple(values.get("assumptions") or ())
    return ThermochemistryReceiptV1(**values, receipt_sha256=receipt_sha256)


__all__ = [
    "ANGLE",
    "DIMENSIONLESS",
    "DIPOLE_MOMENT",
    "ENERGY",
    "ENTROPY",
    "FREQUENCY",
    "LENGTH",
    "AREA",
    "VOLUME",
    "CHARGE",
    "ELECTRIC_POTENTIAL",
    "MASS",
    "MOMENT_OF_INERTIA",
    "PRESSURE",
    "SUPPORTED_PYSCF_SELECTORS",
    "TEMPERATURE",
    "Dimension",
    "QuantityContractError",
    "QuantityExtractionError",
    "QuantityExtractionReceiptV1",
    "QuantitySelectorV1",
    "QuantityValueV1",
    "ResultQuantityExtractionRequestV1",
    "ThermochemistryReceiptV1",
    "ThermochemistryRequestV1",
    "canonical_extraction_receipt_body",
    "canonical_quantity_sha256",
    "derive_pyscf_thermochemistry",
    "derive_result_thermochemistry",
    "extract_pyscf_quantities",
    "make_quantity_value",
    "quantity_map",
    "quantity_value_from_record",
    "quantity_extraction_receipt_from_record",
    "result_file_sha256",
    "thermochemistry_receipt_from_record",
    "validate_pyscf_analysis_artifact",
]
