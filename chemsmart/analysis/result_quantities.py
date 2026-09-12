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

CHARGE: Dimension = (0, 0, 0, 0, 0, 0, 0, 0, 1)
ELECTRIC_POTENTIAL: Dimension = (1, 0, 0, 0, 0, 0, 0, 0, -1)

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
        # The second column of the same population block, read for years
        # and discarded one index from where a session needed it
        # (NOVEL-3 ino3, 2026-09-05); open shells only, sum 2S.
        "mulliken_atomic_spin_populations",
        "loewdin_atomic_spin_populations",
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
) -> tuple[str, ...]:
    """Return the quantity IDs a receipt will carry for this entropy method."""

    names = set(DERIVABLE_THERMOCHEMISTRY_QUANTITIES)
    if str(entropy_method or "rrho").strip().lower() != "rrho":
        names.update(QUASI_HARMONIC_THERMOCHEMISTRY_QUANTITIES)
    return tuple(sorted(names))


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
    derivable.update(QUASI_HARMONIC_THERMOCHEMISTRY_QUANTITIES)
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
        "frequency-bearing result under an explicitly bound temperature, "
        "pressure, and standard state; plan one on the producing "
        "calculation instead of selecting these from the log."
    )


_IDENTIFIER = re.compile(r"^[A-Za-z][A-Za-z0-9_.:-]{0,127}$")
#: Executed-evidence contracts the analysis plane admits, owned by the
#: writer: a previous supported contract is a subset of the current one.
_SUPPORTED_PYSCF_RESULT_CONTRACTS = SUPPORTED_RESULT_CONTRACT_VERSIONS


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
    "energy": ("results/energies",),
    "energies": ("results/energies",),
    "excitation_energies": ("results/excitation_energies",),
    "oscillator_strengths": ("results/oscillator_strengths",),
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
    "spin_square": ("results/spin_square",),
    "spin_square_deviation": ("results/spin_square",),
    "effective_multiplicity": ("results/spin_square_effective_multiplicity",),
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
        if self.selector not in SUPPORTED_SELECTORS:
            elsewhere = QUANTITIES_FROM_ANOTHER_TOOL.get(self.selector)
            detail = (
                f"; that quantity is produced by {elsewhere}, not by result "
                "extraction"
                if elsewhere
                else f"; supported selectors: {sorted(SUPPORTED_SELECTORS)}"
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
        if len(self.dimension) not in {6, 7, 8, 9} or not all(
            isinstance(exponent, int) for exponent in self.dimension
        ):
            raise QuantityContractError(
                "dimension must contain six legacy, seven dipole-extended, "
                "eight mass-extended, or nine charge-extended integers"
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
) -> tuple[str, ...]:
    # PySCF receipts used to keep a shorter "legacy" assumption list at
    # default settings that never named the standard state; every program's
    # receipt now says the same things, and a PySCF receipt additionally
    # says which mass table produced its frequencies (below).
    assumptions = [
        "rigid-rotor harmonic-oscillator thermochemistry for harmonic quantities",
        "ground-state electronic degeneracy equals spin multiplicity",
        (
            "natural-abundance weighted isotopic masses"
            if request.use_weighted_mass
            else "most-abundant isotopic masses"
        ),
        "rotational symmetry derived by the shared ChemSmart engine",
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
    if request.program == "pyscf":
        # The frequencies came from PySCF's harmonic analysis under its
        # isotope-averaged masses while the rotational and translational
        # terms above use the table this engine names; both conventions
        # are on the receipt because the artifact records the first.
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


def derive_result_thermochemistry(
    *,
    request: ThermochemistryRequestV1,
    artifact_path: str | os.PathLike[str],
) -> ThermochemistryReceiptV1:
    """Derive RRHO or quasi-harmonic thermochemistry from a trusted result.

    The formulas and molecular conventions remain owned by ChemSmart's common
    :class:`Thermochemistry` engine.  This function only binds conditions and
    serializes the resulting values with explicit units and provenance.
    """

    artifact = _verify_artifact(artifact_path, request.artifact_sha256)
    if request.program == "pyscf":
        output = PySCFOutput(artifact)
        _require_analysis_ready_pyscf_result(
            artifact=artifact,
            expected_sha256=request.artifact_sha256,
            output=output,
            required_units={
                "results/energies": "Eh",
                "results/positions": "Angstrom",
                "results/hessian": "Eh/Bohr^2",
                "results/vibrational_frequencies": "cm^-1",
            },
        )
        if not output.freq:
            raise QuantityExtractionError(
                "thermochemistry requires a validated PySCF Hessian result"
            )
        if output.result_sha256 != request.artifact_sha256:
            raise QuantityExtractionError(
                "PySCF parser observed substituted bytes"
            )
    engine = Thermochemistry(
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
        raise QuantityExtractionError(
            f"{request.program} result {request.artifact_id!r} carries no "
            "thermochemistry: " + ", ".join(absent) + ". A run whose "
            "optimisation did not converge never reached its frequency "
            "step, so there is no Hessian to derive it from."
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
    if result_file_sha256(artifact) != request.artifact_sha256:
        raise QuantityExtractionError(
            "result artifact changed during thermochemistry derivation"
        )
    assumptions = _thermochemistry_assumptions(request)
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
