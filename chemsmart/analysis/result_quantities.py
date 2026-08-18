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
from typing import Any, Iterable, Mapping

import numpy as np

from chemsmart.analysis.thermochemistry import Thermochemistry
from chemsmart.io.pyscf.output import PySCFOutput
from chemsmart.jobs.pyscf.environment import canonical_sha256
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

SUPPORTED_PYSCF_SELECTORS = frozenset(
    {
        "energy",
        "energies",
        "positions",
        "connectivity",
        "symbols",
        "vibrational_frequencies",
        "scan_coordinate_values",
        "scan_energies",
        "scan_point_indices",
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


_IDENTIFIER = re.compile(r"^[A-Za-z][A-Za-z0-9_.:-]{0,127}$")
_CURRENT_PYSCF_RESULT_CONTRACT = "chemsmart.pyscf-result-contract.v3"
_SELECTOR_RESULT_UNITS = {
    "energy": {"results/energies": "Eh"},
    "energies": {"results/energies": "Eh"},
    "excitation_energies": {"results/excitation_energies": "Eh"},
    "oscillator_strengths": {"results/oscillator_strengths": "dimensionless"},
    "dipole_moment": {"results/dipole_moment": "Debye"},
    "dipole_moment_magnitude": {"results/dipole_moment": "Debye"},
    "positions": {"results/positions": "Angstrom"},
    "connectivity": {"results/positions": "Angstrom"},
    "vibrational_frequencies": {"results/vibrational_frequencies": "cm^-1"},
    "scan_energies": {"results/scan_energies": "Eh"},
    "scan_coordinate_values": {"results/scan_coordinate_values": ""},
    "scan_point_indices": {"results/scan_point_indices": ""},
    "homo": {"results/mo_energy": "Eh"},
    "lumo": {"results/mo_energy": "Eh"},
    "gap": {"results/mo_energy": "Eh"},
    "spin_square": {"results/spin_square": "dimensionless"},
    "spin_square_deviation": {"results/spin_square": "dimensionless"},
    "effective_multiplicity": {
        "results/spin_square_effective_multiplicity": "dimensionless"
    },
}


class QuantityContractError(ValueError):
    """Raised when a quantity request or result violates its typed contract."""


class QuantityExtractionError(QuantityContractError):
    """Raised when trusted result evidence cannot be extracted safely."""


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
        if self.program != "pyscf":
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
        if len(self.dimension) not in {6, 7, 8} or not all(
            isinstance(exponent, int) for exponent in self.dimension
        ):
            raise QuantityContractError(
                "dimension must contain six legacy, seven dipole-extended, "
                "or eight mass-extended integers"
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

    def __post_init__(self) -> None:
        object.__setattr__(self, "quantities", tuple(self.quantities))
        if self.schema_version != "chemsmart.quantity-extraction-receipt.v1":
            raise QuantityContractError(
                "unsupported extraction receipt schema"
            )
        if self.status != "extracted":
            raise QuantityContractError("invalid extraction receipt status")
        body = {
            "schema_version": self.schema_version,
            "artifact_id": self.artifact_id,
            "artifact_sha256": self.artifact_sha256,
            "program": self.program,
            "parser_id": self.parser_id,
            "quantities": self.quantities,
            "status": self.status,
        }
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

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.thermochemistry-request.v1":
            raise QuantityContractError(
                "unsupported thermochemistry request schema"
            )
        _require_identifier(self.artifact_id, "artifact_id")
        _require_sha256(self.artifact_sha256)
        normalized_program = str(self.program).strip().lower()
        object.__setattr__(self, "program", normalized_program)
        if normalized_program != "pyscf":
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
        extended_body = {
            **legacy_body,
            "concentration_mol_l": self.concentration_mol_l,
            "entropy_method": self.entropy_method,
            "entropy_cutoff_cm1": self.entropy_cutoff_cm1,
            "enthalpy_cutoff_cm1": self.enthalpy_cutoff_cm1,
            "alpha": self.alpha,
            "use_weighted_mass": self.use_weighted_mass,
            "frequency_scale_factor": self.frequency_scale_factor,
        }
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


def _require_analysis_ready_pyscf_result(
    *,
    artifact: Path,
    expected_sha256: str,
    output: PySCFOutput,
    required_units: dict[str, str],
) -> dict[str, Any]:
    """Require current, executed, receipt-bound result evidence.

    Fake-preview HDF5 files intentionally resemble real results so downstream
    readers can be exercised.  They are not numerical evidence.  Admission
    therefore requires both the machine contract inside HDF5 and the sibling
    ChemSmart run receipt that binds deterministic checks to these bytes.  A
    Hessian may remain ``unclassified`` until a minimum/transition-state policy
    is supplied; it is still analysis-ready when every invariant is green.
    """

    if (
        not output.normal_termination
        or not output.engine_complete
        or output.failure is not None
        or output.spec.get("preview_only") is not False
        or output.spec.get("result_contract_version")
        != _CURRENT_PYSCF_RESULT_CONTRACT
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
    receipt_state = receipt.get("state")
    scientific_state = receipt.get("scientific_validation_state")
    if (
        receipt.get("fake") is not False
        or receipt.get("engine_complete") is not True
        or receipt.get("child_returncode") != 0
        or receipt.get("findings") not in ([], ())
        or receipt_state not in {"validated", "engine_complete"}
        or scientific_state not in {"validated", "unclassified"}
        or receipt.get("result_sha256") != expected_sha256
    ):
        raise QuantityExtractionError(
            "PySCF run receipt does not admit this exact result for analysis"
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
        # if it were a broken digest.
        if field == "input_artifact_sha256" and expected in (None, ""):
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


def _extract_selector(
    output: PySCFOutput,
    selector: QuantitySelectorV1,
    evidence_ref: str,
) -> QuantityValueV1:
    name = selector.selector
    value: Any
    if name == "energy":
        if not output.energies:
            raise QuantityExtractionError("PySCF artifact has no energy")
        value = _require_finite_numeric(output.energies[-1], name)
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="hartree",
            value=value,
            unit="hartree",
            dimension=ENERGY,
            evidence_ref=evidence_ref,
        )
    if name == "energies":
        value = _require_finite_numeric(output.energies, name)
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="hartree",
            value=value,
            unit="hartree",
            dimension=ENERGY,
            evidence_ref=evidence_ref,
        )
    if name == "excitation_energies":
        if output.excitation_energies is None:
            raise QuantityExtractionError(
                "PySCF artifact has no excitation energies"
            )
        value = _require_finite_numeric(output.excitation_energies, name)
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="hartree",
            value=value,
            unit="hartree",
            dimension=ENERGY,
            evidence_ref=evidence_ref,
        )
    if name == "oscillator_strengths":
        if output.oscillator_strengths is None:
            raise QuantityExtractionError(
                "PySCF artifact has no oscillator strengths"
            )
        value = _require_finite_numeric(output.oscillator_strengths, name)
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="1",
            value=value,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=evidence_ref,
        )
    if name in {
        "spin_square",
        "spin_square_target",
        "spin_square_deviation",
        "effective_multiplicity",
    }:
        multiplicity = output.multiplicity
        if not isinstance(multiplicity, int) or multiplicity <= 0:
            raise QuantityExtractionError(
                "PySCF artifact does not establish a positive multiplicity"
            )
        target = (float(multiplicity) ** 2 - 1.0) / 4.0
        if name == "spin_square_target":
            value = target
        else:
            spin_value = output.results.get("spin_square")
            effective_value = output.results.get(
                "spin_square_effective_multiplicity"
            )
            if spin_value is None or effective_value is None:
                raise QuantityExtractionError(
                    "PySCF artifact has no complete <S^2> diagnostic"
                )
            spin_square = float(
                np.asarray(
                    _require_finite_numeric(spin_value, "spin_square")
                ).reshape(-1)[0]
            )
            effective = float(
                np.asarray(
                    _require_finite_numeric(
                        effective_value, "effective_multiplicity"
                    )
                ).reshape(-1)[0]
            )
            value = {
                "spin_square": spin_square,
                "spin_square_deviation": spin_square - target,
                "effective_multiplicity": effective,
            }[name]
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="1",
            value=value,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=evidence_ref,
        )
    if name in {"dipole_moment", "dipole_moment_magnitude"}:
        if output.dipole_moment is None:
            raise QuantityExtractionError(
                "PySCF artifact has no dipole moment"
            )
        vector = _require_finite_numeric(output.dipole_moment, name)
        value = (
            float(np.linalg.norm(np.asarray(vector, dtype=float)))
            if name == "dipole_moment_magnitude"
            else vector
        )
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="Debye",
            value=value,
            unit="debye",
            dimension=DIPOLE_MOMENT,
            evidence_ref=evidence_ref,
        )
    if name == "positions":
        if output.positions is None:
            raise QuantityExtractionError("PySCF artifact has no positions")
        value = _require_finite_numeric(output.positions, name)
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="angstrom",
            value=value,
            unit="angstrom",
            dimension=LENGTH,
            evidence_ref=evidence_ref,
        )
    if name == "connectivity":
        if output.positions is None or not output.symbols:
            raise QuantityExtractionError(
                "PySCF artifact has no complete geometry for connectivity"
            )
        from chemsmart.analysis.result_readers import _connectivity_matrix

        value = _connectivity_matrix(output.get_molecule())
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="1",
            value=value,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=evidence_ref,
        )
    if name == "symbols":
        value = tuple(str(symbol) for symbol in output.chemical_symbols)
        if not value:
            raise QuantityExtractionError("PySCF artifact has no atom symbols")
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="1",
            value=value,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=evidence_ref,
            data_kind="text_vector",
        )
    if name == "vibrational_frequencies":
        if output.vibrational_frequencies is None:
            raise QuantityExtractionError(
                "PySCF artifact has no vibrational frequencies"
            )
        value = _require_finite_numeric(output.vibrational_frequencies, name)
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="cm^-1",
            value=value,
            unit="cm^-1",
            dimension=FREQUENCY,
            evidence_ref=evidence_ref,
        )
    if name in {"homo", "lumo", "gap"}:
        attribute = {
            "homo": "homo_energy",
            "lumo": "lumo_energy",
            "gap": "fmo_gap",
        }[name]
        value = getattr(output, attribute)
        if value is None:
            raise QuantityExtractionError(
                f"PySCF artifact cannot define the requested {name} value"
            )
        value = float(_require_finite_numeric(value, name))
        normalized = energy_conversion("eV", "hartree", value)
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="eV",
            value=normalized,
            unit="hartree",
            dimension=ENERGY,
            evidence_ref=evidence_ref,
        )
    if name in {"charge", "multiplicity"}:
        value = getattr(output, name)
        if not isinstance(value, int):
            raise QuantityExtractionError(
                f"PySCF artifact has no integer {name} value"
            )
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="1",
            value=value,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=evidence_ref,
            data_kind="integer",
        )
    if name in {"method", "basis"}:
        value = getattr(output, name)
        if not isinstance(value, str) or not value.strip():
            raise QuantityExtractionError(
                f"PySCF artifact has no {name} metadata"
            )
        return _make_quantity(
            quantity_id=selector.quantity_id,
            source_value=value,
            source_unit="1",
            value=value,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=evidence_ref,
            data_kind="text",
        )
    raise QuantityContractError(f"unsupported selector: {name!r}")


def extract_pyscf_quantities(
    *,
    request: ResultQuantityExtractionRequestV1,
    artifact_path: str | os.PathLike[str],
) -> QuantityExtractionReceiptV1:
    """Extract selected values from a host-resolved structured PySCF result."""

    artifact = _verify_artifact(artifact_path, request.artifact_sha256)
    output = PySCFOutput(artifact)
    required_units: dict[str, str] = {}
    for selector in request.selectors:
        required_units.update(
            _SELECTOR_RESULT_UNITS.get(selector.selector, {})
        )
    _require_analysis_ready_pyscf_result(
        artifact=artifact,
        expected_sha256=request.artifact_sha256,
        output=output,
        required_units=required_units,
    )
    if output.result_sha256 != request.artifact_sha256:
        raise QuantityExtractionError(
            "PySCF parser observed substituted bytes"
        )
    evidence_ref = f"artifact:{request.artifact_id}#{request.artifact_sha256}"
    quantities = tuple(
        _extract_selector(output, selector, evidence_ref)
        for selector in request.selectors
    )
    if result_file_sha256(artifact) != request.artifact_sha256:
        raise QuantityExtractionError(
            "result artifact changed during extraction"
        )
    body = {
        "schema_version": "chemsmart.quantity-extraction-receipt.v1",
        "artifact_id": request.artifact_id,
        "artifact_sha256": request.artifact_sha256,
        "program": request.program,
        "parser_id": "chemsmart.io.pyscf.output.PySCFOutput",
        "quantities": quantities,
        "status": "extracted",
    }
    return QuantityExtractionReceiptV1(
        **body, receipt_sha256=canonical_quantity_sha256(body)
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


def _uses_legacy_pyscf_thermochemistry_contract(
    request: ThermochemistryRequestV1,
) -> bool:
    return (
        request.program == "pyscf"
        and request.concentration_mol_l is None
        and request.entropy_method == "rrho"
        and request.entropy_cutoff_cm1 is None
        and request.enthalpy_cutoff_cm1 is None
        and request.alpha == 4
        and request.use_weighted_mass is False
        and request.frequency_scale_factor == 1.0
    )


def _thermochemistry_assumptions(
    request: ThermochemistryRequestV1,
) -> tuple[str, ...]:
    if _uses_legacy_pyscf_thermochemistry_contract(request):
        return (
            "ideal-gas translational partition function",
            "rigid-rotor harmonic-oscillator thermochemistry",
            "ground-state electronic degeneracy equals spin multiplicity",
            "most-abundant isotopic masses",
            "rotational symmetry derived by the shared ChemSmart engine",
        )

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
        check_imaginary_frequencies=True,
    )
    if engine.program != request.program:
        raise QuantityExtractionError(
            "trusted result program differs from the requested thermochemistry "
            f"program: expected {request.program!r}, observed {engine.program!r}"
        )
    engine.check_frequencies()
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
                    quantity_id=(
                        "quasi_harmonic_thermal_gibbs_correction"
                    ),
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
    if not _uses_legacy_pyscf_thermochemistry_contract(request):
        body.update(
            {
                "concentration_mol_l": request.concentration_mol_l,
                "entropy_method": request.entropy_method,
                "entropy_cutoff_cm1": request.entropy_cutoff_cm1,
                "enthalpy_cutoff_cm1": request.enthalpy_cutoff_cm1,
                "alpha": request.alpha,
                "use_weighted_mass": request.use_weighted_mass,
                "frequency_scale_factor": request.frequency_scale_factor,
            }
        )
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
