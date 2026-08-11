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
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

__all__ = [
    "RESULT_READERS",
    "MissingQuantityError",
    "ResultReaderV1",
    "SELECTOR_UNITS",
    "reader_for",
    "registered_reader_programs",
    "registered_reader_selectors",
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
    "gibbs_free_energy": "Eh",
    "oscillator_strengths": "",
    "vibrational_frequencies": "cm^-1",
    "vpt2_harmonic_frequencies": "cm^-1",
    "vpt2_fundamental_frequencies": "cm^-1",
    "vpt2_zero_point_rovibrational_energy": "cm^-1",
    "positions": "Angstrom",
    "symbols": "",
    "charge": "",
    "multiplicity": "",
    "scf_energy": "Eh",
    "correlation_energy": "Eh",
    "dipole_moment": "Debye",
    "dipole_moment_magnitude": "Debye",
}


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


def _symbols(output: Any) -> list[str]:
    return [str(item) for item in output.molecule.chemical_symbols]


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
    #: Native program outputs must prove normal termination.  A standalone
    #: geometry artifact is data rather than an engine run, so that format can
    #: opt out while retaining the same typed quantity path.
    requires_normal_termination: bool = True

    @property
    def selectors(self) -> frozenset[str]:
        return frozenset(self.accessors)

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
        try:
            value = self.accessors[selector](output)
        except MissingQuantityError:
            raise
        except Exception as exc:
            # The parsers raise IndexError/TypeError when a block the run never
            # wrote is requested. That is an absent quantity, not a defect.
            raise MissingQuantityError(
                f"{self.program} result contains no {selector!r} value "
                f"({type(exc).__name__})"
            ) from exc
        if value is None or (isinstance(value, (list, tuple)) and not value):
            raise MissingQuantityError(
                f"{self.program} result contains no {selector!r} value"
            )
        return value, SELECTOR_UNITS[selector]


def _text_output_accessors(
    *, thermochemistry: bool = True
) -> dict[str, Callable[[Any], Any]]:
    """Accessors shared by the log-parsing programs."""

    accessors: dict[str, Callable[[Any], Any]] = {
        "energy": _last_energy,
        "energies": lambda output: [float(item) for item in output.energies],
        "vibrational_frequencies": lambda output: [
            float(item) for item in output.vibrational_frequencies
        ],
        "positions": _positions,
        "symbols": _symbols,
        "charge": lambda output: int(output.charge),
        "multiplicity": lambda output: int(output.multiplicity),
    }
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

    ``ORCAOutput.final_energy`` historically prefers ``final_scf_energy``.
    That is suitable for DFT but silently drops the correlation contribution
    from MP2 and coupled-cluster results.  The last explicit ORCA final-energy
    record is the program's total result and remains correct for DFT jobs too.
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
    values: list[float] = []
    for line in getattr(output, "contents", ()):
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
    """Return total post-SCF correlation as ``E(total) - E(SCF)``."""

    return _orca_total_energy(output) - _orca_scf_energy(output)


def _orca_positions(output: Any) -> list[list[float]]:
    molecule = output.thermochemistry_molecule
    return [[float(value) for value in row] for row in molecule.positions]


def _orca_symbols(output: Any) -> list[str]:
    return [
        str(item)
        for item in output.thermochemistry_molecule.chemical_symbols
    ]


def _orca_accessors() -> dict[str, Callable[[Any], Any]]:
    accessors = _text_output_accessors()
    accessors.update(
        {
            "absorption_wavelengths": lambda output: [
                float(item) for item in output.absorption_wavelengths
            ],
            "excitation_energies": lambda output: [
                float(item) for item in output.excitation_energies_eV
            ],
            "oscillator_strengths": lambda output: [
                float(item) for item in output.oscillator_strengths
            ],
            "energy": _orca_total_energy,
            "entropy_times_temperature": lambda output: float(
                output.entropy_times_temperature
            ),
            "positions": _orca_positions,
            "symbols": _orca_symbols,
            "charge": lambda output: int(output.thermochemistry_charge),
            "multiplicity": lambda output: int(
                output.thermochemistry_multiplicity
            ),
            "scf_energy": _orca_scf_energy,
            "correlation_energy": _orca_correlation_energy,
            "dipole_moment": lambda output: [
                float(item)
                for item in output.dipole_moment_in_debye.reshape(-1)
            ],
            "dipole_moment_magnitude": lambda output: float(
                output.dipole_moment_magnitude_in_debye
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
        }
    )
    return accessors


def _gaussian_output(path: Path) -> Any:
    from chemsmart.io.gaussian.output import Gaussian16Output

    return Gaussian16Output(filename=str(path))


def _gaussian_accessors() -> dict[str, Callable[[Any], Any]]:
    accessors = _text_output_accessors()
    accessors.update(
        {
            "absorption_wavelengths": lambda output: [
                float(item) for item in output.absorptions_in_nm
            ],
            "excitation_energies": lambda output: [
                float(item) for item in output.excitation_energies_eV
            ],
            "oscillator_strengths": lambda output: [
                float(item) for item in output.oscillatory_strengths
            ],
            "dipole_moment": lambda output: [
                float(item) for item in output.all_dipole_moments[-1]
            ],
            "dipole_moment_magnitude": lambda output: float(
                output.all_dipole_moment_magnitudes[-1]
            ),
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


def _xtb_accessors() -> dict[str, Callable[[Any], Any]]:
    # xTB's reader exposes no charge/multiplicity or energy series, so those
    # selectors stay unregistered rather than being faked from the geometry.
    return {
        "energy": _last_energy,
        "vibrational_frequencies": lambda output: [
            float(item) for item in output.vibrational_frequencies
        ],
        "positions": _positions,
        "symbols": _symbols,
        "dipole_moment": lambda output: [
            float(item) for item in output.molecular_dipole_full
        ],
        "dipole_moment_magnitude": lambda output: float(
            output.total_molecular_dipole_moment
        ),
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
        "symbols": _symbols,
    }


#: Programs whose results can be read into typed quantities.  PySCF keeps its
#: dedicated structured-HDF5 path in ``result_quantities``; the entries here are
#: the log-parsing programs that had none.
RESULT_READERS: dict[str, ResultReaderV1] = {
    "orca": ResultReaderV1(
        program="orca",
        artifact_kind="orca_output",
        parser_id="chemsmart.io.orca.output.ORCAOutput",
        open_output=_orca_output,
        accessors=_orca_accessors(),
    ),
    "gaussian": ResultReaderV1(
        program="gaussian",
        artifact_kind="gaussian_output",
        parser_id="chemsmart.io.gaussian.output.Gaussian16Output",
        open_output=_gaussian_output,
        accessors=_gaussian_accessors(),
    ),
    "xtb": ResultReaderV1(
        program="xtb",
        artifact_kind="xtb_output",
        parser_id="chemsmart.io.xtb.output.XTBOutput",
        open_output=_xtb_output,
        accessors=_xtb_accessors(),
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
    "absorption_wavelengths": "LENGTH",
    "energy": "ENERGY",
    "energies": "ENERGY",
    "entropy_times_temperature": "ENERGY",
    "excitation_energies": "ENERGY",
    "gibbs_free_energy": "ENERGY",
    "oscillator_strengths": "DIMENSIONLESS",
    "vibrational_frequencies": "FREQUENCY",
    "vpt2_harmonic_frequencies": "FREQUENCY",
    "vpt2_fundamental_frequencies": "FREQUENCY",
    "vpt2_zero_point_rovibrational_energy": "FREQUENCY",
    "positions": "LENGTH",
    "symbols": "DIMENSIONLESS",
    "charge": "DIMENSIONLESS",
    "multiplicity": "DIMENSIONLESS",
    "scf_energy": "ENERGY",
    "correlation_energy": "ENERGY",
    "dipole_moment": "DIPOLE_MOMENT",
    "dipole_moment_magnitude": "DIPOLE_MOMENT",
}


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
    evidence_ref = f"artifact:{request.artifact_id}#{request.artifact_sha256}"
    quantities = []
    for selector in request.selectors:
        try:
            source_value, source_unit = reader.read(output, selector.selector)
        except MissingQuantityError as exc:
            raise rq.QuantityExtractionError(str(exc)) from exc
        dimension = getattr(rq, _SELECTOR_DIMENSIONS[selector.selector])
        if selector.selector == "symbols":
            value = source_value
            unit = source_unit
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
        quantities.append(
            rq.make_quantity_value(
                quantity_id=selector.quantity_id,
                source_value=source_value,
                source_unit=source_unit,
                value=value,
                unit=unit,
                dimension=dimension,
                evidence_ref=evidence_ref,
            )
        )
    if rq.result_file_sha256(artifact) != request.artifact_sha256:
        raise rq.QuantityExtractionError(
            "result artifact changed during extraction"
        )
    body = {
        "schema_version": "chemsmart.quantity-extraction-receipt.v1",
        "artifact_id": request.artifact_id,
        "artifact_sha256": request.artifact_sha256,
        "program": request.program,
        "parser_id": reader.parser_id,
        "quantities": tuple(quantities),
        "status": "extracted",
    }
    return rq.QuantityExtractionReceiptV1(
        **body, receipt_sha256=rq.canonical_quantity_sha256(body)
    )


def reader_for(program: str) -> ResultReaderV1 | None:
    """Return the registered reader for ``program``, or ``None``."""

    return RESULT_READERS.get(str(program).strip().lower())


def registered_reader_programs() -> tuple[str, ...]:
    """Return every program with a registered log reader, sorted."""

    return tuple(sorted(RESULT_READERS))


def registered_reader_selectors() -> dict[str, tuple[str, ...]]:
    """Return the model-discoverable selector inventory for each reader."""

    return {
        program: tuple(sorted(reader.selectors))
        for program, reader in sorted(RESULT_READERS.items())
    }
