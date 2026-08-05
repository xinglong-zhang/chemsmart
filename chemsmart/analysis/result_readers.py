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
    "energy": "Eh",
    "energies": "Eh",
    "gibbs_free_energy": "Eh",
    "vibrational_frequencies": "cm^-1",
    "positions": "Angstrom",
    "symbols": "",
    "charge": "",
    "multiplicity": "",
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


def _gaussian_output(path: Path) -> Any:
    from chemsmart.io.gaussian.output import Gaussian16Output

    return Gaussian16Output(filename=str(path))


def _xtb_output(path: Path) -> Any:
    from chemsmart.io.xtb.output import XTBOutput

    return XTBOutput(filename=str(path))


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
        accessors=_text_output_accessors(),
    ),
    "gaussian": ResultReaderV1(
        program="gaussian",
        artifact_kind="gaussian_output",
        parser_id="chemsmart.io.gaussian.output.Gaussian16Output",
        open_output=_gaussian_output,
        accessors=_text_output_accessors(),
    ),
    "xtb": ResultReaderV1(
        program="xtb",
        artifact_kind="xtb_output",
        parser_id="chemsmart.io.xtb.output.XTBOutput",
        open_output=_xtb_output,
        accessors=_xtb_accessors(),
    ),
}


#: Physical dimension of each selector, in the shared quantity vocabulary.
_SELECTOR_DIMENSIONS = {
    "energy": "ENERGY",
    "energies": "ENERGY",
    "gibbs_free_energy": "ENERGY",
    "vibrational_frequencies": "FREQUENCY",
    "positions": "LENGTH",
    "symbols": "DIMENSIONLESS",
    "charge": "DIMENSIONLESS",
    "multiplicity": "DIMENSIONLESS",
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
    evidence_ref = f"artifact:{request.artifact_id}#{request.artifact_sha256}"
    quantities = []
    for selector in request.selectors:
        try:
            value, unit = reader.read(output, selector.selector)
        except MissingQuantityError as exc:
            raise rq.QuantityExtractionError(str(exc)) from exc
        dimension = getattr(rq, _SELECTOR_DIMENSIONS[selector.selector])
        quantities.append(
            rq.make_quantity_value(
                quantity_id=selector.quantity_id,
                source_value=value,
                source_unit=unit,
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
