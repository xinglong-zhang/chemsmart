"""Host-side bridge from trusted agent artifacts to scientific postprocessing.

The model-facing runtime passes semantic selectors and expression nodes.  This
module supplies the path only after the host has resolved a
``TrustedArtifactRefV1``; it never accepts a model-authored filesystem path.
"""

from __future__ import annotations

from chemsmart.agent._contracts import ContractError, TrustedArtifactRefV1
from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionReceiptV1,
    QuantityExpressionRequestV1,
    evaluate_quantity_expression,
)
from chemsmart.analysis.result_quantities import (
    QuantityExtractionReceiptV1,
    QuantitySelectorV1,
    ResultQuantityExtractionRequestV1,
    ThermochemistryReceiptV1,
    ThermochemistryRequestV1,
    derive_pyscf_thermochemistry,
    extract_pyscf_quantities,
)


def extract_trusted_result_quantities(
    *,
    artifact: TrustedArtifactRefV1,
    program: str,
    selectors: tuple[QuantitySelectorV1, ...],
) -> QuantityExtractionReceiptV1:
    """Extract quantities after resolving the artifact through host state."""

    from chemsmart.analysis.result_readers import (
        extract_logged_quantities,
        reader_for,
        registered_reader_programs,
    )

    reader = reader_for(program)
    if program != "pyscf" and reader is None:
        raise ContractError(
            "no quantity reader is registered for program; registered: "
            f"{('pyscf',) + registered_reader_programs()}"
        )
    expected_kind = "pyscf_hdf5" if program == "pyscf" else reader.artifact_kind
    if artifact.kind != expected_kind:
        raise ContractError(
            f"{program} quantity extraction requires a bound "
            f"{expected_kind} artifact, not {artifact.kind!r}"
        )
    request = ResultQuantityExtractionRequestV1(
        schema_version="chemsmart.quantity-extraction-request.v1",
        artifact_id=artifact.artifact_id,
        artifact_sha256=artifact.sha256,
        program=program,
        selectors=selectors,
    )
    if program == "pyscf":
        return extract_pyscf_quantities(
            request=request, artifact_path=artifact.path
        )
    return extract_logged_quantities(
        request=request, artifact_path=artifact.path
    )


def derive_trusted_thermochemistry(
    *,
    artifact: TrustedArtifactRefV1,
    program: str,
    temperature_k: float,
    pressure_atm: float,
) -> ThermochemistryReceiptV1:
    """Evaluate shared ChemSmart RRHO thermochemistry for a trusted result."""

    if program != "pyscf":
        raise ContractError("no structured thermochemistry bridge is registered")
    if artifact.kind != "pyscf_hdf5":
        raise ContractError("PySCF thermochemistry requires a bound pyscf_hdf5 artifact")
    request = ThermochemistryRequestV1(
        schema_version="chemsmart.thermochemistry-request.v1",
        artifact_id=artifact.artifact_id,
        artifact_sha256=artifact.sha256,
        program=program,
        temperature_k=temperature_k,
        pressure_atm=pressure_atm,
    )
    return derive_pyscf_thermochemistry(
        request=request,
        artifact_path=artifact.path,
    )


def evaluate_typed_quantity_expression(
    request: QuantityExpressionRequestV1,
) -> QuantityExpressionReceiptV1:
    """Evaluate the bounded AST; no source code or formula text is accepted."""

    return evaluate_quantity_expression(request)


__all__ = [
    "derive_trusted_thermochemistry",
    "evaluate_typed_quantity_expression",
    "extract_trusted_result_quantities",
]
