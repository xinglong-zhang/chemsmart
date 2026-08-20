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
    derive_result_thermochemistry,
    extract_pyscf_quantities,
)


def typed_result_artifact_kind(program: str) -> str:
    """Return the trusted artifact kind consumed by one program's reader."""

    from chemsmart.analysis.result_readers import (
        reader_for,
        registered_reader_programs,
    )

    normalized = str(program).strip().lower()
    if normalized == "pyscf":
        return "pyscf_hdf5"
    reader = reader_for(normalized)
    if reader is None:
        raise ContractError(
            "no typed result reader is registered for program; registered: "
            f"{('pyscf',) + registered_reader_programs()}"
        )
    return reader.artifact_kind


def extract_trusted_result_quantities(
    *,
    artifact: TrustedArtifactRefV1,
    program: str,
    selectors: tuple[QuantitySelectorV1, ...],
) -> QuantityExtractionReceiptV1:
    """Extract quantities after resolving the artifact through host state."""

    from chemsmart.analysis.result_readers import extract_logged_quantities

    normalized_program = str(program).strip().lower()
    expected_kind = typed_result_artifact_kind(normalized_program)
    if artifact.kind != expected_kind:
        raise ContractError(
            f"{normalized_program} quantity extraction requires a bound "
            f"{expected_kind} artifact, not {artifact.kind!r}"
        )
    request = ResultQuantityExtractionRequestV1(
        schema_version="chemsmart.quantity-extraction-request.v1",
        artifact_id=artifact.artifact_id,
        artifact_sha256=artifact.sha256,
        program=normalized_program,
        selectors=selectors,
    )
    if normalized_program == "pyscf":
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
    concentration_mol_l: float | None = None,
    entropy_method: str = "rrho",
    entropy_cutoff_cm1: float | None = None,
    enthalpy_cutoff_cm1: float | None = None,
    alpha: int = 4,
    use_weighted_mass: bool = False,
    frequency_scale_factor: float = 1.0,
) -> ThermochemistryReceiptV1:
    """Evaluate shared RRHO or quasi-harmonic thermochemistry."""

    normalized_program = str(program).strip().lower()
    expected_kind = typed_result_artifact_kind(normalized_program)
    if artifact.kind != expected_kind:
        raise ContractError(
            f"{normalized_program} thermochemistry requires a bound "
            f"{expected_kind} artifact, not {artifact.kind!r}"
        )
    request = ThermochemistryRequestV1(
        schema_version="chemsmart.thermochemistry-request.v1",
        artifact_id=artifact.artifact_id,
        artifact_sha256=artifact.sha256,
        program=normalized_program,
        temperature_k=temperature_k,
        pressure_atm=pressure_atm,
        concentration_mol_l=concentration_mol_l,
        entropy_method=entropy_method,
        entropy_cutoff_cm1=entropy_cutoff_cm1,
        enthalpy_cutoff_cm1=enthalpy_cutoff_cm1,
        alpha=alpha,
        use_weighted_mass=use_weighted_mass,
        frequency_scale_factor=frequency_scale_factor,
    )
    if normalized_program == "pyscf":
        return derive_pyscf_thermochemistry(
            request=request,
            artifact_path=artifact.path,
        )
    return derive_result_thermochemistry(
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
    "typed_result_artifact_kind",
]
