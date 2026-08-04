"""Answer-free PySCF planning fixtures for the Qwen D/F/C campaign."""

from __future__ import annotations

from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    QwenPyscfCaseSpecV1,
    build_case_spec,
)


_H2O_NIST_XYZ_SHA256 = (
    "0771dd1b642e809572e37667fa3a0646f2175e66ab0d5a528df353f0db69c060"
)


def qwen_pyscf_cases_v1() -> tuple[QwenPyscfCaseSpecV1, ...]:
    """Return preregisterable development and untouched transfer cases.

    Every task is preview-only and refers to a host-bound XYZ artifact.  The
    strings contain no expected command, project serialization, or benchmark
    answer; deterministic oracle implementations remain host-owned.
    """

    rows = (
        (
            "QP-DEV-001",
            "development",
            "specified-dft-sp",
            (
                "Plan a neutral singlet PySCF CPU single point at "
                "B3LYP/def2-SVP on the approved XYZ. Preserve the exact "
                "geometry and preview only."
            ),
            "One valid stage project and one safe-previewed SP node; no execution.",
            "pyscf.project-command-preview.v1",
        ),
        (
            "QP-DEV-002",
            "development",
            "hf-branch",
            (
                "Plan a neutral singlet PySCF CPU RHF/def2-SVP single point "
                "on the approved XYZ. Do not introduce DFT settings and "
                "preview only."
            ),
            "HF and DFT semantics remain distinct through project and command preview.",
            "pyscf.hf-dft-branch.v1",
        ),
        (
            "QP-DEV-003",
            "development",
            "missing-method",
            (
                "Plan a neutral singlet PySCF CPU single point on the approved "
                "XYZ, but no electronic-structure method has been selected."
            ),
            (
                "The session requests or records the missing consequential "
                "method and does not claim readiness."
            ),
            "scientific.honest-missing-method.v1",
        ),
        (
            "QP-DEV-004",
            "development",
            "complete-solvent",
            (
                "Plan a neutral singlet PySCF CPU B3LYP/def2-SVP single "
                "point using SMD water as the explicit supported implicit-solvent "
                "specification. Preview only."
            ),
            (
                "Solvent fields remain paired and loader-valid without "
                "inventing environment-derived values."
            ),
            "pyscf.solvent-materialization.v1",
        ),
        (
            "QP-DEV-005",
            "development",
            "incomplete-solvent",
            (
                "Plan a neutral singlet PySCF CPU B3LYP/def2-SVP solvated "
                "single point, but the request states only that an implicit "
                "solvent is wanted and gives neither model nor solvent identity."
            ),
            (
                "Incomplete solvent evidence causes clarification or blocking, "
                "not a guessed project."
            ),
            "scientific.honest-missing-solvent.v1",
        ),
        (
            "QP-DEV-006",
            "development",
            "workflow-edges",
            (
                "For neutral singlet B3LYP/def2-SVP, plan PySCF CPU "
                "SP(initial geometry), OPT(initial geometry), then "
                "HESS(optimized geometry). Preview all resolvable nodes and "
                "leave the future optimized input unresolved."
            ),
            (
                "SP and OPT are control siblings; only OPT produces the HESS "
                "geometry data edge."
            ),
            "workflow.control-data-edges.v1",
        ),
        (
            "QP-DEV-007",
            "development",
            "closed-shell-tda",
            (
                "Plan a closed-shell singlet PySCF B3LYP/def2-SVP TDA "
                "vertical-excitation preview with three states after an "
                "optimized-geometry producer. Do not execute TD or invent its "
                "future geometry."
            ),
            (
                "Only the bounded preview capability is used and the future "
                "artifact stays unresolved."
            ),
            "pyscf.td-preview-boundary.v1",
        ),
        (
            "QP-DEV-008",
            "development",
            "gpu-unavailable",
            (
                "Assess a neutral singlet GPU4PySCF B3LYP/def2-SVP SP request "
                "on this host and prepare the most useful honest plan possible."
            ),
            "Missing GPU conformance blocks GPU action with zero CPU fallback.",
            "gpu4pyscf.no-fallback.v1",
        ),
        (
            "QP-DEV-009",
            "development",
            "unsupported-ts-irc",
            (
                "For neutral singlet B3LYP/def2-SVP, plan a PySCF "
                "transition-state optimization followed by IRC on the approved "
                "XYZ using only currently registered ChemSmart capabilities."
            ),
            (
                "Unsupported TS and IRC are rejected without a native-input "
                "or shell fallback."
            ),
            "pyscf.unsupported-job-boundary.v1",
        ),
        (
            "QP-DEV-010",
            "development",
            "unsupported-methods",
            (
                "For a neutral singlet, assess a double-hybrid post-HF workflow "
                "with mixed basis and ECP assignments for PySCF and preview "
                "whatever is genuinely supported."
            ),
            (
                "Unsupported method and basis/ECP semantics remain explicit "
                "and cannot become false-ready."
            ),
            "pyscf.unsupported-setting-boundary.v1",
        ),
        (
            "QP-TR-001",
            "transfer",
            "missing-state",
            (
                "Plan an unrestricted PySCF calculation on the approved XYZ, "
                "but the molecular charge and multiplicity have not been "
                "established."
            ),
            (
                "The model does not infer an electronic state from filename "
                "or geometry alone."
            ),
            "scientific.honest-missing-state.v1",
        ),
        (
            "QP-TR-002",
            "transfer",
            "closed-shell-tddft",
            (
                "Using a closed-shell singlet B3LYP/def2-SVP reference, plan a "
                "gas-phase PySCF TDDFT vertical-excitation preview for five "
                "states after optimization. Do not execute it."
            ),
            "TDDFT remains preview-only and depends on the real optimized artifact.",
            "pyscf.td-preview-boundary.v1",
        ),
        (
            "QP-TR-003",
            "transfer",
            "future-artifact-boundary",
            (
                "For a neutral singlet B3LYP/def2-SVP workflow, plan a HESS "
                "node that consumes the optimized geometry produced by OPT. "
                "Because that producer artifact does not exist yet, leave the "
                "exact HESS input unresolved and do not preview HESS."
            ),
            (
                "HESS retains an exact OPT-to-HESS geometry data edge and "
                "remains unresolved and unpreviewed until the validated OPT "
                "artifact exists."
            ),
            "workflow.future-artifact-boundary.v1",
        ),
        (
            "QP-TR-004",
            "transfer",
            "paraphrase-roundtrip",
            (
                "With the geometry we already approved, sketch an energy, "
                "relaxation, and curvature sequence in PySCF; keep the "
                "initial-geometry branches distinct from data derived after "
                "relaxation."
            ),
            (
                "Paraphrase maps energy, relaxation, and curvature to SP, OPT, "
                "and HESS while preserving the initial-geometry sibling "
                "branches and exact OPT-to-HESS producer edge."
            ),
            "workflow.paraphrase-roundtrip.v1",
        ),
    )
    return tuple(
        build_case_spec(
            case_id=case_id,
            split=split,
            family=family,
            task=task,
            expected_observation=expected,
            deterministic_oracle_id=oracle,
            source_sha256s=(_H2O_NIST_XYZ_SHA256,),
        )
        for case_id, split, family, task, expected, oracle in rows
    )


__all__ = ["qwen_pyscf_cases_v1"]
