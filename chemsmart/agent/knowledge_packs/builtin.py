"""Small advisory catalog; validators and the live CLI remain authoritative."""

from __future__ import annotations

from chemsmart.agent.knowledge_packs.contracts import build_pack
from chemsmart.agent.skills.conventions import BUILTIN_CONVENTION_RULE_SHA256S

#: Every program pack surfaces the cross-program skills; program bodies stay
#: resolvable through :mod:`chemsmart.agent.skills`.
#:
#: All three are program-neutral, so they are shared rather than split across
#: packs: activation fires per (program, engine), and a skill reaches the prompt
#: index only through an activated pack, so a chemistry-general skill given its
#: own pack would need its own activation terms and would simply never fire.
#:
#: ``scientific-conventions`` covers how a result is expressed,
#: ``method-adequacy`` whether the method could answer the question at all, and
#: ``typed-analysis-contract`` how to express the analysis so the host can carry
#: the evidence.
_CONVENTIONS_SKILL = (
    "method-adequacy",
    "scientific-conventions",
    "typed-analysis-contract",
)


BUILTIN_PROGRAM_PACKS = (
    build_pack(
        pack_id="gaussian-advisory",
        pack_version="1.0.0",
        target_program="gaussian",
        activation_terms=("gaussian", "g16", "oniom", "genecp"),
        exclusions=("periodic", "solid state"),
        advisory_topics=(
            "basis and ECP mapping",
            "frequency and stationary-point checks",
            "route ownership remains in project YAML",
        ),
        source_ids=("gaussian-16-reference",),
        skill_ids=_CONVENTIONS_SKILL,
        convention_rule_sha256s=BUILTIN_CONVENTION_RULE_SHA256S,
    ),
    build_pack(
        pack_id="orca-advisory",
        pack_version="1.0.0",
        target_program="orca",
        activation_terms=("orca", "neb", "dlpno", "cpcm"),
        exclusions=("periodic", "solid state"),
        advisory_topics=(
            "auxiliary basis evidence",
            "NEB endpoint requirements",
            (
                "frequency properties are requested with freq: true in project "
                "YAML rather than a hess CLI job; use opt or ts only when a "
                "geometry search is intended, and use sp for a fixed-geometry "
                "frequency or VPT2 calculation"
            ),
            "structured convergence checks",
        ),
        source_ids=("orca-manual-6.1",),
        skill_ids=_CONVENTIONS_SKILL,
        convention_rule_sha256s=BUILTIN_CONVENTION_RULE_SHA256S,
    ),
    build_pack(
        pack_id="xtb-reference-advisory",
        pack_version="1.0.0",
        target_program="xtb",
        activation_terms=("xtb", "gfn1", "gfn2", "alpb", "gbsa"),
        exclusions=("ab initio replacement", "periodic"),
        advisory_topics=(
            "semiempirical method identity",
            "solvent model and solvent identifier pairing",
            "no silent substitution for DFT or ab initio intent",
        ),
        source_ids=("xtb-docs-6.7.1",),
        skill_ids=_CONVENTIONS_SKILL,
        convention_rule_sha256s=BUILTIN_CONVENTION_RULE_SHA256S,
        reference_only=True,
    ),
    build_pack(
        pack_id="pyscf-cpu-advisory",
        pack_version="1.0.0",
        target_program="pyscf",
        target_engine="cpu",
        activation_terms=("pyscf", "libxc", "density fitting", "hessian"),
        exclusions=("periodic", "double hybrid"),
        advisory_topics=(
            "multiplicity to PySCF spin conversion",
            "target-interpreter dependency evidence",
            "HDF5 provenance verification",
        ),
        source_ids=("pyscf-docs-2.14.0",),
        skill_ids=_CONVENTIONS_SKILL,
        convention_rule_sha256s=BUILTIN_CONVENTION_RULE_SHA256S,
    ),
    build_pack(
        pack_id="gpu4pyscf-advisory",
        pack_version="1.0.0",
        target_program="pyscf",
        target_engine="gpu",
        activation_terms=("gpu4pyscf", "cuda", "cupy", "cutensor", "gpu"),
        exclusions=("double hybrid", "unsupported gpu solvent"),
        advisory_topics=(
            "GPU4PySCF and CuPy version binding",
            "CUDA, cuTENSOR, and device evidence",
            "CPU and GPU results are not assumed bit-identical",
        ),
        source_ids=("gpu4pyscf-docs-1.8.0",),
        skill_ids=_CONVENTIONS_SKILL,
        convention_rule_sha256s=BUILTIN_CONVENTION_RULE_SHA256S,
    ),
)


__all__ = ["BUILTIN_PROGRAM_PACKS"]
