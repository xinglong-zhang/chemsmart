"""Deterministic host fixtures for the 14 program-management cases.

The factory binds repository-owned coordinates and project candidates to the
real :class:`CommandCompiledToolHostV1` constructor.  It does not execute a
chemistry engine, probe a provider, or manufacture conformance evidence.
Caller-supplied component conformance and target-interpreter receipts are the
only way a case can advance from semantic planning to compilation.

H0 and H1 receive the same immutable ``host_inputs`` and the same public
context packet.  Their only experimental difference is the tool surface that
the campaign runner selects.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
from typing import Any, Iterable, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.capabilities import (
    CapabilityQueryReceiptV1,
    CapabilityQueryStatus,
    CapabilityQueryV1,
    ProgramComponentConformanceReceiptV1,
    TrustedComputeEnvironmentReceiptV1,
    build_command_compiled_preview_overlay,
    load_program_capabilities,
    pyscf_compute_environment_target,
    query_capability,
    query_environment,
    resolve_engine_binding,
    resolve_program_binding,
)
from chemsmart.agent.cli_schema import build_live_click_schema
from chemsmart.agent.commands import build_scientific_identity_binding
from chemsmart.agent.experiments.program_management_context import (
    CampaignArtifactSlotV1,
    CampaignHostFixtureV1,
    CampaignPublicReceiptRefV1,
    CampaignToolInputV1,
    build_campaign_public_context,
    host_state_sha256,
    immutable_host_inputs,
    public_artifact,
)
from chemsmart.agent.inspection import inspect_generated_artifact
from chemsmart.agent.projects import validate_project_yaml
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.jobs.pyscf.writer import (
    APPLIED_SPEC_FIELDS,
    PySCFScriptWriter,
    write_pyscf_h5,
)
from chemsmart.settings.pyscf import PySCFProjectSettings


_WATER_RELATIVE_PATH = Path(
    "tests/data/agent/deepseek_program_management/h2o.xyz"
)
_WATER_SHA256 = (
    "42f720e0b1ae9883ad99e814488bf46093068bb386f007358841562629957045"
)
_OPEN_SHELL_RELATIVE_PATH = Path(
    "tests/data/StructuresTests/xyz/crest_best.xyz"
)
_OPEN_SHELL_SHA256 = (
    "551ee5c7eb07ba5927a6241b1c3d154c85f7af41ae668f88519a19fde68046d5"
)


_CPU_HF_PROJECT = b"""sp:
  ab_initio: hf
  basis: def2-svp
  density_fit: false
  engine: cpu
  freq: false
opt:
  ab_initio: hf
  basis: def2-svp
  density_fit: false
  engine: cpu
  freq: false
  opt_maxsteps: 100
  opt_solver: geometric
hess:
  ab_initio: hf
  basis: def2-svp
  density_fit: false
  engine: cpu
  freq: true
"""

_CPU_B3LYP_PROJECT = b"""sp:
  basis: def2-svp
  density_fit: false
  engine: cpu
  freq: false
  functional: b3lyp
opt:
  basis: def2-svp
  density_fit: false
  engine: cpu
  freq: false
  functional: b3lyp
  opt_maxsteps: 100
  opt_solver: geometric
hess:
  basis: def2-svp
  density_fit: false
  engine: cpu
  freq: true
  functional: b3lyp
"""

_GPU_HF_PROJECT = _CPU_HF_PROJECT.replace(b"engine: cpu", b"engine: gpu", 1)

_MIXED_ECP_PROJECT = b"""sp:
  ab_initio: hf
  basis: def2-svp
  density_fit: false
  engine: cpu
  freq: false
  gen_genecp_file: mixed-ecp.gbs
opt:
  ab_initio: hf
  basis: def2-svp
  density_fit: false
  engine: cpu
  freq: false
  gen_genecp_file: mixed-ecp.gbs
  opt_maxsteps: 100
  opt_solver: geometric
hess:
  ab_initio: hf
  basis: def2-svp
  density_fit: false
  engine: cpu
  freq: true
  gen_genecp_file: mixed-ecp.gbs
"""

_MISSING_BASIS_PROJECT = b"""sp:
  ab_initio: hf
  basis: null
  density_fit: false
  engine: cpu
  freq: false
opt:
  ab_initio: hf
  basis: null
  density_fit: false
  engine: cpu
  freq: false
  opt_maxsteps: 100
  opt_solver: geometric
hess:
  ab_initio: hf
  basis: null
  density_fit: false
  engine: cpu
  freq: true
"""


@dataclass(frozen=True)
class _CaseDesign:
    jobtypes: tuple[str, ...]
    engine: str
    project_key: str
    geometry_key: str = "water"
    charge: int = 0
    multiplicity: int = 1
    requested_program: str = "pyscf"
    method_family: str = "hf"
    method_name: str = "hf"
    basis_mode: str = "uniform"
    blocking_rule: str = ""


_CASE_DESIGNS: Mapping[str, _CaseDesign] = {
    "DS-PM-001": _CaseDesign(("sp",), "cpu", "cpu_hf"),
    "DS-PM-002": _CaseDesign(("sp",), "cpu", "cpu_hf"),
    "DS-PM-003": _CaseDesign(("opt", "hess"), "cpu", "cpu_hf"),
    "DS-PM-004": _CaseDesign(
        ("sp",),
        "cpu",
        "cpu_hf",
        requested_program="gaussian",
        blocking_rule="program.substitution.approval_required",
    ),
    "DS-PM-005": _CaseDesign(
        ("ts", "irc"),
        "cpu",
        "",
        requested_program="gaussian",
        blocking_rule="program.substitution.job_family_ineligible",
    ),
    "DS-PM-006": _CaseDesign(
        ("td",),
        "cpu",
        "",
        requested_program="gaussian",
        method_family="dft",
        method_name="b3lyp",
        blocking_rule="program.substitution.job_family_ineligible",
    ),
    "DS-PM-007": _CaseDesign(
        ("sp",), "cpu", "cpu_b3lyp", method_family="dft", method_name="b3lyp"
    ),
    "DS-PM-008": _CaseDesign(
        ("sp",),
        "cpu",
        "mixed_ecp",
        requested_program="gaussian",
        basis_mode="mixed_ecp",
        blocking_rule="pyscf.settings.unsupported",
    ),
    "DS-PM-009": _CaseDesign(
        ("sp",),
        "cpu",
        "cpu_hf",
        geometry_key="open_shell",
        multiplicity=2,
    ),
    "DS-PM-010": _CaseDesign(
        ("sp",),
        "gpu",
        "gpu_hf",
        blocking_rule="environment.compute_receipt.required",
    ),
    "DS-PM-011": _CaseDesign(
        ("sp",),
        "cpu",
        "cpu_hf",
        blocking_rule="command.shell_authority.rejected",
    ),
    "DS-PM-012": _CaseDesign(
        ("sp",),
        "cpu",
        "missing_basis",
        blocking_rule="project.loader.rejected",
    ),
    "DS-PM-013": _CaseDesign(("sp",), "cpu", "cpu_hf"),
    "DS-PM-014": _CaseDesign(("sp",), "cpu", "cpu_hf"),
}

_PROJECT_BYTES = {
    "cpu_hf": _CPU_HF_PROJECT,
    "cpu_b3lyp": _CPU_B3LYP_PROJECT,
    "gpu_hf": _GPU_HF_PROJECT,
    "mixed_ecp": _MIXED_ECP_PROJECT,
    "missing_basis": _MISSING_BASIS_PROJECT,
}


class ProgramManagementHostFixtureFactoryV1:
    """Materialize one reusable, arm-neutral host fixture per campaign case."""

    def __init__(
        self,
        *,
        source_tree_root: str | Path,
        materialization_root: str | Path,
        component_conformance_receipts: Iterable[
            ProgramComponentConformanceReceiptV1
        ] = (),
        compute_environment_receipts: Iterable[
            TrustedComputeEnvironmentReceiptV1
        ] = (),
    ) -> None:
        self.source_tree_root = Path(source_tree_root).resolve()
        self.materialization_root = Path(materialization_root).resolve()
        if not self.source_tree_root.is_dir():
            raise ContractError("campaign source-tree root is not a directory")
        self.materialization_root.mkdir(parents=True, exist_ok=True)
        if self.materialization_root.is_symlink():
            raise ContractError("campaign materialization root cannot be a symlink")
        self.registry = load_program_capabilities()
        self.live_schema = build_live_click_schema()
        self.conformance_receipts = tuple(component_conformance_receipts)
        self.compute_receipts = tuple(compute_environment_receipts)
        self.overlay = build_command_compiled_preview_overlay(
            self.registry,
            conformance_receipts=self.conformance_receipts,
            live_schema=self.live_schema,
        )
        self.environment_targets = (
            pyscf_compute_environment_target("cpu"),
            pyscf_compute_environment_target("gpu"),
        )
        self._validate_external_evidence()
        self._geometries = self._load_geometry_bindings()
        self._projects = self._materialize_project_bindings()
        self._cache: dict[str, CampaignHostFixtureV1] = {}

    def __call__(self, case: Any) -> CampaignHostFixtureV1:
        case_id = str(getattr(case, "case_id", ""))
        if case_id not in _CASE_DESIGNS:
            raise ContractError("campaign fixture case is not preregistered")
        cached = self._cache.get(case_id)
        if cached is None:
            cached = self._build_fixture(case_id)
            self._cache[case_id] = cached
        return cached

    def _validate_external_evidence(self) -> None:
        for receipt in self.conformance_receipts:
            if (
                receipt.registry_sha256 != self.registry.registry_sha256
                or receipt.live_cli_schema_sha256
                != self.live_schema.schema_sha256
            ):
                raise ContractError("external conformance receipt is stale")
        for receipt in self.compute_receipts:
            if receipt.program != "pyscf" or receipt.engine not in {"cpu", "gpu"}:
                raise ContractError("campaign accepts only PySCF compute receipts")

    def _load_geometry_bindings(self) -> dict[str, TrustedArtifactRefV1]:
        return {
            "water": _trusted_source_artifact(
                self.source_tree_root / _WATER_RELATIVE_PATH,
                expected_sha256=_WATER_SHA256,
                artifact_id="geometry.h2o.approved",
                kind="xyz",
            ),
            "open_shell": _trusted_source_artifact(
                self.source_tree_root / _OPEN_SHELL_RELATIVE_PATH,
                expected_sha256=_OPEN_SHELL_SHA256,
                artifact_id="geometry.open-shell.approved",
                kind="xyz",
            ),
        }

    def _materialize_project_bindings(self) -> dict[str, TrustedArtifactRefV1]:
        project_root = self.materialization_root / "projects"
        project_root.mkdir(parents=True, exist_ok=True)
        bindings = {}
        for key, content in sorted(_PROJECT_BYTES.items()):
            path = project_root / f"{key}.yaml"
            _write_exact_once(path, content)
            bindings[key] = _trusted_source_artifact(
                path,
                expected_sha256=file_sha256(path),
                artifact_id=f"project.pyscf.{key.replace('_', '-')}",
                kind="project_yaml",
            )
        return bindings

    def _build_fixture(self, case_id: str) -> CampaignHostFixtureV1:
        design = _CASE_DESIGNS[case_id]
        geometry = self._geometries[design.geometry_key]
        task_spec_sha256 = canonical_sha256(
            {
                "schema_version": "chemsmart.campaign-scientific-task.v1",
                "case_id": case_id,
                "geometry_sha256": geometry.sha256,
                "charge": design.charge,
                "multiplicity": design.multiplicity,
                "requested_program": design.requested_program,
                "selected_program_candidate": "pyscf",
                "engine": design.engine,
                "jobtypes": design.jobtypes,
                "method_family": design.method_family,
                "method_name": design.method_name,
                "basis": "def2-svp",
                "basis_mode": design.basis_mode,
            }
        )
        identity = build_scientific_identity_binding(
            task_spec_sha256=task_spec_sha256,
            geometry_artifact=geometry,
            charge=design.charge,
            multiplicity=design.multiplicity,
        )
        capabilities: list[CapabilityQueryReceiptV1] = []
        environments = []
        program_bindings = []
        engine_bindings = []
        project_validations = []
        project = self._projects.get(design.project_key)
        for jobtype in design.jobtypes:
            capability = query_capability(
                CapabilityQueryV1(
                    program="pyscf", jobtype=jobtype, engine=design.engine
                ),
                registry=self.registry,
                overlay=self.overlay,
                live_schema=self.live_schema,
            )
            environment = query_environment(
                capability,
                targets=self.environment_targets,
                compute_receipts=self.compute_receipts,
            )
            program_binding = resolve_program_binding(capability)
            engine_binding = resolve_engine_binding(
                program_binding, environment
            )
            capabilities.append(capability)
            environments.append(environment)
            program_bindings.append(program_binding)
            engine_bindings.append(engine_binding)
            if project is not None:
                project_validations.append(
                    validate_project_yaml(project, capability=capability)
                )

        artifacts = {geometry.artifact_id: geometry}
        if project is not None:
            artifacts[project.artifact_id] = project
        settings_objects: dict[str, Any] = {}
        run_receipts: dict[str, Mapping[str, Any]] = {}
        result_precheck = None
        if case_id in {"DS-PM-013", "DS-PM-014"}:
            result = self._materialize_seeded_result(
                case_id=case_id,
                geometry=geometry,
                project=project,
            )
            artifacts[result["artifact"].artifact_id] = result["artifact"]
            settings_objects[result["settings_id"]] = result["settings"]
            run_receipts[result["run_receipt_id"]] = result["run_receipt"]
            result_precheck = inspect_generated_artifact(
                program="pyscf",
                settings=result["settings"],
                artifact=result["artifact"],
                project_artifact=project,
                expected_receipt=result["run_receipt"],
            )

        host_values = {
            "artifacts": artifacts,
            "scientific_identities": {identity.binding_sha256: identity},
            "environment_targets": self.environment_targets,
            "compute_environment_receipts": self.compute_receipts,
            "component_conformance_receipts": self.conformance_receipts,
            "settings_objects": settings_objects,
            "run_receipts": run_receipts,
            "scientific_claim_evidence": {},
            "functional_equivalence_receipts": {},
            "substitution_approvals": {},
            "capability_receipts": {
                item.receipt_sha256: item for item in capabilities
            },
            "environment_receipts": {
                item.receipt_sha256: item for item in environments
            },
            "program_binding_receipts": {
                item.binding_sha256: item for item in program_bindings
            },
            "engine_binding_receipts": {
                item.binding_sha256: item for item in engine_bindings
            },
            "project_validation_receipts": {
                item.receipt_sha256: item for item in project_validations
            },
            "registry": self.registry,
            "live_schema": self.live_schema,
        }
        immutable = immutable_host_inputs(host_values)
        host_digest = host_state_sha256(immutable)
        public = self._build_public_context(
            case_id=case_id,
            design=design,
            task_spec_sha256=task_spec_sha256,
            identity=identity,
            geometry=geometry,
            project=project,
            capabilities=tuple(capabilities),
            environments=tuple(environments),
            program_bindings=tuple(program_bindings),
            engine_bindings=tuple(engine_bindings),
            project_validations=tuple(project_validations),
            result_precheck=result_precheck,
            result_fixture=(result if case_id in {"DS-PM-013", "DS-PM-014"} else None),
            host_state_digest=host_digest,
        )
        return CampaignHostFixtureV1(
            schema_version="chemsmart.campaign-host-fixture.v1",
            case_id=case_id,
            host_inputs=immutable,
            public_context=public,
            host_state_sha256=host_digest,
        )

    def _build_public_context(
        self,
        *,
        case_id: str,
        design: _CaseDesign,
        task_spec_sha256: str,
        identity: Any,
        geometry: TrustedArtifactRefV1,
        project: TrustedArtifactRefV1 | None,
        capabilities: tuple[Any, ...],
        environments: tuple[Any, ...],
        program_bindings: tuple[Any, ...],
        engine_bindings: tuple[Any, ...],
        project_validations: tuple[Any, ...],
        result_precheck: Any,
        result_fixture: Mapping[str, Any] | None,
        host_state_digest: str,
    ):
        public_artifacts = [
            public_artifact(
                geometry,
                purpose=(
                    "checked-in exact non-generated coordinate frame"
                    if design.geometry_key == "open_shell"
                    else "user-approved exact coordinate frame"
                ),
                provenance_status=(
                    "checked_in_non_generated_source"
                    if design.geometry_key == "open_shell"
                    else "approved_exact_source"
                ),
            )
        ]
        if project is not None:
            validation_green = bool(project_validations) and all(
                item.status == "valid" for item in project_validations
            )
            project_is_seeded_red = design.project_key in {
                "mixed_ecp",
                "missing_basis",
            }
            public_artifacts.append(
                public_artifact(
                    project,
                    purpose="typed PySCF project fixture",
                    provenance_status=(
                        "seeded_negative_fixture"
                        if project_is_seeded_red
                        else "validated_project_fixture"
                        if validation_green
                        else "host_bound_project_candidate"
                    ),
                )
            )
        if result_fixture is not None:
            public_artifacts.append(
                public_artifact(
                    result_fixture["artifact"],
                    purpose="synthetic structured result for a seeded defect",
                    provenance_status="seeded_negative_fixture",
                )
            )

        receipt_refs = []
        for ordinal, receipt in enumerate(capabilities):
            receipt_refs.append(
                _receipt_ref(
                    role=f"capability.{ordinal}",
                    digest=receipt.receipt_sha256,
                    status=receipt.status.value,
                    rule_ids=receipt.rule_ids,
                    semantic={
                        "engine": receipt.query.engine,
                        "jobtype": receipt.query.jobtype,
                        "program": receipt.query.program,
                    },
                )
            )
        for ordinal, receipt in enumerate(environments):
            receipt_refs.append(
                _receipt_ref(
                    role=f"environment.{ordinal}",
                    digest=receipt.receipt_sha256,
                    status=receipt.status.value,
                    rule_ids=receipt.rule_ids,
                    semantic={"engine": receipt.engine, "program": receipt.program},
                )
            )
        for ordinal, receipt in enumerate(program_bindings):
            receipt_refs.append(
                _receipt_ref(
                    role=f"program-binding.{ordinal}",
                    digest=receipt.binding_sha256,
                    status=receipt.state,
                    rule_ids=receipt.rule_ids,
                    semantic={
                        "requested_program": receipt.requested_program,
                        "selected_program": receipt.selected_program,
                    },
                )
            )
        for ordinal, receipt in enumerate(engine_bindings):
            receipt_refs.append(
                _receipt_ref(
                    role=f"engine-binding.{ordinal}",
                    digest=receipt.binding_sha256,
                    status=receipt.state,
                    rule_ids=receipt.rule_ids,
                    semantic={
                        "engine": receipt.engine,
                        "execution_ready": receipt.execution_ready,
                    },
                )
            )
        for ordinal, receipt in enumerate(project_validations):
            receipt_refs.append(
                _receipt_ref(
                    role=f"project-validation.{ordinal}",
                    digest=receipt.receipt_sha256,
                    status=receipt.status,
                    rule_ids=receipt.rule_ids,
                    semantic={"jobtype": receipt.jobtype, "program": receipt.program},
                )
            )
        if result_precheck is not None:
            receipt_refs.append(
                _receipt_ref(
                    role="result-precheck",
                    digest=result_precheck.receipt_sha256,
                    status=result_precheck.status,
                    rule_ids=tuple(
                        sorted(set(item.rule_id for item in result_precheck.findings))
                    ),
                    semantic={"artifact_id": result_precheck.artifact_id},
                )
            )

        slots = [
            CampaignArtifactSlotV1(
                slot_id="geometry.initial",
                artifact_class="xyz",
                state="external_bound",
                artifact_id=geometry.artifact_id,
                artifact_sha256=geometry.sha256,
                producer_node_id="",
                producer_output_id="",
            )
        ]
        if case_id == "DS-PM-003":
            slots.append(
                CampaignArtifactSlotV1(
                    slot_id="geometry.optimized",
                    artifact_class="xyz",
                    state="future_unresolved",
                    artifact_id="",
                    artifact_sha256="",
                    producer_node_id="node.opt",
                    producer_output_id="optimized_geometry",
                )
            )

        actions = self._common_actions(
            case_id=case_id,
            design=design,
            task_spec_sha256=task_spec_sha256,
            identity=identity,
            geometry=geometry,
            project=project,
            capabilities=capabilities,
            engine_bindings=engine_bindings,
            project_validations=project_validations,
            result_fixture=result_fixture,
        )
        missing = []
        if any(
            item.status is CapabilityQueryStatus.REFERENCE_ONLY
            for item in capabilities
        ):
            missing.append("current component-conformance receipt")
        if any(item.status.value != "available" for item in environments):
            missing.append(
                f"complete observed PySCF {design.engine} environment receipt"
            )
        if design.blocking_rule:
            missing.append(design.blocking_rule)
        if case_id == "DS-PM-003":
            missing.append("optimized geometry artifact and producer receipt")
        if case_id == "DS-PM-014":
            missing.append("deterministic required-property result validator")

        joined = capabilities[0].joined_capability_sha256
        facts = {
            "assurance_ceiling": "semantic_plan_until_all_receipts_are_green",
            "basis": "def2-svp",
            "basis_mode": design.basis_mode,
            "charge": design.charge,
            "coordinate_units": "angstrom",
            "engine": design.engine,
            "jobtypes": design.jobtypes,
            "method_family": design.method_family,
            "method_name": design.method_name,
            "multiplicity": design.multiplicity,
            "requested_program": design.requested_program,
            "selected_program_candidate": "pyscf",
            "scientific_identity_sha256": identity.binding_sha256,
        }
        if case_id == "DS-PM-009":
            facts["pyscf_spin"] = 1
        if case_id == "DS-PM-014":
            facts["required_property"] = "dipole_moment"
            facts["seeded_property_state"] = "missing"
        return build_campaign_public_context(
            case_id=case_id,
            task_spec_sha256=task_spec_sha256,
            registry_sha256=self.registry.registry_sha256,
            live_cli_schema_sha256=self.live_schema.schema_sha256,
            joined_capability_sha256=joined,
            support_overlay_sha256=self.overlay.overlay_sha256,
            scientific_facts=facts,
            artifacts=tuple(public_artifacts),
            artifact_slots=tuple(slots),
            receipts=tuple(receipt_refs),
            next_actions=actions,
            expected_artifact_classes=_expected_artifact_classes(design),
            missing_evidence=tuple(missing),
            host_state_sha256=host_state_digest,
        )

    def _common_actions(
        self,
        *,
        case_id: str,
        design: _CaseDesign,
        task_spec_sha256: str,
        identity: Any,
        geometry: TrustedArtifactRefV1,
        project: TrustedArtifactRefV1 | None,
        capabilities: tuple[Any, ...],
        engine_bindings: tuple[Any, ...],
        project_validations: tuple[Any, ...],
        result_fixture: Mapping[str, Any] | None,
    ) -> tuple[CampaignToolInputV1, ...]:
        if result_fixture is not None:
            return (
                _action(
                    "inspect.seeded-result",
                    "inspect_calculation_artifact",
                    {
                        "artifact_id": result_fixture["artifact"].artifact_id,
                        "program": "pyscf",
                        "project_artifact_id": project.artifact_id,
                        "run_receipt_id": result_fixture["run_receipt_id"],
                        "settings_id": result_fixture["settings_id"],
                    },
                    "verification",
                ),
            )

        workflow_nodes = _workflow_nodes(
            design=design,
            geometry_artifact_id=geometry.artifact_id,
        )
        actions = [
            _action(
                "plan.workflow",
                "plan_command_workflow",
                {
                    "nodes": workflow_nodes,
                    "task_spec_id": task_spec_sha256,
                    "workflow_id": f"workflow.{case_id.lower()}",
                },
                "planning",
            )
        ]
        if case_id == "DS-PM-001":
            return ()
        compile_excluded = {
            "DS-PM-003",
            "DS-PM-004",
            "DS-PM-005",
            "DS-PM-006",
            "DS-PM-008",
            "DS-PM-010",
            "DS-PM-011",
            "DS-PM-012",
        }
        if case_id in compile_excluded or project is None:
            return tuple(actions)
        capability = capabilities[0]
        binding = engine_bindings[0]
        validation = project_validations[0]
        if (
            capability.status
            not in {CapabilityQueryStatus.SUPPORTED, CapabilityQueryStatus.PREVIEW_ONLY}
            or binding.state == "blocked"
            or validation.status != "valid"
        ):
            return tuple(actions)
        fields = {
            "capability_receipt_sha256": capability.receipt_sha256,
            "charge": design.charge,
            "engine_binding_sha256": binding.binding_sha256,
            "execution_target": "run",
            "input_artifact_id": geometry.artifact_id,
            "jobtype": design.jobtypes[0],
            "multiplicity": design.multiplicity,
            "node_id": f"node.{design.jobtypes[0]}",
            "program": "pyscf",
            "project_artifact_id": project.artifact_id,
            "project_validation_receipt_sha256": validation.receipt_sha256,
            "scientific_identity_sha256": identity.binding_sha256,
        }
        actions.append(
            _action(
                "compile.first-node",
                "synthesize_command",
                fields,
                "compile_bound",
            )
        )
        return tuple(actions)

    def _materialize_seeded_result(
        self,
        *,
        case_id: str,
        geometry: TrustedArtifactRefV1,
        project: TrustedArtifactRefV1 | None,
    ) -> Mapping[str, Any]:
        if project is None:
            raise ContractError("seeded PySCF result requires a project")
        settings: PySCFJobSettings = PySCFProjectSettings.from_yaml(
            project.path
        ).sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.title = f"{case_id} seeded verifier fixture"
        settings.validate()
        symbols, positions = _read_xyz(geometry)
        run_id = f"fixture-{case_id.lower()}"
        run_nonce = canonical_sha256({"case_id": case_id, "role": "run_nonce"})
        input_geometry_sha256 = canonical_sha256(
            {"symbols": symbols, "positions": positions, "unit": "Angstrom"}
        )
        spec = {
            "run_id": run_id,
            "run_nonce": run_nonce,
            "label": case_id.lower(),
            "title": settings.title,
            "program": "pyscf",
            "jobtype": "sp",
            "engine": "cpu",
            "stages": ["scf"],
            "symbols": symbols,
            "positions": positions,
            "unit": "Angstrom",
            "xc": None,
            "ab_initio": "hf",
            "method": "hf",
            "basis": "def2-svp",
            "charge": 0,
            "spin": 0,
            "multiplicity": 1,
            "dispersion": None,
            "density_fit": False,
            "aux_basis": None,
            "defgrid": None,
            "atom_grid": None,
            "scf_tol": None,
            "scf_maxiter": None,
            "solvent_model": None,
            "solvent_call": None,
            "solvent_method": None,
            "solvent_id": None,
            "solvent_eps": None,
            "solvent_lebedev_order": None,
            "opt_solver": "geometric",
            "opt_maxsteps": 100,
            "num_threads": 1,
            "max_memory_mb": 1024,
            "input_geometry_sha256": input_geometry_sha256,
            "input_artifact_kind": "xyz",
            "input_artifact_sha256": geometry.sha256,
            "requested_settings_sha256": "",
            "applied_settings_sha256": None,
            "settings_digest": "",
        }
        requested_digest = PySCFScriptWriter.settings_digest(spec)
        spec["requested_settings_sha256"] = requested_digest
        spec["settings_digest"] = requested_digest
        applied_payload = {key: spec.get(key) for key in APPLIED_SPEC_FIELDS}
        applied_digest = PySCFScriptWriter.settings_digest(applied_payload)
        spec["applied_settings_sha256"] = applied_digest
        script_sha256 = canonical_sha256(
            {"case_id": case_id, "fixture": "synthetic-driver-source-v1"}
        )
        input_receipt_sha256 = canonical_sha256(
            {"case_id": case_id, "fixture": "synthetic-input-receipt-v1"}
        )
        environment_receipt_sha256 = canonical_sha256(
            {"case_id": case_id, "fixture": "synthetic-environment-receipt-v1"}
        )
        provenance_applied = (
            canonical_sha256({"case_id": case_id, "defect": "applied-settings"})
            if case_id == "DS-PM-013"
            else applied_digest
        )
        provenance = {
            "run_id": run_id,
            "run_nonce": run_nonce,
            "engine": "cpu",
            "settings_digest": requested_digest,
            "requested_settings_sha256": requested_digest,
            "applied_settings_sha256": provenance_applied,
            "input_geometry_sha256": input_geometry_sha256,
            "input_artifact_kind": "xyz",
            "input_artifact_sha256": geometry.sha256,
            "project_yaml_digest": project.sha256,
            "script_sha256": script_sha256,
            "input_receipt_sha256": input_receipt_sha256,
            "environment_receipt_sha256": environment_receipt_sha256,
            "fixture_only": True,
        }
        properties = {
            "dipole_moment": (
                {
                    "status": "unavailable",
                    "failure": {
                        "type": "seeded_missing_property",
                        "message": "required dipole dataset intentionally absent",
                    },
                }
                if case_id == "DS-PM-014"
                else {"status": "ok"}
            )
        }
        status = {
            "normal_termination": True,
            "failure": None,
            "fixture_only": True,
            "properties": properties,
        }
        results: dict[str, Any] = {
            "energies": [-75.0],
            "positions": positions,
        }
        if case_id != "DS-PM-014":
            results["dipole_moment"] = [0.0, 0.0, 0.0]
        result_root = self.materialization_root / "results"
        result_root.mkdir(parents=True, exist_ok=True)
        result_path = result_root / f"{case_id.lower()}.h5"
        _write_h5_once(
            result_path,
            spec=spec,
            provenance=provenance,
            status=status,
            results=results,
        )
        artifact = _trusted_source_artifact(
            result_path,
            expected_sha256=file_sha256(result_path),
            artifact_id=f"result.{case_id.lower()}",
            kind="pyscf_hdf5",
        )
        expected = {
            "script_sha256": script_sha256,
            "input_receipt_sha256": input_receipt_sha256,
            "environment_receipt_sha256": environment_receipt_sha256,
            "input_geometry_sha256": input_geometry_sha256,
            "input_artifact_kind": "xyz",
            "input_artifact_sha256": geometry.sha256,
            "requested_settings_sha256": requested_digest,
            "project_yaml_digest": project.sha256,
            "require_applied_settings_sha256": True,
            "fixture_only": True,
        }
        if case_id == "DS-PM-014":
            expected["required_properties"] = ("dipole_moment",)
        return {
            "artifact": artifact,
            "settings": settings,
            "settings_id": f"settings.{case_id.lower()}",
            "run_receipt": expected,
            "run_receipt_id": f"run-receipt.{case_id.lower()}",
        }


def _workflow_nodes(
    *, design: _CaseDesign, geometry_artifact_id: str
) -> list[dict[str, Any]]:
    nodes = []
    for ordinal, jobtype in enumerate(design.jobtypes):
        node_id = f"node.{jobtype}"
        dependencies: list[str] = []
        unresolved = []
        inputs: list[dict[str, str]]
        if ordinal == 0:
            inputs = [
                {
                    "binding_id": "geometry.initial",
                    "artifact_id": geometry_artifact_id,
                    "artifact_class": "xyz",
                    "producer_node_id": "",
                    "producer_output_id": "",
                }
            ]
        else:
            producer = f"node.{design.jobtypes[ordinal - 1]}"
            dependencies = [producer]
            inputs = [
                {
                    "binding_id": "geometry.optimized",
                    "artifact_class": "xyz",
                    "producer_node_id": producer,
                    "producer_output_id": "optimized_geometry",
                }
            ]
            unresolved.append("input.optimized_geometry_hash")
        outputs = [
            {"output_id": "structured_result", "artifact_class": "pyscf_hdf5"}
        ]
        if jobtype == "opt" or ordinal < len(design.jobtypes) - 1:
            outputs.insert(
                0,
                {"output_id": "optimized_geometry", "artifact_class": "xyz"},
            )
        if design.blocking_rule:
            unresolved.append(design.blocking_rule)
        nodes.append(
            {
                "node_id": node_id,
                "program": "pyscf",
                "jobtype": jobtype,
                "project_role": "primary",
                "dependencies": dependencies,
                "inputs": inputs,
                "expected_outputs": outputs,
                "unresolved_fields": sorted(set(unresolved)),
            }
        )
    return nodes


def _expected_artifact_classes(design: _CaseDesign) -> tuple[str, ...]:
    values = {"pyscf_hdf5"}
    if "opt" in design.jobtypes:
        values.add("xyz")
    return tuple(sorted(values))


def _receipt_ref(
    *,
    role: str,
    digest: str,
    status: str,
    rule_ids: Iterable[str],
    semantic: Mapping[str, Any],
) -> CampaignPublicReceiptRefV1:
    return CampaignPublicReceiptRefV1(
        role=role,
        receipt_sha256=digest,
        status=str(status),
        rule_ids=tuple(sorted(set(str(item) for item in rule_ids))),
        semantic_fields=tuple(sorted(semantic.items())),
    )


def _action(
    action_id: str,
    tool_name: str,
    fields: Mapping[str, Any],
    assurance_level: str,
) -> CampaignToolInputV1:
    return CampaignToolInputV1(
        action_id=action_id,
        tool_name=tool_name,
        fields=tuple(sorted(fields.items())),
        assurance_level=assurance_level,
    )


def _trusted_source_artifact(
    path: Path, *, expected_sha256: str, artifact_id: str, kind: str
) -> TrustedArtifactRefV1:
    if path.is_symlink():
        raise ContractError("campaign fixture source cannot be a symlink")
    resolved = path.resolve()
    if not resolved.is_file():
        raise ContractError("campaign fixture source must be a regular file")
    observed = file_sha256(resolved)
    if observed != expected_sha256:
        raise ContractError("campaign fixture source digest mismatch")
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=observed,
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )


def _write_exact_once(path: Path, content: bytes) -> None:
    if path.exists():
        if path.is_symlink() or path.read_bytes() != content:
            raise ContractError("existing campaign fixture differs from source")
        return
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            descriptor = -1
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _write_h5_once(
    path: Path,
    *,
    spec: Mapping[str, Any],
    provenance: Mapping[str, Any],
    status: Mapping[str, Any],
    results: Mapping[str, Any],
) -> None:
    if path.exists():
        raise ContractError(
            "seeded HDF5 already exists; use a fresh materialization root"
        )
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    if temporary.exists():
        raise ContractError("campaign HDF5 temporary already exists")
    try:
        write_pyscf_h5(
            temporary,
            spec=spec,
            provenance=provenance,
            status=status,
            results=results,
        )
        if path.exists():
            raise ContractError("campaign HDF5 destination appeared concurrently")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _read_xyz(
    artifact: TrustedArtifactRefV1,
) -> tuple[list[str], list[list[float]]]:
    if file_sha256(artifact.path) != artifact.sha256:
        raise ContractError("XYZ source changed after binding")
    lines = Path(artifact.path).read_text(encoding="utf-8").splitlines()
    try:
        count = int(lines[0].strip())
    except (IndexError, ValueError) as exc:
        raise ContractError("campaign XYZ atom count is invalid") from exc
    rows = lines[2 : 2 + count]
    if len(rows) != count:
        raise ContractError("campaign XYZ is truncated")
    symbols = []
    positions = []
    for row in rows:
        fields = row.split()
        if len(fields) != 4:
            raise ContractError("campaign XYZ row is malformed")
        symbols.append(fields[0])
        positions.append([float(value) for value in fields[1:]])
    return symbols, positions


__all__ = ["ProgramManagementHostFixtureFactoryV1"]
