"""Plan-first live ChemSmart Agent sessions.

This module is the small composition root behind ``chemsmart agent plan`` and
``chemsmart agent run``.  It binds a provider session to exact, pre-existing
XYZ artifacts, the live Click schema, observed program conformance, and a
private Runtime V2 event stream.  It never generates coordinates, engine
inputs, shell commands, or scientific readiness decisions.

The execution profile is deliberately progressive.  Until the tool host has
an approval-bound execution composition API, ``agent run`` uses the complete
planning/preview path and reports ``waiting_for_approval`` rather than
pretending that an engine was invoked.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import inspect
import json
import math
import os
from pathlib import Path
import shutil
import sys
import uuid
from typing import Any, Iterable, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_json,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.adaptive_api_campaign import (
    AdaptiveHypothesisV1,
    AdaptiveNetworkBudgetV1,
)
from chemsmart.agent.api_access import load_secret_lease
from chemsmart.agent.bootstrap import (
    bootstrap_program_conformance,
    probe_python_compute_environment,
)
from chemsmart.agent.capabilities import (
    EnvironmentTargetV1,
    ProgramComponentConformanceReceiptV1,
    TrustedComputeEnvironmentReceiptV1,
    build_program_component_conformance_receipt,
    load_program_capabilities,
)
from chemsmart.agent.cli_schema import build_live_click_schema
from chemsmart.agent.projects import project_document, render_project_yaml
from chemsmart.agent.execution import (
    ApprovedNodeBindingV1,
    ExecutionResourceSpecV1,
    ProducerEdgeRuleV1,
    WorkflowExecutionApprovalV1,
    build_execution_resource_spec,
)
from chemsmart.agent.runtime.contracts import (
    ResourceBudgetV1,
    TaskEnvelopeV1,
    TaskPhase,
)
from chemsmart.agent.runtime.deepseek import (
    DEEPSEEK_OFFICIAL_ENDPOINT,
    DEEPSEEK_V4_FLASH_CONTEXT_TOKENS,
    DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS,
    DeepSeekV4FlashConfigV1,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.services.unified_session import UnifiedSessionRunner
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.tool_specs import (
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
)


_SESSION_WALL_TIME_SECONDS = 90 * 60
_MAX_TOOL_CALLS = 256
_PYSCF_INTERPRETER = Path(
    os.environ.get("CHEMSMART_PYSCF_INTERPRETER", sys.executable)
).expanduser().resolve()
_PRIVATE_ROOT_NAME = ".chemsmart-agent"


@dataclass(frozen=True)
class _XyzObservation:
    artifact: TrustedArtifactRefV1
    atom_count: int
    symbols: tuple[str, ...]

    def public_record(self) -> dict[str, Any]:
        return {
            "artifact_id": self.artifact.artifact_id,
            "artifact_class": self.artifact.kind,
            "sha256": self.artifact.sha256,
            "size_bytes": self.artifact.size_bytes,
            "atom_count": self.atom_count,
            "symbols": self.symbols,
            "provenance_status": "workspace_exact_user_approved",
        }


@dataclass(frozen=True)
class LiveAgentSessionResultV1:
    """Path-free public projection of one local live session."""

    schema_version: str
    session_id: str
    task_spec_sha256: str
    terminal_state: str
    execution_requested: bool
    execution_profile_status: str
    final_text: str
    artifact_records: tuple[dict[str, Any], ...]
    conformance_records: tuple[dict[str, Any], ...]
    public_transcript: tuple[dict[str, Any], ...]
    successful_tool_calls: int
    failed_tool_calls: int
    event_stream_head_sha256: str
    result_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.live-agent-session-result.v1":
            raise ContractError("unsupported live session result schema")
        if self.terminal_state not in {
            "complete",
            "planned",
            "failed",
            "blocked",
            "waiting_for_approval",
        }:
            raise ContractError("invalid live session terminal state")
        body = self._body()
        if self.result_sha256 != canonical_sha256(body):
            raise ContractError("live session result digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            key: value
            for key, value in asdict(self).items()
            if key != "result_sha256"
        }

    def public_summary_json(self) -> str:
        """Return the exact visible result; no local path is included."""

        return canonical_json({**self._body(), "result_sha256": self.result_sha256})


def run_live_agent_session(
    *,
    task: str,
    provider: str,
    secret_file: str | Path,
    workspace: str | Path,
    execution_enabled: bool,
    approval_file: str | Path | None,
) -> LiveAgentSessionResultV1:
    """Run one DeepSeek planning session over exact workspace artifacts.

    A real provider request is made only after the workspace, coordinate
    artifacts, live schema, program previews, and provider budget have been
    bound.  The private run directory contains the append-only event stream
    and sanitized public transcript; it is intentionally outside the source
    tree when the CLI is used as documented.
    """

    normalized_provider = str(provider).strip().lower()
    if normalized_provider != "deepseek":
        raise ContractError("the live v1 driver supports only DeepSeek")
    task = str(task).strip()
    if not task:
        raise ContractError("live agent task must not be empty")
    workspace_path = _validated_workspace(workspace)
    observations = _scan_xyz_artifacts(workspace_path)
    task_spec_sha256 = _task_spec_sha256(task, observations)
    session_id = _session_id(task_spec_sha256)
    run_directory = _private_run_directory(workspace_path, session_id)

    if not observations:
        return _local_result(
            session_id=session_id,
            task_spec_sha256=task_spec_sha256,
            terminal_state="blocked",
            execution_requested=execution_enabled,
            execution_profile_status="not_started",
            final_text=(
                "No exact XYZ artifact is present in the approved workspace. "
                "Add a user-approved coordinate file; coordinates were not generated."
            ),
        )

    registry = load_program_capabilities()
    live_schema = build_live_click_schema()
    conformance, conformance_records = _bootstrap_conformance(
        run_directory=run_directory,
        input_artifact=observations[0].artifact,
        registry_sha256=registry.registry_sha256,
        live_schema=live_schema,
    )
    environment_targets, compute_receipts, environment_records = (
        _observe_environments()
    )
    conformance_records = tuple(
        sorted((*conformance_records, *environment_records), key=_record_sort_key)
    )

    execution_profile_ready = _execution_composition_available()
    use_execution_surface = bool(
        execution_enabled and approval_file and execution_profile_ready
    )
    approved_project_records: tuple[dict[str, Any], ...] = ()
    approved_workflow_record: dict[str, Any] = {}
    surface = (
        build_approved_execution_tool_surface(registry)
        if use_execution_surface
        else build_command_compiled_tool_surface(registry)
    )

    event_store = RuntimeEventStore(
        run_directory / "events.jsonl", session_id=session_id
    )
    host_kwargs: dict[str, Any] = {
        "event_store": event_store,
        "artifacts": {
            item.artifact.artifact_id: item.artifact for item in observations
        },
        "environment_targets": environment_targets,
        "compute_environment_receipts": compute_receipts,
        "component_conformance_receipts": conformance,
        "tool_surface": surface,
        "registry": registry,
        "live_schema": live_schema,
        "task_spec_sha256s": (task_spec_sha256,),
        # Preview candidates stay inside this session's private workspace.
        # The execution composition below replaces this with the exact
        # user-approved workflow workspace.
        "approved_workspace": run_directory,
    }
    if use_execution_surface:
        execution_inputs = _execution_composition_inputs(
            host_type=CommandCompiledToolHostV1,
            workspace=workspace_path,
            run_directory=run_directory,
            approval_file=Path(approval_file),
            task_spec_sha256=task_spec_sha256,
        )
        approved_workflow_record = _public_workflow_approval(
            execution_inputs["workflow_execution_approval"]
        )
        approved_projects = execution_inputs.pop("approved_project_artifacts")
        host_kwargs["artifacts"].update(
            {item.artifact_id: item for item in approved_projects}
        )
        approved_project_records = tuple(
            {
                "artifact_id": item.artifact_id,
                "artifact_class": item.kind,
                "sha256": item.sha256,
                "size_bytes": item.size_bytes,
                "approval_status": "workflow_bound",
            }
            for item in approved_projects
        )
        host_kwargs.update(execution_inputs)
    host = CommandCompiledToolHostV1(**host_kwargs)

    context = _public_context(
        task=task,
        task_spec_sha256=task_spec_sha256,
        observations=observations,
        conformance_records=conformance_records,
        registry_sha256=registry.registry_sha256,
        live_schema_sha256=live_schema.schema_sha256,
        execution_requested=execution_enabled,
        execution_available=use_execution_surface,
        approved_project_records=approved_project_records,
        approved_workflow_record=approved_workflow_record,
    )
    messages = [
        {
            "role": "system",
            "content": _system_prompt(approved_workflow_record),
        },
        {"role": "user", "content": canonical_json(context)},
    ]
    network_budget = _network_budget()
    hypothesis = _hypothesis(
        session_id=session_id,
        messages=messages,
        task_spec_sha256=task_spec_sha256,
        artifact_sha256s=tuple(item.artifact.sha256 for item in observations),
        tool_schema_sha256=surface.tool_schema_sha256,
        network_budget=network_budget,
        execution_requested=execution_enabled,
    )
    envelope = _task_envelope(
        session_id=session_id,
        messages=messages,
        task_spec_sha256=task_spec_sha256,
        workspace_records=tuple(item.public_record() for item in observations),
        tool_schema_sha256=surface.tool_schema_sha256,
        execution_enabled=use_execution_surface,
    )
    lease = load_secret_lease(
        provider=normalized_provider,
        path=secret_file,
        ttl_seconds=_SESSION_WALL_TIME_SECONDS + 60,
    )
    loop_result = UnifiedSessionRunner(
        host=host,
        event_store=event_store,
        credential_lease=lease,
        provider_config=DeepSeekV4FlashConfigV1(
            max_output_tokens=DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS
        ),
    ).run(
        messages=messages,
        envelope=envelope,
        hypothesis=hypothesis,
        network_budget=network_budget,
    )

    execution_status = (
        "approved_profile_active"
        if use_execution_surface
        else "not_requested"
        if not execution_enabled
        else "planning_complete_execution_profile_unavailable"
    )
    terminal_state = loop_result.terminal_state
    final_text = loop_result.final_text
    if execution_enabled and not use_execution_surface:
        terminal_state = "waiting_for_approval"
        suffix = (
            " Planning and preview evidence were retained, but the current "
            "host does not expose a complete approval-bound execution composition. "
            "No chemistry engine was launched."
        )
        final_text = (final_text.rstrip() + suffix).strip()
    body = {
        "schema_version": "chemsmart.live-agent-session-result.v1",
        "session_id": session_id,
        "task_spec_sha256": task_spec_sha256,
        "terminal_state": terminal_state,
        "execution_requested": bool(execution_enabled),
        "execution_profile_status": execution_status,
        "final_text": final_text,
        "artifact_records": tuple(item.public_record() for item in observations),
        "conformance_records": conformance_records,
        "public_transcript": loop_result.public_transcript,
        "successful_tool_calls": loop_result.successful_tool_calls,
        "failed_tool_calls": loop_result.failed_tool_calls,
        "event_stream_head_sha256": loop_result.event_stream_head_sha256,
    }
    return LiveAgentSessionResultV1(
        **body, result_sha256=canonical_sha256(body)
    )


def _validated_workspace(value: str | Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise ContractError("agent workspace must be an absolute path")
    if not path.is_dir() or path.is_symlink():
        raise ContractError("agent workspace must be a regular directory")
    return path.resolve()


def _private_run_directory(workspace: Path, session_id: str) -> Path:
    private_root = workspace / _PRIVATE_ROOT_NAME
    if private_root.is_symlink():
        raise ContractError("private agent directory cannot be a symbolic link")
    private_root.mkdir(exist_ok=True, mode=0o700)
    if not private_root.is_dir():
        raise ContractError("private agent path is not a directory")
    private_root.chmod(0o700)
    root = private_root / "runs"
    if root.is_symlink():
        raise ContractError("private run root cannot be a symbolic link")
    root.mkdir(exist_ok=True, mode=0o700)
    if not root.is_dir():
        raise ContractError("private run root is not a directory")
    root.chmod(0o700)
    target = root / session_id
    target.mkdir(mode=0o700)
    target.chmod(0o700)
    return target


def _scan_xyz_artifacts(workspace: Path) -> tuple[_XyzObservation, ...]:
    observations: dict[str, _XyzObservation] = {}
    private_root = workspace / _PRIVATE_ROOT_NAME
    host_artifact_root = workspace / "artifacts"
    host_node_root = workspace / "nodes"
    for candidate in sorted(workspace.rglob("*.xyz")):
        if any(
            root in candidate.parents
            for root in (private_root, host_artifact_root, host_node_root)
        ):
            continue
        if candidate.is_symlink() or not candidate.is_file():
            continue
        resolved = candidate.resolve()
        try:
            resolved.relative_to(workspace)
        except ValueError as exc:
            raise ContractError("XYZ artifact escapes the approved workspace") from exc
        digest = file_sha256(resolved)
        atom_count, symbols = _inspect_xyz(resolved)
        artifact_id = f"geometry-{digest[:16]}"
        artifact = TrustedArtifactRefV1(
            artifact_id=artifact_id,
            kind="geometry_xyz",
            sha256=digest,
            size_bytes=resolved.stat().st_size,
            path=str(resolved),
            cli_value=str(resolved),
        )
        observations.setdefault(
            digest,
            _XyzObservation(
                artifact=artifact,
                atom_count=atom_count,
                symbols=symbols,
            ),
        )
    return tuple(
        sorted(observations.values(), key=lambda item: item.artifact.artifact_id)
    )


def _inspect_xyz(path: Path) -> tuple[int, tuple[str, ...]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
        count = int(lines[0].strip())
    except (IndexError, UnicodeDecodeError, ValueError) as exc:
        raise ContractError("workspace XYZ has a malformed atom count") from exc
    if count < 1 or len(lines) < count + 2:
        raise ContractError("workspace XYZ is truncated")
    symbols = []
    for row in lines[2 : count + 2]:
        fields = row.split()
        if len(fields) < 4 or not fields[0].isalpha():
            raise ContractError("workspace XYZ has a malformed atom row")
        try:
            coordinates = tuple(float(item) for item in fields[1:4])
        except ValueError as exc:
            raise ContractError("workspace XYZ coordinates are not numeric") from exc
        if not all(math.isfinite(item) for item in coordinates):
            raise ContractError("workspace XYZ coordinates must be finite")
        symbols.append(fields[0])
    return count, tuple(symbols)


def _task_spec_sha256(
    task: str, observations: Iterable[_XyzObservation]
) -> str:
    return canonical_sha256(
        {
            "schema_version": "chemsmart.live-scientific-task.v1",
            "task": task,
            "coordinate_artifacts": tuple(
                {
                    "artifact_id": item.artifact.artifact_id,
                    "sha256": item.artifact.sha256,
                    "atom_count": item.atom_count,
                    "symbols": item.symbols,
                }
                for item in observations
            ),
        }
    )


def _session_id(task_spec_sha256: str) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return f"live-{timestamp}-{task_spec_sha256[:10]}-{uuid.uuid4().hex[:8]}"


def _bootstrap_conformance(
    *,
    run_directory: Path,
    input_artifact: TrustedArtifactRefV1,
    registry_sha256: str,
    live_schema: Any,
) -> tuple[
    tuple[ProgramComponentConformanceReceiptV1, ...],
    tuple[dict[str, Any], ...],
]:
    bootstrap_directory = run_directory / "bootstrap"
    bootstrap_directory.mkdir(mode=0o700)
    project_path = bootstrap_directory / "pyscf-water.yaml"
    document = project_document(
        program="pyscf", sections=_locked_pyscf_sections()
    )
    rendered = render_project_yaml(document)
    _write_private_exact(project_path, rendered.rendered_yaml.encode("utf-8"))
    server_path = bootstrap_directory / "preview-server.yaml"
    _write_private_exact(
        server_path,
        _preview_server_profile().encode("utf-8"),
    )

    receipts: list[ProgramComponentConformanceReceiptV1] = []
    records: list[dict[str, Any]] = []
    pyscf_parts = []
    for engine in ("cpu", "gpu"):
        try:
            receipt = bootstrap_program_conformance(
                program="pyscf",
                engine=engine,
                jobtypes=("sp", "opt", "hess"),
                input_path=input_artifact.path,
                project_path=project_path,
                server_path=server_path,
                charge=0,
                multiplicity=1,
                registry_sha256=registry_sha256,
                live_schema=live_schema,
            )
            pyscf_parts.append(receipt)
            records.append(_conformance_record(receipt, engine=engine))
        except Exception as exc:  # Conformance remains observable and fail-closed.
            records.append(
                {
                    "record_kind": "program_conformance",
                    "program": "pyscf",
                    "engine": engine,
                    "status": "failed",
                    "error_class": type(exc).__name__,
                    "covered_jobtypes": (),
                }
            )
    if pyscf_parts:
        combined = _combine_program_conformance(pyscf_parts)
        receipts.append(combined)
    try:
        xtb = bootstrap_program_conformance(
            program="xtb",
            engine="cpu",
            jobtypes=("sp", "opt", "hess"),
            input_path=input_artifact.path,
            project_path=None,
            server_path=server_path,
            charge=0,
            multiplicity=1,
            registry_sha256=registry_sha256,
            live_schema=live_schema,
        )
        receipts.append(xtb)
        records.append(_conformance_record(xtb, engine="cpu"))
    except Exception as exc:
        records.append(
            {
                "record_kind": "program_conformance",
                "program": "xtb",
                "engine": "cpu",
                "status": "failed",
                "error_class": type(exc).__name__,
                "covered_jobtypes": (),
            }
        )
    return tuple(receipts), tuple(records)


def _combine_program_conformance(
    receipts: Iterable[ProgramComponentConformanceReceiptV1],
) -> ProgramComponentConformanceReceiptV1:
    rows = tuple(receipts)
    if not rows or len({item.program for item in rows}) != 1:
        raise ContractError("combined conformance requires one program")
    passed = tuple(
        item
        for item in rows
        if all(
            value == "passed"
            for value in (
                item.compiler_status,
                item.preview_status,
                item.preflight_status,
                item.verifier_status,
            )
        )
    )
    covered_engines = tuple(
        sorted({engine for item in passed for engine in item.covered_engines})
    )
    job_sets = [set(item.covered_jobtypes) for item in passed]
    covered_jobtypes = tuple(sorted(set.intersection(*job_sets))) if job_sets else ()
    status = "passed" if covered_engines and covered_jobtypes else "failed"
    def aggregate(field: str) -> str:
        return canonical_sha256(tuple(getattr(item, field) for item in rows))
    return build_program_component_conformance_receipt(
        program=rows[0].program,
        registry_sha256=rows[0].registry_sha256,
        live_cli_schema_sha256=rows[0].live_cli_schema_sha256,
        fixture_bundle_sha256=canonical_sha256(
            tuple(item.fixture_bundle_sha256 for item in rows)
        ),
        covered_jobtypes=covered_jobtypes,
        covered_engines=covered_engines,
        compiler_receipt_sha256=aggregate("compiler_receipt_sha256"),
        preview_receipt_sha256=aggregate("preview_receipt_sha256"),
        preflight_receipt_sha256=aggregate("preflight_receipt_sha256"),
        verifier_receipt_sha256=aggregate("verifier_receipt_sha256"),
        compiler_status=status,
        preview_status=status,
        preflight_status=status,
        verifier_status=status,
    )


def _locked_pyscf_sections() -> dict[str, dict[str, Any]]:
    common = {
        "basis": "def2-svp",
        "defgrid": "defgrid2",
        "density_fit": False,
        "functional": "b3lyp",
        "scf_maxiter": 100,
        "scf_tol": 1e-9,
    }
    return {
        "sp": dict(common),
        "opt": {**common, "opt_maxsteps": 100, "opt_solver": "geometric"},
        "hess": dict(common),
    }


def _preview_server_profile() -> str:
    """Scheduler-shaped profile for non-submitting run/sub conformance."""

    xtb = os.environ.get("CHEMSMART_XTB_EXECUTABLE") or shutil.which("xtb")
    xtb_folder = Path(xtb).expanduser().parent if xtb else Path(".")
    return (
        "SERVER:\n"
        "  SCHEDULER: PBS\n"
        "  QUEUE_NAME: preview\n"
        "  NUM_HOURS: 1\n"
        "  MEM_GB: 4\n"
        "  NUM_CORES: 4\n"
        "  NUM_GPUS: 0\n"
        "  NUM_THREADS: 4\n"
        "  SUBMIT_COMMAND: true\n"
        "  SCRATCH_DIR: null\n"
        "  PROJECT: preview\n"
        "  USE_HOSTS: false\n"
        "PYSCF:\n"
        f"  EXEFOLDER: {str(_PYSCF_INTERPRETER.parent)!r}\n"
        "  LOCAL_RUN: true\n"
        "  SCRATCH: false\n"
        "XTB:\n"
        f"  EXEFOLDER: {str(xtb_folder)!r}\n"
        "  LOCAL_RUN: true\n"
        "  SCRATCH: false\n"
    )


def _write_private_exact(path: Path, payload: bytes) -> None:
    if path.exists():
        if path.is_symlink() or path.read_bytes() != payload:
            raise ContractError(
                "private bootstrap artifact conflicts with existing bytes"
            )
        return
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        pending = memoryview(payload)
        while pending:
            written = os.write(descriptor, pending)
            pending = pending[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _observe_environments() -> tuple[
    tuple[EnvironmentTargetV1, ...],
    tuple[TrustedComputeEnvironmentReceiptV1, ...],
    tuple[dict[str, Any], ...],
]:
    targets = [
        EnvironmentTargetV1(
            program="pyscf",
            engine="cpu",
            target_kind="compute_receipt",
            locator="registered-pyscf-cpu-interpreter",
            required_dependencies=("h5py", "numpy", "pyscf"),
            required_dependency_versions=(("pyscf", "2.14.0"),),
        ),
        EnvironmentTargetV1(
            program="pyscf",
            engine="gpu",
            target_kind="compute_receipt",
            locator="registered-pyscf-gpu-interpreter",
            required_dependencies=(
                "cupy",
                "gpu4pyscf",
                "h5py",
                "numpy",
                "pyscf",
            ),
            required_dependency_versions=(
                ("gpu4pyscf", "1.8.0"),
                ("pyscf", "2.14.0"),
            ),
            required_gpu_facts=("device_available", "gpu4pyscf_distribution"),
        ),
        EnvironmentTargetV1(
            program="xtb",
            engine="cpu",
            target_kind="executable",
            locator="xtb",
        ),
    ]
    receipts = []
    records = []
    for engine in ("cpu", "gpu"):
        try:
            receipt = probe_python_compute_environment(
                _PYSCF_INTERPRETER, engine=engine
            )
            receipts.append(receipt)
            gpu = dict(receipt.gpu_evidence)
            records.append(
                {
                    "record_kind": "program_environment",
                    "program": "pyscf",
                    "engine": engine,
                    "status": (
                        "available"
                        if engine == "cpu"
                        else "available"
                        if gpu.get("device_available") is True
                        and gpu.get("gpu4pyscf_distribution") is True
                        else "missing"
                    ),
                    "dependency_versions": receipt.dependency_versions,
                    "gpu_available": bool(gpu.get("device_available", False)),
                    "receipt_sha256": receipt.evidence_sha256,
                }
            )
        except Exception as exc:
            records.append(
                {
                    "record_kind": "program_environment",
                    "program": "pyscf",
                    "engine": engine,
                    "status": "missing",
                    "error_class": type(exc).__name__,
                }
            )
    xtb_location = shutil.which("xtb")
    records.append(
        {
            "record_kind": "program_environment",
            "program": "xtb",
            "engine": "cpu",
            "status": "available" if xtb_location else "missing",
            "observation_method": "shutil.which",
        }
    )
    return tuple(targets), tuple(receipts), tuple(records)


def _conformance_record(
    receipt: ProgramComponentConformanceReceiptV1, *, engine: str
) -> dict[str, Any]:
    status = (
        "passed"
        if all(
            item == "passed"
            for item in (
                receipt.compiler_status,
                receipt.preview_status,
                receipt.preflight_status,
                receipt.verifier_status,
            )
        )
        else "failed"
    )
    return {
        "record_kind": "program_conformance",
        "program": receipt.program,
        "engine": engine,
        "status": status,
        "covered_jobtypes": receipt.covered_jobtypes,
        "receipt_sha256": receipt.receipt_sha256,
    }


def _record_sort_key(value: Mapping[str, Any]) -> tuple[str, str, str]:
    return (
        str(value.get("record_kind", "")),
        str(value.get("program", "")),
        str(value.get("engine", "")),
    )


def _public_context(
    *,
    task: str,
    task_spec_sha256: str,
    observations: tuple[_XyzObservation, ...],
    conformance_records: tuple[dict[str, Any], ...],
    registry_sha256: str,
    live_schema_sha256: str,
    execution_requested: bool,
    execution_available: bool,
    approved_project_records: tuple[dict[str, Any], ...] = (),
    approved_workflow_record: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    public_workflow = dict(approved_workflow_record or {})
    return {
        "schema_version": "chemsmart.live-agent-public-context.v1",
        "task": task,
        "task_spec_sha256": task_spec_sha256,
        "registry_sha256": registry_sha256,
        "live_cli_schema_sha256": live_schema_sha256,
        "artifacts": tuple(item.public_record() for item in observations),
        "approved_project_artifacts": approved_project_records,
        "approved_workflow": public_workflow,
        "approved_execution_contract": _approved_execution_context(
            public_workflow
        ),
        "conformance_observations": conformance_records,
        "execution_requested": bool(execution_requested),
        "execution_tool_available": bool(execution_available),
        "authority": (
            "Use artifact IDs and typed tool fields only. The host owns paths, "
            "CLI argv, project materialization, validation, approval, and execution."
        ),
    }


def _approved_execution_context(
    public_workflow: Mapping[str, Any],
) -> dict[str, Any]:
    """Describe execution mechanics without inventing scientific settings."""

    nodes = tuple(public_workflow.get("nodes", ()))
    if not nodes:
        return {}
    return {
        "programs": tuple(
            sorted({str(item.get("program", "")) for item in nodes})
        ),
        "node_order": tuple(public_workflow["node_order"]),
        "scientific_settings_authority": (
            "task_and_approved_project_artifacts"
        ),
        "producer_input_policy": "validated_approved_edge_only",
        "resources": {
            "cores": 4,
            "memory_gb": 4,
            "gpu_count": 0,
            "scratch": "none",
        },
    }


def _system_prompt(approved_workflow: Mapping[str, Any] | None) -> str:
    execution_available = bool(approved_workflow)
    execution_sentence = (
        "Execution is exposed only as execute_approved_program_node(node_id). "
        "Use the exact approved_workflow node IDs and listed order, and stop on "
        "any red result. When an approved producer edge exists, execute its "
        "producer first; after validation, consume the returned produced_handoffs "
        "artifact and scientific identity before compiling, previewing, "
        "preflighting, and executing that edge's consumer. "
        "Reuse the listed approved project artifact rather than promoting a new "
        "copy. Treat the execution receipt's deterministic validators as the live "
        "result gate; do not call inspect_calculation_artifact for newly executed "
        "nodes because this composition has no separate legacy settings/run IDs."
        if execution_available
        else (
            "Execution is not exposed. Finish with a useful planned or "
            "previewed state."
        )
    )
    return (
        "You are a professional computational-chemistry planning agent operating "
        "ChemSmart 3.1.4. Work plan-first through typed tools. Inspect program "
        "capability and environment, bind exact artifact identity, render and promote "
        "stage-specific project YAML, validate it, build a command DAG, compile safe "
        "commands, and preview every currently resolvable node. Keep every future "
        "producer input unresolved until its validated upstream artifact exists. "
        "Never author native "
        "Gaussian, ORCA, xTB, or PySCF input/script text. Never invent coordinates, "
        "paths, shell syntax, evidence, readiness, or terminal state. Explain method "
        "rationale, alternatives, uncertainty, and diagnostics in concise public "
        "English. Preserve every user-declared stage order as command-DAG control "
        "dependencies, including stages that independently reuse the initial "
        "geometry. Present a method as runnable only when the current project "
        "schema and loader support it; label other scientifically relevant methods "
        "as unsupported architecture alternatives rather than promising preview or "
        "execution. PySCF and xTB project stage keys are exactly sp, opt, and hess; "
        "workflow node IDs may separately express initial or optimized geometry. "
        "For each job, pass the exact receipt_sha256 returned by that job's "
        "inspect_program_capability call into environment inspection and project "
        "validation, then use the engine binding returned by environment inspection. "
        "Do not substitute conformance, joined-capability, or environment receipt "
        "digests for those typed fields. Keep project artifact IDs distinct from "
        "geometry artifact IDs. In workflow inputs, represent an initial artifact "
        "with empty producer_node_id and producer_output_id strings; represent a "
        "future optimized input with its producer IDs and no invented artifact ID. "
        "Omit absent optional settings instead of encoding them as the string none. "
        "When an approved project artifact is supplied, read and validate that exact "
        "artifact instead of rerendering an equivalent project. "
        "If critical evidence is missing, identify it and block honestly. "
        + execution_sentence
    )


def _network_budget() -> AdaptiveNetworkBudgetV1:
    body = {
        "schema_version": "chemsmart.adaptive-network-budget.v1",
        "allowed_provider": "deepseek",
        "endpoint_origin": DEEPSEEK_OFFICIAL_ENDPOINT,
        "purpose": "live-command-compiled-computational-chemistry-planning",
        "max_concurrency": 1,
        "max_input_tokens_per_request": DEEPSEEK_V4_FLASH_CONTEXT_TOKENS,
        "max_output_tokens_per_request": DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS,
        "task_wall_time_seconds": float(_SESSION_WALL_TIME_SECONDS),
        "quota_scope": "current_user_account",
        "top_up_allowed": False,
        "engine_calls": 0,
        "hpc_calls": 0,
    }
    return AdaptiveNetworkBudgetV1(
        **body, budget_sha256=canonical_sha256(body)
    )


def _hypothesis(
    *,
    session_id: str,
    messages: list[dict[str, str]],
    task_spec_sha256: str,
    artifact_sha256s: tuple[str, ...],
    tool_schema_sha256: str,
    network_budget: AdaptiveNetworkBudgetV1,
    execution_requested: bool,
) -> AdaptiveHypothesisV1:
    body = {
        "schema_version": "chemsmart.adaptive-hypothesis.v1",
        "hypothesis_id": session_id,
        "changed_factor": (
            "approval_bound_execution_profile"
            if execution_requested
            else "command_compiled_preview_profile"
        ),
        "comparator_id": "single-agent-natural-language-baseline",
        "expected_outcome": (
            "The model uses typed ChemSmart tools, preserves the exact task, "
            "and leaves every future producer artifact unresolved until host "
            "evidence exists."
        ),
        "deterministic_oracle_id": "live-project-command-preview-gates",
        "source_sha256s": tuple(sorted({task_spec_sha256, *artifact_sha256s})),
        "prompt_sha256": canonical_sha256(messages),
        "tool_schema_sha256": tool_schema_sha256,
        "configuration_sha256": canonical_sha256(
            {
                "network_budget": network_budget,
                "task_spec_sha256": task_spec_sha256,
                "execution_requested": bool(execution_requested),
            }
        ),
        "distinct_from_prior": (
            "Unique live session over an exact task and coordinate-artifact set."
        ),
    }
    return AdaptiveHypothesisV1(
        **body, hypothesis_sha256=canonical_sha256(body)
    )


def _task_envelope(
    *,
    session_id: str,
    messages: list[dict[str, str]],
    task_spec_sha256: str,
    workspace_records: tuple[dict[str, Any], ...],
    tool_schema_sha256: str,
    execution_enabled: bool,
) -> TaskEnvelopeV1:
    resource = ResourceBudgetV1(
        max_input_tokens_per_request=DEEPSEEK_V4_FLASH_CONTEXT_TOKENS,
        max_output_tokens_per_request=DEEPSEEK_V4_FLASH_MAX_OUTPUT_TOKENS,
        max_tool_calls=_MAX_TOOL_CALLS,
        wall_time_seconds=float(_SESSION_WALL_TIME_SECONDS),
        chemistry_engine_calls=4 if execution_enabled else 0,
        hpc_calls=0,
    )
    body = {
        "schema_version": "chemsmart.task-envelope.v1",
        "task_id": f"live-task-{task_spec_sha256[:16]}",
        "session_id": session_id,
        "turn_id": session_id + ".turn-1",
        "request_sha256": canonical_sha256(messages),
        "cwd_sha256": canonical_sha256(workspace_records),
        "phase": TaskPhase.ROUTE,
        "budget": resource,
        "tool_schema_sha256": tool_schema_sha256,
    }
    return TaskEnvelopeV1(
        **body, envelope_sha256=canonical_sha256(body)
    )


def _execution_composition_available() -> bool:
    parameters = inspect.signature(CommandCompiledToolHostV1).parameters
    required = {
        "approved_workspace",
        "workflow_execution_approval",
        "execution_resources",
        "execution_server",
        "execution_environment",
    }
    return required.issubset(parameters)


def _execution_composition_inputs(
    *,
    host_type: type[CommandCompiledToolHostV1],
    workspace: Path,
    run_directory: Path,
    approval_file: Path,
    task_spec_sha256: str,
) -> dict[str, Any]:
    """Load one exact workflow approval and fixed local resource profile."""

    if not _execution_composition_available():
        raise ContractError("approval-bound execution host is unavailable")
    expected_parameters = {
        "approved_workspace",
        "execution_resources",
        "workflow_execution_approval",
        "execution_server",
        "execution_environment",
    }
    if not expected_parameters.issubset(inspect.signature(host_type).parameters):
        raise ContractError("execution host composition API is incomplete")
    if not approval_file.is_file() or approval_file.is_symlink():
        raise ContractError("approval file must be a current regular file")
    try:
        payload = json.loads(approval_file.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError("approval file is not valid JSON") from exc
    if not isinstance(payload, Mapping):
        raise ContractError("approval file root must be an object")
    raw_approval = payload.get("workflow_approval", payload)
    if not isinstance(raw_approval, Mapping):
        raise ContractError("workflow approval must be an object")
    resources = _locked_execution_resources()
    approval = _parse_workflow_approval(raw_approval)
    if Path(approval.workspace).resolve() != workspace:
        raise ContractError("workflow approval targets another workspace")
    if approval.task_spec_sha256 != task_spec_sha256:
        raise ContractError("workflow approval targets another task spec")
    if approval.resource_sha256 != resources.resource_sha256:
        raise ContractError("workflow approval resources differ from locked resources")
    server_profile = _write_execution_server_profile(run_directory)
    path_value = os.environ.get("PATH", "")
    xtb_executable = (
        os.environ.get("CHEMSMART_XTB_EXECUTABLE") or shutil.which("xtb")
    )
    executable_directory = (
        str(Path(xtb_executable).expanduser().parent)
        if xtb_executable
        else ""
    )
    environment = {
        "PATH": (
            path_value
            if not executable_directory
            else executable_directory
            if not path_value
            else executable_directory + os.pathsep + path_value
        ),
        "PYTHONNOUSERSITE": "1",
    }
    approved_projects = _approved_project_artifacts(workspace, approval)
    return {
        "approved_workspace": workspace,
        "execution_resources": resources,
        "workflow_execution_approval": approval,
        "execution_server": str(server_profile),
        "execution_environment": environment,
        "approved_project_artifacts": approved_projects,
    }


def _approved_project_artifacts(
    workspace: Path, approval: WorkflowExecutionApprovalV1
) -> tuple[TrustedArtifactRefV1, ...]:
    """Resolve workflow-bound project bytes without model rematerialization."""

    required = {
        item.project_artifact_sha256 for item in approval.node_bindings
    }
    matches: dict[str, Path] = {}
    project_root = workspace / "projects"
    if project_root.is_dir() and not project_root.is_symlink():
        for candidate in sorted(project_root.glob("*.yaml")):
            if candidate.is_symlink() or not candidate.is_file():
                continue
            digest = file_sha256(candidate)
            if digest in required:
                if digest in matches:
                    raise ContractError(
                        "approved project bytes have multiple workspace identities"
                    )
                matches[digest] = candidate.resolve()
    missing = sorted(required.difference(matches))
    if missing:
        raise ContractError("workflow approval project artifact is unavailable")
    return tuple(
        TrustedArtifactRefV1(
            artifact_id=f"project-approved-{digest[:16]}",
            kind="project_yaml",
            sha256=digest,
            size_bytes=matches[digest].stat().st_size,
            path=str(matches[digest]),
            cli_value=str(matches[digest]),
        )
        for digest in sorted(required)
    )


def _public_workflow_approval(
    approval: WorkflowExecutionApprovalV1,
) -> dict[str, Any]:
    """Expose only the typed execution sequence needed by the model."""

    return {
        "schema_version": "chemsmart.public-approved-workflow.v1",
        "workflow_id": approval.workflow_id,
        "node_order": tuple(item.node_id for item in approval.node_bindings),
        "nodes": tuple(
            {
                "node_id": item.node_id,
                "program": item.program,
                "engine": item.engine,
                "jobtype": item.jobtype,
                "input_mode": item.input_mode,
            }
            for item in approval.node_bindings
        ),
        "producer_edges": tuple(
            {
                "producer_node_id": item.producer_node_id,
                "consumer_node_id": item.consumer_node_id,
                "artifact_kind": item.artifact_kind,
                "selection_rule": item.selection_rule,
            }
            for item in approval.producer_edges
        ),
        "status": approval.status,
    }


def _locked_execution_resources() -> ExecutionResourceSpecV1:
    return build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=4,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=_SESSION_WALL_TIME_SECONDS,
    )


def _parse_workflow_approval(
    value: Mapping[str, Any],
) -> WorkflowExecutionApprovalV1:
    raw = dict(value)
    try:
        raw["node_bindings"] = tuple(
            ApprovedNodeBindingV1(**dict(item))
            for item in raw.get("node_bindings", ())
        )
        raw["producer_edges"] = tuple(
            ProducerEdgeRuleV1(**dict(item))
            for item in raw.get("producer_edges", ())
        )
        return WorkflowExecutionApprovalV1(**raw)
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError("workflow approval does not match the v1 schema") from exc


def _write_execution_server_profile(run_directory: Path) -> Path:
    """Write the deterministic local CPU profile used by approved nodes."""

    xtb = shutil.which("xtb")
    profile = run_directory / "execution-server.yaml"
    text = (
        "SERVER:\n"
        "  SCHEDULER: null\n"
        "  NUM_HOURS: 2\n"
        "  MEM_GB: 4\n"
        "  NUM_CORES: 4\n"
        "  NUM_GPUS: 0\n"
        "  NUM_THREADS: 4\n"
        "  SCRATCH_DIR: null\n"
        "PYSCF:\n"
        f"  EXEFOLDER: {str(_PYSCF_INTERPRETER.parent)!r}\n"
        "  LOCAL_RUN: true\n"
        "  SCRATCH: false\n"
    )
    if xtb:
        text += (
            "XTB:\n"
            f"  EXEFOLDER: {str(Path(xtb).resolve().parent)!r}\n"
            "  LOCAL_RUN: true\n"
            "  SCRATCH: false\n"
        )
    _write_private_exact(profile, text.encode("utf-8"))
    return profile


def _local_result(
    *,
    session_id: str,
    task_spec_sha256: str,
    terminal_state: str,
    execution_requested: bool,
    execution_profile_status: str,
    final_text: str,
) -> LiveAgentSessionResultV1:
    body = {
        "schema_version": "chemsmart.live-agent-session-result.v1",
        "session_id": session_id,
        "task_spec_sha256": task_spec_sha256,
        "terminal_state": terminal_state,
        "execution_requested": bool(execution_requested),
        "execution_profile_status": execution_profile_status,
        "final_text": final_text,
        "artifact_records": (),
        "conformance_records": (),
        "public_transcript": (),
        "successful_tool_calls": 0,
        "failed_tool_calls": 0,
        "event_stream_head_sha256": "",
    }
    return LiveAgentSessionResultV1(
        **body, result_sha256=canonical_sha256(body)
    )


__all__ = ["LiveAgentSessionResultV1", "run_live_agent_session"]
