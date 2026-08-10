"""Typed-to-argv command compilation grounded in the live Click tree."""

from __future__ import annotations

import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import click

from chemsmart.agent._contracts import (
    AuxiliaryArtifactBindingV1,
    ContractError,
    TrustedArtifactRefV1,
    canonical_sha256,
    file_sha256,
    require_identifier,
    require_auxiliary_artifact_bindings,
    require_sha256,
)
from chemsmart.agent.capabilities import (
    CapabilityQueryReceiptV1,
    CapabilityQueryStatus,
    ResolvedEngineBindingV1,
)
from chemsmart.agent.cli_schema import LiveClickSchemaV1, build_live_click_schema
from chemsmart.agent.projects import ProjectValidationReceiptV1


@dataclass(frozen=True)
class ScientificIdentityBindingV1:
    """Host-owned molecular and electronic identity for one geometry frame."""

    schema_version: str
    task_spec_sha256: str
    geometry_artifact_id: str
    geometry_artifact_sha256: str
    charge: int
    multiplicity: int
    binding_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.scientific-identity-binding.v1":
            raise ContractError("unsupported scientific identity binding schema")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        if not str(self.geometry_artifact_id).strip():
            raise ContractError("geometry_artifact_id must not be empty")
        require_sha256(
            self.geometry_artifact_sha256, "geometry_artifact_sha256"
        )
        if self.multiplicity < 1:
            raise ContractError("multiplicity must be positive")
        body = {
            "schema_version": self.schema_version,
            "task_spec_sha256": self.task_spec_sha256,
            "geometry_artifact_id": self.geometry_artifact_id,
            "geometry_artifact_sha256": self.geometry_artifact_sha256,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
        }
        if self.binding_sha256 != canonical_sha256(body):
            raise ContractError("scientific identity binding digest mismatch")


def build_scientific_identity_binding(
    *,
    task_spec_sha256: str,
    geometry_artifact: TrustedArtifactRefV1,
    charge: int,
    multiplicity: int,
) -> ScientificIdentityBindingV1:
    """Bind explicit task state to the exact current geometry bytes."""

    _require_current_artifact(geometry_artifact, "geometry")
    if geometry_artifact.kind not in {"geometry_xyz", "xyz"}:
        raise ContractError(
            "scientific identity requires an exact geometry artifact"
        )
    body = {
        "schema_version": "chemsmart.scientific-identity-binding.v1",
        "task_spec_sha256": require_sha256(
            task_spec_sha256, "task_spec_sha256"
        ),
        "geometry_artifact_id": geometry_artifact.artifact_id,
        "geometry_artifact_sha256": geometry_artifact.sha256,
        "charge": int(charge),
        "multiplicity": int(multiplicity),
    }
    return ScientificIdentityBindingV1(
        **body, binding_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class CommandProposalV1:
    """Model-visible intent with references, never paths, flags, or shell."""

    node_id: str
    execution_target: str
    program: str
    jobtype: str
    project_artifact_id: str
    input_artifact_id: str
    scientific_identity_sha256: str
    charge: int
    multiplicity: int

    def __post_init__(self) -> None:
        require_identifier(self.node_id, "node_id")
        if self.execution_target not in {"run", "sub"}:
            raise ContractError("execution_target must be run or sub")
        require_identifier(self.program, "program")
        require_identifier(self.jobtype, "jobtype")
        require_sha256(
            self.scientific_identity_sha256, "scientific_identity_sha256"
        )
        if self.multiplicity < 1:
            raise ContractError("multiplicity must be positive")


@dataclass(frozen=True)
class ScopedCommandOptionV1:
    scope_path: tuple[str, ...]
    parameter_name: str
    flag: str
    values: tuple[str, ...] = ()


@dataclass(frozen=True)
class CanonicalCommandInvocationV1:
    schema_version: str
    node_id: str
    command_path: tuple[str, ...]
    scoped_options: tuple[ScopedCommandOptionV1, ...]
    argv: tuple[str, ...]
    display_command: str
    live_cli_schema_sha256: str
    joined_capability_sha256: str
    project_receipt_sha256: str
    project_sha256: str
    input_sha256: str
    scientific_identity_sha256: str
    program_engine_binding_sha256: str
    repair_parent_sha256: str
    counterexample_sha256: str
    repair_attempt: int
    status: str
    rule_ids: tuple[str, ...]
    invocation_sha256: str
    auxiliary_input_bindings: tuple[AuxiliaryArtifactBindingV1, ...] = ()

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.canonical-command-invocation.v1":
            raise ContractError("unsupported canonical invocation schema")
        if self.status != "compiled":
            raise ContractError("compilation does not establish preview")
        require_auxiliary_artifact_bindings(self.auxiliary_input_bindings)
        body = {
            "schema_version": self.schema_version,
            "node_id": self.node_id,
            "command_path": self.command_path,
            "scoped_options": self.scoped_options,
            "argv": self.argv,
            "display_command": self.display_command,
            "live_cli_schema_sha256": self.live_cli_schema_sha256,
            "joined_capability_sha256": self.joined_capability_sha256,
            "project_receipt_sha256": self.project_receipt_sha256,
            "project_sha256": self.project_sha256,
            "input_sha256": self.input_sha256,
            "scientific_identity_sha256": self.scientific_identity_sha256,
            "program_engine_binding_sha256": (
                self.program_engine_binding_sha256
            ),
            "repair_parent_sha256": self.repair_parent_sha256,
            "counterexample_sha256": self.counterexample_sha256,
            "repair_attempt": self.repair_attempt,
            "status": self.status,
            "rule_ids": self.rule_ids,
        }
        if self.auxiliary_input_bindings:
            body["auxiliary_input_bindings"] = self.auxiliary_input_bindings
        expected = canonical_sha256(body)
        if self.invocation_sha256 != expected:
            raise ContractError("canonical invocation digest mismatch")


@dataclass(frozen=True)
class CommandInspectionReceiptV1:
    schema_version: str
    invocation_sha256: str
    observed_command_path: tuple[str, ...]
    observed_option_names: tuple[str, ...]
    live_cli_schema_sha256: str
    parser_observation_sha256: str
    status: str
    rule_ids: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.command-inspection-receipt.v1":
            raise ContractError("unsupported command inspection schema")
        if self.status not in {"valid", "invalid"}:
            raise ContractError("invalid command inspection state")
        if self.status == "valid" and not self.parser_observation_sha256:
            raise ContractError("valid inspection requires a Click observation")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "invocation_sha256": self.invocation_sha256,
                "observed_command_path": self.observed_command_path,
                "observed_option_names": self.observed_option_names,
                "live_cli_schema_sha256": self.live_cli_schema_sha256,
                "parser_observation_sha256": (
                    self.parser_observation_sha256
                ),
                "status": self.status,
                "rule_ids": self.rule_ids,
            }
        )
        if self.receipt_sha256 != expected:
            raise ContractError("command inspection receipt digest mismatch")


@dataclass(frozen=True)
class CommandCounterexampleV1:
    schema_version: str
    counterexample_id: str
    invocation_sha256: str
    task_spec_sha256: str
    preflight_receipt_sha256: str
    rule_id: str
    failed_field: str
    expected: Any
    observed: Any
    counterexample_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.command-counterexample.v1":
            raise ContractError("unsupported command counterexample schema")
        require_identifier(self.counterexample_id, "counterexample_id")
        require_identifier(self.rule_id, "rule_id")
        require_identifier(self.failed_field, "failed_field")
        body = {
            "schema_version": self.schema_version,
            "counterexample_id": self.counterexample_id,
            "invocation_sha256": self.invocation_sha256,
            "task_spec_sha256": self.task_spec_sha256,
            "preflight_receipt_sha256": self.preflight_receipt_sha256,
            "rule_id": self.rule_id,
            "failed_field": self.failed_field,
            "expected": self.expected,
            "observed": self.observed,
        }
        if self.counterexample_sha256 != canonical_sha256(body):
            raise ContractError("command counterexample digest mismatch")


def compile_command(
    proposal: CommandProposalV1,
    *,
    capability: CapabilityQueryReceiptV1,
    binding: ResolvedEngineBindingV1,
    project: TrustedArtifactRefV1 | None,
    project_validation: ProjectValidationReceiptV1 | None,
    input_artifact: TrustedArtifactRefV1,
    scientific_identity: ScientificIdentityBindingV1,
    job_artifact_options: Mapping[str, TrustedArtifactRefV1] | None = None,
    live_schema: LiveClickSchemaV1 | None = None,
    server: str = "",
    repair_parent_sha256: str = "",
    counterexample_sha256: str = "",
    repair_attempt: int = 0,
) -> CanonicalCommandInvocationV1:
    """Compile safe-preview argv; never invoke Click or a chemistry engine."""

    live_schema = live_schema or build_live_click_schema()
    if repair_attempt not in {0, 1, 2}:
        raise ContractError("command repair attempt exceeds bounded limit")
    if repair_attempt and not (repair_parent_sha256 and counterexample_sha256):
        raise ContractError("command repair requires parent and counterexample")
    if not repair_attempt and (repair_parent_sha256 or counterexample_sha256):
        raise ContractError("initial compilation cannot carry repair ancestry")
    _validate_compiler_bindings(
        proposal,
        capability=capability,
        binding=binding,
        project=project,
        project_validation=project_validation,
        input_artifact=input_artifact,
        scientific_identity=scientific_identity,
        live_schema=live_schema,
    )
    command_path = (
        proposal.execution_target,
        proposal.program,
        proposal.jobtype,
    )
    execution_scope = (proposal.execution_target,)
    program_scope = execution_scope + (proposal.program,)
    job_scope = program_scope + (proposal.jobtype,)
    options: list[ScopedCommandOptionV1] = []

    if server:
        options.append(
            _scoped_option(live_schema, execution_scope, "server", server)
        )
    options.append(
        _scoped_option(live_schema, execution_scope, "fake", True)
    )
    scratch = live_schema.command(execution_scope).option("scratch")
    if scratch is not None:
        options.append(
            _scoped_option(live_schema, execution_scope, "scratch", False)
        )
    if proposal.execution_target == "sub":
        options.append(
            _scoped_option(live_schema, execution_scope, "test", True)
        )

    if project is not None:
        options.append(
            _scoped_option(
                live_schema, program_scope, "project", project.cli_value
            )
        )
    options.append(
        _scoped_option(
                live_schema,
                program_scope,
                "filename",
                input_artifact.cli_value,
            )
    )
    project_settings = (
        dict(project_validation.settings)
        if project_validation is not None
        else {}
    )
    project_owned = set(
        capability.capability.project_owned_parameters
        if capability.capability is not None
        else ()
    )
    if _emit_effective_override(
        "charge", proposal.charge, project_settings, project_owned
    ):
        options.append(
            _scoped_option(
                live_schema, program_scope, "charge", proposal.charge
            )
        )
    if _emit_effective_override(
        "multiplicity",
        proposal.multiplicity,
        project_settings,
        project_owned,
    ):
        options.append(
            _scoped_option(
                live_schema,
                program_scope,
                "multiplicity",
                proposal.multiplicity,
            )
        )
    gpu_option = live_schema.command(program_scope).option("gpu")
    if gpu_option is not None and binding.engine:
        options.append(
            _scoped_option(
                live_schema,
                program_scope,
                "gpu",
                binding.engine == "gpu",
            )
        )

    # Some canonical ChemSmart jobs consume more than one registered file.
    # Keep the main molecular geometry on the long-standing program-level
    # ``filename`` option, and bind every additional artifact to an option on
    # the exact live job command (for example ORCA NEB's
    # ``ending_xyzfile``).  The model supplies only semantic artifact IDs in
    # its workflow; the host resolves paths here.
    auxiliary_input_bindings = []
    for parameter_name, artifact in sorted(
        (job_artifact_options or {}).items()
    ):
        if parameter_name == "filename":
            raise ContractError(
                "filename is the primary program input, not a job artifact option"
            )
        _require_current_artifact(artifact, parameter_name)
        options.append(
            _scoped_option(
                live_schema,
                job_scope,
                parameter_name,
                artifact.cli_value,
            )
        )
        auxiliary_input_bindings.append(
            AuxiliaryArtifactBindingV1(
                parameter_name=parameter_name,
                artifact_id=artifact.artifact_id,
                artifact_sha256=artifact.sha256,
            )
        )

    ordered_options = tuple(options)
    argv = ["chemsmart", proposal.execution_target]
    argv.extend(_render_scope_options(ordered_options, execution_scope))
    argv.append(proposal.program)
    argv.extend(_render_scope_options(ordered_options, program_scope))
    argv.append(proposal.jobtype)
    argv.extend(_render_scope_options(ordered_options, job_scope))
    argv_tuple = tuple(argv)
    body = {
        "schema_version": "chemsmart.canonical-command-invocation.v1",
        "node_id": proposal.node_id,
        "command_path": command_path,
        "scoped_options": ordered_options,
        "argv": argv_tuple,
        "display_command": shlex.join(argv_tuple),
        "live_cli_schema_sha256": live_schema.schema_sha256,
        "joined_capability_sha256": capability.joined_capability_sha256,
        "project_receipt_sha256": (
            project_validation.receipt_sha256
            if project_validation is not None
            else ""
        ),
        "project_sha256": project.sha256 if project is not None else "",
        "input_sha256": input_artifact.sha256,
        "scientific_identity_sha256": scientific_identity.binding_sha256,
        "program_engine_binding_sha256": binding.binding_sha256,
        "repair_parent_sha256": repair_parent_sha256,
        "counterexample_sha256": counterexample_sha256,
        "repair_attempt": repair_attempt,
        "status": "compiled",
        "rule_ids": ("command.compiler.live_click", "command.preview.flags"),
    }
    if auxiliary_input_bindings:
        body["auxiliary_input_bindings"] = tuple(auxiliary_input_bindings)
    return CanonicalCommandInvocationV1(
        **body, invocation_sha256=canonical_sha256(body)
    )


def inspect_command(
    invocation: CanonicalCommandInvocationV1,
    *,
    live_schema: LiveClickSchemaV1 | None = None,
    root: click.Command | None = None,
) -> CommandInspectionReceiptV1:
    """Independently compare argv/scopes with the currently observed schema."""

    live_schema = live_schema or build_live_click_schema()
    valid = True
    rules = []
    if live_schema.schema_sha256 != invocation.live_cli_schema_sha256:
        valid = False
        rules.append("command.inspect.schema_drift")
    if not live_schema.has_path(invocation.command_path):
        valid = False
        rules.append("command.inspect.path_missing")
    observed_names = []
    for item in invocation.scoped_options:
        command = live_schema.command(item.scope_path)
        option = command.option(item.parameter_name) if command else None
        if option is None or item.flag not in (
            option.flags + option.secondary_flags
        ):
            valid = False
            rules.append("command.inspect.option_scope_mismatch")
        observed_names.append(item.parameter_name)
    reconstructed = ["chemsmart", invocation.command_path[0]]
    reconstructed.extend(
        _render_scope_options(
            invocation.scoped_options, (invocation.command_path[0],)
        )
    )
    reconstructed.append(invocation.command_path[1])
    reconstructed.extend(
        _render_scope_options(
            invocation.scoped_options, invocation.command_path[:2]
        )
    )
    reconstructed.append(invocation.command_path[2])
    reconstructed.extend(
        _render_scope_options(
            invocation.scoped_options, invocation.command_path
        )
    )
    if tuple(reconstructed) != invocation.argv:
        valid = False
        rules.append("command.inspect.argv_roundtrip_mismatch")
    try:
        observed_path, observed_params = _observe_click_parse(
            tuple(invocation.argv[1:]), root=root
        )
        parser_observation_sha256 = canonical_sha256(
            {"path": observed_path, "params": observed_params}
        )
        if observed_path != invocation.command_path:
            valid = False
            rules.append("command.inspect.click_path_mismatch")
        for item in invocation.scoped_options:
            key = "/".join(item.scope_path) + ":" + item.parameter_name
            if key not in observed_params:
                valid = False
                rules.append("command.inspect.click_option_unobserved")
                continue
            observed_value = observed_params[key]
            if item.values:
                expected_values = item.values
                if isinstance(observed_value, (tuple, list)):
                    observed_values = tuple(str(value) for value in observed_value)
                else:
                    observed_values = (str(observed_value),)
                if observed_values != expected_values:
                    valid = False
                    rules.append("command.inspect.click_value_mismatch")
            else:
                expected_flag_value = not item.flag.startswith("--no-")
                if bool(observed_value) is not expected_flag_value:
                    valid = False
                    rules.append("command.inspect.click_flag_mismatch")
    except (click.ClickException, click.UsageError, ContractError, TypeError):
        parser_observation_sha256 = ""
        valid = False
        rules.append("command.inspect.click_parse_failed")
    if valid:
        rules.append("command.inspect.roundtrip_valid")
    body = {
        "schema_version": "chemsmart.command-inspection-receipt.v1",
        "invocation_sha256": invocation.invocation_sha256,
        "observed_command_path": invocation.command_path,
        "observed_option_names": tuple(observed_names),
        "live_cli_schema_sha256": live_schema.schema_sha256,
        "parser_observation_sha256": parser_observation_sha256,
        "status": "valid" if valid else "invalid",
        "rule_ids": tuple(sorted(set(rules))),
    }
    return CommandInspectionReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _validate_compiler_bindings(
    proposal: CommandProposalV1,
    *,
    capability: CapabilityQueryReceiptV1,
    binding: ResolvedEngineBindingV1,
    project: TrustedArtifactRefV1 | None,
    project_validation: ProjectValidationReceiptV1 | None,
    input_artifact: TrustedArtifactRefV1,
    scientific_identity: ScientificIdentityBindingV1,
    live_schema: LiveClickSchemaV1,
) -> None:
    if capability.status not in {
        CapabilityQueryStatus.SUPPORTED,
        CapabilityQueryStatus.PREVIEW_ONLY,
    }:
        raise ContractError("command capability is not supported")
    if (proposal.program, proposal.jobtype) != (
        capability.query.program,
        capability.query.jobtype,
    ):
        raise ContractError("proposal differs from capability query")
    if capability.live_cli_schema_sha256 != live_schema.schema_sha256:
        raise ContractError("capability query uses a stale Click schema")
    if binding.program != proposal.program or binding.state == "blocked":
        raise ContractError("program-engine binding is not usable for preview")
    if binding.capability_receipt_sha256 != capability.receipt_sha256:
        raise ContractError("program-engine binding uses another capability")
    if binding.engine != capability.query.engine:
        raise ContractError("program-engine binding differs from capability")
    if len(binding.program_binding_sha256) != 64:
        raise ContractError("program-engine binding is incomplete")
    declared = capability.capability
    requires_project = bool(
        declared is not None and declared.requires_project_configuration
    )
    supports_project = bool(
        declared is not None and declared.supports_project_configuration
    )
    if requires_project and (project is None or project_validation is None):
        # These are two different mistakes with two different repairs, and one
        # message for both left a live run retrying the same call three times.
        # Promoting a project and validating one are separate steps, so say
        # which is missing and for which program.
        missing = (
            "no project artifact is bound to this node; promote one and "
            "validate it first"
            if project is None
            else (
                f"project {project.artifact_id!r} is bound but has no "
                "validation receipt; call validate_project_yaml on it first"
            )
        )
        raise ContractError(
            f"{proposal.program} requires a validated project: {missing}"
        )
    if not supports_project and project is not None:
        raise ContractError("program does not support project configuration")
    if project is None:
        if proposal.project_artifact_id:
            raise ContractError("proposal names an unbound project artifact")
    elif project.artifact_id != proposal.project_artifact_id:
        raise ContractError("proposal project reference is not bound")
    if input_artifact.artifact_id != proposal.input_artifact_id:
        raise ContractError("proposal input reference is not bound")
    _require_current_artifact(input_artifact, "input")
    if proposal.scientific_identity_sha256 != scientific_identity.binding_sha256:
        raise ContractError("proposal uses another scientific identity binding")
    if scientific_identity.geometry_artifact_id != input_artifact.artifact_id:
        raise ContractError("scientific identity uses another geometry artifact")
    if scientific_identity.geometry_artifact_sha256 != input_artifact.sha256:
        raise ContractError("scientific identity geometry digest mismatch")
    if (proposal.charge, proposal.multiplicity) != (
        scientific_identity.charge,
        scientific_identity.multiplicity,
    ):
        raise ContractError("proposal electronic state differs from task binding")
    if project is not None:
        _require_current_artifact(project, "project")
    if project_validation is not None:
        if project is None:
            raise ContractError("project validation has no artifact binding")
        if project_validation.status != "valid":
            raise ContractError("project loader validation is not green")
        if project_validation.project_sha256 != project.sha256:
            raise ContractError("project validation is stale")
        if (
            project_validation.capability_receipt_sha256
            != capability.receipt_sha256
        ):
            raise ContractError("project validation uses another capability")
    if not live_schema.has_path(
        (proposal.execution_target, proposal.program, proposal.jobtype)
    ):
        raise ContractError("command path is absent from the live Click tree")


def _scoped_option(
    schema: LiveClickSchemaV1,
    scope: tuple[str, ...],
    parameter_name: str,
    value: Any,
) -> ScopedCommandOptionV1:
    command = schema.command(scope)
    option = command.option(parameter_name) if command is not None else None
    if option is None:
        raise ContractError(
            f"live Click scope {'/'.join(scope)} has no {parameter_name} option"
        )
    if option.is_flag and isinstance(value, bool):
        if value:
            flag = option.primary_flag
        else:
            if not option.secondary_flags:
                raise ContractError(
                    f"{parameter_name} has no explicit negative flag"
                )
            flag = sorted(option.secondary_flags)[0]
        values: tuple[str, ...] = ()
    else:
        flag = option.primary_flag
        values = (str(value),)
    return ScopedCommandOptionV1(
        scope_path=scope,
        parameter_name=parameter_name,
        flag=flag,
        values=values,
    )


def _render_scope_options(
    options: tuple[ScopedCommandOptionV1, ...], scope: tuple[str, ...]
) -> list[str]:
    rendered = []
    for item in options:
        if item.scope_path == scope:
            rendered.append(item.flag)
            rendered.extend(item.values)
    return rendered


def _emit_effective_override(
    field: str,
    requested: int,
    project_settings: dict[str, Any],
    project_owned: set[str],
) -> bool:
    if field not in project_owned or field not in project_settings:
        return True
    try:
        return int(project_settings[field]) != int(requested)
    except (TypeError, ValueError):
        return True


def _require_current_artifact(
    binding: TrustedArtifactRefV1, label: str
) -> None:
    path = Path(binding.path)
    if not path.is_file() or path.stat().st_size != binding.size_bytes:
        raise ContractError(f"{label} artifact size differs from its binding")
    if file_sha256(path) != binding.sha256:
        raise ContractError(f"{label} artifact digest differs from its binding")


def _observe_click_parse(
    argv: tuple[str, ...], *, root: click.Command | None
) -> tuple[tuple[str, ...], dict[str, Any]]:
    """Parse every Click scope without invoking any command callback."""

    if root is None:
        from chemsmart.cli.main import entry_point

        root = entry_point
    current = root
    parent: click.Context | None = None
    remaining = list(argv)
    path: list[str] = []
    observed: dict[str, Any] = {}
    while True:
        info_name = current.name or (path[-1] if path else "chemsmart")
        context = current.make_context(
            info_name,
            remaining,
            parent=parent,
            resilient_parsing=False,
        )
        scope = "/".join(path)
        for name, value in sorted(context.params.items()):
            observed[f"{scope}:{name}"] = value
        if not isinstance(current, click.Group):
            remaining = list(context.args)
            if remaining:
                raise ContractError("Click parser left unconsumed arguments")
            break
        unresolved = list(getattr(context, "protected_args", ())) + list(
            context.args
        )
        # Some scientifically complete ChemSmart leaves are Click groups so
        # they can optionally expose a nested variant (for example ORCA
        # ``opt qmmm``).  Direct ``orca opt`` is nevertheless a valid terminal
        # command because that group declares ``invoke_without_command``.
        # Treat that Click-native shape as a leaf instead of asking
        # ``resolve_command`` to index an empty argument list.
        if not unresolved:
            if current.invoke_without_command:
                break
            raise ContractError("Click parser expected another command segment")
        command_name, command, remaining = current.resolve_command(
            context, unresolved
        )
        if command is None or command_name is None:
            raise ContractError("Click parser could not resolve command path")
        path.append(str(command_name))
        parent = context
        current = command
    return tuple(path), observed


__all__ = [
    "CanonicalCommandInvocationV1",
    "CommandInspectionReceiptV1",
    "CommandProposalV1",
    "CommandCounterexampleV1",
    "ScientificIdentityBindingV1",
    "ScopedCommandOptionV1",
    "build_scientific_identity_binding",
    "compile_command",
    "inspect_command",
]
