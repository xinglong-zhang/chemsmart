from __future__ import annotations

import importlib
import inspect
import logging
import os
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any, Literal, get_args, get_origin, get_type_hints

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    TypeAdapter,
    create_model,
)
from pydantic.errors import (
    PydanticInvalidForJsonSchema,
    PydanticSchemaGenerationError,
)

from chemsmart.agent.permissions import RuntimePermissionMode
from chemsmart.agent.tool_protocol import (
    RuntimeToolMetadata,
    is_allowed_in_mode,
)

logger = logging.getLogger(__name__)

# Module-type tool groups. Each registered tool belongs to exactly one group so
# the tool surface can be enabled incrementally — e.g. training/serving a local
# model on "synthesis" only, then adding "project_yaml", then "execution".
# Select groups via ToolRegistry.default(groups=[...]) or the
# CHEMSMART_AGENT_TOOL_GROUPS env var (comma-separated group names).
TOOL_GROUPS: dict[str, frozenset[str]] = {
    "synthesis": frozenset(
        {
            "synthesize_command",
            "repair_command",
        }
    ),
    "project_yaml": frozenset(
        {
            "extract_project_protocol",
            "render_project_yaml",
            "validate_project_yaml",
            "critic_project_yaml",
            "write_project_yaml",
            "read_project_yaml",
            "update_project_yaml",
            "search_basis_sets",
        }
    ),
    "harness_jobs": frozenset(
        {
            "build_molecule",
            "recommend_method",
            "build_gaussian_settings",
            "build_orca_settings",
            "build_job",
            "dry_run_input",
            "validate_runtime",
            "extract_optimized_geometry",
        }
    ),
    "execution": frozenset(
        {
            "execute_chemsmart_command",
            "run_local",
            "submit_hpc",
        }
    ),
    "wizard": frozenset(
        {
            "wizard_probe",
            "wizard_refresh",
            "wizard_verify",
            "wizard_write",
        }
    ),
    "diagnostics": frozenset(
        {
            "read",
            "ssh_probe",
            "scheduler_query",
            "log_tail",
        }
    ),
}

_TOOL_GROUP_BY_NAME: dict[str, str] = {
    tool_name: group_name
    for group_name, tool_names in TOOL_GROUPS.items()
    for tool_name in tool_names
}

_TOOL_GROUPS_ENV_VAR = "CHEMSMART_AGENT_TOOL_GROUPS"


def tool_group(tool_name: str) -> str | None:
    """Return the module group a registered tool belongs to."""

    return _TOOL_GROUP_BY_NAME.get(tool_name)


def resolve_tool_groups(
    groups: Iterable[str] | None = None,
) -> frozenset[str] | None:
    """Resolve the enabled tool-name set from groups or the env override.

    Returns None when no restriction applies (all tools enabled).
    """

    if groups is None:
        raw = os.environ.get(_TOOL_GROUPS_ENV_VAR, "").strip()
        if not raw:
            return None
        groups = [part.strip() for part in raw.split(",") if part.strip()]
    selected = list(groups)
    if not selected:
        return None
    unknown = sorted(set(selected) - set(TOOL_GROUPS))
    if unknown:
        known = ", ".join(sorted(TOOL_GROUPS))
        raise ValueError(
            f"Unknown tool group(s) {unknown!r}. Known groups: {known}"
        )
    enabled: set[str] = set()
    for group_name in selected:
        enabled |= TOOL_GROUPS[group_name]
    return frozenset(enabled)


class ToolInputModel(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")


@dataclass(frozen=True)
class ToolSpec:
    name: str
    func: Any
    input_schema: type[ToolInputModel]
    description: str | None = None
    accepts_kwargs: bool = False
    schema_overrides: dict[str, dict[str, Any]] = field(default_factory=dict)
    metadata: RuntimeToolMetadata = field(default_factory=RuntimeToolMetadata)

    def openai_tool_def(self) -> dict[str, Any]:
        schema = self._schema_with_overrides()
        return {
            "type": "function",
            "function": {
                "name": self.name,
                "description": self.description
                or inspect.getdoc(self.func)
                or self.name,
                "parameters": schema,
            },
        }

    def anthropic_tool_def(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description
            or inspect.getdoc(self.func)
            or self.name,
            "input_schema": self._schema_with_overrides(),
        }

    def _schema_with_overrides(self) -> dict[str, Any]:
        schema = self.input_schema.model_json_schema()
        schema.pop("title", None)
        schema.pop("$defs", None)
        properties = schema.get("properties", {})
        for field_name, override in self.schema_overrides.items():
            if field_name not in properties:
                continue
            properties[field_name] = {
                key: value
                for key, value in properties[field_name].items()
                if key != "title"
            }
            properties[field_name].update(override)
        return schema


class ToolRegistry:
    def __init__(self, tools: list[ToolSpec]):
        self._tools = {tool.name: tool for tool in tools}

    @classmethod
    def default(
        cls,
        groups: Iterable[str] | None = None,
    ) -> "ToolRegistry":
        enabled_tools = resolve_tool_groups(groups)
        tool_sources = [
            ("build_molecule", "chemsmart.agent.tools", None, None),
            ("recommend_method", "chemsmart.agent.tools", None, None),
            (
                "extract_project_protocol",
                "chemsmart.agent.project_yaml",
                "Extract chemsmart project-YAML method facts from a literature protocol or natural-language method description.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template="Extract project protocol facts",
                ),
            ),
            (
                "render_project_yaml",
                "chemsmart.agent.project_yaml",
                "Render a chemsmart Gaussian/ORCA project YAML candidate from extracted method facts.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template="Render project YAML {project_name}",
                ),
            ),
            (
                "validate_project_yaml",
                "chemsmart.agent.project_yaml",
                "Validate project YAML by loading it through chemsmart project settings.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template="Validate project YAML {project_name}",
                ),
            ),
            (
                "critic_project_yaml",
                "chemsmart.agent.project_yaml",
                "Critique whether project YAML matches a literature protocol and chemsmart runtime semantics.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template="Critique project YAML {project_name}",
                ),
            ),
            (
                "write_project_yaml",
                "chemsmart.agent.project_yaml",
                "Write a validated project YAML into the current workspace .chemsmart/<program> after explicit approval.",
                RuntimeToolMetadata(
                    read_only=False,
                    ui_summary_template="Write project YAML {project_name}",
                    side_effect="writes a user project YAML file",
                ),
            ),
            (
                "read_project_yaml",
                "chemsmart.agent.project_yaml",
                "Read the active workspace project YAML and summarize the chemsmart runtime settings it loads.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template="Read project YAML {project_name}",
                ),
            ),
            (
                "update_project_yaml",
                "chemsmart.agent.project_yaml",
                "Patch an existing workspace project YAML by dotted path, validate it, and write only after approval.",
                RuntimeToolMetadata(
                    read_only=False,
                    ui_summary_template="Update project YAML {project_name}",
                    side_effect="writes a workspace project YAML file",
                ),
            ),
            (
                "search_basis_sets",
                "chemsmart.agent.harness.basis_sets.catalog",
                "Search BSE-backed basis-set names for a short user phrase; returns top candidates only, never the full catalog.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template="Search basis sets {query}",
                ),
            ),
            (
                "synthesize_command",
                "chemsmart.agent.tools_command",
                "Synthesize one grounded chemsmart CLI command using the CLI schema, project YAML check, adapter, and runtime semantic gate.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template="Synthesize chemsmart command",
                ),
            ),
            (
                "repair_command",
                "chemsmart.agent.tools_command",
                "Repair a failed chemsmart CLI command and re-run runtime semantic validation.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template="Repair chemsmart command",
                ),
            ),
            (
                "execute_chemsmart_command",
                "chemsmart.agent.tools_command",
                "Execute an already semantic-gated chemsmart CLI command after explicit approval.",
                RuntimeToolMetadata(
                    read_only=False,
                    ui_summary_template="Execute chemsmart command",
                    side_effect="runs a local chemsmart command",
                ),
            ),
            ("build_gaussian_settings", "chemsmart.agent.tools", None, None),
            ("build_orca_settings", "chemsmart.agent.tools", None, None),
            ("build_job", "chemsmart.agent.tools", None, None),
            ("dry_run_input", "chemsmart.agent.tools", None, None),
            ("validate_runtime", "chemsmart.agent.tools", None, None),
            ("run_local", "chemsmart.agent.tools", None, None),
            (
                "read",
                "chemsmart.agent.tools_fs",
                "Read a local text file with 1-based line numbers. Use start_line/limit to page large files.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template="Read {path} L{start_line}-{end_line}",
                    side_effect=None,
                ),
            ),
            (
                "ssh_probe",
                "chemsmart.agent.tools_hpc",
                "Run a predefined read-only probe on a remote HPC server. probe_name must be one of the catalog entries; free-form commands are not allowed.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template=("SSH probe {probe_name} on {server}"),
                ),
            ),
            (
                "scheduler_query",
                "chemsmart.agent.tools_hpc",
                "Inspect HPC scheduler state — queue/partition aggregates (job_id omitted) or per-job status (with job_id). slurm/pbs/sge/lsf. Read-only.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template=(
                        "Scheduler query {scheduler} on {server}"
                    ),
                ),
            ),
            (
                "log_tail",
                "chemsmart.agent.tools_hpc",
                "Tail a remote log file with optional grep filter. Returns last N lines + summary of detected error signatures (OOM, walltime, missing module, node failure, scheduler reject, segfault). Read-only.",
                RuntimeToolMetadata(
                    read_only=True,
                    ui_summary_template=("Tail {path} on {server} ({lines}L)"),
                ),
            ),
            (
                "extract_optimized_geometry",
                "chemsmart.agent.tools",
                None,
                None,
            ),
            ("submit_hpc", "chemsmart.agent.tools", None, None),
            ("wizard_probe", "chemsmart.agent.wizard.tools", None, None),
            ("wizard_refresh", "chemsmart.agent.wizard.tools", None, None),
            ("wizard_verify", "chemsmart.agent.wizard.tools", None, None),
            ("wizard_write", "chemsmart.agent.wizard.tools", None, None),
        ]
        return cls(
            [
                _build_tool_spec(
                    _load_agent_tool(name, module_name),
                    registered_name=name,
                    description=description,
                    metadata=metadata,
                )
                for name, module_name, description, metadata in tool_sources
                if enabled_tools is None or name in enabled_tools
            ]
        )

    def list_tools(self) -> list[ToolSpec]:
        return list(self._tools.values())

    def get_tool(self, name: str) -> ToolSpec | None:
        return self._tools.get(name)

    def openai_tool_defs(
        self,
        tools: list[ToolSpec] | None = None,
    ) -> list[dict[str, Any]]:
        return [tool.openai_tool_def() for tool in tools or self.list_tools()]

    def anthropic_tool_defs(
        self,
        tools: list[ToolSpec] | None = None,
    ) -> list[dict[str, Any]]:
        return [
            tool.anthropic_tool_def() for tool in tools or self.list_tools()
        ]

    def tool_defs_for_provider(
        self,
        provider_name: str,
        tools: list[ToolSpec] | None = None,
    ) -> list[dict[str, Any]]:
        if provider_name == "anthropic":
            return self.anthropic_tool_defs(tools)
        return self.openai_tool_defs(tools)

    def assemble_tool_pool(
        self,
        mode: RuntimePermissionMode,
        profile: Any = None,
    ) -> list[ToolSpec]:
        del profile
        return [
            tool
            for tool in self.list_tools()
            if is_allowed_in_mode(tool, mode)
        ]

    def normalize_args(
        self,
        name: str,
        args: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        if name not in self._tools:
            return dict(args or {})

        tool = self._tools[name]
        payload = dict(args or {})
        try:
            validated = tool.input_schema.model_validate(payload)
        except Exception:
            return payload

        normalized = dict(validated.model_dump(exclude_defaults=True))
        if tool.accepts_kwargs and validated.model_extra:
            normalized.update(validated.model_extra)
        return normalized

    def describe_tool(self, name: str) -> str:
        tool = self.get_tool(name)
        if tool is None:
            return name
        doc = tool.description or inspect.getdoc(tool.func) or name
        return doc.splitlines()[0].strip()

    def call(self, name: str, args: dict[str, Any] | None = None) -> Any:
        if name not in self._tools:
            known = ", ".join(sorted(self._tools))
            raise ValueError(
                f"Unknown tool {name!r}. Registered tools: {known}"
            )

        tool = self._tools[name]
        args = args or {}
        try:
            validated = tool.input_schema.model_validate(args)
        except Exception as exc:
            return {
                "ok": False,
                "error": {
                    "type": exc.__class__.__name__,
                    "message": str(exc),
                    "tool": name,
                },
            }

        payload = dict(validated.model_dump())
        if tool.accepts_kwargs and validated.model_extra:
            payload.update(validated.model_extra)

        try:
            return tool.func(**payload)
        except Exception as exc:
            return {
                "ok": False,
                "error": {
                    "type": exc.__class__.__name__,
                    "message": str(exc),
                    "tool": name,
                },
            }


def _build_tool_spec(
    func: Any,
    registered_name: str | None = None,
    description: str | None = None,
    metadata: RuntimeToolMetadata | None = None,
) -> ToolSpec:
    fields: dict[str, Any] = {}
    schema_overrides: dict[str, dict[str, Any]] = {}
    accepts_kwargs = False
    signature = inspect.signature(func)
    resolved_hints = get_type_hints(func)
    for param in signature.parameters.values():
        if param.kind is inspect.Parameter.VAR_KEYWORD:
            accepts_kwargs = True
            continue
        if param.kind is inspect.Parameter.VAR_POSITIONAL:
            continue
        annotation = resolved_hints.get(param.name, param.annotation)
        if annotation is inspect.Signature.empty:
            annotation = Any
        schema_override = _annotation_to_schema(annotation)
        if schema_override is not None:
            schema_overrides[param.name] = schema_override
        annotation = _schema_friendly_annotation(
            annotation,
            tool_name=registered_name or func.__name__,
            field_name=param.name,
        )
        default = param.default
        if default is inspect.Signature.empty:
            fields[param.name] = (annotation, Field(...))
        else:
            fields[param.name] = (annotation, Field(default=default))

    schema_name = "".join(
        part.capitalize() for part in func.__name__.split("_")
    )
    model = create_model(
        f"{schema_name}Input",
        __base__=ToolInputModel,
        **fields,
    )
    model.model_rebuild()
    return ToolSpec(
        name=registered_name or func.__name__,
        func=func,
        input_schema=model,
        description=description,
        accepts_kwargs=accepts_kwargs,
        schema_overrides=schema_overrides,
        metadata=metadata or RuntimeToolMetadata(),
    )


def _annotation_to_schema(annotation: Any) -> dict[str, Any] | None:
    if get_origin(annotation) is not Literal:
        return None

    values = list(get_args(annotation))
    if not values:
        return None

    value_types = {type(value) for value in values}
    if value_types == {str}:
        return {"type": "string", "enum": values}
    if value_types == {int}:
        return {"type": "integer", "enum": values}
    if value_types == {float}:
        return {"type": "number", "enum": values}
    if value_types == {bool}:
        return {"type": "boolean", "enum": values}
    return {"enum": values}


def _schema_friendly_annotation(
    annotation: Any,
    *,
    tool_name: str,
    field_name: str,
) -> Any:
    try:
        TypeAdapter(annotation).json_schema()
    except (
        PydanticInvalidForJsonSchema,
        PydanticSchemaGenerationError,
        TypeError,
    ):
        logger.warning(
            "Falling back to Any for tool schema field %s.%s with annotation %r",
            tool_name,
            field_name,
            annotation,
        )
        return Any
    return annotation


def _load_agent_tool(name: str, module_name: str) -> Any:
    module = importlib.import_module(module_name)
    return getattr(module, name)
