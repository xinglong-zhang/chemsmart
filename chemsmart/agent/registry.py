from __future__ import annotations

import importlib
import inspect
import logging
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
    def default(cls) -> "ToolRegistry":
        tool_sources = [
            ("build_molecule", "chemsmart.agent.tools", None, None),
            ("recommend_method", "chemsmart.agent.tools", None, None),
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
