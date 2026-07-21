from __future__ import annotations

import pytest
from pydantic import create_model

from chemsmart.agent.permissions import RuntimePermissionMode
from chemsmart.agent.registry import ToolInputModel, ToolRegistry, ToolSpec
from chemsmart.agent.tool_protocol import RuntimeToolMetadata
from chemsmart.jobs.gaussian.settings import GaussianJobSettings


def test_registry_round_trips_build_gaussian_settings_call():
    registry = ToolRegistry.default()

    result = registry.call(
        "build_gaussian_settings",
        {"functional": "B3LYP", "basis": "6-31G*", "charge": 0},
    )

    assert isinstance(result, GaussianJobSettings)
    assert result.functional == "B3LYP"
    assert result.basis == "6-31G*"

    tool_defs = registry.openai_tool_defs()
    names = [tool_def["function"]["name"] for tool_def in tool_defs]
    assert "build_gaussian_settings" in names
    assert "wizard_probe" in names
    assert "wizard_refresh" in names
    assert "wizard_verify" in names
    assert "wizard_write" in names
    build_job_def = next(
        tool_def
        for tool_def in tool_defs
        if tool_def["function"]["name"] == "build_job"
    )
    gaussian_settings_def = next(
        tool_def
        for tool_def in tool_defs
        if tool_def["function"]["name"] == "build_gaussian_settings"
    )
    assert build_job_def["function"]["parameters"]["properties"]["kind"] == {
        "type": "string",
        "enum": [
            "gaussian.opt",
            "gaussian.ts",
            "gaussian.freq",
            "gaussian.sp",
            "gaussian.irc",
            "gaussian.scan",
            "orca.opt",
            "orca.ts",
            "orca.freq",
            "orca.sp",
            "orca.irc",
            "orca.scan",
            "xtb.opt",
            "xtb.sp",
            "xtb.hess",
        ],
    }
    gaussian_settings_props = gaussian_settings_def["function"]["parameters"][
        "properties"
    ]
    assert "title" in gaussian_settings_props
    assert "freq" in gaussian_settings_props
    assert "numfreq" in gaussian_settings_props
    assert "additional_opt_options_in_route" in gaussian_settings_props
    assert "additional_route_parameters" in gaussian_settings_props


def test_registry_exposes_typed_handles_and_hides_runtime_owned_fields():
    registry = ToolRegistry.default()
    definitions = {
        item["function"]["name"]: item["function"]["parameters"]
        for item in registry.openai_tool_defs()
    }

    build_job = definitions["build_job"]
    # extract_optimized_geometry results (geom handles) feed chained jobs,
    # e.g. an xTB pre-optimization consumed by a Gaussian/ORCA build_job.
    assert build_job["properties"]["molecule"]["pattern"].startswith(
        "^(?:mol|geom)"
    )
    save_geometry = definitions["save_geometry"]
    assert save_geometry["properties"]["molecule"]["pattern"].startswith(
        "^(?:mol|geom)"
    )
    assert build_job["properties"]["settings"]["pattern"].startswith(
        "^(?:gset|oset|xset)"
    )
    assert "jobrunner" not in build_job["properties"]

    for tool_name in {
        "dry_run_input",
        "validate_runtime",
        "run_local",
        "extract_optimized_geometry",
        "submit_hpc",
    }:
        job_schema = definitions[tool_name]["properties"]["job"]
        assert job_schema["type"] == "string"
        assert job_schema["pattern"].startswith("^job_")

    assert "transport" not in definitions["submit_hpc"]["properties"]


def test_registry_model_facing_properties_never_use_empty_any_schema():
    registry = ToolRegistry.default()

    empty_properties = []
    for definition in registry.openai_tool_defs():
        function = definition["function"]
        for name, schema in function["parameters"]["properties"].items():
            if not schema:
                empty_properties.append(f"{function['name']}.{name}")

    assert empty_properties == []


def test_registry_unknown_tool_name_raises_clearly():
    registry = ToolRegistry.default()

    with pytest.raises(ValueError) as excinfo:
        registry.call("does_not_exist", {})

    message = str(excinfo.value)
    assert "does_not_exist" in message
    assert "build_molecule" in message


def test_registry_forbids_unexpected_tool_arguments():
    registry = ToolRegistry.default()

    result = registry.call(
        "build_gaussian_settings",
        {
            "functional": "B3LYP",
            "basis": "6-31G*",
            "unexpected": "boom",
        },
    )

    assert result["ok"] is False
    assert result["error"]["type"] == "ValidationError"
    assert "unexpected" in result["error"]["message"]


def test_registry_exposes_machine_readable_tool_descriptions():
    registry = ToolRegistry.default()

    build_molecule_def = next(
        tool_def
        for tool_def in registry.openai_tool_defs()
        if tool_def["function"]["name"] == "build_molecule"
    )

    assert "Load one molecule" in build_molecule_def["function"]["description"]


def test_registry_default_registration_sets_read_tool_metadata():
    registry = ToolRegistry.default()

    read_tool = registry.get_tool("read")
    ssh_probe_tool = registry.get_tool("ssh_probe")
    scheduler_query_tool = registry.get_tool("scheduler_query")
    log_tail_tool = registry.get_tool("log_tail")
    project_tools = {
        name: registry.get_tool(name)
        for name in {
            "extract_project_protocol",
            "render_project_yaml",
            "validate_project_yaml",
            "critic_project_yaml",
            "write_project_yaml",
            "read_project_yaml",
            "update_project_yaml",
            "search_basis_sets",
            "synthesize_command",
            "repair_command",
            "execute_chemsmart_command",
            "inspect_calculation",
        }
    }

    assert read_tool is not None
    assert ssh_probe_tool is not None
    assert scheduler_query_tool is not None
    assert log_tail_tool is not None
    assert all(tool is not None for tool in project_tools.values())
    assert read_tool.metadata == RuntimeToolMetadata(
        read_only=True,
        ui_summary_template="Read {path} L{start_line}-{end_line}",
        side_effect=None,
    )
    assert ssh_probe_tool.metadata == RuntimeToolMetadata(
        read_only=True,
        ui_summary_template="SSH probe {probe_name} on {server}",
    )
    assert scheduler_query_tool.metadata == RuntimeToolMetadata(
        read_only=True,
        ui_summary_template="Scheduler query {scheduler} on {server}",
    )
    assert log_tail_tool.metadata == RuntimeToolMetadata(
        read_only=True,
        ui_summary_template="Tail {path} on {server} ({lines}L)",
    )
    assert project_tools["extract_project_protocol"].metadata.read_only is True
    assert project_tools["render_project_yaml"].metadata.read_only is True
    assert project_tools["validate_project_yaml"].metadata.read_only is True
    assert project_tools["critic_project_yaml"].metadata.read_only is True
    assert project_tools["write_project_yaml"].metadata.read_only is False
    assert project_tools["read_project_yaml"].metadata.read_only is True
    assert project_tools["update_project_yaml"].metadata.read_only is False
    assert project_tools["synthesize_command"].metadata.read_only is True
    assert project_tools["repair_command"].metadata.read_only is True
    assert (
        project_tools["execute_chemsmart_command"].metadata.read_only is False
    )
    assert project_tools["search_basis_sets"].metadata == RuntimeToolMetadata(
        read_only=True,
        ui_summary_template="Search basis sets {query}",
    )
    assert project_tools[
        "inspect_calculation"
    ].metadata == RuntimeToolMetadata(
        read_only=True,
        ui_summary_template="Inspect calculation {run_id}",
    )
    list_workspace_tool = registry.get_tool("list_workspace")
    assert list_workspace_tool is not None
    assert list_workspace_tool.metadata == RuntimeToolMetadata(
        read_only=True,
        ui_summary_template="List workspace {subdir}",
    )
    save_geometry_tool = registry.get_tool("save_geometry")
    assert save_geometry_tool is not None
    assert save_geometry_tool.metadata == RuntimeToolMetadata(
        read_only=False,
        ui_summary_template="Save geometry {filename}",
        side_effect="writes a geometry file in the workspace",
    )
    assert all(
        tool.metadata == RuntimeToolMetadata()
        for tool in registry.list_tools()
        if tool.name
        not in {
            "read",
            "list_workspace",
            "save_geometry",
            "ssh_probe",
            "scheduler_query",
            "log_tail",
            "extract_project_protocol",
            "render_project_yaml",
            "validate_project_yaml",
            "critic_project_yaml",
            "write_project_yaml",
            "read_project_yaml",
            "update_project_yaml",
            "search_basis_sets",
            "synthesize_command",
            "repair_command",
            "execute_chemsmart_command",
            "inspect_calculation",
        }
    )


def test_assemble_tool_pool_plan_mode_returns_empty_list():
    registry = ToolRegistry([_make_tool_spec("legacy_tool")])

    pool = registry.assemble_tool_pool(RuntimePermissionMode.PLAN)

    assert pool == []


@pytest.mark.parametrize(
    ("mode", "tool_name", "expected"),
    [
        (RuntimePermissionMode.READ_ONLY, "legacy_tool", True),
        (RuntimePermissionMode.READ_ONLY, "read_only_tool", True),
        (RuntimePermissionMode.READ_ONLY, "edit_safe_tool", False),
        (RuntimePermissionMode.READ_ONLY, "bypass_only_tool", False),
        (RuntimePermissionMode.ACCEPT_EDITS, "legacy_tool", True),
        (RuntimePermissionMode.ACCEPT_EDITS, "read_only_tool", True),
        (RuntimePermissionMode.ACCEPT_EDITS, "edit_safe_tool", True),
        (RuntimePermissionMode.ACCEPT_EDITS, "bypass_only_tool", False),
        (RuntimePermissionMode.BYPASS, "legacy_tool", True),
        (RuntimePermissionMode.BYPASS, "read_only_tool", True),
        (RuntimePermissionMode.BYPASS, "edit_safe_tool", True),
        (RuntimePermissionMode.BYPASS, "bypass_only_tool", True),
        (RuntimePermissionMode.PLAN, "legacy_tool", False),
        (RuntimePermissionMode.PLAN, "read_only_tool", False),
        (RuntimePermissionMode.PLAN, "edit_safe_tool", False),
        (RuntimePermissionMode.PLAN, "bypass_only_tool", False),
    ],
)
def test_assemble_tool_pool_runtime_mode_matrix(mode, tool_name, expected):
    registry = ToolRegistry(
        [
            _make_tool_spec("legacy_tool"),
            _make_tool_spec(
                "read_only_tool",
                metadata=RuntimeToolMetadata(read_only=True),
            ),
            _make_tool_spec(
                "edit_safe_tool",
                metadata=RuntimeToolMetadata(edit_safe=True),
            ),
            _make_tool_spec(
                "bypass_only_tool",
                metadata=RuntimeToolMetadata(
                    requires_mode=frozenset({RuntimePermissionMode.BYPASS})
                ),
            ),
        ]
    )

    names = {tool.name for tool in registry.assemble_tool_pool(mode)}

    assert (tool_name in names) is expected


@pytest.mark.parametrize(
    ("mode", "expected_present"),
    [
        (RuntimePermissionMode.READ_ONLY, True),
        (RuntimePermissionMode.ACCEPT_EDITS, True),
        (RuntimePermissionMode.BYPASS, True),
        (RuntimePermissionMode.PLAN, False),
    ],
)
def test_default_registry_exposes_ssh_probe_in_runtime_modes(
    mode,
    expected_present,
):
    registry = ToolRegistry.default()

    names = {tool.name for tool in registry.assemble_tool_pool(mode)}

    assert ("ssh_probe" in names) is expected_present


def _make_tool_spec(
    name: str,
    *,
    metadata: RuntimeToolMetadata | None = None,
) -> ToolSpec:
    def tool() -> dict[str, str]:
        """Test tool."""

        return {"tool": name}

    schema = create_model(
        f"{name.title().replace('_', '')}Input",
        __base__=ToolInputModel,
    )
    return ToolSpec(
        name=name,
        func=tool,
        input_schema=schema,
        metadata=metadata or RuntimeToolMetadata(),
    )


def test_tool_groups_cover_every_registered_tool_exactly_once():
    from chemsmart.agent.registry import TOOL_GROUPS

    registry = ToolRegistry.default()
    all_names = {tool.name for tool in registry.list_tools()}
    grouped: list[str] = [
        name for names in TOOL_GROUPS.values() for name in names
    ]

    assert set(grouped) == all_names  # no ungrouped, no phantom entries
    assert len(grouped) == len(set(grouped))  # each tool in exactly one group


def test_default_registry_restricts_to_selected_groups():
    registry = ToolRegistry.default(groups=["synthesis", "project_yaml"])
    names = {tool.name for tool in registry.list_tools()}

    assert "synthesize_command" in names
    assert "read_project_yaml" in names
    assert "run_local" not in names
    assert "execute_chemsmart_command" not in names


def test_default_registry_honors_tool_groups_env(monkeypatch):
    monkeypatch.setenv("CHEMSMART_AGENT_TOOL_GROUPS", "synthesis")

    registry = ToolRegistry.default()

    assert {tool.name for tool in registry.list_tools()} == {
        "synthesize_command",
        "repair_command",
    }


def test_default_registry_rejects_unknown_group():
    with pytest.raises(ValueError, match="Unknown tool group"):
        ToolRegistry.default(groups=["bogus"])


def test_tool_group_lookup():
    from chemsmart.agent.registry import tool_group

    assert tool_group("synthesize_command") == "synthesis"
    assert tool_group("execute_chemsmart_command") == "execution"
    assert tool_group("write_project_yaml") == "project_yaml"
    assert tool_group("nonexistent_tool") is None
