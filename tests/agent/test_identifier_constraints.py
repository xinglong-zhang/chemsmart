"""An identifier rule must be visible before the call, not only after it.

Measured: 36 of 37 arguments the runtime passes through its public-identifier
validator were exposed to the model as unconstrained strings.  The rule existed
only in the rejection.

Two live sessions wrote whole sentences into identifier fields and learned the
rule the hard way -- a dependency phrase copied from the workflow frontier into
``stage_order``, and

    'gas-phase B3LYP/aug-cc-pVTZ geometry optimization of neutral closed-shell
     NH3 with conventional four-index integrals (no density fitting)'

into ``project_role``.  Both are reasonable things to write if nothing says the
field wants a bare identifier.
"""

import re

import pytest

from chemsmart.agent.tool_specs import (
    IDENTIFIER_ARGUMENTS,
    OPTIONAL_IDENTIFIER_ARGUMENTS,
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
)

_PATTERN = "^[a-z][a-z0-9_.-]*$"
_OPTIONAL = "^$|^[a-z][a-z0-9_.-]*$"


def _identifier_schemas(surface):
    def walk(path, schema):
        if not isinstance(schema, dict):
            return
        for key, value in (schema.get("properties") or {}).items():
            if key in (
                IDENTIFIER_ARGUMENTS | OPTIONAL_IDENTIFIER_ARGUMENTS
            ) and isinstance(value, dict):
                yield f"{path}.{key}", value
            yield from walk(f"{path}.{key}", value)
        items = schema.get("items")
        if isinstance(items, dict):
            yield from walk(path + "[]", items)

    for item in surface.tool_definitions:
        function = item["function"]
        yield from walk(function["name"], function.get("parameters") or {})


@pytest.mark.parametrize(
    "build",
    [build_command_compiled_tool_surface, build_approved_execution_tool_surface],
)
def test_every_identifier_argument_declares_its_rule(build):
    loose = [
        path
        for path, schema in _identifier_schemas(build())
        if schema.get("type") == "string"
        and schema.get("pattern") not in (_PATTERN, _OPTIONAL)
        and not schema.get("enum")
    ]
    assert not loose, (
        f"{len(loose)} identifier arguments are exposed as unconstrained "
        f"strings: {loose[:10]}"
    )


def test_the_constraint_reaches_inside_arrays_of_objects():
    """Workflow nodes live in an array; that is where the failures happened."""

    surface = build_command_compiled_tool_surface()
    paths = {path for path, _ in _identifier_schemas(surface)}
    assert any("[]" in path for path in paths), paths
    node_fields = {
        path.rsplit(".", 1)[-1]
        for path in paths
        if "calculation_nodes[]" in path or "nodes[]" in path
    }
    assert {"node_id", "program", "jobtype"} <= node_fields, node_fields


def test_the_declared_pattern_is_the_one_the_runtime_enforces():
    """A schema stricter or looser than the validator would mislead.

    ScientificWorkflowNodeV2 is the node type whose validator produced the
    live rejection; CommandNodeIntentV1 does not guard project_role, which is
    itself worth knowing -- the rule is not uniform across node types even
    though the argument name is.
    """

    from chemsmart.agent._contracts import ContractError
    from chemsmart.agent.workflows import ScientificWorkflowNodeV2

    compiled = re.compile(_PATTERN)
    assert compiled.match("sp-ccpvqz")
    assert compiled.match("opt_anti")
    assert not compiled.match("Gas-phase B3LYP optimization")
    assert not compiled.match("2-stage")

    def _node(project_role):
        return ScientificWorkflowNodeV2(
            node_id="opt",
            stage="opt",
            requested_program="pyscf",
            program="pyscf",
            engine="cpu",
            project_role=project_role,
            unresolved_fields=(),
        )

    _node("nh3-b3lyp-avtz-opt")
    with pytest.raises(ContractError, match="project_role"):
        _node(
            "gas-phase B3LYP/aug-cc-pVTZ geometry optimization of neutral "
            "closed-shell NH3 with conventional four-index integrals"
        )
