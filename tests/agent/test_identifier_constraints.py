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

from chemsmart.agent._contracts import ContractError

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


def test_the_optional_class_is_what_the_validators_accept_not_a_guess():
    """artifact_id was misclassed on the first attempt and a live run hit it.

    The two halves of a workflow edge are mutually exclusive: an input bound to
    an initial artifact names artifact_id and leaves the producer fields empty;
    an input fed by an upstream node names the producer fields and leaves
    artifact_id empty.  Constraining either half forbids the other, and the
    live session that hit it was chaining an optimization into a single point --
    the shape of every composite protocol.

    Membership is therefore probed against the runtime rather than asserted.
    """

    from chemsmart.agent._contracts import ContractError
    from chemsmart.agent.tool_specs import (
        IDENTIFIER_ARGUMENTS,
        OPTIONAL_IDENTIFIER_ARGUMENTS,
    )
    from chemsmart.agent.workflows import (
        ArtifactInputIntentV1,
        ArtifactOutputIntentV1,
    )

    def _accepts_empty(build):
        try:
            build()
        except ContractError:
            return False
        return True

    probes = {
        "artifact_id": lambda: ArtifactInputIntentV1(
            binding_id="g",
            artifact_class="geometry_xyz",
            artifact_id="",
            producer_node_id="opt",
            producer_output_id="geom",
        ),
        "producer_node_id": lambda: ArtifactInputIntentV1(
            binding_id="g",
            artifact_class="geometry_xyz",
            artifact_id="start.xyz",
            producer_node_id="",
            producer_output_id="",
        ),
        "artifact_class": lambda: ArtifactInputIntentV1(
            binding_id="g",
            artifact_class="",
            artifact_id="start.xyz",
            producer_node_id="",
            producer_output_id="",
        ),
        "output_id": lambda: ArtifactOutputIntentV1(
            output_id="", artifact_class="energy"
        ),
    }
    for name, build in probes.items():
        permissive = _accepts_empty(build)
        classed_optional = name in OPTIONAL_IDENTIFIER_ARGUMENTS
        assert permissive == classed_optional, (
            f"{name}: the runtime "
            f"{'accepts' if permissive else 'refuses'} an empty value but the "
            f"schema classes it as {'optional' if classed_optional else 'strict'}"
        )
        assert name in (IDENTIFIER_ARGUMENTS | OPTIONAL_IDENTIFIER_ARGUMENTS)


def test_a_producer_fed_input_passes_the_schema_it_previously_failed():
    """The exact argument that a live composite-protocol plan was rejected on."""

    from chemsmart.agent.tool_runtime import _validate_json_value

    surface = build_command_compiled_tool_surface()
    definition = next(
        item["function"]
        for item in surface.tool_definitions
        if item["function"]["name"] == "plan_scientific_workflow"
    )
    schema = definition["parameters"]["properties"]["calculation_nodes"][
        "items"
    ]["properties"]["inputs"]["items"]["properties"]["artifact_id"]
    _validate_json_value("calculation_nodes[1].inputs[0].artifact_id", "", schema)
    _validate_json_value(
        "calculation_nodes[0].inputs[0].artifact_id", "start.xyz", schema
    )
    with pytest.raises(ContractError):
        _validate_json_value("x", "Start.XYZ", schema)
