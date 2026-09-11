"""A collection field is read the way its owner declares it.

``node_observations`` is declared ``tuple[dict[str, Any], ...]`` because
it is displayed as rows and rides a canonical body. The first consumer
that wanted one node's lines called ``.get`` on it -- a mapping's method
on a tuple -- and the ``AttributeError`` surfaced inside the executor's
pre-launch check, which killed a goal unsettled on the first live case
after the repair shipped (``sm1-formaldehyde``, 2026-09-11: 23 minutes,
0 engine calls, ledger ending at ``run_started``).

Three things were wrong at once and each gets an instrument here:

1. the shape was re-derived at the consumer instead of being asked for,
   so the owner could not keep its own promise;
2. no test pinned the launch refusal at all -- the charter sentence and
   the commit message both named the test that holds the *probe's*
   receipt, which is a different thing; and
3. the defect is mechanically detectable, which is what the lint below
   is for: a tuple-typed field is never read with ``.get``.
"""

import ast
import pathlib

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.execution import WorkflowExecutionApprovalBundleV1
from chemsmart.agent.executor import ApprovedWorkflowExecutor

PACKAGE = pathlib.Path("chemsmart")

#: The lines the host writes for an aborted probe, in the shape
#: ``_build_execution_review`` emits them.
ABORTED = (
    "input-check probe: aborted in 0.1 s of 20 s, not charged -- ORCA "
    "aborted the run at its input check",
    "  orca: Use one of the keywords",
    "  orca: RIJDX : treat the HF exchange exactly (equals RIJONX)",
    "  orca: [file orca_main/main_input_check.cpp, line 3594]: Error "
    "(ORCA_MAIN): ... aborting the run",
)


def _tuple_typed_field_names() -> frozenset[str]:
    """Every dataclass field in the execution plane annotated as a tuple.

    Computed from the annotations rather than listed, so a field added
    later is covered on the day it is added.
    """

    import chemsmart.agent.execution as execution

    names: set[str] = set()
    for value in vars(execution).values():
        fields = getattr(value, "__dataclass_fields__", None)
        if not fields:
            continue
        for name, field in fields.items():
            if str(field.type).replace(" ", "").startswith("tuple["):
                names.add(name)
    return frozenset(names)


def _names_that_are_ever_a_mapping() -> frozenset[str]:
    """Names the package itself declares or assigns as a mapping.

    A name-only scan over-captures, which is the same mistake a marker
    regex without a left boundary made earlier in this round:
    ``requested_observable_declarations`` is a tuple on the approval
    bundle *and* a genuine ``dict[str, dict[str, Any]]`` on the tool
    host, so reading it with ``.get`` there is correct. A name is only in
    scope for this lint when nothing in the package ever gives it a
    mapping, which makes the rule exactly "a name that means a tuple
    everywhere is never read as a mapping".
    """

    mapping_names: set[str] = set()
    for path in sorted(PACKAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            target = None
            annotation = None
            if isinstance(node, ast.AnnAssign):
                target, annotation = node.target, node.annotation
            elif isinstance(node, ast.Assign) and len(node.targets) == 1:
                target = node.targets[0]
            if target is None:
                continue
            name = (
                target.attr
                if isinstance(target, ast.Attribute)
                else target.id if isinstance(target, ast.Name) else ""
            )
            if not name:
                continue
            spelled = ast.unparse(annotation) if annotation else ""
            value = node.value
            assigns_mapping = isinstance(value, ast.Dict) or (
                isinstance(value, ast.Call)
                and isinstance(value.func, ast.Name)
                and value.func.id == "dict"
            )
            if (
                spelled.replace(" ", "").startswith(
                    ("dict[", "Mapping[", "MutableMapping[")
                )
                or assigns_mapping
            ):
                mapping_names.add(name)
    return frozenset(mapping_names)


def _mapping_reads_of_tuple_fields() -> list[str]:
    """``something.<tuple field>.get(...)`` anywhere in the package."""

    wanted = _tuple_typed_field_names() - _names_that_are_ever_a_mapping()
    offenders: list[str] = []
    for path in sorted(PACKAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            method = node.func
            if not isinstance(method, ast.Attribute) or method.attr != "get":
                continue
            owner = method.value
            # ``x.field.get(...)`` and ``getattr(x, "field", {}).get(...)``
            name = ""
            if isinstance(owner, ast.Attribute):
                name = owner.attr
            elif (
                isinstance(owner, ast.Call)
                and isinstance(owner.func, ast.Name)
                and owner.func.id == "getattr"
                and len(owner.args) >= 2
                and isinstance(owner.args[1], ast.Constant)
                and isinstance(owner.args[1].value, str)
            ):
                name = owner.args[1].value
            if name in wanted:
                offenders.append(f"{path}:{node.lineno}: {name}.get(...)")
    return offenders


@pytest.mark.capability("gate:executor.launch_refuses_a_refused_input")
def test_no_tuple_field_is_read_with_get():
    """The lint that would have caught it before it shipped."""

    scoped = _tuple_typed_field_names() - _names_that_are_ever_a_mapping()
    assert "node_observations" in scoped, (
        "the field whose mapping read killed a live goal is out of this "
        "lint's scope, so the lint no longer covers the defect it exists "
        f"for: {sorted(scoped)[:8]}"
    )
    offenders = _mapping_reads_of_tuple_fields()
    assert not offenders, (
        "these read a tuple-annotated field as if it were a mapping, "
        "which raises AttributeError at the call site rather than at the "
        f"declaration: {offenders}"
    )


@pytest.mark.capability("gate:executor.launch_refuses_a_refused_input")
def test_the_lint_catches_a_planted_offender(tmp_path):
    """A lint that cannot fail is not a lint."""

    wanted = sorted(
        _tuple_typed_field_names() - _names_that_are_ever_a_mapping()
    )
    planted = tmp_path / "offender.py"
    planted.write_text(
        f"def read(bundle):\n" f'    return bundle.{wanted[0]}.get("node")\n',
        encoding="utf-8",
    )
    tree = ast.parse(planted.read_text(encoding="utf-8"))
    hits = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "get"
        and isinstance(node.func.value, ast.Attribute)
        and node.func.value.attr in set(wanted)
    ]
    assert hits, "the planted mapping read was not detected"


@pytest.mark.capability("gate:executor.launch_refuses_a_refused_input")
def test_the_owner_answers_for_one_node():
    """The reduction lives on the class that declares the shape."""

    bundle = WorkflowExecutionApprovalBundleV1.__new__(
        WorkflowExecutionApprovalBundleV1
    )
    object.__setattr__(
        bundle,
        "node_observations",
        (
            {"node_id": "scan-c4-r2", "observations": ABORTED},
            {"node_id": "opt-clean", "observations": ("geom cap: 3N = 63",)},
        ),
    )
    assert bundle.node_observation_lines("scan-c4-r2") == ABORTED
    assert bundle.node_observation_lines("opt-clean") == ("geom cap: 3N = 63",)
    assert bundle.node_observation_lines("absent") == ()


@pytest.mark.capability("gate:executor.launch_refuses_a_refused_input")
def test_an_aborted_check_refuses_the_launch_and_a_clean_one_does_not():
    """The gate itself, which nothing pinned when it shipped.

    A crash is not a refusal: the assertion is on ``ContractError``
    specifically, because the shipped defect raised ``AttributeError``
    here and an exception-agnostic check would have called that a pass.
    """

    bundle = WorkflowExecutionApprovalBundleV1.__new__(
        WorkflowExecutionApprovalBundleV1
    )
    object.__setattr__(
        bundle,
        "node_observations",
        (
            {"node_id": "scan-c4-r2", "observations": ABORTED},
            {"node_id": "opt-clean", "observations": ("geom cap: 3N = 63",)},
        ),
    )
    executor = ApprovedWorkflowExecutor.__new__(ApprovedWorkflowExecutor)
    object.__setattr__(executor, "execution_bundle", bundle)
    check = executor._refuse_launch_the_program_already_refused

    with pytest.raises(ContractError) as refused:
        check("scan-c4-r2")
    spoken = str(refused.value)
    assert "is not launched" in spoken
    assert "RIJDX" in spoken, (
        "the refusal does not quote the program's own lines, so the "
        f"session cannot see which field to repair: {spoken}"
    )
    assert (
        "compile the node again" in spoken
    ), f"the refusal names no route out of itself: {spoken}"

    # A node whose check passed, and a node with no observations at all,
    # are both admitted in silence.
    check("opt-clean")
    check("never-probed")
