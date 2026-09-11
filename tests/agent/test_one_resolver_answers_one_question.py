"""Two host organs that answer one question call one function.

``CONDUCT.md`` states this rule and nothing enforced it, which is why it
has now been broken three times by three different mechanisms:

- the frontier admitted what the review refused, because each computed
  readiness its own way;
- ``scientific_workflow_plans`` is keyed by digest and three readers
  resolved "the current plan" as ``tuple(values())[-1]``. A dict
  re-insertion keeps its key in place, so after plan A, plan B, and an
  amendment restoring A, all three stood on B -- the plan the session
  had abandoned. ``dd29df23`` added ``_current_scientific_plan`` and
  repaired the three readers spelled ``tuple(...)[-1]``;
- **it missed a fourth**, spelled as a list comprehension with an extra
  ``if plan.nodes`` filter, inside ``prepare_execution_review`` -- the
  gate on the single human decision. Measured 2026-09-11 through the
  public planning tool: after A -> B -> A the eligibility loop answered
  B while the reviewed packet answered A, so a node ineligible in one
  and eligible in the other crosses into an approved packet unchecked.
  That survived a targeted repair, a 63-witness bank and 2,146 passing
  tests, because the repair matched a *spelling*.

So the enforcement cannot be a spelling either. This lint names the
owner of a question and forbids anyone else from reading its state
directly, whatever the expression looks like.
"""

import ast
import pathlib

import pytest

AGENT = pathlib.Path(__file__).resolve().parents[2] / "chemsmart"

#: attribute -> (owning function, why the question has one answer).
#: A reader outside the owner must call the owner instead.
OWNED_STATE = {
    "scientific_workflow_plans": (
        "_current_scientific_plan",
        "which plan the session is standing on; insertion order answers "
        "about the plan a restoring amendment abandoned",
    ),
}


def _order_dependent_readers(
    root: pathlib.Path, attribute: str, owner: str
) -> list[str]:
    """Functions that pick one entry out of ``self.<attribute>`` by position.

    The forbidden act is narrow and exact: reducing a **digest-keyed**
    container to a single entry by its *position* -- ``[-1]``, ``[0]``,
    ``tuple(...)[-1]``, a list comprehension then ``[-1]``. A dict keeps
    a re-inserted key in place, so position answers about whatever was
    inserted first, which after a restoring amendment is the plan the
    session abandoned.

    Deliberately allowed, because they are different questions:
    aggregating over every entry (``any(plan.nodes for plan in ...)``),
    collecting from all of them (guide activation reads every plan's
    jobtypes), fetching one by an id the caller already holds, writing,
    and ``next(... if <criterion>)`` -- searching by a property is not
    resolving by order.
    """

    offenders: list[str] = []
    for path in sorted(root.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except (OSError, SyntaxError):
            continue
        for function in [
            node
            for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        ]:
            if function.name == owner:
                continue
            reads_container = any(
                isinstance(node, ast.Attribute)
                and node.attr == attribute
                and isinstance(node.value, ast.Name)
                and node.value.id == "self"
                for node in ast.walk(function)
            )
            if not reads_container:
                continue
            positional = False
            for node in ast.walk(function):
                if not isinstance(node, ast.Subscript):
                    continue
                index = node.slice
                if isinstance(index, ast.UnaryOp) and isinstance(
                    index.op, ast.USub
                ):
                    index = index.operand
                if not isinstance(index, ast.Constant):
                    continue
                if not isinstance(index.value, int):
                    continue
                # the subscripted thing must derive from the container
                if any(
                    isinstance(inner, ast.Attribute)
                    and inner.attr == attribute
                    for inner in ast.walk(node.value)
                ):
                    positional = True
                    break
                # ... or from a local bound to a comprehension over it
                if isinstance(node.value, ast.Name):
                    for assign in ast.walk(function):
                        if not isinstance(assign, ast.Assign):
                            continue
                        if not any(
                            isinstance(t, ast.Name) and t.id == node.value.id
                            for t in assign.targets
                        ):
                            continue
                        if any(
                            isinstance(inner, ast.Attribute)
                            and inner.attr == attribute
                            for inner in ast.walk(assign.value)
                        ):
                            positional = True
                            break
                if positional:
                    break
            if positional:
                offenders.append(
                    f"{path.relative_to(root.parent)}::{function.name}"
                )
    return sorted(set(offenders))


@pytest.mark.capability("gate:resolver.one_answer_per_question")
def test_owned_state_is_read_only_by_its_owner():
    """Every reader of an ambiguous question calls the one resolver."""

    findings = []
    for attribute, (owner, why) in OWNED_STATE.items():
        offenders = _order_dependent_readers(AGENT, attribute, owner)
        if offenders:
            findings.append(
                f"self.{attribute} ({why}) is resolved directly by: "
                + ", ".join(offenders)
                + f"; call {owner}() instead"
            )
    assert not findings, "\n".join(findings)


@pytest.mark.capability("gate:resolver.one_answer_per_question")
def test_the_lint_catches_the_shape_it_forbids(tmp_path):
    """Falsifiable: plant the offender and require detection.

    Without this, a lint whose matcher silently found nothing would pass
    over a package full of offenders -- the class of defect this
    laboratory has paid for before, most recently when my own
    53-parameter audit reported `basis` and `functional` as broken while
    they were demonstrably reaching the program.
    """

    package = tmp_path / "pkg"
    package.mkdir()
    (package / "offender.py").write_text(
        "class Host:\n"
        "    def _current_scientific_plan(self):\n"
        "        return tuple(self.scientific_workflow_plans.values())[-1]\n"
        "\n"
        "    def some_other_reader(self):\n"
        "        plans = [\n"
        "            plan\n"
        "            for plan in self.scientific_workflow_plans.values()\n"
        "            if plan.nodes\n"
        "        ]\n"
        "        return plans[-1]\n",
        encoding="utf-8",
    )
    found = _order_dependent_readers(
        package, "scientific_workflow_plans", "_current_scientific_plan"
    )
    assert any("some_other_reader" in item for item in found), found
    assert not any(
        "_current_scientific_plan" in item.split("::")[-1] for item in found
    ), f"the owner must be exempt: {found}"


@pytest.mark.capability("gate:resolver.one_answer_per_question")
def test_a_lookup_by_a_held_id_is_not_the_ambiguous_question(tmp_path):
    """Reading one plan by a digest you already hold is not resolving.

    The forbidden act is asking the container *which* plan is current.
    Fetching a named entry is ordinary, and forbidding it would push
    callers into worse shapes.
    """

    package = tmp_path / "pkg"
    package.mkdir()
    (package / "fine.py").write_text(
        "class Host:\n"
        "    def reader(self, digest):\n"
        "        return self.scientific_workflow_plans.get(digest)\n"
        "\n"
        "    def writer(self, plan):\n"
        "        self.scientific_workflow_plans[plan.sha] = plan\n"
        "\n"
        "    def membership(self, digest):\n"
        "        return digest in self.scientific_workflow_plans\n",
        encoding="utf-8",
    )
    assert not _order_dependent_readers(
        package, "scientific_workflow_plans", "_current_scientific_plan"
    )
