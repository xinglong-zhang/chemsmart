"""One calculation, one folder a chemist can read unaided.

Evidence for a calculation was scattered: the engine's own files landed in
``<workspace>/nodes/<node_id>/`` while the project YAML it was run from,
the exact CHEMSMART command that ran, and the receipts lived elsewhere or
only inside an event stream. A human asking "what was actually computed
here, and how" had to reconstruct it from three places and a JSONL file.

Under a wave cohort that gets worse, not better: several independent
calculations land in one cycle, and the cycle is the reasoning unit while
the *calculation* is the evidence unit. So the branch is the folder, the
cycle is its parent, and everything one calculation produced or consumed
is inside it (owner ruling, 2026-09-16).

The layout is one authority because four call sites answered it
separately: the launch site, the plan-time occupied-id guard, and two
workspace scans that must not mistake a node's outputs for a user's own
input geometry.
"""

from __future__ import annotations

from pathlib import Path

from chemsmart.agent.execution import (
    existing_node_branches,
    node_branch_directory,
    node_branch_roots,
)


def test_a_branch_lives_under_its_cycle(tmp_path):
    branch = node_branch_directory(
        tmp_path, "conformer-a-opt", cycle_label="cycle-0003"
    )
    assert branch == tmp_path / "cycle-0003" / "conformer-a-opt"


def test_without_a_cycle_the_legacy_layout_is_unchanged(tmp_path):
    """Every recorded run keeps the path it was written at."""

    branch = node_branch_directory(tmp_path, "water-opt", cycle_label=None)
    assert branch == tmp_path / "nodes" / "water-opt"


def test_every_branch_root_is_barred_from_a_workspace_scan(tmp_path):
    """A node's own output is not a user-supplied geometry.

    The scan bars ``nodes/``; a cycle folder is the same thing under a
    new name, and a scan that missed it would offer an engine's output
    back to the model as though a human had put it there.
    """

    (tmp_path / "nodes").mkdir()
    (tmp_path / "cycle-0001").mkdir()
    (tmp_path / "cycle-0002").mkdir()
    (tmp_path / "inputs").mkdir()
    roots = node_branch_roots(tmp_path)
    assert tmp_path / "nodes" in roots
    assert tmp_path / "cycle-0001" in roots
    assert tmp_path / "cycle-0002" in roots
    assert tmp_path / "inputs" not in roots


def test_a_node_id_is_used_once_per_workspace_across_cycles(tmp_path):
    """The occupied-id refusal keeps its meaning under the new layout.

    'A re-run takes a fresh id; the earlier directory is evidence' was
    enforced by looking in one place. Per-cycle folders would have made
    the same id in a later cycle look unused, quietly weakening a
    recorded refusal, so the guard looks in every cycle.
    """

    occupied = tmp_path / "cycle-0001" / "water-opt"
    occupied.mkdir(parents=True)
    (occupied / "water.log").write_text("output", encoding="utf-8")

    assert existing_node_branches(tmp_path, "water-opt") == (occupied,)
    assert existing_node_branches(tmp_path, "water-sp") == ()

    # An empty directory is not evidence and does not occupy the id.
    (tmp_path / "cycle-0002" / "water-sp").mkdir(parents=True)
    assert existing_node_branches(tmp_path, "water-sp") == ()


def test_the_launch_site_and_the_guard_read_one_authority():
    """Four call sites answered 'where does a node's evidence live'."""

    import ast

    source = Path("chemsmart/agent/tool_runtime.py")
    tree = ast.parse(source.read_text(encoding="utf-8"))
    # Only path construction: `something / "nodes"`. The same string is
    # also an ordinary mapping key for a workflow's nodes, which is a
    # different question and not this one.
    literals = [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.BinOp)
        and isinstance(node.op, ast.Div)
        and isinstance(node.right, ast.Constant)
        and node.right.value == "nodes"
    ]
    assert not literals, (
        "tool_runtime still builds a node evidence path from the literal "
        f'"nodes" at lines {literals}; the layout has one author'
    )


def test_a_branch_says_what_was_asked_before_the_engine_answers(tmp_path):
    """The five questions, from the folder alone.

    Written before launch on purpose: a node killed mid-engine still
    carries what was asked of it, and a chemist reading the folder is not
    left with an empty directory and a terminal state word.
    """

    from chemsmart.agent.tool_runtime import _write_branch_request

    project = tmp_path / "projects" / "b3lyp.yaml"
    project.parent.mkdir(parents=True)
    project.write_text("gas:\n  functional: b3lyp\n", encoding="utf-8")

    branch = tmp_path / "cycle-0001" / "conformer-a-opt"
    branch.mkdir(parents=True)
    command = [
        "/env/bin/python",
        "-m",
        "chemsmart",
        "run",
        "--no-fake",
        "-p",
        str(project),
        "-f",
        "conformer-a.xyz",
        "orca",
        "opt",
    ]
    _write_branch_request(branch, command)

    # What CHEMSMART command actually ran -- the argv itself, not a
    # reconstruction of it.
    recorded = (branch / "command.txt").read_text(encoding="utf-8").strip()
    assert "chemsmart run --no-fake" in recorded
    assert str(project) in recorded
    assert recorded.endswith("orca opt")

    # What exact project configuration it used, copied so the branch does
    # not depend on that file still existing or still saying this.
    copied = (branch / "project.yaml").read_text(encoding="utf-8")
    assert "functional: b3lyp" in copied
    assert str(project) in copied.splitlines()[0]
    project.unlink()
    assert "functional: b3lyp" in (branch / "project.yaml").read_text(
        encoding="utf-8"
    )


def test_a_branch_records_an_unreadable_project_rather_than_failing(tmp_path):
    """The host does not fail a calculation over its own bookkeeping."""

    from chemsmart.agent.tool_runtime import _write_branch_request

    branch = tmp_path / "cycle-0001" / "n"
    branch.mkdir(parents=True)
    _write_branch_request(
        branch,
        [
            "python",
            "-m",
            "chemsmart",
            "run",
            "-p",
            str(tmp_path / "gone.yaml"),
        ],
    )
    assert (branch / "command.txt").is_file()
    assert not (branch / "project.yaml").exists()
    assert "gone.yaml" in (branch / "project.yaml.missing").read_text(
        encoding="utf-8"
    )
