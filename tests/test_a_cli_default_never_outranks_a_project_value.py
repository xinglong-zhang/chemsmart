"""A CLI option the user never typed must not overwrite project YAML.

FUNDAMENTAL: project YAML is where the computational-chemistry rationale
lives, and the CLI is how a human or an agent drives it. Every ``jobtype``
command therefore reads its project's settings first and applies an option
only when the caller actually supplied one -- the pattern its own comment
states as "only update value if user explicitly specifies a value for the
attribute to preserve project defaults".

That pattern is expressed as ``if option is not None:`` over a Click option
whose default is ``None``. A non-``None`` default silently defeats it: the
guard is always true and the flag nobody passed wins.

po3-r17 (2026-09-11) lost a whole 12-hour window to exactly that. The
session declared ``recalc_hess: 999`` in project YAML so a saddle search
could not become an unbounded chain of numerical second derivatives; the
written ORCA input carried ``5``; the preview validator correctly reported
``expected 999, observed 5`` and the workflow was refused. ``recalc_hess``
was the only option in the whole CLI shaped that way, and one line fixed
it -- so this test pins the invariant rather than the case, because the
next such option is what it exists to catch.
"""

import ast
import pathlib
import re

import pytest

CLI_ROOT = pathlib.Path(__file__).resolve().parents[1] / "chemsmart" / "cli"


def _signature_defaults(function: ast.FunctionDef) -> dict[str, object]:
    """Literal defaults in a command callback's signature, by name."""

    args = function.args
    defaults: dict[str, object] = {}
    pairs = list(zip(args.kwonlyargs, args.kw_defaults))
    if args.defaults:
        offset = len(args.args) - len(args.defaults)
        pairs += list(zip(args.args[offset:], args.defaults))
    for argument, node in pairs:
        if node is None:
            continue
        try:
            defaults[argument.arg] = ast.literal_eval(node)
        except (ValueError, SyntaxError):
            continue
    return defaults


def _offenders() -> list[str]:
    """Every option whose own default defeats its own guard."""

    found: list[str] = []
    for path in sorted(CLI_ROOT.rglob("*.py")):
        source = path.read_text(encoding="utf-8")
        for function in [
            node
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.FunctionDef)
        ]:
            body = ast.get_source_segment(source, function) or ""
            for name, value in _signature_defaults(function).items():
                if value is None or value is False:
                    continue
                assigns = re.search(
                    rf"\.{re.escape(name)}\s*=\s*{re.escape(name)}\b", body
                )
                if not assigns:
                    continue
                guarded = re.search(
                    rf"if\s+{re.escape(name)}\s+is\s+not\s+None\s*:", body
                ) or re.search(rf"if\s+{re.escape(name)}\s*:", body)
                if guarded:
                    found.append(
                        f"{path.relative_to(CLI_ROOT.parents[1])}"
                        f"::{function.name} option {name!r} defaults to "
                        f"{value!r} under a 'did the user type it' guard"
                    )
    return found


@pytest.mark.capability("setting:orca:recalc_hess")
def test_no_cli_option_default_defeats_its_own_guard():
    """The guard means "the caller typed this", so the default is None."""

    offenders = _offenders()
    assert not offenders, (
        "a CLI option would overwrite the project's own value with a "
        "default nobody passed:\n  " + "\n  ".join(offenders)
    )


@pytest.mark.capability("setting:orca:recalc_hess")
def test_the_scan_finds_a_planted_offender():
    """The scan is falsifiable: it must catch the shape it forbids.

    Without this, a scan that silently matched nothing would pass over a
    repository full of offenders -- the class of defect this laboratory
    has paid for before ("a cheap proxy is not a measurement").
    """

    planted = ast.parse(
        "def cmd(ctx, recalc_hess=5):\n"
        "    settings = ctx.obj['s']\n"
        "    if recalc_hess is not None:\n"
        "        settings.recalc_hess = recalc_hess\n"
    )
    function = planted.body[0]
    source = (
        "def cmd(ctx, recalc_hess=5):\n"
        "    settings = ctx.obj['s']\n"
        "    if recalc_hess is not None:\n"
        "        settings.recalc_hess = recalc_hess\n"
    )
    defaults = _signature_defaults(function)
    assert defaults.get("recalc_hess") == 5
    body = ast.get_source_segment(source, function) or ""
    assert re.search(r"if\s+recalc_hess\s+is\s+not\s+None\s*:", body)
    assert re.search(r"\.recalc_hess\s*=\s*recalc_hess\b", body)


@pytest.mark.capability("program_jobtype:orca:cpu:ts")
def test_a_project_recalc_hess_reaches_the_job(tmp_path):
    """The live loss, driven through the real command.

    Both halves matter: the project's declared value must survive, and a
    project that declares nothing must still get the settings class's own
    5, so restoring project authority changes no existing behaviour.
    """

    from unittest.mock import MagicMock, patch

    import yaml
    from click.testing import CliRunner

    from chemsmart.cli.orca.orca import orca as orca_cli

    molecule = tmp_path / "h2.xyz"
    molecule.write_text("2\nh2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n")
    projects = tmp_path / "project_yaml"
    projects.mkdir()

    def run(name, declared):
        body = {"ts": {"functional": "b3lyp", "basis": "def2-svp"}}
        body["ts"].update(declared)
        (projects / f"{name}.yaml").write_text(
            yaml.safe_dump(body), encoding="utf-8"
        )
        runner = CliRunner()
        with patch("chemsmart.jobs.orca.ts.ORCATSJob") as job:
            job.return_value = MagicMock()
            result = runner.invoke(
                orca_cli,
                [
                    "-p",
                    str(projects / name),
                    "-f",
                    str(molecule),
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "ts",
                ],
                obj={},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        return job.call_args.kwargs["settings"]

    declared = run("declares", {"recalc_hess": 999, "geom_maxiter": 300})
    assert declared.recalc_hess == 999, (
        "the project declared recalc_hess 999 and an option nobody typed "
        f"overwrote it with {declared.recalc_hess}"
    )
    assert declared.geom_maxiter == 300

    silent = run("silent", {})
    assert silent.recalc_hess == 5, (
        "a project that declares nothing must still receive the settings "
        f"class's own default, not {silent.recalc_hess}"
    )
