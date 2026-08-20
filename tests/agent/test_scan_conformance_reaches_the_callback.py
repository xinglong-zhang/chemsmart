"""A declared family must reach the same Click callback a user reaches.

Conformance drives real `--fake` invocations, so a family whose callback needs
an argument the probe never supplies reports `failed` and never becomes
previewable. `scan` and `modred` are exactly that case: both refuse to build a
job without the coordinate they drive or hold.

The probe renders that coordinate through the *same* translation the production
path uses rather than a hand-written argv, so this test also covers the case
where the two drift apart. `opt` is carried as a control: it needs no
coordinate, and if it ever fails here the fixture is broken rather than the
family.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from chemsmart.agent.bootstrap import _CONFORMANCE_COORDINATES
from chemsmart.agent.cli_schema import build_live_click_schema
from chemsmart.agent.commands import native_coordinate_options
from chemsmart.agent.projects import project_document, render_project_yaml
from chemsmart.cli.main import entry_point

_SERVER = Path(
    "chemsmart/settings/templates/.chemsmart/server/SLURM.yaml"
).resolve()


@pytest.fixture(scope="module")
def fixtures(tmp_path_factory):
    root = tmp_path_factory.mktemp("conformance")
    geometry = root / "probe.xyz"
    geometry.write_text(
        "3\nprobe\nO 0.0 0.0 0.0\nH 0.0 0.0 0.96\nH 0.93 0.0 -0.24\n",
        encoding="utf-8",
    )
    project = root / "probe_project.yaml"
    project.write_text(
        render_project_yaml(
            project_document(
                program="orca",
                sections={"gas": {"functional": "b3lyp", "basis": "def2-svp"}},
            )
        ).rendered_yaml,
        encoding="utf-8",
    )
    return geometry, project


@pytest.mark.parametrize("jobtype", ("opt", "scan", "modred"))
def test_the_probe_argv_reaches_the_callback(jobtype, fixtures):
    geometry, project = fixtures
    argv = [
        "run",
        "--server",
        str(_SERVER),
        "--fake",
        "--no-scratch",
        "orca",
        "--project",
        str(project),
        "--filename",
        str(geometry),
        "--charge",
        "0",
        "--multiplicity",
        "1",
        jobtype,
    ]
    for name, value in sorted(
        native_coordinate_options(
            "orca", _CONFORMANCE_COORDINATES.get(jobtype)
        ).items()
    ):
        argv.extend((f"--{name.replace('_', '-')}", str(value)))

    with CliRunner().isolated_filesystem():
        result = CliRunner().invoke(entry_point, argv)

    assert (
        result.exit_code == 0
    ), f"{jobtype} did not reach a clean fake preview: {result.exception}"


@pytest.mark.parametrize("jobtype", ("scan", "modred"))
def test_every_rendered_option_exists_on_the_live_command(jobtype):
    """The probe must not invent a flag the program does not accept."""

    schema = build_live_click_schema()
    for program in ("orca", "gaussian"):
        command = schema.command(("run", program, jobtype))
        assert command is not None, f"{program} {jobtype} is not in the schema"
        for name in native_coordinate_options(
            program, _CONFORMANCE_COORDINATES[jobtype]
        ):
            assert (
                command.option(name) is not None
            ), f"{program} {jobtype} has no option named {name!r}"


def test_the_probe_coordinate_is_legal_for_the_fixture_geometry():
    """Atoms 1 and 2 must exist in the probe molecule, or nothing is proven."""

    for spec in _CONFORMANCE_COORDINATES.values():
        entries = list(spec.get("constrained") or ())
        if spec.get("scan"):
            entries.append(spec["scan"])
        for entry in entries:
            assert max(entry["atoms"]) <= 3, (
                "the conformance geometry has three atoms; a coordinate "
                "naming a higher index would never reach the callback"
            )
