"""Gaussian reaction-path stages retain their program-relative CLI form."""

from pathlib import Path
from types import SimpleNamespace

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.live_session import _conformance_jobtypes
from chemsmart.agent.program_verifiers import _validate_native_input
from chemsmart.io.gaussian.route import GaussianRoute


def _write_geometry(path: Path) -> TrustedArtifactRefV1:
    path.write_text(
        "2\nH2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n",
        encoding="utf-8",
    )
    return TrustedArtifactRefV1(
        artifact_id="reaction-path-geometry",
        kind="geometry_xyz",
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )


def _write_input(path: Path, route: str) -> None:
    path.write_text(
        f"{route}\n\nReaction path preview\n\n"
        "0 1\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n\n",
        encoding="utf-8",
    )


def test_gaussian_declares_canonical_ts_and_irc_preview_stages():
    assert _conformance_jobtypes("gaussian", "cpu") == (
        "irc",
        "link",
        "opt",
        "sp",
        "td",
        "ts",
    )


def test_irc_maxpoints_is_not_misread_as_transition_state():
    forward = GaussianRoute("# B3LYP 6-31G* irc(calcfc,forward,maxpoints=32)")
    reverse = GaussianRoute("# B3LYP 6-31G* irc(calcfc,reverse,maxpoints=32)")
    assert forward.jobtype == "ircf"
    assert reverse.jobtype == "ircr"


def test_one_irc_node_validates_as_forward_and_reverse_bundle(tmp_path):
    geometry = _write_geometry(tmp_path / "input.xyz")
    forward = tmp_path / "path_ircf.com"
    reverse = tmp_path / "path_ircr.com"
    _write_input(
        forward,
        "# B3LYP 6-31G* irc(calcfc,recalc=6,forward,maxpoints=32)",
    )
    _write_input(
        reverse,
        "# B3LYP 6-31G* irc(calcfc,recalc=6,reverse,maxpoints=32)",
    )
    expectation = SimpleNamespace(
        program="gaussian",
        jobtype="irc",
        settings=(
            ("basis", "6-31G*"),
            ("freq", False),
            ("functional", "B3LYP"),
            ("jobtype", "irc"),
        ),
        charge=0,
        multiplicity=1,
        input_artifact=geometry,
    )

    assert _validate_native_input(expectation, (forward, reverse)) == []
