"""Bootstrap conformance binds a state the supplied molecule can have.

The conformance probe fake-previews every declared program on the
workspace's own first geometry, and it bound charge 0, multiplicity 1 to
whatever that geometry was.  A neutral singlet is a state no
odd-electron molecule can have, and PySCF's preflight is the one that
says so (``pyscf.electrons.spin_parity``): on a workspace whose supplied
structure was the allyl radical, PySCF's conformance went red at
bootstrap, every PySCF job type became reference-only, and a correctly
planned single point returned the goal to the human at zero engine calls
(PySCF round, g3, 2026-09-12).  The probe now derives the neutral
multiplicity the electron count permits, from the input's own atoms.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.capabilities import load_program_capabilities
from chemsmart.agent.cli_schema import build_live_click_schema
from chemsmart.agent.live_session import _bootstrap_conformance

pytestmark = pytest.mark.capability("program_jobtype:pyscf:cpu:sp")

RADICAL = "4\nmethyl radical, 9 electrons\nC 0 0 0\nH 1.08 0 0\nH -0.54 0.935 0\nH -0.54 -0.935 0\n"
CLOSED = "2\nH2, 2 electrons\nH 0 0 0\nH 0 0 0.74\n"


def _artifact(path: Path) -> TrustedArtifactRefV1:
    return TrustedArtifactRefV1(
        artifact_id=path.stem,
        kind="geometry_xyz",
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )


@pytest.mark.parametrize(
    ("name", "text"), [("methyl.xyz", RADICAL), ("h2.xyz", CLOSED)]
)
def test_every_program_bootstraps_over_the_supplied_molecule(
    tmp_path, name, text
):
    xyz = tmp_path / name
    xyz.write_text(text, encoding="utf-8")
    registry = load_program_capabilities()
    receipts, records = _bootstrap_conformance(
        run_directory=tmp_path,
        input_artifact=_artifact(xyz),
        registry_sha256=registry.registry_sha256,
        live_schema=build_live_click_schema(),
    )
    by_program = {receipt.program: receipt for receipt in receipts}
    assert "pyscf" in by_program, [
        (r.get("program"), r.get("status"), r.get("error_class"))
        for r in records
    ]
    pyscf = by_program["pyscf"]
    assert pyscf.compiler_status == "passed"
    assert pyscf.preview_status == "passed"
    assert pyscf.preflight_status == "passed"
    assert "sp" in pyscf.covered_jobtypes
    # Coverage is per stage, and the response stage names its manifold
    # after the state the probe binds: a doublet workspace previews td on
    # the unrestricted manifold instead of reporting a gap for a singlet
    # it cannot have, so every CPU stage is covered on both molecules.
    assert ("cpu", "sp") in pyscf.effective_engine_job_pairs
    assert ("cpu", "opt") in pyscf.effective_engine_job_pairs
    assert ("cpu", "hess") in pyscf.effective_engine_job_pairs
    assert ("cpu", "td") in pyscf.effective_engine_job_pairs
