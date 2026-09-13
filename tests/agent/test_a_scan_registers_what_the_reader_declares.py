"""A workspace scan registers every job type the reader declares.

The PySCF scan admitted results by a hand-written set, ``{"sp", "opt",
"hess"}``, beside a reader that declares its job types by selector
table. When the response stage became executable the reader declared
``td`` and the scan did not follow, so the first archived spectra sat in
a workspace and registered nothing -- every other organ for them was in
place, and the witness bank found it on that tree (PySCF round 2,
2026-09-13). One authority per question: the scan now asks the reader.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from chemsmart.agent.live_session import discover_registered_result_artifacts
from chemsmart.analysis.result_readers import reader_for

FIXTURES = (
    Path(__file__).resolve().parents[1] / "data" / "PySCFTests" / "outputs"
)

#: One green archived run per job type the reader declares.
CASES = {
    "sp": "water_ccsdt_sp",
    "opt": "formaldehyde_s1_opt",
    "hess": "water_hess",
    "td": "water_td_singlet",
}


@pytest.mark.capability("gate:resolver.one_answer_per_question")
def test_the_scan_registers_every_job_type_the_reader_declares(tmp_path):
    reader = reader_for("pyscf")
    declared = {name for name, _selectors in reader.jobtype_selectors}
    assert declared == set(CASES), (
        "a job type the reader declares has no archived run here; add one "
        "so the scan is checked against it"
    )
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    for case in CASES.values():
        for path in (FIXTURES / case).iterdir():
            if path.suffix == ".h5" or path.name.endswith(
                (".receipt.json", ".input.json", ".environment.json")
            ):
                shutil.copy2(path, workspace / path.name)
    registered = discover_registered_result_artifacts(workspace)
    found = {}
    for artifact in registered:
        output = reader.open_output(Path(artifact.path))
        found[str(output.jobtype)] = artifact.artifact_id
    assert set(found) == declared, (
        f"the scan registered {sorted(found)} of the declared "
        f"{sorted(declared)}"
    )
