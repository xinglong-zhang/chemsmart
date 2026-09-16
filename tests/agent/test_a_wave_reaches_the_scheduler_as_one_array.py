"""A wave is submitted as one throttled array and one dependent wake.

Two builders existed and nothing called them: `dispatch_run_to_scheduler`
always built the single-job script, so a cohort of seven would have been
submitted as one job running seven calculations in series. The array
machinery was reachable from its own tests and from nothing else, which
is how its off-by-two survived.

What the dispatcher owes the cohort, and what is asserted here over the
real submitter and the real scheduler command:

- the manifest is written before anything is submitted, because the
  index-to-node mapping is what an element resolves itself through and a
  running element cannot wait for its author;
- one array job carries every member, throttled by the host's own
  concurrency number, never by the Agent's scientific width;
- no element wakes the goal -- N tails would wake the model N times and a
  wave is one reasoning point;
- the wake is its own job, gated `afterany` on the array, asking for one
  core rather than inheriting the cohort's allocation; and
- the receipt names both jobs, because a goal parked on the array alone
  could not say what would re-enter it.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from chemsmart.agent.cohort import COHORT_MANIFEST_FILE, read_cohort_manifest
from chemsmart.agent.dispatch import (
    DISPATCH_RECEIPT_FILE,
    dispatch_run_to_scheduler,
)
from chemsmart.settings.server import Server

_NODES = ("conformer-a-opt", "conformer-b-opt", "conformer-c-opt")
#: What a real bundle declares. The manifest must carry *this*, not a
#: hash of the file: admission compares the bundle's own digest.
_BUNDLE_DIGEST = "7" * 64


def _server():
    return Server(
        "canned-slurm",
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=8,
        MEM_GB=60,
        NUM_HOURS=24,
        QUEUE_NAME="compute",
    )


@pytest.fixture
def submitted(monkeypatch):
    """Every sbatch this dispatch makes, with the script it submitted."""

    seen: list[dict] = []

    def fake_run(argv, **kwargs):
        script = Path(kwargs.get("cwd") or ".") / Path(argv[-1]).name
        seen.append(
            {
                "argv": list(argv),
                "cwd": kwargs.get("cwd"),
                "script": script.read_text(encoding="utf-8"),
                "path": script,
            }
        )
        return subprocess.CompletedProcess(
            argv, 0, stdout=f"Submitted batch job 4{len(seen)}2\n", stderr=""
        )

    monkeypatch.setattr(subprocess, "run", fake_run)
    return seen


def _dispatch(tmp_path, **kwargs):
    run_directory = tmp_path / "run"
    bundle = tmp_path / "bundle.json"
    bundle.write_text(
        json.dumps({"bundle_sha256": _BUNDLE_DIGEST}), encoding="utf-8"
    )
    return dispatch_run_to_scheduler(
        approval_file=bundle,
        workspace=tmp_path / "ws",
        run_directory=run_directory,
        goal_id="g1",
        cycle=1,
        server="canned-slurm",
        python="/opt/env/bin/python",
        **kwargs,
    )


@pytest.fixture(autouse=True)
def _profile(monkeypatch):
    monkeypatch.setattr(
        Server, "from_servername", classmethod(lambda cls, name: _server())
    )


def test_a_cohort_is_one_array_job_and_one_wake_job(tmp_path, submitted):
    receipt = _dispatch(tmp_path, cohort_node_ids=_NODES)

    assert len(submitted) == 2, (
        "a wave must reach the scheduler as an array and a dependent "
        f"wake; {len(submitted)} job(s) were submitted"
    )
    array, wake = submitted

    array_directives = [
        line
        for line in array["script"].splitlines()
        if line.startswith("#SBATCH")
    ]
    assert any(
        line.startswith("#SBATCH --array=0-2%") for line in array_directives
    ), array_directives
    assert "--cohort-element" in array["script"]
    assert "agent wake" not in array["script"], (
        "an element woke the goal: N elements would wake the model N "
        "times and a wave is one reasoning point"
    )

    wake_directives = [
        line
        for line in wake["script"].splitlines()
        if line.startswith("#SBATCH")
    ]
    assert "#SBATCH --dependency=afterany:412" in wake_directives
    assert "#SBATCH --kill-on-invalid-dep=yes" in wake_directives
    assert "agent wake" in wake["script"]
    assert "--cohort-element" not in wake["script"]

    assert receipt.job_id == "412"
    assert receipt.wake_job_id == "422"
    assert tuple(receipt.cohort_node_ids) == _NODES


def test_the_manifest_is_written_before_anything_is_submitted(
    tmp_path, submitted
):
    """An element resolves itself through the manifest as it starts."""

    manifests: list[bool] = []
    run_directory = tmp_path / "run"

    import subprocess as _sp

    original = _sp.run

    def watching(argv, **kwargs):
        manifests.append((run_directory / COHORT_MANIFEST_FILE).is_file())
        return original(argv, **kwargs)

    _sp.run = watching
    try:
        _dispatch(tmp_path, cohort_node_ids=_NODES)
    finally:
        _sp.run = original

    assert manifests and all(manifests), (
        "the array was submitted before its manifest existed, so an "
        "element that starts promptly cannot say which node it is"
    )
    manifest = read_cohort_manifest(run_directory)
    assert manifest is not None
    assert tuple(manifest.node_ids) == _NODES
    assert manifest.node_for_element(1) == "conformer-b-opt"


def test_the_host_throttles_the_wave_the_agent_sized(tmp_path, submitted):
    """Scientific width is the Agent's; concurrency is the host's."""

    from chemsmart.settings.submitters import DEFAULT_MAX_CONCURRENT_TASKS

    wide = tuple(f"node-{index}" for index in range(7))
    _dispatch(tmp_path, cohort_node_ids=wide)

    array = submitted[0]["script"]
    directive = next(
        line
        for line in array.splitlines()
        if line.startswith("#SBATCH --array=")
    )
    assert directive == (
        f"#SBATCH --array=0-6%{DEFAULT_MAX_CONCURRENT_TASKS}"
    ), (
        "seven independent calculations are one cohort and one turn; the "
        "host renders the concurrency and the Agent never splits the wave"
    )


def test_a_run_with_no_cohort_still_dispatches_as_one_job(tmp_path, submitted):
    """The single-node path is unchanged: one job, its own tail waking."""

    receipt = _dispatch(tmp_path)
    assert len(submitted) == 1
    assert "--array=" not in submitted[0]["script"]
    assert "agent wake" in submitted[0]["script"]
    assert receipt.wake_job_id == ""
    assert receipt.cohort_node_ids == ()


def test_the_receipt_on_disk_names_both_jobs(tmp_path, submitted):
    _dispatch(tmp_path, cohort_node_ids=_NODES)
    record = json.loads(
        (tmp_path / "run" / DISPATCH_RECEIPT_FILE).read_text(encoding="utf-8")
    )
    assert record["job_id"] == "412"
    assert record["wake_job_id"] == "422"
    assert record["cohort_node_ids"] == list(_NODES)


def test_the_element_redirect_survives_a_path_with_a_space(tmp_path):
    """Every other path in the script is quoted; this one was not.

    `bash` redirects to the first word: with a workspace under "My
    Drive", `> /My Drive/ws/run/execution-result.${i}.json` redirects to
    `/My` and hands the rest to `chemsmart agent run` as stray
    positional arguments. Every element fails before doing anything, and
    the wave then cannot satisfy its own barrier.
    """

    from types import SimpleNamespace

    from chemsmart.agent.dispatch import build_cohort_dispatch_script

    spaced = tmp_path / "My Drive" / "run dir"
    spaced.mkdir(parents=True)
    submitter = _server().get_submitter(
        SimpleNamespace(label="goal-g1-cycle-1", PROGRAM=None, folder=".")
    )
    script = build_cohort_dispatch_script(
        submitter=submitter,
        python="/env/bin/python",
        approval_file=tmp_path / "bundle.json",
        workspace=tmp_path / "ws",
        run_directory=spaced,
        cohort_size=2,
    )
    line = next(item for item in script.splitlines() if "agent run" in item)
    redirect = line.split(">", 1)[1].strip()
    assert redirect.startswith("'") or redirect.startswith(
        '"'
    ), f"the redirect target is unquoted: {redirect}"
    # And the variable must still expand: quoting the whole thing with
    # single quotes would redirect to a literal '${SLURM_ARRAY_TASK_ID}'.
    assert "${SLURM_ARRAY_TASK_ID}" in redirect
    assert "My Drive" in redirect


def test_the_manifest_names_the_digest_admission_compares(tmp_path, submitted):
    """Composition, not spelling: dispatch, then admit through it.

    The manifest was bound to a hash of the approval *file* while
    `authorise_cohort_element` compares the bundle's own declared digest
    -- `canonical_sha256` over its content -- so the two never agree.
    Measured live on CUHK: the one element of a four-member wave that
    got as far as running was refused against its own approval, and the
    wave was blocked rather than wrong, which is the quiet kind.
    """

    from chemsmart.agent.cohort import authorise_cohort_element

    _dispatch(tmp_path, cohort_node_ids=_NODES)
    run_directory = tmp_path / "run"

    for element, expected in enumerate(_NODES):
        assert (
            authorise_cohort_element(
                run_directory,
                element=element,
                bundle_sha256=_BUNDLE_DIGEST,
            )
            == expected
        ), (
            "the dispatcher wrote a manifest this approval's own "
            "elements cannot be admitted through"
        )


def test_a_bundle_that_names_no_digest_is_refused_before_sbatch(tmp_path, submitted):
    """A cohort that cannot be bound is not submitted and then discovered."""

    from chemsmart.agent._contracts import ContractError

    (tmp_path / "bundle.json").write_text("{}", encoding="utf-8")
    run_directory = tmp_path / "run"
    with pytest.raises(ContractError, match="bundle_sha256"):
        dispatch_run_to_scheduler(
            approval_file=tmp_path / "bundle.json",
            workspace=tmp_path / "ws",
            run_directory=run_directory,
            goal_id="g1",
            cycle=1,
            server="canned-slurm",
            python="/opt/env/bin/python",
            cohort_node_ids=_NODES,
        )
    assert submitted == [], "the scheduler was called before the check"
