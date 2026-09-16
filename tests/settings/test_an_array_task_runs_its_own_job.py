"""An array task runs the job it was submitted for.

The array authority answers one question in three places -- the runscript
filenames, the ``--array`` directive, and the shell that maps a task id to
a filename -- and on ``b1a2a13d`` the three disagreed: files were written
``0..N-1`` by ``enumerate``, the directive said ``1-N``, and the shell
added one again. Every task therefore ran a different job than the one it
was submitted for, and two tasks ran nothing at all.

Nothing caught it because nothing ever drove the consumer. The producer
was tested by reading the string it wrote; this module executes the shell
block the producer emits, for every task id the directive declares, and
asks which file the scheduler would actually run. That is the seam
production crosses.
"""

from __future__ import annotations

import os
import re
import subprocess
from types import SimpleNamespace

import pytest

from chemsmart.settings.server import Server
from chemsmart.settings.user import CHEMSMARTUserSettings

#: The directive is the scheduler's own word for which task ids exist.
_ARRAY_DIRECTIVE = re.compile(
    r"^#SBATCH --array=(\d+)-(\d+)(?:%(\d+))?\s*$", re.MULTILINE
)

_SERVER_NAME = "canned-slurm"

#: The array writer emits program specifics, which resolve an Executable
#: from a server file on disk. The witness therefore gives it a real one
#: rather than stubbing the production path it is trying to drive.
_SERVER_YAML = """\
SERVER:
    SCHEDULER: SLURM
    QUEUE_NAME: compute
    NUM_HOURS: 24
    MEM_GB: 52
    NUM_CORES: 6
    NUM_GPUS: 0
    NUM_THREADS: 6
    SUBMIT_COMMAND: sbatch
    SCRATCH_DIR: null
    USE_HOSTS: false
XTB:
    EXEFOLDER: null
    LOCAL_RUN: true
    SCRATCH: false
    CONDA_ENV: |
        # no environment needed for this witness
    MODULES: |
    SCRIPTS: |
    ENVARS: |
"""


@pytest.fixture()
def server(tmp_path, monkeypatch):
    """A real SLURM server profile, resolvable by name."""

    config = tmp_path / "server-config"
    config.mkdir()
    (config / f"{_SERVER_NAME}.yaml").write_text(_SERVER_YAML)
    monkeypatch.setattr(
        CHEMSMARTUserSettings,
        "user_server_dir",
        property(lambda self: str(config)),
    )
    return Server(
        _SERVER_NAME,
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=6,
        MEM_GB=52,
        NUM_HOURS=24,
        QUEUE_NAME="compute",
    )


def _write_array(server, folder, count, num_nodes=None):
    """Write an array of ``count`` distinguishable jobs into ``folder``."""

    jobs = [
        SimpleNamespace(label=f"mol{index}", folder=str(folder), PROGRAM="XTB")
        for index in range(count)
    ]
    # Each job carries an argument no other job carries, so a resolved
    # runscript names the job it belongs to rather than merely existing.
    cli_args = [["run", "--label", f"mol{index}"] for index in range(count)]
    submitter = server.get_submitter(jobs[0])
    submitter.write_array_job(jobs=jobs, num_nodes=num_nodes, cli_args=cli_args)
    return submitter, jobs


def _task_ids(script_text):
    match = _ARRAY_DIRECTIVE.search(script_text)
    assert match is not None, (
        "the array submit script declares no --array directive, so no task "
        f"id exists to drive:\n{script_text}"
    )
    first, last, throttle = match.groups()
    return (
        list(range(int(first), int(last) + 1)),
        None if throttle is None else int(throttle),
    )


def _resolve(script_text, task_id, cwd):
    """What file this task would run, by executing the emitted shell.

    ``python`` is shadowed by a shell function so the mapping is observed
    without running a job, and ``SLURM_SUBMIT_DIR`` is set because the
    script's own ``cd`` would otherwise land in ``$HOME``.
    """

    shim = 'python() { echo "RESOLVED:$1"; }\n'
    completed = subprocess.run(
        ["bash", "-c", shim + script_text],
        cwd=str(cwd),
        env={
            **os.environ,
            "SLURM_ARRAY_TASK_ID": str(task_id),
            "SLURM_SUBMIT_DIR": str(cwd),
        },
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, (
        f"task {task_id}: the dispatch shell failed: {completed.stderr}"
    )
    resolved = [
        line[len("RESOLVED:") :]
        for line in completed.stdout.splitlines()
        if line.startswith("RESOLVED:")
    ]
    assert len(resolved) == 1, (
        f"task {task_id}: the shell chose {len(resolved)} scripts, not one: "
        f"{resolved}"
    )
    return resolved[0]


def test_every_array_task_runs_exactly_its_own_job(server, tmp_path, monkeypatch):
    """The whole contract, driven end to end: task k runs job k, once."""

    folder = tmp_path / "jobs"
    folder.mkdir()
    # The submitter must write into the job's own folder, so the test runs
    # from somewhere else on purpose.
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)

    submitter, jobs = _write_array(server, folder, count=3)
    script_path = folder / submitter.array_submit_script
    assert script_path.is_file(), (
        "the array submit script was not written into the job folder; "
        f"the folder holds {sorted(p.name for p in folder.iterdir())} and "
        f"the process cwd holds {sorted(p.name for p in elsewhere.iterdir())}"
    )

    script_text = script_path.read_text()
    task_ids, _ = _task_ids(script_text)
    assert len(task_ids) == len(jobs), (
        f"{len(jobs)} jobs were written and the directive declares "
        f"{len(task_ids)} tasks"
    )

    # The whole mapping is collected before anything is asserted, so one
    # failure prints the entire task -> file -> job table. Stopping at the
    # first missing file hides the more dangerous half: a task that
    # resolves a file that *does* exist and belongs to another job.
    table = []
    for task_id in task_ids:
        resolved = _resolve(script_text, task_id, folder)
        target = folder / resolved
        owner = None
        if target.is_file():
            body = target.read_text()
            owners = [job.label for job in jobs if f"'{job.label}'" in body]
            owner = owners[0] if len(owners) == 1 else f"ambiguous:{owners}"
        table.append((task_id, resolved, target.is_file(), owner))

    report = "\n".join(
        f"  task {task_id} -> {resolved} "
        f"({'exists' if found else 'MISSING'}) runs {owner or 'NOTHING'}"
        for task_id, resolved, found, owner in table
    )
    on_disk = sorted(p.name for p in folder.glob("chemsmart_run_array_*.py"))

    missing = [row for row in table if not row[2]]
    assert not missing, (
        f"{len(missing)} of {len(table)} array tasks would run a file that "
        f"does not exist.\nfiles on disk: {on_disk}\n{report}"
    )

    ran = {task_id: owner for task_id, _, _, owner in table}
    assert sorted(ran.values()) == sorted(job.label for job in jobs), (
        "the tasks did not run each job exactly once; a task that runs "
        "another task's job is silent in every log.\n"
        f"files on disk: {on_disk}\n{report}"
    )


def test_the_concurrency_throttle_is_the_scheduler_s_own(server, tmp_path):
    """``num_nodes`` is the ``%N`` cap on simultaneously running tasks."""

    folder = tmp_path / "jobs"
    folder.mkdir()
    submitter, _ = _write_array(server, folder, count=6, num_nodes=4)
    script_text = (folder / submitter.array_submit_script).read_text()
    task_ids, throttle = _task_ids(script_text)
    assert len(task_ids) == 6
    assert throttle == 4


def test_an_unthrottled_array_declares_every_task(server, tmp_path):
    folder = tmp_path / "jobs"
    folder.mkdir()
    submitter, _ = _write_array(server, folder, count=2)
    script_text = (folder / submitter.array_submit_script).read_text()
    task_ids, throttle = _task_ids(script_text)
    assert len(task_ids) == 2
    assert throttle is None


def test_a_scheduler_without_an_array_directive_is_refused_by_name(
    tmp_path, monkeypatch
):
    """A scheduler failure must not wear a job's failure word.

    PBS and LSF write no array directive, so an array on them would be one
    job run N times under one id. That is refused where it is asked for,
    naming the scheduler and what wiring it would need.
    """

    from chemsmart.settings.submitters import PBSSubmitter, SLFSubmitter

    folder = tmp_path / "jobs"
    folder.mkdir()
    for submitter_class in (PBSSubmitter, SLFSubmitter):
        assert submitter_class.ARRAY_TASK_ID_VARIABLE is None, (
            f"{submitter_class.__name__} declares an array task-id variable "
            "while writing no array directive: a declaration without wiring "
            "is the defect this module exists to prevent"
        )
        job = SimpleNamespace(label="mol0", folder=str(folder), PROGRAM="XTB")
        submitter = submitter_class(
            name=submitter_class.NAME, job=job, server=None
        )
        with pytest.raises(ValueError, match="declares no array task-id"):
            submitter.write_array_job(jobs=[job], cli_args=[["run"]])


def test_an_array_submission_names_the_job_the_scheduler_created(
    server, tmp_path, monkeypatch
):
    """An array that names no job cannot be waited on, asked after, or
    parked on. The array path used to call ``Popen`` and discard the id
    ``sbatch`` had just printed, so it produced no receipt at all while
    the single-job path beside it produced one.
    """

    import subprocess as _subprocess

    from chemsmart.settings import server as server_module

    folder = tmp_path / "jobs"
    folder.mkdir()
    jobs = [
        SimpleNamespace(label=f"mol{index}", folder=str(folder), PROGRAM="XTB")
        for index in range(3)
    ]
    seen = {}

    def fake_run(argv, **kwargs):
        seen["argv"] = argv
        seen["cwd"] = kwargs.get("cwd")
        return _subprocess.CompletedProcess(
            argv, 0, stdout="Submitted batch job 4242\n", stderr=""
        )

    monkeypatch.setattr(server_module.subprocess, "run", fake_run)
    monkeypatch.setattr(
        server_module.Server, "_check_running_jobs", lambda self, job: None
    )

    receipt = server.submit_array_job(
        jobs=jobs,
        num_nodes=4,
        cli_args=[["run", "--label", f"mol{i}"] for i in range(3)],
    )

    assert receipt is not None, "an array submission returned no receipt"
    assert receipt.job_id == "4242"
    assert receipt.scheduler == "SLURM"
    # The command names the array script, and runs in the job's folder --
    # the same folder the scripts were written into.
    assert receipt.submit_command.endswith("chemsmart_sub_array_mol0.sh")
    assert seen["cwd"] == str(folder)
    assert receipt.submit_script == str(folder / "chemsmart_sub_array_mol0.sh")


def test_a_test_submission_writes_the_scripts_and_submits_nothing(
    server, tmp_path, monkeypatch
):
    from chemsmart.settings import server as server_module

    folder = tmp_path / "jobs"
    folder.mkdir()
    jobs = [
        SimpleNamespace(label=f"mol{index}", folder=str(folder), PROGRAM="XTB")
        for index in range(2)
    ]

    def refuse(*args, **kwargs):  # pragma: no cover - must never run
        raise AssertionError("a test submission reached the scheduler")

    monkeypatch.setattr(server_module.subprocess, "run", refuse)
    monkeypatch.setattr(
        server_module.Server, "_check_running_jobs", lambda self, job: None
    )

    assert (
        server.submit_array_job(
            jobs=jobs, test=True, cli_args=[["run"], ["run"]]
        )
        is None
    )
    assert (folder / "chemsmart_sub_array_mol0.sh").is_file()
    assert sorted(p.name for p in folder.glob("chemsmart_run_array_*.py")) == [
        "chemsmart_run_array_0.py",
        "chemsmart_run_array_1.py",
    ]
