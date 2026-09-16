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
    submitter.write_array_job(
        jobs=jobs, num_nodes=num_nodes, cli_args=cli_args
    )
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
    assert (
        completed.returncode == 0
    ), f"task {task_id}: the dispatch shell failed: {completed.stderr}"
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


def test_every_array_task_runs_exactly_its_own_job(
    server, tmp_path, monkeypatch
):
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


def test_an_array_with_no_explicit_cap_takes_the_host_default(
    server, tmp_path
):
    """There is no unthrottled cohort.

    This asserted `throttle is None`, which was the old opt-in
    behaviour: a cohort submitted with no explicit argument ran with no
    cap at all. Physical concurrency is the host's and it is on by
    default (owner ruling, 2026-09-16), so the assertion moves with the
    contract rather than the contract moving with the assertion.
    """

    from chemsmart.settings.submitters import DEFAULT_MAX_CONCURRENT_TASKS

    folder = tmp_path / "jobs"
    folder.mkdir()
    submitter, _ = _write_array(server, folder, count=2)
    script_text = (folder / submitter.array_submit_script).read_text()
    task_ids, throttle = _task_ids(script_text)
    assert len(task_ids) == 2
    assert throttle == DEFAULT_MAX_CONCURRENT_TASKS


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


def test_every_registered_scheduler_can_write_its_own_directives(
    tmp_path, monkeypatch
):
    """A submitter that cannot write a header is a scheduler CHEMSMART
    advertises and cannot use.

    ``SLFSubmitter`` read ``self.server.num_nodes``, which ``Server`` does
    not define, so every LSF submission raised AttributeError -- on a
    class the registry hands out by name. Nothing drove it.
    """

    import io

    from chemsmart.settings import submitters as submitters_module
    from chemsmart.settings.submitters import Submitter

    # Keys a scheduler legitimately requires of the user's own settings.
    # Supplying them means anything still failing is a code defect rather
    # than an unconfigured host.
    if submitters_module.user_settings is not None:
        data = dict(submitters_module.user_settings.data)
        data.setdefault("RSCGRP", "small")
        data.setdefault("PROJECT", "probe-project")
        monkeypatch.setattr(
            submitters_module.user_settings, "data", data, raising=False
        )

    job = SimpleNamespace(label="probe", folder=str(tmp_path), PROGRAM="XTB")
    failures = {}
    for submitter_class in Submitter.subclasses():
        name = submitter_class.NAME
        if not name:
            continue
        profile = Server(
            f"canned-{name}",
            SCHEDULER=name,
            SUBMIT_COMMAND="true",
            NUM_CORES=6,
            MEM_GB=52,
            NUM_HOURS=24,
            NUM_GPUS=0,
            QUEUE_NAME="compute",
        )
        submitter = submitter_class(name=name, job=job, server=profile)
        buffer = io.StringIO()
        try:
            submitter._write_scheduler_options(buffer)
        except Exception as error:  # noqa: BLE001 - the point of the test
            failures[name] = f"{type(error).__name__}: {error}"
            continue
        written = buffer.getvalue()
        assert written.strip(), f"{name} wrote no directives at all"

    assert not failures, (
        "these registered schedulers cannot write their own directives: "
        f"{failures}"
    )


def test_an_array_carries_the_operator_s_own_scheduler_directives(
    tmp_path, monkeypatch
):
    """Reservation, QoS and every site directive survive the array path.

    The single-job writer calls _write_extra_scheduler_directives; the
    array writer did not. Switching the Agent to an array would therefore
    have dropped exactly the directives the plan lists as a *fairness*
    criterion -- and on CUHK the group's reservation is the only route to
    the nodes a large ceiling needs.
    """

    config = tmp_path / "server-config"
    config.mkdir()
    (config / f"{_SERVER_NAME}.yaml").write_text(_SERVER_YAML)
    monkeypatch.setattr(
        CHEMSMARTUserSettings,
        "user_server_dir",
        property(lambda self: str(config)),
    )
    server = Server(
        _SERVER_NAME,
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=6,
        MEM_GB=52,
        NUM_HOURS=24,
        QUEUE_NAME="compute",
        EXTRA_SCHEDULER_DIRECTIVES=(
            "#SBATCH --reservation=xlzhang_1\n#SBATCH --qos=high\n"
        ),
    )
    folder = tmp_path / "jobs"
    folder.mkdir()
    jobs = [
        SimpleNamespace(label=f"mol{i}", folder=str(folder), PROGRAM="XTB")
        for i in range(2)
    ]
    submitter = server.get_submitter(jobs[0])

    import io

    single = io.StringIO()
    submitter._write_scheduler_options(single)
    array = io.StringIO()
    submitter.jobs = jobs
    submitter._write_array_scheduler_options(array, None)

    for directive in ("--reservation=xlzhang_1", "--qos=high"):
        assert directive in single.getvalue(), "fixture is wrong"
        assert directive in array.getvalue(), (
            f"the array path drops {directive!r}, which the single-job "
            "path carries"
        )


def test_each_array_element_runs_in_its_own_job_s_directory(server, tmp_path):
    """An array of N molecules in N directories must not run N times in
    the first molecule's directory.

    The single-job runscript passes execution_cwd; the array runscript did
    not, so RunScript emitted `pass` instead of os.chdir and every element
    inherited $SLURM_SUBMIT_DIR -- jobs[0].folder. N sets of outputs then
    collide on program-default filenames, in one directory, silently.
    """

    submit_folder = tmp_path / "jobs"
    submit_folder.mkdir()
    folders = []
    jobs = []
    for index in range(3):
        own = tmp_path / f"mol{index}-dir"
        own.mkdir()
        folders.append(own)
        jobs.append(
            SimpleNamespace(
                label=f"mol{index}",
                folder=str(submit_folder) if index == 0 else str(own),
                PROGRAM="XTB",
            )
        )
    # The array is submitted from the first job's folder, as the submitter
    # already assumes; each element still belongs somewhere of its own.
    jobs[0].submission_execution_cwd = str(folders[0])
    for index in range(1, 3):
        jobs[index].submission_execution_cwd = str(folders[index])

    submitter = server.get_submitter(jobs[0])
    submitter.write_array_job(
        jobs=jobs,
        num_nodes=None,
        cli_args=[["run", "--label", f"mol{i}"] for i in range(3)],
    )

    task_ids = submitter.array_task_ids(len(jobs))
    for index, task_id in enumerate(task_ids):
        body = (
            submit_folder / submitter.array_run_script(task_id)
        ).read_text()
        assert str(folders[index]) in body, (
            f"element {task_id} does not enter its own job's directory; "
            f"it would run wherever the array was submitted from:\n{body}"
        )
