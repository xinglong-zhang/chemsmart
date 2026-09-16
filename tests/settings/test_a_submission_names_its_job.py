"""A submission returns the job it created, and the scheduler's later
answers about that job parse from canned real outputs (the aiida method).

Recorded on this host (slurm 26.05.2, accounting storage disabled): a
finished job stays visible to scontrol and squeue only for MinJobAge, so
"unknown" is a first-class answer the parsers must give rather than an
error they must hide.
"""

from __future__ import annotations

import subprocess

import pytest

from chemsmart.settings.probe import ProbeUnitError
from chemsmart.settings.probe.scheduler_job import (
    TERMINAL_JOB_STATES,
    parse_elapsed_seconds,
    parse_scontrol_job,
    parse_submission,
)
from chemsmart.settings.server import Server
from chemsmart.settings.submitters import SubmissionReceiptV1

_SCONTROL = """\
JobId=191 JobName=probe-wake
   UserId=chemsmart(1000) GroupId=chemsmart(1000) MCS_label=N/A
   Priority=4294901569 Nice=0 Account=(null) QOS=(null)
   JobState=COMPLETED Reason=None Dependency=(null)
   Requeue=1 Restarts=0 BatchFlag=1 Reboot=0 ExitCode=0:0
   RunTime=00:00:02 TimeLimit=00:01:00 TimeMin=N/A
   SubmitTime=2026-09-02T12:39:50 EligibleTime=2026-09-02T12:39:50
   StartTime=2026-09-02T12:39:50 EndTime=2026-09-02T12:39:52 Deadline=N/A
   Partition=compute AllocNode:Sid=localhost.localdomain:4140442
   NodeList=chemsmart-hpc
   Command=(null)
   SubmitLine=/opt/slurm/current/bin/sbatch -p compute -t 1 -J probe-wake --wrap hostname; sleep 2
   WorkDir=/tmp/scratch
"""
_SQUEUE = "191|COMPLETED|0:02|2026-09-02T12:39:50|2026-09-02T12:39:50|2026-09-02T12:39:52\n"
_UNKNOWN_STDERR = "slurm_load_jobs error: Invalid job id specified\n"


def test_sbatch_qsub_and_parsable_lines_name_the_job():
    assert parse_submission(0, "Submitted batch job 191\n", "") == "191"
    assert parse_submission(0, "191\n", "") == "191"
    assert parse_submission(0, "191;cluster\n", "") == "191"
    assert parse_submission(0, "4242.pbsserver\n", "") == "4242.pbsserver"


def test_a_failed_or_silent_submission_is_refused_with_its_own_words():
    with pytest.raises(ProbeUnitError, match="submission exited 1"):
        parse_submission(1, "", "sbatch: error: invalid partition")
    with pytest.raises(ProbeUnitError, match="printed no job id"):
        parse_submission(0, "Job submitted.\n", "")


def test_scontrol_reports_a_finished_job_inside_the_min_job_age_window():
    state = parse_scontrol_job(0, _SCONTROL, "", job_id="191")
    assert state.known
    assert state.state == "COMPLETED"
    assert state.terminal
    assert state.exit_code == "0:0"
    assert state.run_seconds == 2
    assert (state.submit_time, state.end_time) == (
        "2026-09-02T12:39:50",
        "2026-09-02T12:39:52",
    )


def test_a_forgotten_job_is_unknown_not_an_error():
    state = parse_scontrol_job(1, "", _UNKNOWN_STDERR, job_id="999999")
    assert not state.known
    assert not state.terminal


def test_elapsed_fields_in_every_slurm_shortening():
    assert parse_elapsed_seconds("0:02") == 2
    assert parse_elapsed_seconds("12:34") == 12 * 60 + 34
    assert parse_elapsed_seconds("1:02:03") == 3723
    assert parse_elapsed_seconds("2-01:00:00") == 2 * 86400 + 3600
    assert parse_elapsed_seconds("N/A") is None


def test_cancelled_by_a_user_is_still_terminal():
    text = _SCONTROL.replace(
        "JobState=COMPLETED", "JobState=CANCELLED by 1000"
    )
    assert parse_scontrol_job(0, text, "", job_id="191").terminal
    assert "CANCELLED" in TERMINAL_JOB_STATES


class _Job:
    label = "agent-run"
    PROGRAM = None

    def __init__(self, folder):
        self.folder = str(folder)

    def is_complete(self):
        return False


def test_server_submit_returns_the_receipt_with_the_job_id(
    tmp_path, monkeypatch
):
    server = Server(
        "canned",
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=1,
        MEM_GB=1,
        NUM_HOURS=1,
        QUEUE_NAME="compute",
    )
    seen = {}

    def fake_run(argv, **kwargs):
        seen["argv"] = argv
        seen["cwd"] = kwargs.get("cwd")
        return subprocess.CompletedProcess(
            argv, 0, stdout="Submitted batch job 191\n", stderr=""
        )

    monkeypatch.setattr("chemsmart.settings.server.subprocess.run", fake_run)
    receipt = server._submit_job(_Job(tmp_path))
    assert isinstance(receipt, SubmissionReceiptV1)
    assert receipt.job_id == "191"
    assert receipt.scheduler == "SLURM"
    assert seen["argv"] == ["sbatch", "chemsmart_sub_agent-run.sh"]
    assert seen["cwd"] == str(tmp_path)
    assert receipt.submit_script.endswith("chemsmart_sub_agent-run.sh")
    assert receipt.submitted_at.endswith("+00:00")
