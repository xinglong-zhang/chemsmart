"""`chemsmart wizard` is a new user's first command; it must be honest.

Detection and probing are monkeypatched to canned facts (the real probes
have their own live and canned tests); these tests pin the flow: announced
evidence, queue-derived defaults, the scheduler-named file, the honest
local downgrade for unqualified schedulers, and --yes non-interactivity.
"""

from __future__ import annotations

import yaml
from click.testing import CliRunner

import chemsmart.cli.wizard as wizard_module
from chemsmart.cli.wizard import wizard
from chemsmart.settings.probe import (
    DetectionV1,
    HostFactsV1,
    QueueFactsV1,
    SchedulerFactsV1,
)

_HOST = HostFactsV1(
    hostname="ax41",
    os_name="Linux",
    os_release="5.15",
    machine="x86_64",
    cpu_count=12,
    mem_kb=65761448,
    scratch_candidates=("/tmp",),
    has_module_command=True,
)
_QUEUE = QueueFactsV1(
    name="compute",
    is_default=True,
    available=True,
    max_time_seconds=172800,
    cores_per_node=6,
    mem_kb_per_node=61440000,
    node_count=1,
)
_SLURM = SchedulerFactsV1(
    scheduler="SLURM",
    version="slurm 26.05.2",
    submit_path="/opt/slurm/current/bin/sbatch",
    queues=(_QUEUE,),
)


def _wire(monkeypatch, tmp_path, *, detection, scheduler):
    monkeypatch.setattr(wizard_module, "gather_host_facts", lambda: _HOST)
    monkeypatch.setattr(wizard_module, "detect_scheduler", lambda: detection)
    monkeypatch.setattr(
        wizard_module, "_probe_scheduler", lambda _d: scheduler
    )
    monkeypatch.setattr(
        wizard_module, "_user_server_dir", lambda: tmp_path / "server"
    )
    monkeypatch.setattr(
        wizard_module,
        "run_verification",
        lambda _choices: (),
    )


def test_yes_mode_writes_the_scheduler_named_file(tmp_path, monkeypatch):
    _wire(
        monkeypatch,
        tmp_path,
        detection=DetectionV1("SLURM", ("sinfo responds",)),
        scheduler=_SLURM,
    )

    result = CliRunner().invoke(wizard, ["--server", "--yes"])

    assert result.exit_code == 0, result.output
    target = tmp_path / "server" / "SLURM.yaml"
    assert target.exists()
    payload = yaml.safe_load(target.read_text(encoding="utf-8"))
    assert payload["SERVER"]["SCHEDULER"] == "SLURM"
    assert payload["SERVER"]["QUEUE_NAME"] == "compute"
    assert payload["SERVER"]["NUM_CORES"] == 6, "queue, not nproc"
    assert payload["SERVER"]["MEM_GB"] == 52
    assert "XTB" in payload, "template program skeletons present"
    assert "sinfo responds" in result.output
    assert "chemsmart sub --test --fake -s SLURM" in result.output


def test_prompts_show_the_queue_ceiling_and_accept_edits(
    tmp_path, monkeypatch
):
    _wire(
        monkeypatch,
        tmp_path,
        detection=DetectionV1("SLURM", ("sinfo responds",)),
        scheduler=_SLURM,
    )

    # cores=4, mem=40, hours=2, scratch none, no program overwrite,
    # skip project/email.
    result = CliRunner().invoke(
        wizard,
        ["--server"],
        input="4\n40\n2\nnone\n\n\n",
    )

    assert result.exit_code == 0, result.output
    assert "queue advertises 6" in result.output
    assert "queue advertises 60000 MB" in result.output
    assert "queue allows up to 48 h" in result.output
    payload = yaml.safe_load(
        (tmp_path / "server" / "SLURM.yaml").read_text(encoding="utf-8")
    )
    assert payload["SERVER"]["NUM_CORES"] == 4
    assert payload["SERVER"]["MEM_GB"] == 40
    assert payload["SERVER"]["NUM_HOURS"] == 2
    assert payload["SERVER"]["SCRATCH_DIR"] is None


def test_an_unqualified_scheduler_downgrades_honestly(tmp_path, monkeypatch):
    _wire(
        monkeypatch,
        tmp_path,
        detection=DetectionV1("SGE", ("SGE_ROOT is set",)),
        scheduler=SchedulerFactsV1(scheduler="SGE", version="8.1.9"),
    )

    result = CliRunner().invoke(wizard, ["--server", "--yes"])

    assert result.exit_code == 0, result.output
    assert "not release-qualified" in result.output
    target = tmp_path / "server" / "local.yaml"
    assert target.exists()
    payload = yaml.safe_load(target.read_text(encoding="utf-8"))
    assert payload["SERVER"]["SCHEDULER"] is None
    assert payload["SERVER"]["SUBMIT_COMMAND"] is None
    assert payload["SERVER"]["NUM_CORES"] == 12, "host facts for local runs"


def test_an_existing_file_is_backed_up_and_programs_survive(
    tmp_path, monkeypatch
):
    _wire(
        monkeypatch,
        tmp_path,
        detection=DetectionV1("SLURM", ("sinfo responds",)),
        scheduler=_SLURM,
    )
    server_dir = tmp_path / "server"
    server_dir.mkdir(parents=True)
    (server_dir / "SLURM.yaml").write_text(
        "SERVER:\n    SCHEDULER: SLURM\n    NUM_CORES: 4\n"
        "GAUSSIAN:\n    EXEFOLDER: /opt/chem/gaussian/g16  # licensed\n",
        encoding="utf-8",
    )

    result = CliRunner().invoke(wizard, ["--server", "--yes"])

    assert result.exit_code == 0, result.output
    assert "Backed up the previous file" in result.output
    backups = list(server_dir.glob("SLURM.yaml.bak-*"))
    assert len(backups) == 1
    merged = (server_dir / "SLURM.yaml").read_text(encoding="utf-8")
    assert "# licensed" in merged
    assert yaml.safe_load(merged)["SERVER"]["NUM_CORES"] == 6
