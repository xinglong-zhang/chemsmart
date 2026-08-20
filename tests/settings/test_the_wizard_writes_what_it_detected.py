"""The wizard writes detected truth and preserves hand-tuned work.

Defaults come from the chosen queue, never from the login node (this host
advertises 6 CPUs on its partition while nproc says 12); an existing file
is backed up and only its SERVER block is regenerated -- program blocks
survive byte-for-byte, comments and all.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from chemsmart.settings.probe import (
    HostFactsV1,
    QueueFactsV1,
    SchedulerFactsV1,
)
from chemsmart.settings.wizard import (
    derive_choices,
    extract_extra_commands,
    render_server_block,
    splice_server_block,
    write_server_yaml,
)

_HOST = HostFactsV1(
    hostname="ax41",
    os_name="Linux",
    cpu_count=12,
    mem_kb=65761448,
    scratch_candidates=("/scratch/chemsmart", "/tmp"),
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


def _choices():
    return derive_choices(
        scheduler=_SLURM,
        queue=_QUEUE,
        host=_HOST,
        scratch_dir="/scratch/chemsmart",
    )


def test_defaults_come_from_the_queue_not_the_login_node():
    choices = _choices()

    assert choices.num_cores == 6, "sinfo is the authority, not nproc"
    assert choices.mem_gb == int(0.9 * 61440000 / 1024 / 1024)  # 52
    assert choices.num_hours == 24  # min(48h queue limit, 24)
    assert choices.scheduler == "SLURM"
    assert choices.queue_name == "compute"
    assert choices.submit_command == "/opt/slurm/current/bin/sbatch"


def test_an_unsubmittable_scheduler_is_written_as_local_with_the_reason():
    sge = SchedulerFactsV1(scheduler="SGE", version="8.1.9")
    choices = derive_choices(
        scheduler=sge, queue=None, host=_HOST, scratch_dir=None
    )

    assert choices.scheduler is None
    assert choices.submit_command is None
    assert choices.num_cores == 12  # host fallback for local runs
    assert any("not release-qualified" in note for note in choices.provenance)


def test_the_rendered_block_parses_and_states_its_provenance():
    block = render_server_block(_choices(), host=_HOST)
    payload = yaml.safe_load(block)

    assert payload["SERVER"]["SCHEDULER"] == "SLURM"
    assert payload["SERVER"]["QUEUE_NAME"] == "compute"
    assert payload["SERVER"]["NUM_CORES"] == 6
    assert payload["SERVER"]["MEM_GB"] == 52
    assert "# detected SLURM (slurm 26.05.2)" in block
    assert "queue 'compute': 6 cores, 60000 MB per node" in block


_EXISTING = (
    "# my own note above the server block\n"
    "SERVER:\n"
    "    SCHEDULER: SLURM\n"
    "    NUM_CORES: 4\n"
    "GAUSSIAN:\n"
    "    EXEFOLDER: /opt/chem/gaussian/g16  # licensed install\n"
    "    MODULES: |\n"
    "        module purge\n"
    "ORCA:\n"
    "    EXEFOLDER: /opt/chem/orca/6.1.1\n"
)


def test_splice_replaces_server_and_preserves_programs_byte_for_byte():
    block = render_server_block(_choices(), host=_HOST)
    merged = splice_server_block(_EXISTING, block)

    program_tail = _EXISTING[_EXISTING.index("GAUSSIAN:") :]
    assert merged.endswith(program_tail)
    assert "NUM_CORES: 6" in merged
    assert "NUM_CORES: 4" not in merged
    assert (
        "# my own note above the server block" not in merged
    ), "comments describing the OLD server block go with it"
    assert "# licensed install" in merged


def test_splice_prepends_when_no_server_block_exists():
    block = render_server_block(_choices(), host=_HOST)
    merged = splice_server_block("GAUSSIAN:\n    EXEFOLDER: null\n", block)

    assert merged.startswith("# Written by `chemsmart wizard`")
    assert merged.endswith("GAUSSIAN:\n    EXEFOLDER: null\n")


def test_write_backs_up_then_splices(tmp_path: Path):
    target = tmp_path / "SLURM.yaml"
    target.write_text(_EXISTING, encoding="utf-8")
    block = render_server_block(_choices(), host=_HOST)

    backup = write_server_yaml(target, block, template_text="")

    assert backup is not None and backup.exists()
    assert backup.read_text(encoding="utf-8") == _EXISTING
    assert backup.name.startswith("SLURM.yaml.bak-")
    merged = target.read_text(encoding="utf-8")
    assert "# licensed install" in merged
    assert yaml.safe_load(merged)["SERVER"]["NUM_CORES"] == 6
    assert (target.stat().st_mode & 0o777) == 0o600


def test_a_new_file_grows_from_the_template(tmp_path: Path):
    template = "SERVER:\n    SCHEDULER: null\nXTB:\n    LOCAL_RUN: true\n"
    target = tmp_path / "SLURM.yaml"
    block = render_server_block(_choices(), host=_HOST)

    backup = write_server_yaml(target, block, template_text=template)

    assert backup is None
    merged = target.read_text(encoding="utf-8")
    assert yaml.safe_load(merged)["SERVER"]["SCHEDULER"] == "SLURM"
    assert merged.endswith("XTB:\n    LOCAL_RUN: true\n")


def test_hand_tuned_extra_commands_are_carried_through_a_refresh():
    """`ulimit -s unlimited` is operational knowledge the host cannot
    re-derive; a refresh must not silently drop it (live-host finding)."""

    previous = (
        "SERVER:\n"
        "  SCHEDULER: SLURM\n"
        "  EXTRA_COMMANDS: |\n"
        "    ulimit -s unlimited\n"
        "    module load mkl\n"
        "\n"
        "ORCA:\n"
        "    EXEFOLDER: /opt/chem/orca\n"
    )
    carried = extract_extra_commands(previous)
    assert carried == "ulimit -s unlimited\nmodule load mkl"

    block = render_server_block(
        _choices(), host=_HOST, carried_extra_commands=carried
    )
    payload = yaml.safe_load(block)
    assert payload["SERVER"]["EXTRA_COMMANDS"] == (
        "ulimit -s unlimited\nmodule load mkl\n"
    )
    assert "# carried from the previous configuration" in block


def test_the_placeholder_and_absence_carry_nothing():
    placeholder = (
        "SERVER:\n"
        "    EXTRA_COMMANDS: |\n"
        "        # Host commands required before execution.\n"
    )
    assert extract_extra_commands(placeholder) is None
    assert extract_extra_commands("SERVER:\n    SCHEDULER: null\n") is None
    assert extract_extra_commands("") is None

    inline = "SERVER:\n    EXTRA_COMMANDS: ulimit -s unlimited\n"
    assert extract_extra_commands(inline) == "ulimit -s unlimited"
