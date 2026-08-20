"""``chemsmart wizard`` — the first command a new ChemSmart user runs.

It analyzes the host it runs on (scheduler, queues, cores, memory,
scratch), confirms every choice with the detected value as the default,
writes the scheduler-named server YAML the existing ``-s``-omitted lookup
already finds, verifies the result with tabulated non-fatal checks, and
ends with the documented qualification step. Detected facts are always
announced; nothing is silent, and the scheduler — never the login node —
is the resource authority.
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path

import click
import yaml

from chemsmart.settings.probe import (
    HostFactsV1,
    QueueFactsV1,
    SchedulerFactsV1,
    detect_scheduler,
    gather_host_facts,
)
from chemsmart.settings.wizard import (
    SUBMITTABLE_SCHEDULERS,
    ServerChoicesV1,
    derive_choices,
    extract_top_level_block,
    extract_extra_commands,
    render_server_block,
    run_verification,
    splice_top_level_block,
    write_server_yaml,
)

_PROGRAM_BLOCKS = ("GAUSSIAN", "ORCA", "PYSCF", "XTB", "NCIPLOT")


def _template_text() -> str:
    template = resources.files("chemsmart.settings").joinpath(
        "templates/.chemsmart/server/local.yaml"
    )
    return template.read_text(encoding="utf-8")


def _user_server_dir() -> Path:
    from chemsmart.settings.user import CHEMSMARTUserSettings

    return Path(CHEMSMARTUserSettings().user_server_dir)


def _announce(lines) -> None:
    for line in lines:
        click.echo(f"  · {line}")


def _probe_scheduler(detection) -> SchedulerFactsV1 | None:
    if detection.scheduler == "SLURM":
        from chemsmart.settings.probe.slurm import probe_slurm

        return probe_slurm()
    if detection.scheduler == "PBS":
        from chemsmart.settings.probe.pbs import probe_pbs

        return probe_pbs()
    if detection.scheduler is not None:
        return SchedulerFactsV1(
            scheduler=detection.scheduler, evidence=detection.evidence
        )
    return None


def _choose_queue(
    scheduler: SchedulerFactsV1 | None, *, assume_yes: bool
) -> QueueFactsV1 | None:
    if scheduler is None or not scheduler.queues:
        return None
    usable = [queue for queue in scheduler.queues if queue.available]
    if not usable:
        usable = list(scheduler.queues)
    default = next(
        (queue for queue in usable if queue.is_default), usable[0]
    )
    if assume_yes or len(usable) == 1:
        return default
    names = [queue.name for queue in usable]
    chosen = click.prompt(
        "Queue/partition",
        type=click.Choice(names),
        default=default.name,
        show_choices=True,
    )
    return next(queue for queue in usable if queue.name == chosen)


def _confirm_choices(
    derived: ServerChoicesV1,
    queue: QueueFactsV1 | None,
    host: HostFactsV1,
    *,
    assume_yes: bool,
) -> ServerChoicesV1:
    if assume_yes:
        return derived
    ceiling = ""
    if queue is not None and queue.cores_per_node:
        ceiling = f" (queue advertises {queue.cores_per_node})"
    cores = click.prompt(
        f"Cores per job{ceiling}", default=derived.num_cores, type=int
    )
    mem_note = ""
    if queue is not None and queue.mem_kb_per_node:
        mem_note = (
            f" (queue advertises {queue.mem_kb_per_node // 1024} MB per node)"
        )
    elif host.mem_kb:
        mem_note = f" (host has {host.mem_kb // 1024} MB)"
    mem_gb = click.prompt(
        f"Memory per job in GB{mem_note}", default=derived.mem_gb, type=int
    )
    hours_note = ""
    if queue is not None and queue.max_time_seconds:
        hours_note = (
            f" (queue allows up to {queue.max_time_seconds // 3600} h)"
        )
    hours = click.prompt(
        f"Wall time in hours{hours_note}", default=derived.num_hours, type=int
    )
    scratch_default = derived.scratch_dir or "none"
    scratch = click.prompt(
        "Scratch directory ('none' disables scratch)",
        default=scratch_default,
    )
    return ServerChoicesV1(
        scheduler=derived.scheduler,
        queue_name=derived.queue_name,
        num_hours=hours,
        mem_gb=mem_gb,
        num_cores=cores,
        num_gpus=derived.num_gpus,
        num_threads=cores,
        submit_command=derived.submit_command,
        scratch_dir=None if scratch.strip().lower() == "none" else scratch,
        provenance=derived.provenance,
    )


def _offer_program_overwrites(
    target: Path, template_text: str, *, assume_yes: bool
) -> None:
    if assume_yes or not target.exists():
        return
    existing = target.read_text(encoding="utf-8")
    present = [
        name for name in _PROGRAM_BLOCKS if f"\n{name}:" in "\n" + existing
    ]
    if not present:
        return
    if not click.confirm(
        "Reconfigure any program blocks from the template? (kept as-is "
        "otherwise)",
        default=False,
    ):
        return
    for name in present:
        if not click.confirm(f"  overwrite {name}?", default=False):
            continue
        block = extract_top_level_block(template_text, name)
        if not block:
            continue
        existing = splice_top_level_block(existing, name, block)
        folder = click.prompt(
            f"  {name} executable folder (Enter to skip)",
            default="",
            show_default=False,
        ).strip()
        target.write_text(existing, encoding="utf-8")
        if folder:
            from chemsmart.cli.config import set_program_exefolder

            set_program_exefolder(target.parent, name, folder)
            existing = target.read_text(encoding="utf-8")
    target.write_text(existing, encoding="utf-8")


def _offer_usersettings(*, assume_yes: bool) -> None:
    if assume_yes:
        return
    project = click.prompt(
        "Scheduler account/project for --account directives (Enter to skip)",
        default="",
        show_default=False,
    ).strip()
    email = click.prompt(
        "Notification email (Enter to skip)", default="", show_default=False
    ).strip()
    if not project and not email:
        return
    from chemsmart.settings.user import CHEMSMARTUserSettings

    settings_path = (
        Path(CHEMSMARTUserSettings().config_dir) / "usersettings.yaml"
    )
    payload = {}
    if settings_path.exists():
        payload = yaml.safe_load(
            settings_path.read_text(encoding="utf-8")
        ) or {}
    if project:
        payload["PROJECT"] = project
    if email:
        payload["EMAIL"] = email
    settings_path.parent.mkdir(parents=True, exist_ok=True)
    settings_path.write_text(
        yaml.safe_dump(payload, sort_keys=True), encoding="utf-8"
    )
    click.echo(f"  · wrote {settings_path}")


def _print_verification(choices: ServerChoicesV1) -> int:
    click.echo("\nVerification:")
    failures = 0
    for check in run_verification(choices):
        mark = {"ok": "ok  ", "warn": "warn", "fail": "FAIL", "skipped": "skip"}[
            check.status
        ]
        detail = f" — {check.detail}" if check.detail else ""
        click.echo(f"  [{mark}] {check.name}{detail}")
        if check.status == "fail":
            failures += 1
    if failures:
        click.echo(
            f"{failures} check(s) failed; the configuration was still "
            "written — fix the environment and re-run the wizard."
        )
    return failures


@click.command(name="wizard")
@click.option(
    "--server",
    "server_only",
    is_flag=True,
    help="Run the server setup only (bare `chemsmart wizard` also offers "
    "the agent-layer setup afterwards).",
)
@click.option(
    "--yes",
    "assume_yes",
    is_flag=True,
    help="Accept every detected default without prompting.",
)
@click.option(
    "--no-verify",
    is_flag=True,
    help="Skip the post-write verification checks.",
)
def wizard(server_only, assume_yes, no_verify):
    """Analyze this host and write a working server configuration."""

    host = gather_host_facts()
    click.echo(
        f"Host: {host.hostname} · {host.os_name} {host.os_release} · "
        f"{host.machine} · {host.cpu_count} cores · "
        f"{(host.mem_kb or 0) // 1024} MB"
    )
    detection = detect_scheduler()
    _announce(detection.evidence)
    scheduler = _probe_scheduler(detection)
    if scheduler is not None and scheduler.version:
        click.echo(f"Scheduler: {scheduler.scheduler} ({scheduler.version})")
    elif scheduler is not None:
        click.echo(f"Scheduler: {scheduler.scheduler}")
    if scheduler is not None and scheduler.queues:
        click.echo("Queues:")
        for queue in scheduler.queues:
            star = "*" if queue.is_default else " "
            state = "up" if queue.available else "down"
            hours = (
                f"{queue.max_time_seconds // 3600}h"
                if queue.max_time_seconds
                else "no limit"
            )
            click.echo(
                f"  {star} {queue.name}: {state}, "
                f"{queue.cores_per_node or '?'} cores, "
                f"{(queue.mem_kb_per_node or 0) // 1024} MB, "
                f"{queue.node_count or '?'} node(s), max {hours}"
            )
    queue = _choose_queue(scheduler, assume_yes=assume_yes)
    scratch = host.scratch_candidates[0] if host.scratch_candidates else None
    derived = derive_choices(
        scheduler=scheduler, queue=queue, host=host, scratch_dir=scratch
    )
    if (
        scheduler is not None
        and scheduler.scheduler not in SUBMITTABLE_SCHEDULERS
    ):
        click.echo(
            f"{scheduler.scheduler} was detected, but ChemSmart submission "
            "is not release-qualified for it; writing a local-run "
            "configuration that records the detected facts."
        )
    choices = _confirm_choices(
        derived, queue, host, assume_yes=assume_yes
    )

    name = choices.scheduler or "local"
    target = _user_server_dir() / f"{name}.yaml"
    template_text = _template_text()
    carried = (
        extract_extra_commands(target.read_text(encoding="utf-8"))
        if target.exists()
        else None
    )
    if carried is not None:
        click.echo(
            "Carrying the previous EXTRA_COMMANDS forward (host commands "
            "are not detectable)."
        )
    block = render_server_block(
        choices, host=host, carried_extra_commands=carried
    )
    backup = write_server_yaml(target, block, template_text=template_text)
    if backup is not None:
        click.echo(f"Backed up the previous file to {backup}")
    click.echo(f"Wrote {target}")
    _offer_program_overwrites(
        target, template_text, assume_yes=assume_yes
    )
    from chemsmart.cli.config import ensure_program_blocks

    ensure_program_blocks(target.parent)
    _offer_usersettings(assume_yes=assume_yes)

    if not no_verify:
        _print_verification(choices)

    click.echo(
        "\nQualify the configuration without submitting anything:\n"
        f"  chemsmart sub --test --fake -s {name} xtb -p test "
        "-f your_molecule.xyz sp\n"
        "Inspect the generated script (modules, scratch, cores, memory) "
        "before a real submission."
    )
    if not server_only:
        click.echo(
            "\nNext: configure the agent layer (provider, model, key) with\n"
            "  chemsmart config agent"
        )


__all__ = ["wizard"]
