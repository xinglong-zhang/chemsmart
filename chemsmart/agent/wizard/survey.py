"""Scheduler survey orchestration for the wizard."""

from __future__ import annotations

from dataclasses import dataclass

from chemsmart.agent.wizard.normalize import SCHEDULER_SUBMIT, choose_queue
from chemsmart.agent.wizard.parsers import (
    PbsNodeFacts,
    PbsServerFacts,
    QueueFacts,
    parse_lsf_bqueues_l,
    parse_pbs_pbsnodes_av,
    parse_pbs_qmgr_list_server,
    parse_pbs_qstat_qf_json,
    parse_pbs_qstat_qf_text,
    parse_sge_qconf_sq,
    parse_sge_qconf_sql,
    parse_sge_qhost,
    parse_sge_qstat_gc,
    parse_slurm_scontrol_partition_oneliner,
    parse_slurm_sinfo_json,
)
from chemsmart.agent.wizard.probe import (
    ALL_PROBE_SPECS,
    ProbeSpec,
    run_local_probe,
    run_ssh_probe,
)
from chemsmart.agent.wizard.topology import Topology


@dataclass(frozen=True)
class ScheduleSurvey:
    scheduler: str
    submit_command: str
    queues: list[QueueFacts]
    chosen_queue: str | None
    evidence: dict[str, str]


class AmbiguousSchedulerError(Exception):
    """Raised when more than one scheduler family appears valid."""


def run_schedule_survey(runner, topology: Topology) -> ScheduleSurvey:
    """Probe scheduler commands and normalize queue facts."""

    successes: list[tuple[str, list[QueueFacts], dict[str, str]]] = []

    slurm = _probe_slurm(runner, topology)
    if slurm is not None:
        successes.append(slurm)

    pbs = _probe_pbs(runner, topology)
    if pbs is not None:
        successes.append(pbs)

    lsf = _probe_lsf(runner, topology)
    if lsf is not None:
        successes.append(lsf)

    sge = _probe_sge(runner, topology)
    if sge is not None:
        successes.append(sge)

    if len(successes) > 1:
        families = ", ".join(name for name, _, _ in successes)
        raise AmbiguousSchedulerError(
            f"Multiple scheduler families detected: {families}"
        )
    if not successes:
        raise RuntimeError("No scheduler family produced parseable output.")

    scheduler, queues, evidence = successes[0]
    return ScheduleSurvey(
        scheduler=scheduler,
        submit_command=SCHEDULER_SUBMIT[scheduler],
        queues=queues,
        chosen_queue=choose_queue(queues),
        evidence=evidence,
    )


def _probe_slurm(runner, topology: Topology):
    evidence: dict[str, str] = {}
    sinfo = _run_probe(
        runner, topology, ALL_PROBE_SPECS["survey.slurm.sinfo_json"]
    )
    queues = []
    if sinfo.returncode == 0 and sinfo.stdout.strip():
        parsed = _safe_parse(parse_slurm_sinfo_json, sinfo.stdout)
        if parsed:
            queues = parsed
            evidence["sinfo --json"] = "parsed"
    scontrol = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["survey.slurm.scontrol_partition"],
    )
    if scontrol.returncode == 0 and scontrol.stdout.strip():
        parsed = _safe_parse(
            parse_slurm_scontrol_partition_oneliner,
            scontrol.stdout,
        )
        if parsed:
            queues = parsed
            evidence["scontrol show partition --oneliner"] = "parsed"
    if queues:
        return "SLURM", queues, evidence
    return None


def _probe_pbs(runner, topology: Topology):
    evidence: dict[str, str] = {}
    qstat_json = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["survey.pbs.qstat_json"],
    )
    queues = []
    if qstat_json.returncode == 0 and qstat_json.stdout.strip():
        parsed = _safe_parse(parse_pbs_qstat_qf_json, qstat_json.stdout)
        if parsed:
            queues = _merge_queue_lists(queues, parsed)
            evidence["qstat -Q -f -F json"] = "parsed"
    qstat_text = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["survey.pbs.qstat_text"],
    )
    if qstat_text.returncode == 0 and qstat_text.stdout.strip():
        parsed = _safe_parse(parse_pbs_qstat_qf_text, qstat_text.stdout)
        if parsed:
            queues = _merge_queue_lists(queues, parsed)
            evidence["qstat -Q -f"] = "parsed"
    server_facts = _probe_pbs_server_facts(runner, topology, evidence)
    node_facts = _probe_pbs_node_facts(runner, topology, evidence)
    if queues:
        return (
            "PBS",
            _enrich_pbs_queues(queues, server_facts, node_facts),
            evidence,
        )
    return None


def _probe_lsf(runner, topology: Topology):
    evidence: dict[str, str] = {}
    bqueues_json = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["survey.lsf.bqueues_json"],
    )
    if bqueues_json.returncode == 0:
        evidence["bqueues -o ... -json"] = "ok"
    bqueues_text = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["survey.lsf.bqueues_text"],
    )
    if bqueues_text.returncode == 0 and bqueues_text.stdout.strip():
        parsed = _safe_parse(parse_lsf_bqueues_l, bqueues_text.stdout)
        if parsed:
            evidence["bqueues -l"] = "parsed"
            return "LSF", parsed, evidence
    return None


def _probe_sge(runner, topology: Topology):
    evidence: dict[str, str] = {}
    qconf_sql = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["survey.sge.qconf_sql"],
    )
    if qconf_sql.returncode != 0 or not qconf_sql.stdout.strip():
        return None
    names = _safe_parse(parse_sge_qconf_sql, qconf_sql.stdout) or []
    if not names:
        return None
    evidence["qconf -sql"] = "parsed"

    queue_slots = _probe_sge_qstat_gc(runner, topology, evidence)

    # Probe qhost once — fallback for mem (h_vmem=INFINITY) and cores (slots=1)
    qhost_mem_gb: int | None = None
    qhost_ncpu: int | None = None
    qhost_result = _run_probe(
        runner, topology, ALL_PROBE_SPECS["survey.sge.qhost"]
    )
    if qhost_result.returncode == 0 and qhost_result.stdout.strip():
        parsed = _safe_parse(parse_sge_qhost, qhost_result.stdout)
        if parsed is not None:
            qhost_mem_gb, qhost_ncpu = parsed
            if qhost_mem_gb is not None:
                evidence["qhost mem_total"] = f"{qhost_mem_gb}G"
            if qhost_ncpu is not None:
                evidence["qhost ncpu"] = str(qhost_ncpu)

    queues: list[QueueFacts] = []
    for name in names:
        result = _run_probe(
            runner,
            topology,
            ALL_PROBE_SPECS["survey.sge.qconf_sq"],
            queue_name=name,
        )
        if result.returncode != 0 or not result.stdout.strip():
            continue
        queue = _safe_parse(parse_sge_qconf_sq, result.stdout)
        if queue is None:
            continue
        # SGE slots=1 means "1 slot per serial job", not physical core count.
        # Fall back to qhost NCPU when slots <= 1 so NUM_CORES reflects reality.
        needs_mem = queue.default_mem_gb is None and qhost_mem_gb is not None
        needs_cpu = (queue.default_cores or 0) <= 1 and qhost_ncpu is not None
        slots_total = queue_slots.get(queue.name, queue.slots_total)
        if needs_mem or needs_cpu or slots_total != queue.slots_total:
            queue = QueueFacts(
                name=queue.name,
                default=queue.default,
                max_walltime_hours=queue.max_walltime_hours,
                default_walltime_hours=queue.default_walltime_hours,
                default_mem_gb=(
                    qhost_mem_gb if needs_mem else queue.default_mem_gb
                ),
                default_cores=qhost_ncpu if needs_cpu else queue.default_cores,
                gpus_per_node=queue.gpus_per_node,
                enabled=queue.enabled,
                started=queue.started,
                slots_total=slots_total,
            )
        queues.append(queue)
        evidence[f"qconf -sq {name}"] = "parsed"
    if queues:
        return "SGE", queues, evidence
    return None


def _run_probe(runner, topology: Topology, spec: ProbeSpec, **slots: str):
    if topology.mode == "A":
        return run_local_probe(runner, spec, **slots)
    if topology.mode == "B" and topology.host:
        return run_ssh_probe(runner, topology.host, spec, **slots)
    raise ValueError(f"Unsupported topology: {topology}")


def _safe_parse(parser, payload):
    try:
        return parser(payload)
    except (TypeError, ValueError, KeyError):
        return None


def _probe_sge_qstat_gc(
    runner,
    topology: Topology,
    evidence: dict[str, str],
) -> dict[str, int]:
    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["survey.sge.qstat_gc"],
    )
    if result.returncode != 0 or not result.stdout.strip():
        return {}

    parsed = _safe_parse(parse_sge_qstat_gc, result.stdout)
    if not parsed:
        return {}

    evidence["qstat -g c"] = "parsed"
    return parsed


def _probe_pbs_server_facts(
    runner,
    topology: Topology,
    evidence: dict[str, str],
) -> PbsServerFacts | None:
    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["survey.pbs.qmgr_server"],
    )
    if result.returncode != 0 or not result.stdout.strip():
        return None

    parsed = _safe_parse(parse_pbs_qmgr_list_server, result.stdout)
    if parsed is None:
        return None

    evidence['qmgr -c "list server"'] = "parsed"
    return parsed


def _probe_pbs_node_facts(
    runner,
    topology: Topology,
    evidence: dict[str, str],
) -> PbsNodeFacts | None:
    result = _run_probe(
        runner,
        topology,
        ALL_PROBE_SPECS["survey.pbs.pbsnodes_text"],
    )
    if result.returncode != 0 or not result.stdout.strip():
        return None

    parsed = _safe_parse(parse_pbs_pbsnodes_av, result.stdout)
    if parsed is None:
        return None

    evidence["pbsnodes -av"] = "parsed"
    return parsed


def _merge_queue_lists(
    existing: list[QueueFacts],
    incoming: list[QueueFacts],
) -> list[QueueFacts]:
    if not existing:
        return list(incoming)

    merged = {queue.name: queue for queue in existing}
    order = [queue.name for queue in existing]
    for queue in incoming:
        if queue.name in merged:
            merged[queue.name] = _merge_queue(merged[queue.name], queue)
            continue
        merged[queue.name] = queue
        order.append(queue.name)
    return [merged[name] for name in order]


def _merge_queue(primary: QueueFacts, secondary: QueueFacts) -> QueueFacts:
    return QueueFacts(
        name=primary.name,
        default=primary.default or secondary.default,
        max_walltime_hours=_prefer_value(
            primary.max_walltime_hours,
            secondary.max_walltime_hours,
        ),
        default_walltime_hours=_prefer_value(
            primary.default_walltime_hours,
            secondary.default_walltime_hours,
        ),
        default_mem_gb=_prefer_value(
            primary.default_mem_gb,
            secondary.default_mem_gb,
        ),
        default_cores=_prefer_value(
            primary.default_cores,
            secondary.default_cores,
        ),
        gpus_per_node=_prefer_value(
            primary.gpus_per_node,
            secondary.gpus_per_node,
        ),
        enabled=primary.enabled and secondary.enabled,
        started=primary.started and secondary.started,
        slots_total=_prefer_value(
            primary.slots_total,
            secondary.slots_total,
        ),
    )


def _enrich_pbs_queues(
    queues: list[QueueFacts],
    server_facts: PbsServerFacts | None,
    node_facts: PbsNodeFacts | None,
) -> list[QueueFacts]:
    return [
        QueueFacts(
            name=queue.name,
            default=queue.default
            or (
                server_facts is not None
                and server_facts.default_queue == queue.name
            ),
            max_walltime_hours=queue.max_walltime_hours,
            default_walltime_hours=queue.default_walltime_hours,
            default_mem_gb=_enrich_pbs_mem_gb(queue, server_facts, node_facts),
            default_cores=_enrich_pbs_cores(queue, server_facts, node_facts),
            gpus_per_node=queue.gpus_per_node,
            enabled=queue.enabled,
            started=queue.started,
            slots_total=queue.slots_total,
        )
        for queue in queues
    ]


def _enrich_pbs_mem_gb(
    queue: QueueFacts,
    server_facts: PbsServerFacts | None,
    node_facts: PbsNodeFacts | None,
) -> int | None:
    if queue.default_mem_gb is not None:
        return queue.default_mem_gb
    if server_facts is not None and server_facts.default_mem_gb is not None:
        return server_facts.default_mem_gb
    if node_facts is not None:
        return node_facts.min_mem_gb
    return None


def _enrich_pbs_cores(
    queue: QueueFacts,
    server_facts: PbsServerFacts | None,
    node_facts: PbsNodeFacts | None,
) -> int | None:
    cores = queue.default_cores
    if cores is None and server_facts is not None:
        cores = server_facts.default_cores
    if (
        node_facts is not None
        and node_facts.min_ncpus is not None
        and (cores is None or cores <= 1)
    ):
        return node_facts.min_ncpus
    return cores


def _prefer_value(primary, secondary):
    return primary if primary is not None else secondary
