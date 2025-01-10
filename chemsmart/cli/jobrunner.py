import functools
import logging
import click


logger = logging.getLogger(__name__)


def jobrunner_options(f):
    @click.option(
        "-E",
        "--exclude-localhost/--no-exclude-localhost",
        default=False,
        help="if using queue manager, this excludes local host where python script will be running",
    )
    @click.option(
        "--hosts",
        "--host",
        type=str,
        default=None,
        multiple=True,
        help="Multiple hosts can be specified",
    )
    @click.option(
        "-n", "--num-nodes", type=int, default=1, help="number of nodes"
    )
    @click.option(
        "-c",
        "--num-cpus",
        type=int,
        default=None,
    )
    @click.option(
        "-g",
        "--num-gpus",
        "--num-gpus-per-node",
        type=int,
        default=None,
        help="Number of gpus per node. Defaults to number of GPUs on specified server if None.",
    )
    @click.option(
        "-m", "--memory-gigs", type=int, default=None, help="Memory in gigs"
    )
    @click.option(
        "-G",
        "--num-gpus-per-node-per-job",
        type=int,
        default=None,
        help="Number of gpus per node for a single job. Defaults to number of GPUs on specified server if None.",
    )
    @click.option(
        "-s",
        "--servername",
        "--server",
        type=str,
        default=None,
        help="Server. If not specified, will try to automatically determine and use the current server.",
    )
    @click.option(
        "--use-host-queues/--no-use-host-queues",
        default=True,
        help="Use host queues",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def create_jobrunner(
    servername,
    num_nodes,
    num_cpus,
    num_gpus,
    num_gpus_per_node_per_job,
    memory_gigs,
    exclude_localhost,
    ship=False,
    fake=False,
    hosts=None,
    use_host_queues=True,
):
    from chemsmart.settings.server import Server
    from chemsmart.jobs.runner import JobRunner

    server = Server.from_server(servername)

    if ship:
        use_host_queues = False

    host_queue = (
        None
        if not use_host_queues
        else server.create_host_queue(
            exclude_localhost=exclude_localhost, hostnames=hosts
        )
    )
    gpu_queue = (
        None
        if not use_host_queues
        else server.create_gpu_queue(exclude_localhost=exclude_localhost)
    )

    if num_gpus_per_node_per_job is None:
        num_gpus_per_node_per_job = num_gpus

    return JobRunner(
        memory_gigs=memory_gigs,
        num_nodes=num_nodes,
        num_gpus_per_node=num_gpus_per_node_per_job,
        num_tasks_per_node=num_cpus,
        server=server,
        fake=fake,
        exclude_localhost=exclude_localhost,
        gpu_queue=gpu_queue,
        host_queue=host_queue,
    )
