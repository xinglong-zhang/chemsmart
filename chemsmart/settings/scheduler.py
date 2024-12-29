import ast
import logging
import os
import shlex
import socket
import subprocess
import sys
from abc import abstractmethod
from subprocess import check_output
from tempfile import TemporaryDirectory

from pyatoms.cli.submitters import (
    Bridges2SubmitscriptWriter,
    FugakuSubmitscriptWriter,
    NullSubmitscriptWriter,
    PBSSubmitscriptWriter,
    SlurmSubmitscriptWriter,
)
from pyatoms.io.monty import MSONable
from pyatoms.jobs.host import Host, HostQueue
from pyatoms.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)


class SchedulerFactory:
    def __init__(self, style, **kwargs):
        self.style = style
        self.kwargs = kwargs

    def __call__(self):
        scheduler_cls = [
            s for s in Scheduler.subclasses() if self.style == s.STYLE
        ]
        if len(scheduler_cls) == 0:
            raise ValueError(f"No such server of style: {self.style}")

        if len(scheduler_cls) > 1:
            raise RuntimeError(
                f"Multiple servers of style {self.style}: {scheduler_cls}"
            )

        scheduler_cls = scheduler_cls[0]
        return scheduler_cls(**self.kwargs)

    @classmethod
    def from_dict(cls, d):
        return cls(**d)


class Scheduler(RegistryMixin, MSONable):
    """Shares methods for getting information about hosts from scheduler.

    Args:
        kwargs (dict): kwargs contains all the keyword arguments to be passed to the SubmitscriptWriter
    """

    SUBMITSCRIPT_WRITER = NotImplemented

    def __init__(
        self,
        default_queue=None,
        default_num_hours=None,
        extra_commands=None,
        fqdn_hostnames=True,
        use_localhost=False,
        **kwargs,
    ):
        self.default_queue = default_queue
        self.default_num_hours = default_num_hours
        self.extra_commands = extra_commands
        self.fqdn_hostnames = fqdn_hostnames
        self.use_localhost = use_localhost
        self.kwargs = kwargs

    def __repr__(self):
        return (
            f"{self.__class__.__qualname__}<default_queue={self.default_queue}, "
            f"default_num_hours={self.default_num_hours}>"
        )

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return NotImplemented
        return self.__dict__ == other.__dict__

    @abstractmethod
    def hostnames(self):
        raise NotImplementedError

    @abstractmethod
    def _submit_command(self, *args, **kwargs):
        raise NotImplementedError

    def submitscript_kwargs(self, *args, **kwargs):
        return {
            "default_queue": self.default_queue,
            "default_num_hours": self.default_num_hours,
            "extra_commands": self.extra_commands,
        } | self.kwargs

    def submitscript_writer(
        self,
        job,
        server,
        num_gpus=None,
        num_cpus=None,
        memory_gigs=None,
        **kwargs,
    ):
        num_gpus = server.NUM_GPUS if num_gpus is None else num_gpus
        num_cpus = server.NUM_CORES if num_cpus is None else num_cpus

        # Scale memory to number of CPUs used
        memory_gigs = (
            int(server.MAX_MEM_GIGS / server.NUM_CORES * num_cpus)
            if memory_gigs is None
            else memory_gigs
        )

        # Scale number of CPU with number of GPUs used
        if num_gpus > 0:
            num_cpus = int(num_gpus / server.NUM_GPUS * num_cpus)

        return self.SUBMITSCRIPT_WRITER(
            job=job,
            num_cpus=num_cpus,
            memory_gigs=memory_gigs,
            num_gpus=num_gpus,
            servername=server.servername,
            **kwargs,
            **self.submitscript_kwargs(),
        )

    def submit_job(self, job, server, test=False, **kwargs):
        submitscript = self._write_submitscript(
            job=job, server=server, **kwargs
        )
        if not test:
            self._check_existing_jobs(job)
            self.submit_file(submitscript)

    def submit_file(self, submitscript):
        self._run_submit_command(submitscript)

    def _write_submitscript(self, job, server, **kwargs):
        writer = self.submitscript_writer(job=job, server=server, **kwargs)
        return writer.write()

    def _run_submit_command(self, submitscript):
        submit_folder = os.path.dirname(submitscript)
        submit_file = os.path.basename(submitscript)
        command = self._submit_command(submit_file)

        if command is None:
            raise ValueError(f"Cannot submit on {self}")

        p = subprocess.Popen(shlex.split(command), cwd=submit_folder)
        return p.wait()

    def _check_existing_jobs(self, job):
        from pyatoms.io.gaussian.utils import ClusterHelper
        from pyatoms.jobs.gaussian import GaussianJob

        if not isinstance(job, GaussianJob) or job.label is None:
            return

        cluster_helper = ClusterHelper()
        running_job_ids, running_job_names = (
            cluster_helper.get_gaussian_running_jobs()
        )

        if job.label in running_job_names:
            logger.info(
                f"Warning: submitting job with duplicate name: {job.label}"
            )
            sys.exit(f"Duplicate job NOT submitted: {job.label}")

    def create_gpu_queue(self, exclude_localhost=False):
        if (hostnames := self.hostnames()) is None and self.use_localhost:
            # Use localhost if no hostnames detected
            hostnames = [self._localhost]

        if hostnames is None:
            return None

        if exclude_localhost:
            hostnames = self._remove_localhost(hostnames)

        hosts = [
            Host(
                name=name,
                gpus=self._get_gpus_from_cuda_visible_devices(hostname=name),
                priority=0,
            )
            for name in hostnames
        ]
        return HostQueue.from_hosts(hosts)

    def _get_gpus_from_cuda_visible_devices(self, hostname):
        """Read CUDA_VISIBLE_DEVICES from localhost environment."""
        cuda_visible_devices = os.environ.get("CUDA_VISIBLE_DEVICES")
        logger.info(f"CUDA_VISIBLE_DEVICES: {cuda_visible_devices}")
        if cuda_visible_devices is None:
            return []

        gpus = cuda_visible_devices.split(",")
        logger.info(f"Got gpus: {gpus} from {hostname}")
        return gpus

    def _get_least_busy_gpus(self, hostname, num_gpus):
        exe = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "executables",
            "get_least_busy_gpus.py",
        )
        cmd = f'ssh {hostname} "{exe} {num_gpus}"'
        out = check_output(
            shlex.split(cmd), encoding="UTF-8"
        )  # if no encoding will give bytes
        gpus = ast.literal_eval(out)
        logger.info(f"Got gpus: {gpus} from {hostname}")
        return gpus

    def create_host_queue(self, hostnames=None, exclude_localhost=False):
        hostnames = hostnames or self.hostnames()

        if hostnames is None and self.use_localhost:
            logger.info("No hostnames detected. Using localhost")
            hostnames = [self._localhost]

        if hostnames is None:
            raise ValueError(f"Could not get any hostnames for {self}")

        if exclude_localhost:
            hostnames = self._remove_localhost(hostnames)

        return HostQueue.from_hostnames(
            hostnames=hostnames,
            localhost=self._localhost,
            priority=0,
            localhost_priority=10,
        )

    @property
    def _localhost(self):
        local_hostname = (
            socket.getfqdn()
        )  # more specific than socket.gethostname()
        if not self.fqdn_hostnames:
            local_hostname = local_hostname.split(".")[0]
        return local_hostname

    def _remove_localhost(self, hostnames):
        hostnames = hostnames.copy()
        local_host_index = [
            i for i, host in enumerate(hostnames) if self._localhost in host
        ]
        if len(local_host_index) > 1:
            raise ValueError("Found more than one match for localhost")

        del hostnames[local_host_index[0]]
        logger.info(f"Hosts after excluding localhost: {hostnames}")

        if len(hostnames) == 0:
            raise ValueError("No more hosts after removing localhost!")

        return hostnames

    def _read_node_names(self, nodefile):
        with open(nodefile) as f:
            hosts = list(set(f.readlines()))
        hosts = [host.strip() for host in hosts]

        if not self.fqdn_hostnames:
            # NSCC OpenMPI 4.1.2 does not recognize --host specification with FQDN, so
            # we will need to remove the domain name e.g.
            # x1000c0s2b0n0.hostmgmt2000.cm.asp2a.nscc.sg -> x1000c0s2b0n0
            hosts = [h.split(".")[0] for h in hosts]

        logger.info(f"Hosts: {hosts}")
        return hosts

    @classmethod
    def from_name(cls, name, **kwargs):
        scheduler_cls = cls.class_from_name(name)
        return scheduler_cls(**kwargs)

    @classmethod
    def class_from_name(cls, name):
        for c in cls.subclasses():
            if name == c.STYLE:
                return c
        raise ValueError(f"No such scheduler mixin: {name}")


class FugakuScheduler(Scheduler):
    STYLE = "fugaku"
    SUBMITSCRIPT_WRITER = FugakuSubmitscriptWriter

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def hostnames(self):
        nodefile = os.environ.get("PBS_NODEFILE")
        logger.debug(f"Nodefile: {nodefile}")
        if nodefile is None:
            return None

        return self._read_node_names(nodefile)

    def _submit_command(self, submitscript):
        return f"pjsub {submitscript}"


class NullScheduler(Scheduler):
    """Scheduler that does nothing for servers that do not have a scheduler."""

    STYLE = "null"
    SUBMITSCRIPT_WRITER = NullSubmitscriptWriter

    def hostnames(self):
        # return ['127.0.0.1']
        import socket

        return [socket.gethostname()]

    def _submit_command(self, server, *args, **kwargs):
        raise ValueError(f"No scheduler for {server}.")

    def _write_submitscript(self, job, server, *args, **kwargs):
        from pyatoms.cli.submitters import NullSubmitscriptWriter

        return NullSubmitscriptWriter(job=job).write()

    def submitscript_kwargs(self, *args, **kwargs):
        return {}


class PBSScheduler(Scheduler):
    """Default PBS scheduler."""

    STYLE = "pbs"
    SUBMITSCRIPT_WRITER = PBSSubmitscriptWriter

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def hostnames(self):
        nodefile = os.environ.get("PBS_NODEFILE")
        logger.debug(f"Nodefile: {nodefile}")
        if nodefile is None:
            return None

        return self._read_node_names(nodefile)

    def _submit_command(self, submitscript):
        return f"qsub {submitscript}"

    def create_gpu_queue(self, hostnames=None, exclude_localhost=False):
        hostnames = self.hostnames() if hostnames is None else hostnames

        if hostnames is None and self.use_localhost:
            hostnames = [self._localhost]

        if hostnames is None:
            raise ValueError(f"Could not get any hostnames for {self}")

        if exclude_localhost:
            hostnames = self._remove_localhost(hostnames)
            self._set_localhost_cuda_visible_devices()

        hosts = [
            Host(name=name, gpus=self._get_gpus_from_host(name), priority=0)
            for name in hostnames
        ]
        logger.info(f"GPU hosts: {hosts}")
        return HostQueue.from_hosts(hosts)

    def _set_localhost_cuda_visible_devices(self):
        gpus = self._get_gpus_from_host(hostname=self._localhost)
        gpus = ",".join(gpus)
        os.environ["CUDA_VISIBLE_DEVICES"] = gpus
        logger.info(f"Set localhost CUDA_VISIBLE_DEVICES to {gpus}")

    def _get_gpus_from_host(self, hostname):
        """Return GPUs in terms of 0, 1, etc.

        Needed to facilitate multi node training where mpirun -x CUDA_VISIBLE_DEVICES has to be specified
        and the same one needs to be specified for all nodes. This is not possible if using the
        detailed GPU ID.
        """
        gpus = self._get_gpus_from_host_specific(hostname)
        return [str(i) for i in range(len(gpus))]

    def _get_gpus_from_host_specific(self, hostname):
        """Get specific GPU IDs.

        Works for ASPIRE2A, not sure about other supercomputing resources using PBS.

        Note: Fails for asp2acpu jobs as the ${PBS_JOBID}.env file seems to be created only
        for asp2agpu jobs.
        """
        exe = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "executables",
            "get_host_gpus.sh",
        )
        pbs_jobid = os.environ.get("PBS_JOBID")

        if pbs_jobid is None:
            return self._get_gpus_from_cuda_visible_devices(hostname)

        logger.info(f"Job ID: {pbs_jobid}")
        cmd = f'ssh  {hostname} "{exe} {pbs_jobid}"'

        try:
            out = check_output(
                shlex.split(cmd), encoding="UTF-8", env=os.environ
            )  # if no encoding will give bytes
        except subprocess.CalledProcessError:
            # Error in getting host gpus.
            logger.exception("Found exception when getting host GPUs")
            out = ""

        gpus = out.strip().split(",")

        logger.info(f"Got gpus: {gpus} from {hostname}")
        return gpus


class SlurmScheduler(Scheduler):
    """Default SLURM scheduler."""

    STYLE = "slurm"
    SUBMITSCRIPT_WRITER = SlurmSubmitscriptWriter

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def hostnames(self):
        with TemporaryDirectory() as tempdir:
            cmd = "scontrol show hostname"
            outfile = os.path.join(tempdir, "hostname.out")
            errfile = os.path.join(tempdir, "hostname.err")
            with open(outfile, "w") as out, open(errfile, "w") as err:
                p = subprocess.Popen(shlex.split(cmd), stdout=out, stderr=err)
                p.communicate()
            return self._read_node_names(outfile)

    def _submit_command(self, submitscript):
        return f"sbatch {submitscript}"


class Bridges2SlurmScheduler(SlurmScheduler):
    """Slurm scheduler with submit script writer customized for Bridges."""

    STYLE = "bridges2-slurm"
    SUBMITSCRIPT_WRITER = Bridges2SubmitscriptWriter
