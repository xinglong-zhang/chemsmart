import copy
import logging
import os
import re
import socket
from abc import abstractmethod
from functools import lru_cache
from itertools import chain

from pyatoms.io.monty import MSONable
from pyatoms.io.yaml import YAMLCreatableMixin
from pyatoms.settings.server.execution import (
    ExecutionMixin,
    LocalExecutionMixin,
    MpirunExecutionMixin,
    SlurmExecutionMixin,
)
from pyatoms.settings.server.scheduler import (
    NullScheduler,
    PBSScheduler,
    Scheduler,
    SchedulerFactory,
    SlurmScheduler,
)
from pyatoms.settings.user import UserServerSettingsMixin
from pyatoms.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)


class ServerFactory(YAMLCreatableMixin):
    """Callable factory for creating `Server`s.

    The appropriate servers are created based on `style`, which reference
    to useful predefined `Servers`. They have the 'servername` attribute,
    which can be any string that labels the created `Server` for easy reference
    later on.

    Usage:

    Can be initialized from yaml file, if so `servername` argument
    will be the name of the yaml file without the suffix.

        factory = ServerFactory.from_yaml('myserver.yaml')
        server = factory()

    Can also be initialized simply by specifying the desired `style`.
    The below initializes a 'pluto' style `Server` with a servername
    'earth'.

        factory = ServerFactory(style='pluto', servername='earth')
        server = factory()
    """

    YAML_KEY = "SERVER"

    def __init__(self, servername, style, **kwargs):
        self.servername = servername
        self.style = style
        self.kwargs = kwargs

    def __call__(self):
        server_cls = [s for s in Server.subclasses() if self.style == s.STYLE]
        if len(server_cls) == 0:
            raise ValueError(f"No such server of style: {self.style}")

        if len(server_cls) > 1:
            raise RuntimeError(
                f"Multiple servers of style {self.style}: {server_cls}"
            )

        server_cls = server_cls[0]
        kwargs = self.kwargs.copy()

        scheduler = self._create_scheduler(kwargs)
        executor = self._create_covalent_executor(kwargs)

        try:
            return server_cls(
                servername=self.servername,
                scheduler=scheduler,
                covalent_executor=executor,
                **kwargs,
            )
        except TypeError as e:
            if "unexpected keyword argument" in str(e):
                raise TypeError(
                    f"Cannot create server {self.servername} due to unexpected keywords"
                ) from e

    def _create_scheduler(self, kwargs):
        if "scheduler" not in kwargs:
            return NullScheduler()

        factory = SchedulerFactory.from_dict(kwargs["scheduler"])
        scheduler = factory()
        del kwargs["scheduler"]
        return scheduler

    def _create_covalent_executor(self, kwargs):
        if "covalent_executor" not in kwargs:
            return None

        from pyatoms.settings.server.covalent import CovalentExecutorFactory

        factory = CovalentExecutorFactory.from_dict(
            kwargs["covalent_executor"]
        )
        executor = factory()
        del kwargs["covalent_executor"]
        return executor

    @classmethod
    def _other_yaml_kwargs(cls, filename):
        name = os.path.splitext(os.path.basename(filename))[0]
        return {"servername": name}


class Server(RegistryMixin, UserServerSettingsMixin, MSONable):
    """Stores information and settings for a Server.

    For example, knows about the correct commands for run ning executables.
    """

    SCRATCH_FOLDER = None
    MAX_MEM_GIGS = NotImplemented
    TMP_DIR = NotImplemented
    NUM_GPUS = NotImplemented
    NUM_CORES = NotImplemented
    NUM_THREADS = NotImplemented
    DEFAULT_NUM_HOURS = NotImplemented
    DEFAULT_QUEUE_NAME = NotImplemented
    USE_HOSTS = NotImplemented
    SUITABLE_SERVERS = NotImplemented
    FACTORY = ServerFactory

    CURRENT_SERVER_SENTINEL = "current"

    RESERVED_NAMES = [CURRENT_SERVER_SENTINEL]

    CUSTOM_REGISTERED_SERVERS = []

    def __init__(
        self,
        servername,
        scheduler,
        covalent_executor=None,
        hostnames=None,
        set_env=None,
    ):
        if hostnames is None:
            hostnames = []

        if set_env is None:
            set_env = {}

        if not isinstance(hostnames, list | tuple):
            raise ValueError(
                f"hostnames must be an iterable of associated hostnames. "
                f'Instead was: "{hostnames}". (Error for server={servername})'
            )

        if isinstance(scheduler, str):
            scheduler = Scheduler.from_name(scheduler)

        if servername in self.RESERVED_NAMES:
            raise ValueError(
                f"Server name {servername} cannot be in the reserved names ({self.RESERVED_NAMES})"
            )

        hostnames = tuple(hostnames)

        self.set_env = set_env
        self.servername = servername
        self.scheduler = scheduler
        self.hostnames = hostnames
        self.covalent_executor = covalent_executor

    def __repr__(self):
        return (
            f"{self.__class__.__qualname__}<servername={self.servername}, num_gpus={self.NUM_GPUS}, "
            f"num_cores={self.NUM_CORES}, scheduler={self.scheduler}, use_hosts={self.USE_HOSTS}>"
        )

    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented

        return self.__dict__ == other.__dict__

    @property
    def num_cores_running(self):
        # Subclasses can implement server specific overrides of NUM_CORES
        return self.NUM_CORES

    def matches_hostname(self):
        hostname = socket.gethostname()
        return any(re.match(h, hostname) for h in self.hostnames)

    def submitscript_writer(self, job, **kwargs):
        return self.scheduler.submitscript_writer(
            server=self, job=job, **kwargs
        )

    @property
    def scheduler_style(self):
        return self.scheduler.STYLE

    def get_hostnames(self):
        return self.scheduler.hostnames()

    def create_host_queue(self, **kwargs):
        if not self.USE_HOSTS:
            return None

        return self.scheduler.create_host_queue(**kwargs)

    def create_gpu_queue(self, **kwargs):
        if not self.USE_HOSTS or self.NUM_GPUS == 0:
            return None

        return self.scheduler.create_gpu_queue(**kwargs)

    @abstractmethod
    def parallel_command(
        self,
        folder,
        executable,
        resources,
        host,
        env_list,
        env_dict=None,
        mpirun=None,
        *args,
        **kwargs,
    ):
        raise NotImplementedError

    def serial_command(
        self, folder, executable, resources, env_list=None, env_dict=None
    ):
        return f"{executable}"

    def _submit_checks(self, **kwargs):
        """Can be subclassed to check for any configuration issues before job submission."""
        return

    def submit(self, job, scheduler=None, **kwargs):
        if scheduler is None:
            scheduler = self.scheduler

        self._submit_checks()
        scheduler.submit_job(job=job, server=self, **kwargs)

    def submit_job_script(self, submitscript):
        self._submit_checks()
        self.scheduler.submit_file(submitscript)

    def command(
        self,
        folder,
        executable,
        resources,
        host=None,
        env_list=None,
        env_dict=None,
        local=False,
        mpirun=None,
    ):
        """Generate the command for running in a subprocess.

        Args:
            folder (str): Folder to run the executable in.
            executable (str): Executable to run.
            resources (Resources): Resources to use for running the executable.
            host (Host): List of Hosts to run the executable on.
            local (bool): If True, will run the executable locally. Defaults to False.
            env_list (list): List of environment variables to export from the current environment.
            env_dict (dict): Dictionary of environment variables to set and export.
            mpirun (str): Command to use for running executables in parallel. Defaults to None, which will use
                the simple "mpirun" command.
        """
        if env_dict is None:
            env_dict = {}

        if set(self.set_env).intersection(env_dict):
            raise RuntimeError(
                f"Overlapping env vars in self.set_env ({self.set_env}) and env_dict ({env_dict})"
            )

        env_dict |= self.set_env

        if local or (resources.num_tasks == 0 and resources.num_gpus == 0):
            command = self.serial_command(
                folder=folder,
                executable=executable,
                resources=resources,
                env_list=env_list,
                env_dict=env_dict,
            )
        else:
            command = self.parallel_command(
                folder=folder,
                executable=executable,
                resources=resources,
                env_dict=env_dict,
                env_list=env_list,
                host=host,
                mpirun=mpirun,
            )
        return command

    @classmethod
    def current_server(cls):
        return cls.from_name(current_server())

    @classmethod
    def from_servername(cls, servername, **kwargs):
        return cls.from_name(servername, **kwargs)

    @classmethod
    def from_yaml(cls, filename):
        factory = cls.FACTORY.from_yaml(filename)
        return factory()

    @classmethod
    def from_registered_servers(cls, servername):
        for server in cls.CUSTOM_REGISTERED_SERVERS:
            if server.servername == servername:
                logger.info(
                    f"Found server {servername} in custom registered servers"
                )
                return server
        raise ValueError(f"Could not find server with name {servername}")

    @classmethod
    def from_name(cls, servername, return_none=True):
        """Create a server from the name.

        Args:
            servername (str): Name of server to create.
            return_none (bool): If True, returns None if server cannot be created from the name, else RuntimeError
                 will be raised. Defaults to True.
        """
        # return_none should be True here as we don't want to fail if there is no user-defined server settings,
        # there may be still prebuilt user settingsa
        if servername is cls.CURRENT_SERVER_SENTINEL:
            servername = current_server()

        server = (
            cls.from_user_settings(servername, return_none=True)
            or cls.from_prebuilt_classes(servername=servername)
            or cls.from_registered_servers(servername=servername)
        )

        if server is not None:
            logger.info(f"Got server: {server}")
            return server

        # Could not find project
        if return_none:
            return None

        all_fixed_servers = list(
            chain.from_iterable([r.SUITABLE_SERVERS for r in cls.subclasses()])
        )
        all_user_yaml = cls.user_servers()
        available_servers = all_fixed_servers + all_user_yaml
        raise ValueError(
            f"No settings implemented for {servername}. Configure pyatoms by:\n\n"
            f'    1.  Running "pyatoms config server" (RECOMMENDED) \n'
            f"    2.  Manually creating a new server.yaml file in {cls.handler_path()}.\n\n"
            f" Currently available servers: {available_servers}"
        )

    @classmethod
    def from_prebuilt_classes(cls, servername):
        # Hard coded classes
        for r in cls.subclasses():
            if servername in r.SUITABLE_SERVERS:
                return r(servername=servername)
        return None


class MpirunGenericServer(MpirunExecutionMixin, Server):
    SUITABLE_SERVERS = ["mpirun_generic"]
    STYLE = "mpirun_generic"

    MAX_MEM_GIGS = 32  # give conservative defaults, user can override
    NUM_THREADS = 4
    NUM_CORES = 4
    NUM_GPUS = 0
    TMP_DIR = None

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = NullScheduler()
        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )


class SlurmGenericServer(SlurmExecutionMixin, Server):
    SUITABLE_SERVERS = ["slurm_generic"]
    STYLE = "slurm_generic"

    MAX_MEM_GIGS = 32  # give conservative defaults, user can override
    NUM_THREADS = 4
    NUM_CORES = 4
    NUM_GPUS = 0
    TMP_DIR = None

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = SlurmScheduler()
        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )


class LocalServer(LocalExecutionMixin, Server):
    SUITABLE_SERVERS = ["local"]
    STYLE = "local"

    NUM_CORES = 4
    NUM_THREADS = 4
    NUM_GPUS = 0
    USE_HOSTS = False
    MAX_MEM_GIGS = 115
    TMP_DIR = None

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = NullScheduler()

        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )


class InpackageServer(LocalServer):
    """A separate server for using executables in pyatoms package.

    The executables for individual programs can be specified in
    jobs/<program>/execution.py under <program>Executable class.
    Separate from the one defined in <program>Executable for local server.
    """

    SUITABLE_SERVERS = ["inpackage"]
    STYLE = "inpackage"


class FermiServer(LocalServer):
    SUITABLE_SERVERS = ["fermi"]
    STYLE = "fermi"

    NUM_CORES = 16
    NUM_THREADS = 16
    NUM_GPUS = 0
    USE_HOSTS = False
    MAX_MEM_GIGS = 20
    TMP_DIR = os.path.expanduser("~/scratch")

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = NullScheduler()

        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )


class PlutoServer(Server):
    SUITABLE_SERVERS = ["pluto"]
    STYLE = "pluto"

    NUM_CORES = 48
    NUM_THREADS = 48
    USE_HOSTS = False
    NUM_GPUS = 0
    MAX_MEM_GIGS = 182
    TMP_DIR = None

    def __init__(
        self,
        servername=STYLE,
        hostnames=("hpc-login.*", "hpc-com.*", "hpc-accel.*"),
        scheduler=None,
        **kwargs,
    ):
        if scheduler is None:
            scheduler = SlurmScheduler(
                default_queue="normal", default_num_hours=72
            )

        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )

    @property
    def num_cores_running(self):
        # SLURM_NTASKS is only set if user specifically adds --ntasks-per-node (-n) to their sbatch script
        # If user does not add --ntasks-per-node, then use default server.NUM_CORES
        env = os.environ
        if "SLURM_NTASKS" in env and "SLURM_NNODES" in env:
            return int(int(env["SLURM_NTASKS"]) / int(env["SLURM_NNODES"]))
        return self.NUM_CORES

    def parallel_command(
        self,
        folder,
        executable,
        resources,
        host,
        env_list=None,
        env_dict=None,
        mpirun=None,
    ):
        if mpirun is None:
            mpirun = "srun"

        return f"{mpirun}  -D {folder} -N {resources.num_nodes} -n {resources.num_tasks} {executable}"


class PlutoMpirunServer(PlutoServer):
    SUITABLE_SERVERS = ["pluto-mpirun"]
    STYLE = "pluto-mpirun"

    def parallel_command(
        self,
        folder,
        executable,
        resources,
        host,
        num_gpus=None,
        env_list=None,
        env_dict=None,
        mpirun=None,
    ):
        # --bind-to none to allow application to use all cores if needed: https://stackoverflow.com/a/61010441
        if host:
            host_string = ",".join(host)
            host_string = f"-host {host_string}"
        else:
            host_string = ""

        if mpirun is None:
            mpirun = "mpirun"

        return f"{mpirun} -launcher ssh --bind-to none -wdir {folder} -np {resources.num_tasks} {host_string} {executable}"


class CoronaJobServer(MpirunExecutionMixin, Server):
    SUITABLE_SERVERS = ["corona"]
    STYLE = "corona"

    NUM_CORES = 40
    NUM_THREADS = 40
    NUM_GPUS = 0
    USE_HOSTS = True
    MAX_MEM_GIGS = 192
    TMP_DIR = "/gpfs1/users/astar/ihpc/chenwjb/tmp"

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = PBSScheduler(
                default_queue="normal", default_num_hours=72, software="vasp"
            )

        if hostnames is None:
            hostnames = ["corona"]
        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )


class StratusJobServer(MpirunExecutionMixin, Server):
    SUITABLE_SERVERS = ["stratus"]
    STYLE = "stratus"

    NUM_CORES = 64
    NUM_THREADS = 64
    NUM_GPUS = 0
    USE_HOSTS = True
    MAX_MEM_GIGS = 160
    TMP_DIR = "/gpfs1/users/astar/ihpc/chenwjb/tmp"

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = PBSScheduler(
                default_queue="normal", default_num_hours=72, software="vasp"
            )

        if hostnames is None:
            hostnames = ["stratus", "r0.*", "r1.*"]
        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )


class EuclidJobServer(Server):
    SUITABLE_SERVERS = ["euclid"]
    STYLE = "euclid"

    NUM_CORES = 16
    NUM_THREADS = 16
    NUM_GPUS = 0
    USE_HOSTS = True
    MAX_MEM_GIGS = 60
    TMP_DIR = None

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = SlurmScheduler(
                default_queue="long", default_num_hours=168
            )

        if hostnames is None:
            hostnames = ["eu.*"]
        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )

    def parallel_command(
        self,
        folder,
        executable,
        resources,
        host,
        env_list=None,
        env_dict=None,
        mpirun=None,
    ):
        if mpirun is None:
            mpirun = "mpirun"

        if host:
            host_string = ",".join(host)

            return f"{mpirun} -wdir {folder} -host {host_string} {executable}"
        return f"{mpirun}  -wdir {folder} -perhost 16 {executable}"


class ThalesJobServer(Server):
    SUITABLE_SERVERS = ["thales"]
    STYLE = "thales"

    NUM_CORES = 16
    NUM_THREADS = 16
    NUM_GPUS = 0
    USE_HOSTS = True
    MAX_MEM_GIGS = 60
    TMP_DIR = None

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = SlurmScheduler(
                default_queue="long", default_num_hours=168
            )

        if hostnames is None:
            hostnames = ["th.*"]
        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )

    def parallel_command(
        self,
        folder,
        executable,
        resources,
        host,
        env_list=None,
        env_dict=None,
        mpirun=None,
    ):
        if mpirun is None:
            mpirun = "mpirun"

        if host:
            host_string = ",".join(host)
            return f"{mpirun} -wdir {folder} -host {host_string} {executable}"
        return f"{mpirun}  -wdir {folder} -perhost 16 {executable}"


class AndersenJobServer(Server):
    SUITABLE_SERVERS = ["andersen"]
    STYLE = "andersen"

    NUM_CORES = 16
    NUM_THREADS = 32
    NUM_GPUS = 0
    USE_HOSTS = True
    MAX_MEM_GIGS = 256
    TMP_DIR = None

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = SlurmScheduler(
                default_queue="normal", default_num_hours=72
            )
        if hostnames is None:
            hostnames = ["andersen"]
        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )

    def parallel_command(
        self,
        folder,
        executable,
        resources,
        host,
        env_list=None,
        env_dict=None,
        mpirun=None,
    ):
        """Get the command for running the executable in a subprocess with mpirun.

        -mca btl ^openib prevents the below warning:
        WARNING: There was an error initializing an OpenFabrics device.
        """
        if host:
            # Sorting the hosts could prevent the below fatal error:
            # Failed to receive UCX worker address: Not found
            host_string = ",".join(sorted(host))
            host_string = f"-host {host_string}"
        else:
            host_string = ""

        if mpirun is None:
            mpirun = "mpirun"

        return (
            f"{mpirun} -wdir {folder} -mca btl ^openib -npernode 16 --oversubscribe "
            f"-np {resources.num_tasks} {host_string} {executable}"
        )


class HelixServer(MpirunExecutionMixin, Server):
    """Helix server does not have scheduler, but we inject PBSScheduler for testing purposes only."""

    SUITABLE_SERVERS = ["helix"]
    STYLE = "helix"

    NUM_CORES = 12
    NUM_THREADS = 12
    NUM_GPUS = 1
    USE_HOSTS = False
    MAX_MEM_GIGS = 256
    TMP_DIR = None

    def __init__(
        self, servername=STYLE, scheduler=None, hostnames=None, **kwargs
    ):
        if scheduler is None:
            scheduler = PBSScheduler()

        if hostnames is None:
            hostnames = ["helix", "IHP-SPD-E0M273"]
        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )


class HelixIntelServer(LocalServer):
    SUITABLE_SERVERS = ["helix_intel"]
    STYLE = "helix_intel"

    NUM_CORES = 12
    NUM_THREADS = 12
    NUM_GPUS = 0
    USE_HOSTS = False
    MAX_MEM_GIGS = 256
    TMP_DIR = None

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = PBSScheduler()

        if hostnames is None:
            hostnames = ["helix", "IHP-SPD-E0M273"]

        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )

    def parallel_command(
        self,
        folder,
        executable,
        resources,
        host,
        env_list=None,
        env_dict=None,
        mpirun=None,
    ):
        host_string = ""
        if host:
            host_string = f'-host {",".join(host)}'

        # --bind-to none to allow application to use all cores if needed: https://stackoverflow.com/a/61010441
        bind_to_none = "--bind-to none"

        if mpirun is None:
            mpirun = "/home/benjamin/intel/parallel_studio_xe_2020/compilers_and_libraries_2020/linux/mpi/intel64/bin/mpirun"

        return f"{mpirun} -wdir {folder} -np {resources.num_tasks} {bind_to_none} {host_string} {executable}"


class Bridges2Server(SlurmExecutionMixin, Server):
    SUITABLE_SERVERS = ["bridges2"]
    STYLE = "bridges2"

    NUM_CORES = 64
    NUM_THREADS = 64
    NUM_GPUS = 0
    USE_HOSTS = False
    SCRATCH_FOLDER = os.path.abspath("/tmp")
    MAX_MEM_GIGS = 64
    TMP_DIR = None

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if scheduler is None:
            scheduler = SlurmScheduler(
                default_queue="normal", default_num_hours=72
            )

        if hostnames is None:
            hostnames = ["bridges2"]
        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )


class GekkoJobServer(MpirunExecutionMixin, Server):
    SUITABLE_SERVERS = ["gekko"]
    STYLE = "gekko"

    NUM_CORES = 64
    NUM_THREADS = 64
    NUM_GPUS = 0
    USE_HOSTS = False
    MAX_MEM_GIGS = 96
    TMP_DIR = "/tmp"

    def __init__(
        self, servername=STYLE, hostnames=None, scheduler=None, **kwargs
    ):
        if hostnames is None:
            hostnames = ["gekko"]

        if scheduler is None:
            from pyatoms.settings.user import PyatomsUserSettings

            user_settings = PyatomsUserSettings()
            project = user_settings.get("PYATOMS_PROJECT_GROUP", None)
            scheduler = PBSScheduler(
                default_queue="q256_free", default_num_hours=8, project=project
            )

        super().__init__(
            servername=servername,
            hostnames=hostnames,
            scheduler=scheduler,
            **kwargs,
        )

    def _submit_checks(self):
        if self.scheduler.project is None:
            raise ValueError(
                "Specify PROJECT_GROUP in ~/.pyatoms/pyatoms.yaml"
            )


class CustomServer(Server):
    """Custom server for user customization. Usually created from yaml files."""

    SUITABLE_SERVERS = []
    STYLE = "custom"

    def __init__(
        self,
        servername,
        num_cores,
        num_threads,
        max_mem_gigs,
        execution_type,
        scheduler,
        use_hosts=False,
        tmp_dir=None,
        hostnames=None,
        num_gpus=0,
        default_num_hours=12,
        default_queue_name=None,
        covalent_executor=None,
        set_env=None,
    ):
        super().__init__(
            servername=servername,
            scheduler=scheduler,
            hostnames=hostnames,
            covalent_executor=covalent_executor,
            set_env=set_env,
        )

        self.num_cores = num_cores
        self.num_threads = num_threads
        self.num_gpus = num_gpus
        self.default_num_hours = default_num_hours
        self.default_queue_name = default_queue_name
        self.use_hosts = use_hosts
        self.max_mem_gigs = max_mem_gigs
        self.tmp_dir = tmp_dir
        self.execution_type = execution_type

    def register(self):
        """Register server in the list of custom servers. Useful for testing."""
        if self in self.CUSTOM_REGISTERED_SERVERS:
            raise ValueError(f"{self} already registered.")

        return self.CUSTOM_REGISTERED_SERVERS.append(self)

    def unregister(self):
        """Unregister server in the list of custom servers. Useful for testing."""
        if self not in self.CUSTOM_REGISTERED_SERVERS:
            raise ValueError(f"{self} not registered.")

        return self.CUSTOM_REGISTERED_SERVERS.remove(self)

    @property
    def NUM_CORES(self):
        return self.num_cores

    @property
    def NUM_THREADS(self):
        return self.num_threads

    @property
    def NUM_GPUS(self):
        return self.num_gpus

    @property
    def DEFAULT_NUM_HOURS(self):
        return self.default_num_hours

    @property
    def DEFAULT_QUEUE_NAME(self):
        return self.default_queue_name

    @property
    def USE_HOSTS(self):
        return self.use_hosts

    @property
    def MAX_MEM_GIGS(self):
        return self.max_mem_gigs

    @property
    def TMP_DIR(self):
        return self.tmp_dir

    @property
    def EXECUTION_TYPE(self):
        return self.execution_type

    def parallel_command(self, *args, **kwargs):
        new = combine_with_mixin(
            self, ExecutionMixin.class_from_name(self.EXECUTION_TYPE)
        )
        return new.parallel_command(*args, **kwargs)

    def submit_file(self, *args, **kwargs):
        return self.scheduler.submit_file(*args, **kwargs)

    @classmethod
    def _other_yaml_kwargs(cls, filename):
        name = os.path.basename(filename).removesuffix(".yaml")
        return {"name": name}


def combine_with_mixin(obj, mixin):
    """Dynamically adds a mixin to an object. Returns a copy of the object."""
    obj_copy = copy.copy(obj)

    # Don't register this dynamic class
    obj_copy.__class__ = type(
        obj.__class__.__name__, (mixin, type(obj)), {"REGISTERABLE": False}
    )
    return obj_copy


@lru_cache
def current_server():
    """Guesses the current server based on its hostname.

    Looks up the server's hostname in a predefined dictionary of servernames, including
    user defined servers in ~/.pyatoms/server

    Returns the `name` of the current server as to be used with the `Server` class.
    """

    def _detect_server():
        from pyatoms.settings.user import UserServerYamlManager

        manager = UserServerYamlManager()

        # User servers may not have server information (only having VASP information for eg),
        # so we need return_none=True
        user_servers = [
            Server.from_name(name, return_none=True)
            for name in manager.user_servers()
        ]
        logger.debug(f"User servers: {user_servers}")
        hardcoded_servers = [
            s() for s in Server.subclasses() if s is not CustomServer
        ]

        servers = user_servers + hardcoded_servers
        servers = [s for s in servers if s is not None]
        matching_servers = [
            server for server in servers if server.matches_hostname()
        ]

        server = (
            LocalServer()
            if len(matching_servers) == 0
            else matching_servers[0]
        )
        return server.servername

    from pyatoms.settings.user import PyatomsUserSettings

    user_settings = PyatomsUserSettings()
    current_server = user_settings.default_server() or _detect_server()
    logger.info(f"Current server detected as: {current_server}")
    return current_server
