import os
from io import StringIO

import pytest

from chemsmart.jobs.gaussian.batch import GaussianBatchJob
from chemsmart.settings.executable import GaussianExecutable, ORCAExecutable
from chemsmart.settings.server import Server
from chemsmart.settings.submitters import (
    PBSSubmitter,
    SLFSubmitter,
    SLURMSubmitter,
)


class TestServer:
    def test_server_yaml(self, server_yaml_file):
        assert os.path.exists(server_yaml_file)
        assert os.path.isfile(server_yaml_file)
        server = Server.from_yaml(name=server_yaml_file)
        assert server.scheduler.lower() == "pbs"
        assert server.queue_name == "normal"
        assert server.num_hours == 24
        assert server.mem_gb == 375
        assert server.num_cores == 64
        assert server.num_gpus == 0
        assert server.num_threads == 64
        assert server.submit_command == "qsub"
        assert server.scratch_dir is None
        assert server.use_hosts is True
        assert (
            server.extra_commands == """export PATH=$HOME/bin/chemsmart:$PATH
export PATH=$HOME/bin/chemsmart/chemsmart/cli:$PATH
export PATH=$HOME/bin/chemsmart/chemsmart/scripts:$PATH
export PYTHONPATH=$HOME/bin/chemsmart:$PYTHONPATH
"""
        )
        assert server.extra_scheduler_directives == "#PBS -m abe\n"

    def test_gaussian_executable(self, server_yaml_file):
        gaussian_executable = GaussianExecutable.from_servername(
            server_yaml_file
        )
        assert gaussian_executable.executable_folder == os.path.expanduser(
            "~/programs/g16"
        )
        assert gaussian_executable.local_run is True

        gaussian_conda_env = """source ~/anaconda3/etc/profile.d/conda.sh
conda activate ~/anaconda3/envs/chemsmart
"""
        assert gaussian_executable.conda_env == gaussian_conda_env

        gaussian_modules = """module purge
module load craype-x86-rome
module load libfabric/1.11.0.4.125
"""
        assert gaussian_executable.modules == gaussian_modules

        assert (
            gaussian_executable.scripts
            == 'tcsh -c "source ~/programs/g16/bsd/g16.login"\n'
        )

        gassian_envars = """export SCRATCH=~/scratch
export GAUSS_EXEDIR=~/programs/g16
export g16root=~/programs/g16

"""
        assert gaussian_executable.envars == gassian_envars

    def test_orca_executable(self, server_yaml_file):
        orca_executable = ORCAExecutable.from_servername(server_yaml_file)
        assert orca_executable.executable_folder == os.path.expanduser(
            "~/programs/orca_6_0_0"
        )
        assert orca_executable.local_run is False

        assert orca_executable.conda_env is None

        assert orca_executable.modules is None

        assert orca_executable.scripts is None

        orca_envars = """export PATH=~/programs/openmpi-4.1.6/build/bin:$PATH
export LD_LIBRARY_PATH=~/programs/openmpi-4.1.6/build/lib:$LD_LIBRARY_PATH
"""
        assert orca_executable.envars == orca_envars

    def test_slurm_submitter_writes_extra_scheduler_directives(self):
        server = Server(
            "custom-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=24,
            NUM_GPUS=0,
            EXTRA_SCHEDULER_DIRECTIVES="#SBATCH --reservation=xlzhang_1\n",
        )
        job = type("DummyJob", (), {"label": "job1"})()
        submitter = SLURMSubmitter(job=job, server=server)

        buffer = StringIO()
        submitter._write_scheduler_options(buffer)
        assert "#SBATCH --reservation=xlzhang_1\n" in buffer.getvalue()

    def test_pbs_submitter_writes_extra_scheduler_directives(self):
        server = Server(
            "custom-pbs",
            SCHEDULER="PBS",
            NUM_CORES=8,
            MEM_GB=24,
            NUM_GPUS=0,
            EXTRA_SCHEDULER_DIRECTIVES="#PBS -m abe\n",
        )
        job = type("DummyJob", (), {"label": "job1"})()
        submitter = PBSSubmitter(job=job, server=server)

        buffer = StringIO()
        submitter._write_scheduler_options(buffer)
        assert "#PBS -m abe\n" in buffer.getvalue()


class TestArraySubmitInfrastructure:
    """Phase 2.2: array submit script writing (1-based task ids)."""

    @staticmethod
    def _stub_program_sections(monkeypatch):
        monkeypatch.setattr(
            SLURMSubmitter,
            "_write_program_specifics",
            lambda self, f: None,
        )
        monkeypatch.setattr(
            SLURMSubmitter,
            "_write_extra_commands",
            lambda self, f: None,
        )

    def test_slurm_array_submit_script_name_uses_batch_label(self):
        server = Server(
            "array-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=16,
            MEM_GB=32,
            NUM_GPUS=0,
        )
        child = type("Child", (), {"label": "child1", "PROGRAM": "gaussian"})()
        submitter = SLURMSubmitter(job=child, server=server)
        submitter.batch_label = "mols_batch"
        assert (
            submitter.array_submit_script
            == "chemsmart_sub_array_mols_batch.sh"
        )

    def test_slurm_write_array_job_creates_1based_scripts(
        self, tmp_path, monkeypatch
    ):
        self._stub_program_sections(monkeypatch)
        monkeypatch.chdir(tmp_path)
        server = Server(
            "array-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=16,
            MEM_GB=32,
            NUM_GPUS=0,
            NUM_HOURS=12,
            QUEUE_NAME="normal",
        )
        children = [
            type(
                "Child",
                (),
                {
                    "label": f"mol{i}",
                    "PROGRAM": "gaussian",
                    "folder": str(tmp_path),
                },
            )()
            for i in range(1, 5)
        ]
        template = children[0]
        submitter = SLURMSubmitter(job=template, server=server)
        shared_cli = ["gaussian", "-f", "mols.xyz", "-i", "1,2,3,4", "opt"]

        submitter.write_array_job(
            jobs=children,
            array_concurrency=2,
            cli_args=shared_cli,
            batch_label="mols_batch",
        )

        assert (tmp_path / "chemsmart_sub_array_mols_batch.sh").exists()
        for task_id in range(1, 5):
            assert (tmp_path / f"chemsmart_run_array_{task_id}.py").exists()
        assert not (tmp_path / "chemsmart_run_array_0.py").exists()

        submit_text = (
            tmp_path / "chemsmart_sub_array_mols_batch.sh"
        ).read_text()
        assert "#SBATCH --array=1-4%2\n" in submit_text
        assert "--nodes=1 --ntasks-per-node=16 --mem=32G" in submit_text
        assert "TASK_ID=$SLURM_ARRAY_TASK_ID" in submit_text
        assert "SLURM_ARRAY_TASK_ID + 1" not in submit_text
        assert "python chemsmart_run_array_${TASK_ID}.py" in submit_text

        run_text = (tmp_path / "chemsmart_run_array_2.py").read_text()
        assert "gaussian" in run_text
        assert "mols.xyz" in run_text

    def test_submit_array_job_test_mode_writes_without_queueing(
        self, tmp_path, monkeypatch
    ):
        self._stub_program_sections(monkeypatch)
        monkeypatch.chdir(tmp_path)
        server = Server(
            "array-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=16,
            NUM_GPUS=0,
            SUBMIT_COMMAND="sbatch",
        )
        monkeypatch.setattr(
            Server,
            "_check_running_jobs",
            staticmethod(lambda job: None),
        )
        submitted = []

        def _fake_submit(self, job, submitter):
            submitted.append(submitter.array_submit_script)

        monkeypatch.setattr(Server, "_submit_array_job", _fake_submit)

        children = [
            type(
                "Child",
                (),
                {
                    "label": "a",
                    "PROGRAM": "gaussian",
                    "folder": str(tmp_path),
                },
            )(),
            type(
                "Child",
                (),
                {
                    "label": "b",
                    "PROGRAM": "gaussian",
                    "folder": str(tmp_path),
                },
            )(),
        ]
        server.submit_array_job(
            jobs=children,
            array_concurrency=1,
            test=True,
            cli_args=["gaussian", "opt"],
            batch_label="pka_batch",
        )

        assert (tmp_path / "chemsmart_sub_array_pka_batch.sh").exists()
        assert (tmp_path / "chemsmart_run_array_1.py").exists()
        assert (tmp_path / "chemsmart_run_array_2.py").exists()
        assert submitted == []

        submit_text = (
            tmp_path / "chemsmart_sub_array_pka_batch.sh"
        ).read_text()
        assert "#SBATCH --array=1-2%1\n" in submit_text

    @staticmethod
    def _array_children(tmp_path, count=4):
        return [
            type(
                "Child",
                (),
                {
                    "label": f"mol{i}",
                    "PROGRAM": "gaussian",
                    "folder": str(tmp_path),
                },
            )()
            for i in range(1, count + 1)
        ]

    def test_pbs_write_array_job_creates_array_directives(
        self, tmp_path, monkeypatch
    ):
        self._stub_program_sections(monkeypatch)
        monkeypatch.setattr(
            PBSSubmitter,
            "_write_program_specifics",
            lambda self, f: None,
        )
        monkeypatch.setattr(
            PBSSubmitter,
            "_write_extra_commands",
            lambda self, f: None,
        )
        monkeypatch.chdir(tmp_path)
        server = Server(
            "array-pbs",
            SCHEDULER="PBS",
            NUM_CORES=8,
            MEM_GB=16,
            NUM_GPUS=0,
            NUM_HOURS=12,
            QUEUE_NAME="normal",
        )
        children = self._array_children(tmp_path)
        submitter = PBSSubmitter(job=children[0], server=server)

        submitter.write_array_job(
            jobs=children,
            array_concurrency=2,
            cli_args=["gaussian", "opt"],
            batch_label="mols_batch",
        )

        submit_text = (
            tmp_path / "chemsmart_sub_array_mols_batch.sh"
        ).read_text()
        assert "#PBS -J 1-4%2\n" in submit_text
        assert (
            "#PBS -o mols_batch_array_${PBS_ARRAYID}.pbsout\n" in submit_text
        )
        assert "TASK_ID=$PBS_ARRAYID" in submit_text
        assert "python chemsmart_run_array_${TASK_ID}.py" in submit_text

    def test_lsf_write_array_job_creates_array_directives(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.setattr(
            SLFSubmitter,
            "_write_program_specifics",
            lambda self, f: None,
        )
        monkeypatch.setattr(
            SLFSubmitter,
            "_write_extra_commands",
            lambda self, f: None,
        )
        monkeypatch.chdir(tmp_path)
        server = Server(
            "array-lsf",
            SCHEDULER="SLF",
            NUM_CORES=8,
            MEM_GB=16,
            NUM_GPUS=0,
            NUM_HOURS=12,
            NUM_NODES=1,
        )
        children = self._array_children(tmp_path, count=3)
        submitter = SLFSubmitter(job=children[0], server=server)

        submitter.write_array_job(
            jobs=children,
            array_concurrency=1,
            cli_args=["gaussian", "opt"],
            batch_label="pka_batch",
        )

        submit_text = (
            tmp_path / "chemsmart_sub_array_pka_batch.sh"
        ).read_text()
        assert "#BSUB -J pka_batch_array[1-3%1]\n" in submit_text
        assert "#BSUB -o pka_batch_array_%I.bsubout\n" in submit_text
        assert "TASK_ID=$LSB_JOBINDEX" in submit_text


class TestSchedulerArrayPolicy:
    def test_no_run_in_parallel_forces_throttle_one(self):
        from chemsmart.settings.server import SchedulerArrayPolicy

        policy = SchedulerArrayPolicy(
            no_run_in_parallel=True, array_concurrency=4, max_concurrent=8
        )
        assert policy.array_throttle(10) == 1

    def test_array_concurrency_preferred_over_max_concurrent(self):
        from chemsmart.settings.server import SchedulerArrayPolicy

        policy = SchedulerArrayPolicy(
            no_run_in_parallel=False, array_concurrency=3, max_concurrent=8
        )
        assert policy.array_throttle(10) == 3

    def test_max_concurrent_caps_num_jobs(self):
        from chemsmart.settings.server import SchedulerArrayPolicy

        policy = SchedulerArrayPolicy(
            no_run_in_parallel=False, array_concurrency=None, max_concurrent=2
        )
        assert policy.array_throttle(10) == 2

    def test_parallel_without_explicit_n_runs_all_tasks(self):
        from chemsmart.settings.server import SchedulerArrayPolicy

        policy = SchedulerArrayPolicy(
            no_run_in_parallel=False,
            array_concurrency=None,
            max_concurrent=None,
        )
        assert policy.array_throttle(10) == 10

    def test_from_jobrunner_uses_cli_array_concurrency_not_server_default(
        self, pbs_server
    ):
        from chemsmart.jobs.runner import JobRunner
        from chemsmart.settings.server import SchedulerArrayPolicy

        runner = JobRunner(
            server=pbs_server,
            fake=True,
            no_run_in_parallel=False,
            num_cores=64,
        )
        policy = SchedulerArrayPolicy.from_jobrunner(runner)
        assert policy.array_concurrency is None
        assert policy.array_throttle(5) == 5

    def test_from_jobrunner(self, pbs_server):
        from chemsmart.jobs.runner import JobRunner
        from chemsmart.settings.server import SchedulerArrayPolicy

        runner = JobRunner(
            server=pbs_server,
            fake=True,
            no_run_in_parallel=True,
            array_concurrency=4,
            num_cores=8,
        )
        policy = SchedulerArrayPolicy.from_jobrunner(runner)
        assert policy.no_run_in_parallel is True
        assert policy.array_concurrency == 4
        assert policy.array_throttle(10) == 1


class TestSubmitBatch:
    def test_submit_batch_delegates_to_submit_array_job(
        self, tmp_path, monkeypatch
    ):
        from chemsmart.jobs.gaussian.batch import GaussianBatchJob
        from chemsmart.settings.server import (
            SchedulerArrayPolicy,
            Server,
        )

        monkeypatch.chdir(tmp_path)
        server = Server(
            "batch-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=16,
            NUM_GPUS=0,
            SUBMIT_COMMAND="sbatch",
        )
        monkeypatch.setattr(
            Server,
            "_check_running_jobs",
            staticmethod(lambda job: None),
        )
        captured = {}

        def _fake_submit_array(
            self,
            jobs,
            array_concurrency=None,
            test=False,
            cli_args=None,
            batch_label=None,
            **kwargs,
        ):
            captured.update(
                {
                    "jobs": list(jobs),
                    "array_concurrency": array_concurrency,
                    "test": test,
                    "cli_args": cli_args,
                    "batch_label": batch_label,
                }
            )

        monkeypatch.setattr(Server, "submit_array_job", _fake_submit_array)

        children = [
            type(
                "Child",
                (),
                {
                    "label": "a",
                    "PROGRAM": "gaussian",
                    "folder": str(tmp_path),
                },
            )(),
            type(
                "Child",
                (),
                {
                    "label": "b",
                    "PROGRAM": "gaussian",
                    "folder": str(tmp_path),
                },
            )(),
        ]
        batch = GaussianBatchJob(jobs=children, label="mols_batch")
        server.submit_batch(
            batch,
            policy=SchedulerArrayPolicy(no_run_in_parallel=True),
            test=True,
            cli_args=["gaussian", "-f", "mols.xyz", "-i", "1,2", "opt"],
        )

        assert captured["batch_label"] == "mols_batch"
        assert captured["array_concurrency"] == 1
        assert captured["test"] is True
        assert len(captured["jobs"]) == 2
        assert captured["cli_args"] == [
            "gaussian",
            "-f",
            "mols.xyz",
            "-i",
            "1,2",
            "opt",
        ]


class TestCheckRunningJobs:
    class _MockClusterHelper:
        running_job_names = []

        def get_gaussian_running_jobs(self):
            return [], self.running_job_names

    def test_rejects_duplicate_batch_job_label(self, monkeypatch):
        self._MockClusterHelper.running_job_names = ["pka_batch"]
        monkeypatch.setattr(
            "chemsmart.utils.cluster.ClusterHelper",
            self._MockClusterHelper,
        )

        batch = GaussianBatchJob(jobs=[], label="pka_batch")
        with pytest.raises(
            SystemExit, match="Duplicate job NOT submitted: pka_batch"
        ):
            Server._check_running_jobs(batch)

    def test_allows_unique_batch_job_label(self, monkeypatch):
        self._MockClusterHelper.running_job_names = ["other_batch"]
        monkeypatch.setattr(
            "chemsmart.utils.cluster.ClusterHelper",
            self._MockClusterHelper,
        )

        batch = GaussianBatchJob(jobs=[], label="pka_batch")
        Server._check_running_jobs(batch)

    def test_batch_job_checks_container_label_not_children(self, monkeypatch):
        self._MockClusterHelper.running_job_names = ["acid1_pka"]
        monkeypatch.setattr(
            "chemsmart.utils.cluster.ClusterHelper",
            self._MockClusterHelper,
        )

        child = type("ChildJob", (), {"label": "acid1_pka"})()
        batch = GaussianBatchJob(jobs=[child], label="acids_pka_batch")
        Server._check_running_jobs(batch)

    def test_skips_jobs_without_scheduler_label(self, monkeypatch):
        def _fail_if_called(self):
            raise AssertionError("cluster query should not run")

        monkeypatch.setattr(
            self._MockClusterHelper,
            "get_gaussian_running_jobs",
            _fail_if_called,
        )
        monkeypatch.setattr(
            "chemsmart.utils.cluster.ClusterHelper",
            self._MockClusterHelper,
        )

        batch = GaussianBatchJob(jobs=[], label=None)
        Server._check_running_jobs(batch)
