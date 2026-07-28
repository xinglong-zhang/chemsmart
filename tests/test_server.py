import os
from io import StringIO

import pytest

from chemsmart.jobs.gaussian.batch import GaussianBatchJob
from chemsmart.settings.executable import GaussianExecutable, ORCAExecutable
from chemsmart.settings.server import Server
from chemsmart.settings.submitters import (
    FUGAKUSubmitter,
    PBSSubmitter,
    SLFSubmitter,
    SLURMSubmitter,
    _without_mail_directives,
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
        assert (tmp_path / "chemsmart_run_array_mols_batch_array.py").exists()
        for task_id in range(1, 5):
            assert not (
                tmp_path / f"chemsmart_run_array_{task_id}.py"
            ).exists()

        submit_text = (
            tmp_path / "chemsmart_sub_array_mols_batch.sh"
        ).read_text()
        assert "#SBATCH --array=1-4%2\n" in submit_text
        assert "--nodes=1 --ntasks-per-node=16 --mem=32G" in submit_text
        assert "#SBATCH --job-name=mols_batch_array\n" in submit_text
        assert "#SBATCH --output=mols_batch_array.slurmout\n" in submit_text
        assert "#SBATCH --error=mols_batch_array.slurmerr\n" in submit_text
        assert "#SBATCH --open-mode=append\n" in submit_text
        assert "%a.slurmout" not in submit_text
        assert "TASK_ID=$SLURM_ARRAY_TASK_ID" in submit_text
        assert "SLURM_ARRAY_TASK_ID + 1" not in submit_text
        assert "python chemsmart_run_array_mols_batch_array.py" in submit_text
        assert "python chemsmart_run_array_${TASK_ID}.py" not in submit_text
        assert "===== BEGIN array task ${TASK_ID} =====" in submit_text
        assert "===== END array task ${TASK_ID} (exit=${status}) =====" in (
            submit_text
        )
        assert "flock 9 || exit 1" in submit_text
        assert ".mols_batch_array.loglock" in submit_text

        run_text = (
            tmp_path / "chemsmart_run_array_mols_batch_array.py"
        ).read_text()
        assert "TASK_CLI" in run_text
        assert "gaussian" in run_text
        assert "mols.xyz" in run_text
        assert "SLURM_ARRAY_TASK_ID" in run_text
        assert "1:" in run_text and "4:" in run_text

    def test_slurm_write_array_job_omits_throttle_when_concurrency_none(
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
            NUM_HOURS=12,
            QUEUE_NAME="normal",
        )
        children = self._array_children(tmp_path, count=3)
        submitter = SLURMSubmitter(job=children[0], server=server)

        submitter.write_array_job(
            jobs=children,
            array_concurrency=None,
            cli_args=["gaussian", "opt"],
            batch_label="no_throttle",
        )

        submit_text = (
            tmp_path / "chemsmart_sub_array_no_throttle.sh"
        ).read_text()
        assert "#SBATCH --array=1-3\n" in submit_text
        assert "%3" not in submit_text
        assert "#SBATCH --array=1-3%" not in submit_text

    def test_pbs_write_array_job_omits_throttle_when_concurrency_none(
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
        children = self._array_children(tmp_path, count=2)
        submitter = PBSSubmitter(job=children[0], server=server)

        submitter.write_array_job(
            jobs=children,
            array_concurrency=None,
            cli_args=["gaussian", "opt"],
            batch_label="pbs_no_throttle",
        )

        submit_text = (
            tmp_path / "chemsmart_sub_array_pbs_no_throttle.sh"
        ).read_text()
        assert "#PBS -J 1-2\n" in submit_text
        assert "%2" not in submit_text

    def test_write_array_job_empty_jobs_warns_without_scripts(
        self, tmp_path, monkeypatch, caplog
    ):
        import logging

        self._stub_program_sections(monkeypatch)
        monkeypatch.chdir(tmp_path)
        server = Server(
            "array-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=16,
            NUM_GPUS=0,
        )
        submitter = SLURMSubmitter(
            job=type("Child", (), {"label": "x", "PROGRAM": "gaussian"})(),
            server=server,
        )

        with caplog.at_level(logging.WARNING):
            submitter.write_array_job(jobs=[], batch_label="empty_batch")

        assert "No jobs provided for array job" in caplog.text
        assert not (tmp_path / "chemsmart_sub_array_empty_batch.sh").exists()

    def test_write_array_job_rejects_unsupported_scheduler(self):
        server = Server(
            "fugaku",
            SCHEDULER="FUGAKU",
            NUM_CORES=8,
            MEM_GB=16,
            NUM_GPUS=0,
        )
        submitter = FUGAKUSubmitter(
            job=type("Child", (), {"label": "x", "PROGRAM": "gaussian"})(),
            server=server,
        )
        with pytest.raises(
            ValueError, match="Batch array submission is not supported"
        ):
            submitter.write_array_job(
                jobs=[submitter.job],
                batch_label="unsupported",
            )

    def test_without_mail_directives_strips_scheduler_mail_lines(self):
        assert (
            _without_mail_directives(
                "#PBS -m abe\n#PBS -N job\n#PBS -M user@ex.com\n"
            )
            == "#PBS -N job\n"
        )
        assert (
            _without_mail_directives(
                [
                    "#SBATCH --mail-type=ALL",
                    "#SBATCH --job-name=x",
                    "#SBATCH --mail-user=a@b.c",
                ]
            )
            == "#SBATCH --job-name=x\n"
        )
        assert _without_mail_directives("#PBS -m abe\n") is None
        assert _without_mail_directives(None) is None

    def test_slurm_write_array_job_embeds_per_job_cli_args(
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
        children = self._array_children(tmp_path, count=2)
        submitter = SLURMSubmitter(job=children[0], server=server)
        per_job_cli = [
            ["gaussian", "-f", "1a.xyz", "-l", "1a", "pka", "submit"],
            ["gaussian", "-f", "2a.xyz", "-l", "2a", "pka", "submit"],
        ]

        submitter.write_array_job(
            jobs=children,
            array_concurrency=2,
            cli_args=per_job_cli,
            batch_label="pka_batch",
        )

        run_text = (
            tmp_path / "chemsmart_run_array_pka_batch_array.py"
        ).read_text()
        assert "1a.xyz" in run_text
        assert "2a.xyz" in run_text
        assert "TASK_CLI" in run_text
        submit_text = (
            tmp_path / "chemsmart_sub_array_pka_batch.sh"
        ).read_text()
        assert "python chemsmart_run_array_pka_batch_array.py" in submit_text

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
        assert (tmp_path / "chemsmart_run_array_pka_batch_array.py").exists()
        assert not (tmp_path / "chemsmart_run_array_1.py").exists()
        assert not (tmp_path / "chemsmart_run_array_2.py").exists()
        assert submitted == []

        submit_text = (
            tmp_path / "chemsmart_sub_array_pka_batch.sh"
        ).read_text()
        assert "#SBATCH --array=1-2%1\n" in submit_text
        assert "python chemsmart_run_array_pka_batch_array.py" in submit_text

        run_text = (
            tmp_path / "chemsmart_run_array_pka_batch_array.py"
        ).read_text()
        assert "TASK_CLI" in run_text
        assert "1:" in run_text and "2:" in run_text
        assert "gaussian" in run_text
        assert "opt" in run_text

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
            EXTRA_SCHEDULER_DIRECTIVES="#PBS -m abe\n",
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
        assert "#PBS -o mols_batch_array.pbsout\n" in submit_text
        assert "#PBS -e mols_batch_array.pbserr\n" in submit_text
        assert "${PBS_ARRAYID}.pbsout" not in submit_text
        assert "TASK_ID=$PBS_ARRAYID" in submit_text
        assert "python chemsmart_run_array_mols_batch_array.py" in submit_text
        assert "python chemsmart_run_array_${TASK_ID}.py" not in submit_text
        assert "===== BEGIN array task ${TASK_ID} =====" in submit_text
        assert "flock 9 || exit 1" in submit_text
        assert (tmp_path / "chemsmart_run_array_mols_batch_array.py").exists()
        assert "#PBS -m n\n" in submit_text
        assert submit_text.count("#PBS -m abe\n") == 0

    def test_slurm_array_submit_script_mails_once_for_whole_array(
        self, tmp_path, monkeypatch, mocker
    ):
        self._stub_program_sections(monkeypatch)
        monkeypatch.chdir(tmp_path)
        mock_settings = mocker.Mock()
        mock_settings.data = {"EMAIL": "user@example.com"}
        mocker.patch(
            "chemsmart.settings.submitters.user_settings",
            mock_settings,
        )
        server = Server(
            "array-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=16,
            MEM_GB=32,
            NUM_GPUS=0,
            NUM_HOURS=12,
            QUEUE_NAME="normal",
            EXTRA_SCHEDULER_DIRECTIVES=(
                "#SBATCH --mail-type=END,FAIL,ARRAY_TASKS\n"
            ),
        )
        children = self._array_children(tmp_path, count=3)
        submitter = SLURMSubmitter(job=children[0], server=server)
        submitter.write_array_job(
            jobs=children,
            array_concurrency=1,
            cli_args=["gaussian", "opt"],
            batch_label="mols_batch",
        )

        submit_text = (
            tmp_path / "chemsmart_sub_array_mols_batch.sh"
        ).read_text()
        assert "#SBATCH --mail-user=user@example.com\n" in submit_text
        assert "#SBATCH --mail-type=END,FAIL\n" in submit_text
        assert "ARRAY_TASKS" not in submit_text

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
        assert "#BSUB -o pka_batch_array.bsubout\n" in submit_text
        assert "#BSUB -e pka_batch_array.bsuberr\n" in submit_text
        assert "%I.bsubout" not in submit_text
        assert "TASK_ID=$LSB_JOBINDEX" in submit_text
        assert "===== BEGIN array task ${TASK_ID} =====" in submit_text
        assert "flock 9 || exit 1" in submit_text
        assert "python chemsmart_run_array_pka_batch_array.py" in submit_text
        assert (tmp_path / "chemsmart_run_array_pka_batch_array.py").exists()


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
        from chemsmart.jobs.batch import (
            make_batch_cli_rewriter,
            set_job_batch_entry,
        )
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
        set_job_batch_entry(
            children[0], {"filepath": "mols.xyz", "molecule_index": 1}
        )
        set_job_batch_entry(
            children[1], {"filepath": "mols.xyz", "molecule_index": 2}
        )
        batch = GaussianBatchJob(
            jobs=children,
            label="mols_batch",
            rewrite_cli=make_batch_cli_rewriter("opt"),
        )
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
            ["gaussian", "-f", "mols.xyz", "-i", "1", "opt"],
            ["gaussian", "-f", "mols.xyz", "-i", "2", "opt"],
        ]

    def test_submit_batch_rejects_non_batch_job(self):
        from chemsmart.settings.server import Server

        server = Server(
            "batch-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=16,
            NUM_GPUS=0,
        )
        with pytest.raises(TypeError, match="expects a BatchJob"):
            server.submit_batch(type("NotBatch", (), {"label": "x"})())

    def test_submit_batch_rejects_empty_batch(self):
        from chemsmart.jobs.gaussian.batch import GaussianBatchJob
        from chemsmart.settings.server import Server

        server = Server(
            "batch-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=16,
            NUM_GPUS=0,
        )
        with pytest.raises(ValueError, match="no child jobs"):
            server.submit_batch(GaussianBatchJob(jobs=[], label="empty"))

    def test_submit_batch_requires_rewrite_cli_when_batch_entries_present(
        self,
    ):
        from chemsmart.jobs.batch import set_job_batch_entry
        from chemsmart.jobs.gaussian.batch import GaussianBatchJob
        from chemsmart.settings.server import Server

        server = Server(
            "batch-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=16,
            NUM_GPUS=0,
        )
        child = type(
            "Child",
            (),
            {"label": "a", "PROGRAM": "gaussian", "folder": "."},
        )()
        set_job_batch_entry(
            child, {"filepath": "mols.xyz", "molecule_index": 1}
        )
        batch = GaussianBatchJob(jobs=[child], label="needs_rewrite")
        with pytest.raises(ValueError, match="rewrite_cli"):
            server.submit_batch(batch, test=True, cli_args=["gaussian", "opt"])

    def test_submit_array_job_checks_orca_batch_container_label(
        self, tmp_path, monkeypatch
    ):
        from chemsmart.jobs.orca.batch import ORCABatchJob
        from chemsmart.settings.server import Server

        monkeypatch.chdir(tmp_path)
        server = Server(
            "batch-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=16,
            NUM_GPUS=0,
            SUBMIT_COMMAND="sbatch",
        )
        checked = []

        def _capture(job):
            checked.append(job)

        monkeypatch.setattr(
            Server, "_check_running_jobs", staticmethod(_capture)
        )
        monkeypatch.setattr(
            Server,
            "get_submitter",
            lambda self, job, **kwargs: type(
                "Stub",
                (),
                {
                    "write_array_job": staticmethod(
                        lambda **kw: checked.append(("wrote", kw))
                    ),
                    "array_submit_script": "unused.sh",
                },
            )(),
        )

        children = [
            type(
                "Child",
                (),
                {
                    "label": "o1",
                    "PROGRAM": "orca",
                    "folder": str(tmp_path),
                },
            )(),
            type(
                "Child",
                (),
                {
                    "label": "o2",
                    "PROGRAM": "orca",
                    "folder": str(tmp_path),
                },
            )(),
        ]
        server.submit_array_job(
            jobs=children,
            array_concurrency=1,
            test=True,
            cli_args=["orca", "opt"],
            batch_label="orca_batch",
        )

        assert isinstance(checked[0], ORCABatchJob)
        assert checked[0].label == "orca_batch"
        assert checked[0].jobs == children


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
