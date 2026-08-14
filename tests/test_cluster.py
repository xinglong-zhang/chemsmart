"""Tests for chemsmart.utils.cluster."""

from unittest.mock import MagicMock

import pytest
import requests

from chemsmart.utils.cluster import (
    ClusterHelper,
    is_pubchem_api_available,
    is_pubchem_network_available,
)


class FakeCompletedProcess:
    """Minimal stand-in for the Popen object returned by subprocess.Popen."""

    def __init__(self, stdout_bytes, stderr_bytes=b""):
        self._stdout_bytes = stdout_bytes
        self._stderr_bytes = stderr_bytes
        # Piped Popen chains call `.stdout.close()` on the upstream process.
        self.stdout = MagicMock()

    def communicate(self):
        return self._stdout_bytes, self._stderr_bytes


class TestClusterHelperUsername:
    """Tests for ClusterHelper._get_username."""

    def test_get_username_uses_os_getlogin(self, mocker):
        mocker.patch("os.getlogin", return_value="alice")
        helper = ClusterHelper()
        assert helper.username == "alice"

    def test_get_username_falls_back_to_whoami(self, mocker):
        mocker.patch("os.getlogin", side_effect=OSError())
        mock_popen = mocker.patch(
            "chemsmart.utils.cluster.subprocess.Popen",
            return_value=FakeCompletedProcess(b"bob\n"),
        )
        helper = ClusterHelper()
        assert helper.username == "bob"
        mock_popen.assert_called_once()


class TestClusterHelperSlurm:
    """Tests for SLURM-related job discovery."""

    def test_get_gaussian_running_job_ids_on_slurm(self, mocker):
        mocker.patch("os.getlogin", return_value="alice")
        squeue_output = (
            b"  JOBID PARTITION     NAME     USER\n"
            b"  12345      main   job1   alice\n"
            b"  67890      main   job2   alice\n"
        )
        mocker.patch(
            "chemsmart.utils.cluster.subprocess.Popen",
            side_effect=[
                FakeCompletedProcess(b""),
                FakeCompletedProcess(squeue_output),
            ],
        )
        helper = ClusterHelper()
        job_ids = helper._get_gaussian_running_job_ids_on_slurm()
        assert job_ids == [12345, 67890]

    def test_get_gaussian_running_job_ids_on_slurm_no_matches(self, mocker):
        mocker.patch("os.getlogin", return_value="alice")
        mocker.patch(
            "chemsmart.utils.cluster.subprocess.Popen",
            side_effect=[
                FakeCompletedProcess(b""),
                FakeCompletedProcess(b""),
            ],
        )
        helper = ClusterHelper()
        assert helper._get_gaussian_running_job_ids_on_slurm() == []

    def test_get_gaussian_running_job_names_on_slurm(self, mocker):
        mocker.patch("os.getlogin", return_value="alice")
        squeue_output = b"  12345      main   job1   alice\n"
        sacct_output = (
            b"User       JobID       JobName\n"
            b"alice      12345       chemsmart_sub_mymol.sh\n"
        )

        popen_returns = [
            FakeCompletedProcess(b""),
            FakeCompletedProcess(squeue_output),
            FakeCompletedProcess(sacct_output),
        ]
        mocker.patch(
            "chemsmart.utils.cluster.subprocess.Popen",
            side_effect=popen_returns,
        )
        helper = ClusterHelper()
        names = helper._get_gaussian_running_job_names_on_slurm()
        assert names == ["mymol"]

    def test_get_gaussian_running_jobs_slurm_success(self, mocker):
        mocker.patch("os.getlogin", return_value="alice")
        helper = ClusterHelper()
        mocker.patch.object(
            helper,
            "_get_gaussian_running_job_ids_on_slurm",
            return_value=[1, 2],
        )
        mocker.patch.object(
            helper,
            "_get_gaussian_running_job_names_on_slurm",
            return_value=["a", "b"],
        )
        ids, names = helper.get_gaussian_running_jobs()
        assert ids == [1, 2]
        assert names == ["a", "b"]

    def test_get_gaussian_running_jobs_falls_back_to_torque(self, mocker):
        mocker.patch("os.getlogin", return_value="alice")
        helper = ClusterHelper()
        mocker.patch.object(
            helper,
            "_get_gaussian_running_job_ids_on_slurm",
            side_effect=RuntimeError("no slurm"),
        )
        mocker.patch.object(
            helper,
            "_get_gaussian_running_jobs_on_torque",
            return_value=([3], ["c"]),
        )
        ids, names = helper.get_gaussian_running_jobs()
        assert ids == [3]
        assert names == ["c"]

    def test_get_gaussian_running_jobs_no_scheduler_available(self, mocker):
        mocker.patch("os.getlogin", return_value="alice")
        helper = ClusterHelper()
        mocker.patch.object(
            helper,
            "_get_gaussian_running_job_ids_on_slurm",
            side_effect=RuntimeError("no slurm"),
        )
        mocker.patch.object(
            helper,
            "_get_gaussian_running_jobs_on_torque",
            side_effect=FileNotFoundError(),
        )
        ids, names = helper.get_gaussian_running_jobs()
        assert ids == []
        assert names == []


class TestClusterHelperTorque:
    """Tests for PBS/Torque job discovery."""

    def test_get_gaussian_running_jobs_on_torque(self, mocker):
        mocker.patch("os.getlogin", return_value="alice")
        qstat_output = (
            b"Job Id: 555.server\n"
            b"    Job_Owner = alice@server\n"
            b"    Job_Name = chemsmart_sub_mymol.sh\n"
        )
        mocker.patch(
            "chemsmart.utils.cluster.subprocess.Popen",
            side_effect=[
                FakeCompletedProcess(b""),
                FakeCompletedProcess(qstat_output),
            ],
        )
        helper = ClusterHelper()
        job_ids, job_names = helper._get_gaussian_running_jobs_on_torque()
        assert job_ids == [555]
        assert job_names == ["mymol"]

    def test_get_gaussian_running_jobs_on_torque_no_duplicate_ids(
        self, mocker
    ):
        mocker.patch("os.getlogin", return_value="alice")
        qstat_output = (
            b"Job Id: 555.server\n"
            b"Job Id: 555.server\n"
            b"    Job_Name = chemsmart_sub_mymol.sh\n"
        )
        mocker.patch(
            "chemsmart.utils.cluster.subprocess.Popen",
            side_effect=[
                FakeCompletedProcess(b""),
                FakeCompletedProcess(qstat_output),
            ],
        )
        helper = ClusterHelper()
        job_ids, _ = helper._get_gaussian_running_jobs_on_torque()
        assert job_ids == [555]


class TestPubchemNetworkAvailability:
    """Tests for is_pubchem_network_available."""

    def test_returns_true_when_connection_succeeds(self, mocker):
        mocker.patch(
            "chemsmart.utils.cluster.socket.create_connection",
            return_value=None,
        )
        assert is_pubchem_network_available() is True

    def test_returns_false_on_os_error(self, mocker):
        mocker.patch(
            "chemsmart.utils.cluster.socket.create_connection",
            side_effect=OSError(),
        )
        assert is_pubchem_network_available() is False


class TestPubchemApiAvailability:
    """Tests for is_pubchem_api_available."""

    @pytest.mark.parametrize("status_code", [200, 404])
    def test_returns_true_for_success_and_not_found(self, mocker, status_code):
        mock_response = mocker.Mock(status_code=status_code)
        mocker.patch(
            "chemsmart.utils.cluster.requests.get",
            return_value=mock_response,
        )
        assert is_pubchem_api_available() is True

    @pytest.mark.parametrize("status_code", [500, 503])
    def test_returns_false_for_server_errors(self, mocker, status_code):
        mock_response = mocker.Mock(status_code=status_code)
        mocker.patch(
            "chemsmart.utils.cluster.requests.get",
            return_value=mock_response,
        )
        assert is_pubchem_api_available() is False

    def test_returns_false_on_request_exception(self, mocker):
        mocker.patch(
            "chemsmart.utils.cluster.requests.get",
            side_effect=requests.exceptions.Timeout(),
        )
        assert is_pubchem_api_available() is False

    def test_returns_false_on_os_error(self, mocker):
        mocker.patch(
            "chemsmart.utils.cluster.requests.get",
            side_effect=OSError(),
        )
        assert is_pubchem_api_available() is False
