"""
Cluster utilities for job monitoring and connectivity checks.

Helpers for interacting with common HPC schedulers and network services
used by CHEMSMART workflows. Currently supports SLURM (via `squeue`/`sacct`)
and PBS/Torque (via `qstat`) to discover running jobs, along with a
simple PubChem reachability check.

Key functionality:
- Detect current username reliably.
- Query running jobs (IDs and names) on SLURM or PBS/Torque.
- Validate network connectivity to PubChem for external lookups.
"""

import os
import shlex
import socket
import subprocess

import requests


class ClusterHelper:
    """
    Helper for querying running jobs on HPC schedulers.

    Provides methods to obtain IDs and names of running jobs on
    SLURM and PBS/Torque clusters.

    Attributes:
        username (str): Current user's username for filtering scheduler output.
    """

    def __init__(self):
        """
        Initialize cluster helper with current username detection.
        """
        self.username = self._get_username()

    def _get_username(self):
        """
        Determine the current username with fallbacks.

        Tries `os.getlogin()` first and falls back to invoking `whoami`
        if the login name is unavailable (e.g., no controlling TTY).

        Returns:
            str: Current username.
        """
        try:
            username = (
                os.getlogin()
            )  # OSError: [Errno 6] No such device or address
        except OSError:
            cmd = "whoami"
            p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
            out, err = p.communicate()
            username = out.decode("utf-8").strip()
        return username

    def get_running_jobs(self):
        """
        Get IDs and names of currently running jobs.

        Queries SLURM first and, on failure, falls back to PBS/Torque.

        Returns:
            tuple[list[int], list[str]]: (job_ids, job_names). Returns two
            empty lists if no jobs are found or no supported scheduler is
            available.
        """
        try:
            running_job_ids = self._get_running_job_ids_on_slurm()
            running_job_names = self._get_running_job_names_on_slurm()
        except Exception:
            try:
                running_job_ids, running_job_names = (
                    self._get_running_jobs_on_torque()
                )
            except FileNotFoundError:
                running_job_ids, running_job_names = [], []
        return running_job_ids, running_job_names

    def _get_running_job_ids_on_slurm(self):
        """
        Collect job IDs running on SLURM.

        Invokes `squeue` and filters rows by the current username.

        Returns:
            list[int]: SLURM job IDs for running jobs.
        """
        running_job_ids = []
        p1 = subprocess.Popen(shlex.split("squeue"), stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            shlex.split(f"grep {self.username}"),
            stdin=p1.stdout,
            stdout=subprocess.PIPE,
        )
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        out, err = p2.communicate()
        for raw_line in out.decode("utf-8").split("\n"):
            line = raw_line.strip()
            if self.username in line:
                line_elem = line.split()
                job_id = int(line_elem[0])
                running_job_ids.append(job_id)
        return running_job_ids

    def _get_running_job_names_on_slurm(self):
        """
        Collect job names running on SLURM.

        Invokes `sacct` and filters by username, then cross-references IDs
        from `squeue`. Job names are returned as reported by the scheduler.

        Returns:
            list[str]: Job names for currently running jobs.
        """
        running_job_ids = self._get_running_job_ids_on_slurm()
        running_job_names = []
        cmd = 'sacct --format="User,JobID,JobName%100"'
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
        out, err = p.communicate()
        for raw_line in out.decode("utf-8").split("\n"):
            line = raw_line.strip()
            if self.username in line:
                line_elem = line.split()
                job_id = int(line_elem[1])
                # Checks for currently running jobs
                if job_id in running_job_ids:
                    running_job_names.append(line_elem[-1])
        return running_job_names

    def _get_running_jobs_on_torque(self):
        """
        Collect job IDs and names on PBS/Torque.

        Invokes `qstat -f`, filters surrounding lines by username, parses
        job IDs and names as reported by the scheduler.

        Returns:
            tuple[list[int], list[str]]: (job_ids, job_names).
        """
        running_job_ids = []
        running_job_names = []
        p1 = subprocess.Popen(shlex.split("qstat -f"), stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            shlex.split(f"grep -C 2 {self.username}@"),
            stdin=p1.stdout,
            stdout=subprocess.PIPE,
        )
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        out, err = p2.communicate()
        for line in out.decode("utf-8").split("\n"):
            if "Job Id:" in line:
                job_id_name = line.split("Job Id:")[-1].strip()
                job_id = int(job_id_name.split(".")[0])
                if job_id not in running_job_ids:
                    running_job_ids.append(job_id)
            if "Job_Name" in line:
                running_job_names.append(line.split("=")[-1].strip())
        return running_job_ids, running_job_names


def is_pubchem_network_available():
    """
    Check PubChem (NCBI) HTTPS connectivity.

    Attempts to open a TCP connection to pubchem.ncbi.nlm.nih.gov:443 with
    a short timeout to validate reachability.

    Returns:
        bool: True if reachable, False otherwise.
    """
    try:
        socket.create_connection(("pubchem.ncbi.nlm.nih.gov", 443), timeout=5)
        return True
    except OSError:
        return False


def is_pubchem_api_available():
    """
    Check if PubChem REST API is available and responding.

    Makes a lightweight API request to verify the service is operational
    and not returning server errors (503, 500, etc.).

    Returns:
        bool: True if API is available and responding, False if service is
              down, busy (503), or unreachable.
    """

    # Use a simple, lightweight endpoint (aspirin - CID 2244)
    test_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula/JSON"

    try:
        response = requests.get(test_url, timeout=5)
        # Accept 200 (success) or even 404 (not
        # found) - both mean API is responding
        # Reject 503 (busy), 500 (server error), etc.
        return response.status_code in (200, 404)
    except (requests.exceptions.RequestException, OSError):
        # Network error, timeout, or connection failure
        return False
