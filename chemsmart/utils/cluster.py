
"""Utils for cluster-related stuff."""
import os
import shlex
import subprocess

class ClusterHelper:
    """Class for cluster helper to obtain the running jobs. Currently only works
    for SLURM and Torque/PBS queuing systems. Can be expanded later to others."""
    def __init__(self):
        self.username = self._get_username()

    def _get_username(self):
        try:
            username = os.getlogin()  # OSError: [Errno 6] No such device or address
        except OSError:
            cmd = 'whoami'
            p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
            out, err = p.communicate()
            username = out.decode('utf-8').strip()
        return username

    def get_gaussian_running_jobs(self):
        try:
            running_job_ids = self._get_gaussian_running_job_ids_on_slurm()
            running_job_names = self._get_gaussian_running_job_names_on_slurm()
        except Exception:
            try:
                running_job_ids, running_job_names = self._get_gaussian_running_jobs_on_torque()
            except FileNotFoundError:
                running_job_ids, running_job_names = [], []
        return running_job_ids, running_job_names

    def _get_gaussian_running_job_ids_on_slurm(self):
        """Function for getting a list of the job ids of gaussian jobs that are running on slurm."""
        running_job_ids = []
        p1 = subprocess.Popen(shlex.split('squeue'), stdout=subprocess.PIPE)
        p2 = subprocess.Popen(shlex.split(f'grep {self.username}'), stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        out, err = p2.communicate()
        for raw_line in out.decode('utf-8').split('\n'):
            line = raw_line.strip()
            if self.username in line:
                line_elem = line.split()
                job_id = int(line_elem[0])
                running_job_ids.append(job_id)
        return running_job_ids

    def _get_gaussian_running_job_names_on_slurm(self):
        """Function for getting a list of the job names of gaussian jobs that are running on slurm."""
        running_job_ids = self._get_gaussian_running_job_ids_on_slurm()
        running_job_names = []
        cmd = 'sacct --format="User,JobID,JobName%100"'
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
        out, err = p.communicate()
        for raw_line in out.decode('utf-8').split('\n'):
            line = raw_line.strip()
            if self.username in line:
                line_elem = line.split()
                job_id = int(line_elem[1])
                # checks for currently running jobs
                if job_id in running_job_ids:
                    job_name_submitscript = line_elem[-1]
                    job_name_submitscript_no_ext = os.path.splitext(job_name_submitscript)[0]
                    # gaussian submission script names
                    job_name = job_name_submitscript_no_ext.split('chemsmart_sub_')[-1]
                    running_job_names.append(job_name)
        return running_job_names

    def _get_gaussian_running_jobs_on_torque(self):
        """Function for getting a list of the job ids of gaussian jobs that are running on torque."""
        running_job_ids = []
        running_job_names = []
        p1 = subprocess.Popen(shlex.split('qstat -f'), stdout=subprocess.PIPE)
        p2 = subprocess.Popen(shlex.split(f'grep -C 2 {self.username}@'), stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        out, err = p2.communicate()
        for line in out.decode('utf-8').split('\n'):
            if 'Job Id:' in line:
                job_id_name = line.split('Job Id:')[-1].strip()
                job_id = int(job_id_name.split('.')[0])
                if job_id not in running_job_ids:
                    running_job_ids.append(job_id)
            if 'Job_Name' in line:
                job_name_submitscript = line.split('=')[-1].strip()
                job_name_submitscript_no_ext = os.path.splitext(job_name_submitscript)[0]
                # gaussian submission script names
                job_name = job_name_submitscript_no_ext.split('chemsmart_sub_')[-1]
                running_job_names.append(job_name)
        return running_job_ids, running_job_names