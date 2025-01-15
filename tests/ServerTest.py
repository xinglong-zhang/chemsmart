import os
from chemsmart.settings.server import Server
from chemsmart.settings.executable import GaussianExecutable


class TestServer:
    def test_server_yaml(self, server_yaml_file):
        assert os.path.exists(server_yaml_file)
        assert os.path.isfile(server_yaml_file)
        server = Server.from_yaml(name=server_yaml_file)
        assert server.scheduler.lower() == "pbs"
        assert server.queue_name == "normal"
        assert server.num_hours == 24
        assert server.mem_gb == 400
        assert server.num_cores == 64
        assert server.num_gpus == 0
        assert server.num_threads == 64
        assert server.submit_command == "qsub"
        assert server.scratch_dir is None
        assert server.use_hosts is True
        assert server.extra_commands is None

    def test_gaussian_executable(self, server_yaml_file):
        gaussian_executable = GaussianExecutable.from_servername(
            server_yaml_file
        )
        assert gaussian_executable.executable_folder == os.path.expanduser(
            "~/programs/g16"
        )
        assert gaussian_executable.local_run is True

        gaussian_conda_env = """source ~/anaconda3/etc/profile.d/conda.sh
conda activate ~/anaconda3/envs/pyatoms-dev2
"""
        assert gaussian_executable.conda_env == gaussian_conda_env

        gaussian_modules = """module purge
module load craype-x86-rome
module load libfabric/1.11.0.4.125
"""
        assert gaussian_executable.modules == gaussian_modules

        assert (
            gaussian_executable.scripts == "source GAUSS_DIR=~/programs/g16\n"
        )

        gassian_envars = """export SCRATCH=~/scratch
export GAUSS_EXEDIR=~/programs/g16
export g16root=~/programs/g16

"""
        assert gaussian_executable.envars == gassian_envars
        print(gaussian_executable.executable_folder)

        # server_yaml = YAMLFile(filename=server_yaml_file)
        # assert len(server_yaml.yaml_contents_dict) == 3
        # assert server_yaml.yaml_contents_dict["SERVER"]["SCHEDULER"] == "pbs"
        # assert len(server_yaml.yaml_contents_dict["SERVER"].keys()) == 12
        # assert list(server_yaml.yaml_contents_dict.keys())[0] == "SERVER"
        # assert list(server_yaml.yaml_contents_dict.keys())[1] == "GAUSSIAN"
        # assert list(server_yaml.yaml_contents_dict.keys())[2] == "ORCA"
        # assert server_yaml.yaml_contents_dict["GAUSSIAN"]["G16FOLDER"] == "~/programs/g16"
        # assert server_yaml.yaml_contents_dict["ORCA"]["ORCAFOLDER"] == "~/programs/orca_6_0_0"
