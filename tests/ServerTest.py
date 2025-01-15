import os
from chemsmart.settings.server import Server


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
        extra_commands = """source ~/anaconda3/etc/profile.d/conda.sh
conda activate /home/users/astar/bmsi/zhangx5/anaconda3/envs/pyatoms-dev2
module purge
module load craype-x86-rome
module load libfabric/1.11.0.4.125
module load craype-network-ofi
module load perftools-base/22.04.0
module load cce/13.0.2
module load craype/2.7.15
module load cray-dsmml/0.2.2
module load cray-mpich/8.1.15
module load cray-libsci/21.08.1.2
module load cray-pals/1.1.6
module load PrgEnv-cray/8.3.3
module load gcc/11.2.0
module load mkl/2022.0.2
module load intel/2022.0.2
"""
        assert server.extra_commands == extra_commands

        # server_yaml = YAMLFile(filename=server_yaml_file)
        # assert len(server_yaml.yaml_contents_dict) == 3
        # assert server_yaml.yaml_contents_dict["SERVER"]["SCHEDULER"] == "pbs"
        # assert len(server_yaml.yaml_contents_dict["SERVER"].keys()) == 12
        # assert list(server_yaml.yaml_contents_dict.keys())[0] == "SERVER"
        # assert list(server_yaml.yaml_contents_dict.keys())[1] == "GAUSSIAN"
        # assert list(server_yaml.yaml_contents_dict.keys())[2] == "ORCA"
        # assert server_yaml.yaml_contents_dict["GAUSSIAN"]["G16FOLDER"] == "~/programs/g16"
        # assert server_yaml.yaml_contents_dict["ORCA"]["ORCAFOLDER"] == "~/programs/orca_6_0_0"
