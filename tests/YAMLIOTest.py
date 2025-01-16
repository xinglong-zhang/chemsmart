import os
from chemsmart.io.yaml import YAMLFile


class TestYAMLFile:
    def test_server_yaml(self, server_yaml_file):
        assert os.path.exists(server_yaml_file)
        assert os.path.isfile(server_yaml_file)
        server_yaml = YAMLFile(filename=server_yaml_file)
        assert len(server_yaml.yaml_contents_dict) == 3
        assert server_yaml.yaml_contents_dict["SERVER"]["SCHEDULER"] == "PBS"
        assert len(server_yaml.yaml_contents_dict["SERVER"].keys()) == 12
        assert list(server_yaml.yaml_contents_dict.keys())[0] == "SERVER"
        assert list(server_yaml.yaml_contents_dict.keys())[1] == "GAUSSIAN"
        assert list(server_yaml.yaml_contents_dict.keys())[2] == "ORCA"
        assert (
            server_yaml.yaml_contents_dict["GAUSSIAN"]["EXEFOLDER"]
            == "~/programs/g16"
        )
        assert (
            server_yaml.yaml_contents_dict["ORCA"]["EXEFOLDER"]
            == "~/programs/orca_6_0_0"
        )
