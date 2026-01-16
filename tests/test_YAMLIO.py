import os

from chemsmart.io.yaml import YAMLFile


class TestYAMLFile:
    def test_server_yaml(self, server_yaml_file):
        """Test reading and parsing server YAML configuration file."""
        assert os.path.exists(server_yaml_file)
        assert os.path.isfile(server_yaml_file)
        server_yaml = YAMLFile(filename=server_yaml_file)

        assert len(server_yaml.yaml_contents_dict) == 4
        expected_keys = ["SERVER", "GAUSSIAN", "ORCA", "XTB"]
        assert list(server_yaml.yaml_contents_dict.keys()) == expected_keys

        # Test SERVER section
        assert "SERVER" in server_yaml.yaml_contents_dict
        server_config = server_yaml.yaml_contents_dict["SERVER"]
        assert server_config["SCHEDULER"] == "PBS"
        assert server_config["QUEUE_NAME"] == "normal"
        assert server_config["NUM_HOURS"] == 24
        assert server_config["MEM_GB"] == 400
        assert server_config["NUM_CORES"] == 64
        assert server_config["NUM_GPUS"] == 0
        assert server_config["NUM_THREADS"] == 64
        assert server_config["SUBMIT_COMMAND"] == "qsub"
        assert server_config["PROJECT"] == 13003611
        assert server_config["SCRATCH_DIR"] is None
        assert server_config["USE_HOSTS"] is True

        # Test GAUSSIAN section
        assert "GAUSSIAN" in server_yaml.yaml_contents_dict
        gaussian_config = server_yaml.yaml_contents_dict["GAUSSIAN"]
        assert gaussian_config["EXEFOLDER"] == "~/programs/g16"
        assert gaussian_config["LOCAL_RUN"] is True
        assert gaussian_config["SCRATCH"] is True

        # Test ORCA section
        assert "ORCA" in server_yaml.yaml_contents_dict
        orca_config = server_yaml.yaml_contents_dict["ORCA"]
        assert orca_config["EXEFOLDER"] == "~/programs/orca_6_0_0"
        assert orca_config["LOCAL_RUN"] is False

        # Test XTB section
        assert "XTB" in server_yaml.yaml_contents_dict
        xtb_config = server_yaml.yaml_contents_dict["XTB"]
        assert xtb_config["EXEFOLDER"] is None
        assert xtb_config["LOCAL_RUN"] is True
        assert xtb_config["SCRATCH"] is False
