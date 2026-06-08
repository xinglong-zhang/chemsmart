from chemsmart.agent.wizard import validate_server_yaml

VALID_YAML = """
SERVER:
  SCHEDULER: SLURM
  QUEUE_NAME: debug
  NUM_HOURS: 8
  MEM_GB: 128
  NUM_CORES: 32
  NUM_GPUS: 0
  NUM_THREADS: 32
  SUBMIT_COMMAND: sbatch
GAUSSIAN:
  EXEFOLDER: /apps/gaussian
  LOCAL_RUN: true
  SCRATCH: true
  CONDA_ENV: /opt/conda/envs/chemsmart
  MODULES: ''
  ENVARS: ''
"""


def test_validate_server_yaml_accepts_valid_yaml():
    result = validate_server_yaml(VALID_YAML)

    assert result.ok is True
    assert result.errors == []
    assert result.parsed["SERVER"]["SCHEDULER"] == "SLURM"


def test_validate_server_yaml_reports_missing_scheduler():
    result = validate_server_yaml(
        VALID_YAML.replace("  SCHEDULER: SLURM\n", "")
    )

    assert result.ok is False
    assert "Missing required SERVER.SCHEDULER." in result.errors


def test_validate_server_yaml_rejects_nonsense_payload():
    result = validate_server_yaml("this is not server yaml")

    assert result.ok is False
    assert result.parsed is None
    assert result.errors == ["YAML root must be a mapping."]
