"""
GROMACS project settings.

This module defines lightweight settings objects for GROMACS workflows.
The first supported workflow is the prepared workflow, where users provide
existing MDP, structure, and topology files.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from chemsmart.io.yaml import YAMLFile


@dataclass
class GromacsProjectSettings:
    """
    Project-level settings for GROMACS workflows.

    These settings do not automatically decide simulation parameters.
    They record user-defined project information so that GROMACS jobs can be
    created and executed consistently.
    """

    project_name: Optional[str] = None
    project_dir: Optional[Path] = None
    workflow: str = "prepared"

    mdp_file: Optional[Path] = None
    structure_file: Optional[Path] = None
    top_file: Optional[Path] = None
    tpr_file: Optional[Path] = None
    index_file: Optional[Path] = None
    itp_files: list[Path] = field(default_factory=list)

    force_field: Optional[str] = None
    water_model: Optional[str] = None

    timestep: Optional[float] = None
    temperature: Optional[float] = None
    pressure: Optional[float] = None
    thermostat: Optional[str] = None
    barostat: Optional[str] = None
    constraints: Optional[str] = None
    constraint_algorithm: Optional[str] = None

    @classmethod
    def from_yaml(cls, filename):
        """
        Create GromacsProjectSettings from a YAML file.

        Relative input paths are resolved against the YAML file directory.
        """
        filename = Path(filename)
        yaml_file = YAMLFile(filename=filename)
        data = yaml_file.yaml_contents_dict

        settings_data = {
            "project_dir": filename.parent,
        }

        project_data = data.get("project", {})
        inputs_data = data.get("inputs", {})
        gromacs_data = data.get("gromacs_settings", {})

        if "name" in project_data:
            settings_data["project_name"] = project_data["name"]

        if "workflow" in project_data:
            settings_data["workflow"] = project_data["workflow"]
        elif "mode" in project_data:
            settings_data["workflow"] = project_data["mode"]

        settings_data.update(gromacs_data)

        input_key_map = {
            "mdp_file": "mdp_file",
            "structure_file": "structure_file",
            "topology_file": "top_file",
            "top_file": "top_file",
            "tpr_file": "tpr_file",
            "index_file": "index_file",
            "itp_files": "itp_files",
        }

        for yaml_key, settings_key in input_key_map.items():
            if yaml_key in inputs_data:
                settings_data[settings_key] = inputs_data[yaml_key]

        for key, value in data.items():
            if key not in {
                "project",
                "inputs",
                "gromacs_settings",
                "programs",
                "workflow",
            }:
                settings_data.setdefault(key, value)

        return cls.from_dict(settings_data)

    @classmethod
    def from_dict(cls, data):
        """
        Create GromacsProjectSettings from a dictionary.
        """
        data = data.copy()

        project_dir = data.get("project_dir")
        if project_dir is not None:
            project_dir = Path(project_dir)
            data["project_dir"] = project_dir

        for key in [
            "mdp_file",
            "structure_file",
            "top_file",
            "tpr_file",
            "index_file",
        ]:
            if data.get(key) is not None:
                data[key] = Path(data[key])

        if data.get("itp_files") is not None:
            data["itp_files"] =[
                cls._make_path(file, project_dir) for file in data["itp_files"]
            ]

        return cls(**data)

    @staticmethod
    def _make_path(path, project_dir=None):
        """
        Convert a path-like value to Path.

        If project_dir is provided and the path is relative, resolve it against
        the project directory.
        """
        path = Path(path)

        if project_dir is not None and not path.is_absolute():
            return Path(project_dir) / path

        return path

    def validate_prepared_inputs(self):
        """
        Validate that the prepared workflow has the minimum required inputs.
        """
        missing = []

        required_files = {
            "mdp_file": self.mdp_file,
            "structure_file": self.structure_file,
            "top_file": self.top_file,
        }

        for name, path in required_files.items():
            if path is None:
                missing.append(name)

        if missing:
            raise ValueError(
                "Missing required GROMACS prepared input settings: "
                + ", ".join(missing)
            )

    def to_job_kwargs(self):
        """
        Convert project settings into keyword arguments for GromacsJob.

        This allows project-level settings to be passed into job creation
        without duplicating file/path logic.
        """
        return {
            "mdp_file": self.mdp_file,
            "structure_file": self.structure_file,
            "top_file": self.top_file,
            "tpr_file": self.tpr_file,
            "itp_files": self.itp_files,
            "index_file": self.index_file,
            "workflow": self.workflow,
        }