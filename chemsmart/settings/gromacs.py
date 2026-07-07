"""
GROMACS project-level settings.

This module stores reusable project-level configuration for GROMACS jobs.
The main purpose is to support project.yaml -> settings -> job creation.
"""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import List, Optional

import yaml


@dataclass
class GromacsProjectSettings:
    """
    Project-level settings for GROMACS workflows.

    These settings are intentionally separated from the job object. The settings
    object stores reusable project configuration, while the job object stores one
    concrete executable job.
    """

    project_name: Optional[str] = None
    project_dir: Optional[Path] = None

    workflow: str = "prepared"
    job_type: str = "em"

    mdp_file: Optional[Path] = None
    ions_mdp_file: Optional[Path] = None
    structure_file: Optional[Path] = None
    input_pdb: Optional[Path] = None
    top_file: Optional[Path] = None
    tpr_file: Optional[Path] = None
    index_file: Optional[Path] = None
    itp_files: Optional[List[Path]] = None

    processed_structure_file: Optional[Path] = None
    boxed_structure_file: Optional[Path] = None
    solvated_structure_file: Optional[Path] = None
    ions_tpr_file: Optional[Path] = None
    ionized_structure_file: Optional[Path] = None

    force_field: Optional[str] = None
    water_model: Optional[str] = None
    timestep: Optional[float] = None
    temperature: Optional[float] = None
    pressure: Optional[float] = None
    thermostat: Optional[str] = None
    barostat: Optional[str] = None
    constraints: Optional[str] = None
    constraint_algorithm: Optional[str] = None


    box_type: Optional[str] = "cubic"
    box_distance: Optional[float] = 1.0
    solvent_file: Optional[Path] = None

    positive_ion: Optional[str] = "NA"
    negative_ion: Optional[str] = "CL"
    neutral: bool = True
    genion_group: str = "SOL"

    grompp_maxwarn: Optional[int] = None
    mdrun_threads: Optional[int] = None
    mdrun_ntmpi: Optional[int] = None
    mdrun_ntomp: Optional[int] = None
    mdrun_extra_args: Optional[List[str]] = None

    @classmethod
    def from_dict(cls, data):
        """
        Create settings from a dictionary.

        Relative paths are resolved against project_dir if project_dir is given.
        """
        project_dir = data.get("project_dir")
        if project_dir is not None:
            project_dir = Path(project_dir)

        mapped = {
            "project_name": data.get("project_name"),
            "workflow": data.get("workflow", data.get("mode", "prepared")),
            "job_type": data.get("job_type", "em"),
            "mdp_file": data.get("mdp_file"),
            "ions_mdp_file": data.get("ions_mdp_file"),
            "structure_file": data.get("structure_file"),
            "input_pdb": data.get("input_pdb"),
            "top_file": data.get("top_file") or data.get("topology_file"),
            "tpr_file": data.get("tpr_file"),
            "index_file": data.get("index_file"),
            "itp_files": data.get("itp_files"),
            "processed_structure_file": data.get("processed_structure_file"),
            "boxed_structure_file": data.get("boxed_structure_file"),
            "solvated_structure_file": data.get("solvated_structure_file"),
            "ions_tpr_file": data.get("ions_tpr_file"),
            "ionized_structure_file": data.get("ionized_structure_file"),
            "force_field": data.get("force_field"),
            "water_model": data.get("water_model"),
            "timestep": data.get("timestep"),
            "temperature": data.get("temperature"),
            "pressure": data.get("pressure"),
            "thermostat": data.get("thermostat"),
            "barostat": data.get("barostat"),
            "constraints": data.get("constraints"),
            "constraint_algorithm": data.get("constraint_algorithm"),
            "box_type": data.get("box_type", "cubic"),
            "box_distance": data.get("box_distance", 1.0),
            "solvent_file": data.get("solvent_file"),
            "positive_ion": data.get("positive_ion", "NA"),
            "negative_ion": data.get("negative_ion", "CL"),
            "neutral": data.get("neutral", True),
            "genion_group": data.get("genion_group", "SOL"),
            "grompp_maxwarn": data.get("grompp_maxwarn"),
            "mdrun_threads": data.get("mdrun_threads"),
            "mdrun_ntmpi": data.get("mdrun_ntmpi"),
            "mdrun_ntomp": data.get("mdrun_ntomp"),
            "mdrun_extra_args": data.get("mdrun_extra_args"),
        }

        path_fields = {
            "mdp_file",
            "ions_mdp_file",
            "structure_file",
            "input_pdb",
            "top_file",
            "tpr_file",
            "index_file",
            "processed_structure_file",
            "boxed_structure_file",
            "solvated_structure_file",
            "ions_tpr_file",
            "ionized_structure_file",
            "solvent_file",
        }

        resolved = {
            "project_name": mapped["project_name"],
            "project_dir": project_dir,
            "workflow": mapped["workflow"],
            "job_type": mapped["job_type"],
        }

        for key, value in mapped.items():
            if key in {"project_name", "workflow", "job_type"}:
                continue

            if key in path_fields:
                resolved[key] = cls._resolve_path(value, project_dir)
            elif key == "itp_files":
                resolved[key] = cls._resolve_paths(value, project_dir)
            elif key == "mdrun_extra_args":
                resolved[key] = list(value or []) or None
            else:
                resolved[key] = value

        return cls(**resolved)

    @classmethod
    def from_yaml(cls, yaml_file):
        """
        Load settings from a GROMACS project YAML file.
        """
        yaml_path = Path(yaml_file).resolve()
        raw = yaml.safe_load(yaml_path.read_text(encoding="utf-8")) or {}

        project = raw.get("project", {})
        inputs = raw.get("inputs", {})
        gromacs_settings = raw.get("gromacs_settings", {})
        runtime = raw.get("runtime", {})

        flat = {
            "project_name": project.get("name"),
            "workflow": project.get(
                "mode", project.get("workflow", "prepared")
            ),
            "job_type": project.get("job_type", "em"),
            "project_dir": yaml_path.parent,
        }

        flat.update(inputs)
        flat.update(gromacs_settings)
        flat.update(runtime)

        return cls.from_dict(flat)

    @staticmethod
    def _resolve_path(value, project_dir):
        """
        Resolve one path-like value.
        """
        if value in (None, ""):
            return None

        path = Path(value)

        if path.is_absolute() or project_dir is None:
            return path

        return project_dir / path

    @classmethod
    def _resolve_paths(cls, values, project_dir):
        """
        Resolve a list of path-like values.
        """
        if not values:
            return []

        return [cls._resolve_path(value, project_dir) for value in values]

    def validate_prepared_inputs(self):
        """
        Validate the minimum inputs for a prepared workflow.

        A prepared workflow requires structure and topology files, then runs
        grompp -> mdrun. The MDP file is optional because it can be generated
        automatically by GromacsInputWriter.
        """
        missing = []

        if self.structure_file is None:
            missing.append("structure_file")
        if self.top_file is None:
            missing.append("top_file")

        if missing:
            raise ValueError(
                "Missing required GROMACS prepared input settings: "
                + ", ".join(missing)
            )

    def validate_full_setup_inputs(self):
        """
        Validate the minimum inputs for a full setup workflow.

        The first full setup implementation covers the common non-interactive
        happy path only. The MDP file is optional because it can be generated
        automatically by GromacsInputWriter.
        """
        missing = []

        if self.input_pdb is None and self.structure_file is None:
            missing.append("input_pdb|structure_file")
        if self.force_field is None:
            missing.append("force_field")
        if self.water_model is None:
            missing.append("water_model")

        if missing:
            raise ValueError(
                "Missing required GROMACS full_setup settings: "
                + ", ".join(missing)
            )

    def validate(self):
        """
        Validate settings according to the selected workflow.
        """
        if self.workflow == "prepared":
            self.validate_prepared_inputs()
            return

        if self.workflow == "full_setup":
            self.validate_full_setup_inputs()
            return

        raise ValueError(f"Unsupported GROMACS workflow: {self.workflow}")

    def with_overrides(self, **overrides):
        """
        Return a new settings object with selected values overridden.

        None values are ignored so that CLI options can safely override YAML
        only when the user explicitly provides them.
        """
        data = asdict(self)

        for key, value in overrides.items():
            if value is not None:
                data[key] = value

        return self.from_dict(data)

    def to_job_kwargs(self):
        """
        Convert project settings into keyword arguments for GromacsJob.
        """
        data = {
            "mdp_file": self.mdp_file,
            "ions_mdp_file": self.ions_mdp_file,
            "structure_file": self.structure_file,
            "input_pdb": self.input_pdb,
            "top_file": self.top_file,
            "tpr_file": self.tpr_file,
            "index_file": self.index_file,
            "itp_files": self.itp_files or [],
            "workflow": self.workflow,
            "processed_structure_file": self.processed_structure_file,
            "boxed_structure_file": self.boxed_structure_file,
            "solvated_structure_file": self.solvated_structure_file,
            "ions_tpr_file": self.ions_tpr_file,
            "ionized_structure_file": self.ionized_structure_file,
            "force_field": self.force_field,
            "water_model": self.water_model,
            "timestep": self.timestep,
            "temperature": self.temperature,
            "pressure": self.pressure,
            "thermostat": self.thermostat,
            "barostat": self.barostat,
            "constraints": self.constraints,
            "constraint_algorithm": self.constraint_algorithm,
            "box_type": self.box_type,
            "box_distance": self.box_distance,
            "solvent_file": self.solvent_file,
            "positive_ion": self.positive_ion,
            "negative_ion": self.negative_ion,
            "neutral": self.neutral,
            "genion_group": self.genion_group,
            "grompp_maxwarn": self.grompp_maxwarn,
            "mdrun_threads": self.mdrun_threads,
            "mdrun_ntmpi": self.mdrun_ntmpi,
            "mdrun_ntomp": self.mdrun_ntomp,
            "mdrun_extra_args": self.mdrun_extra_args or [],
        }

        return {key: value for key, value in data.items() if value is not None}
