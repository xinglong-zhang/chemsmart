"""Stable data contracts for deterministic ChemSmart command parsing."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class ParsedModelCommand:
    command: str
    workspace: str
    tokens: list[str]
    parse_error: str | None = None
    entrypoint: str | None = None
    action: str | None = None
    program: str | None = None
    job: str | None = None
    server: str | None = None
    dry_run: bool = False
    project: str | None = None
    project_p_flag_meaning: str | None = None
    top_level_program: str | None = None
    filename: str | None = None
    record_index: str | None = None
    record_id: str | None = None
    structure_index: str | None = None
    structure_id: str | None = None
    molecule_id: str | None = None
    label: str | None = None
    charge: str | None = None
    multiplicity: str | None = None
    functional: str | None = None
    ab_initio: str | None = None
    basis: str | None = None
    aux_basis: str | None = None
    extrapolation_basis: str | None = None
    defgrid: str | None = None
    scf_tol: str | None = None
    scf_algorithm: str | None = None
    solvent_model: str | None = None
    solvent_id: str | None = None
    route_parameters: str | None = None
    opt_options: str | None = None
    structural_options: dict[str, str] = field(default_factory=dict)
    resources: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "command": self.command,
            "workspace": self.workspace,
            "tokens": self.tokens,
            "parse_error": self.parse_error,
            "entrypoint": self.entrypoint,
            "action": self.action,
            "program": self.program,
            "job": self.job,
            "server": self.server,
            "dry_run": self.dry_run,
            "project": self.project,
            "project_p_flag_meaning": self.project_p_flag_meaning,
            "top_level_program": self.top_level_program,
            "filename": self.filename,
            "record_index": self.record_index,
            "record_id": self.record_id,
            "structure_index": self.structure_index,
            "structure_id": self.structure_id,
            "molecule_id": self.molecule_id,
            "label": self.label,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "functional": self.functional,
            "ab_initio": self.ab_initio,
            "basis": self.basis,
            "aux_basis": self.aux_basis,
            "extrapolation_basis": self.extrapolation_basis,
            "defgrid": self.defgrid,
            "scf_tol": self.scf_tol,
            "scf_algorithm": self.scf_algorithm,
            "solvent_model": self.solvent_model,
            "solvent_id": self.solvent_id,
            "route_parameters": self.route_parameters,
            "opt_options": self.opt_options,
            "structural_options": dict(self.structural_options),
            "resources": dict(self.resources),
            "warnings": list(self.warnings),
        }


__all__ = ["ParsedModelCommand"]
