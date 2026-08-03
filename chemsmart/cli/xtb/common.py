"""Shared job construction for the xTB CLI leaves."""

import logging

logger = logging.getLogger(__name__)


def build_xtb_jobs(ctx, job_cls, settings, skip_completed, kwargs):
    jobrunner = ctx.obj["jobrunner"]
    molecules = ctx.obj["molecules"]
    molecule_indices = ctx.obj["molecule_indices"]
    label = ctx.obj["label"]
    source_filename = ctx.obj["filename"]
    project_reference = ctx.obj.get("project_reference")
    project_source_file = ctx.obj.get("project_source_file")
    explicit_state_fields = ctx.obj.get("explicit_state_fields", frozenset())

    def settings_for_molecule(molecule):
        resolved = settings.copy()
        if (
            "charge" not in explicit_state_fields
            and molecule.charge is not None
        ):
            resolved.charge = molecule.charge
        if (
            "multiplicity" not in explicit_state_fields
            and molecule.multiplicity is not None
        ):
            resolved.multiplicity = molecule.multiplicity
        resolved.validate_for_molecule(molecule)
        return resolved

    if len(molecules) > 1 and molecule_indices is not None:
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            molecule_label = f"{label}_idx{idx}"
            jobs.append(
                job_cls(
                    molecule=molecule,
                    settings=settings_for_molecule(molecule),
                    label=molecule_label,
                    jobrunner=jobrunner,
                    skip_completed=skip_completed,
                    source_filename=source_filename,
                    project_reference=project_reference,
                    project_source_file=project_source_file,
                    source_index=idx,
                    **kwargs,
                )
            )
        return jobs

    return job_cls(
        molecule=molecules[-1],
        settings=settings_for_molecule(molecules[-1]),
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        source_filename=source_filename,
        project_reference=project_reference,
        project_source_file=project_source_file,
        source_index=(
            molecule_indices[0]
            if molecule_indices is not None and len(molecule_indices) == 1
            else None
        ),
        **kwargs,
    )
