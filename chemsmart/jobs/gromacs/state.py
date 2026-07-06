"""Workflow state for GROMACS multi-step setup jobs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from chemsmart.jobs.gromacs.job import GromacsJob


@dataclass(slots=True)
class GromacsWorkflowState:
    """
    Store intermediate file names for a GROMACS full setup workflow.

    Keeping these files in one state object avoids spreading intermediate file
    names across job, runner and CLI layers.
    """

    working_dir: Path
    label: str

    processed_structure_file: Path
    boxed_structure_file: Path
    solvated_structure_file: Path
    ions_tpr_file: Path
    ionized_structure_file: Path
    em_tpr_file: Path

    @classmethod
    def from_job(cls, job: GromacsJob) -> GromacsWorkflowState:
        """
        Build workflow state from a GROMACS job.
        """
        working_dir = Path(job.folder).resolve()
        label = job.label

        return cls(
            working_dir=working_dir,
            label=label,
            processed_structure_file=Path(
                job.processed_structure_file or working_dir / "processed.gro"
            ),
            boxed_structure_file=Path(
                job.boxed_structure_file or working_dir / "boxed.gro"
            ),
            solvated_structure_file=Path(
                job.solvated_structure_file or working_dir / "solvated.gro"
            ),
            ions_tpr_file=Path(job.ions_tpr_file or working_dir / "ions.tpr"),
            ionized_structure_file=Path(
                job.ionized_structure_file or working_dir / "ionized.gro"
            ),
            em_tpr_file=Path(job.tpr_file or working_dir / f"{label}.tpr"),
        )
