from pathlib import Path

from chemsmart.jobs.gromacs.job import GromacsJob
from chemsmart.jobs.gromacs.workflow import GromacsWorkflow


def test_workflow_passes_stage_outputs(
    monkeypatch,
    tmp_path,
):
    initial_structure = tmp_path / "input.gro"
    topology = tmp_path / "topol.top"

    initial_structure.write_text(
        "structure",
        encoding="utf-8",
    )
    topology.write_text(
        "topology",
        encoding="utf-8",
    )

    calls = []

    def fake_run(job):
        calls.append(
            {
                "type": job.TYPE,
                "structure_file": job.structure_file,
                "checkpoint_file": job.checkpoint_file,
                "top_file": job.top_file,
            }
        )

        stem = Path(job.tpr_file).with_suffix("")

        stem.with_suffix(".tpr").write_text(
            "tpr",
            encoding="utf-8",
        )
        stem.with_suffix(".gro").write_text(
            "gro",
            encoding="utf-8",
        )
        stem.with_suffix(".cpt").write_text(
            "cpt",
            encoding="utf-8",
        )
        stem.with_suffix(".log").write_text(
            "Finished mdrun",
            encoding="utf-8",
        )

    monkeypatch.setattr(
        GromacsJob,
        "run",
        fake_run,
        raising=False,
    )

    workflow = GromacsWorkflow(
        structure_file=initial_structure,
        top_file=topology,
        jobrunner=object(),
        folder=tmp_path / "workflow",
    )

    jobs = workflow.run()

    assert list(jobs) == [
        "em",
        "nvt",
        "npt",
        "md",
    ]

    workflow_dir = tmp_path / "workflow"

    assert calls[0]["structure_file"] == initial_structure
    assert calls[0]["checkpoint_file"] is None

    assert calls[1]["structure_file"] == (
        workflow_dir / "em.gro"
    )
    assert calls[1]["checkpoint_file"] is None

    assert calls[2]["structure_file"] == (
        workflow_dir / "nvt.gro"
    )
    assert calls[2]["checkpoint_file"] == (
        workflow_dir / "nvt.cpt"
    )

    assert calls[3]["structure_file"] == (
        workflow_dir / "npt.gro"
    )
    assert calls[3]["checkpoint_file"] == (
        workflow_dir / "npt.cpt"
    )

    assert all(
        call["top_file"] == topology
        for call in calls
    )


def test_workflow_rejects_missing_initial_structure(
    tmp_path,
):
    topology = tmp_path / "topol.top"
    topology.write_text(
        "topology",
        encoding="utf-8",
    )

    missing_structure = tmp_path / "missing.gro"

    try:
        GromacsWorkflow(
            structure_file=missing_structure,
            top_file=topology,
            jobrunner=object(),
            folder=tmp_path / "workflow",
        )
    except FileNotFoundError as error:
        assert str(missing_structure) in str(error)
    else:
        raise AssertionError(
            "Expected FileNotFoundError"
        )