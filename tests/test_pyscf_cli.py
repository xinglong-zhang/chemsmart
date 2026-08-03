"""Focused tests for the PySCF Click surface and option propagation."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemsmart.cli.pyscf import pyscf
from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.settings.executable import GaussianExecutable, PySCFExecutable
from chemsmart.settings.server import Server

_JOB_CLASS = {
    "sp": "chemsmart.jobs.pyscf.singlepoint.PySCFSinglePointJob",
    "opt": "chemsmart.jobs.pyscf.opt.PySCFOptJob",
    "hess": "chemsmart.jobs.pyscf.hess.PySCFHessJob",
}


def test_only_library_backends_may_omit_exefolder(tmp_path):
    server_yaml = Path(tmp_path) / "server.yaml"
    server_yaml.write_text(
        "GAUSSIAN:\n  LOCAL_RUN: true\n"
        "PYSCF:\n  LOCAL_RUN: true\n",
        encoding="utf-8",
    )

    with pytest.raises(KeyError, match="EXEFOLDER"):
        GaussianExecutable.from_servername(str(server_yaml))
    assert (
        PySCFExecutable.from_servername(str(server_yaml)).executable_folder
        is None
    )


def _run_and_capture_settings(
    cli_args, jobtype="sp", num_gpus=0, preflight_violations=None
):
    runner = CliRunner()
    jobrunner = MagicMock(num_gpus=num_gpus)

    patches = [patch(_JOB_CLASS[jobtype])]
    # Stage B adds this symbol to the shared builder. Keeping the test helper
    # tolerant of both checkpoints lets the Stage A propagation tests remain
    # useful while the typed preflight is integrated.
    try:
        import chemsmart.cli.pyscf.common as common

        if hasattr(common, "preflight"):
            patches.append(
                patch.object(
                    common,
                    "preflight",
                    return_value=preflight_violations or [],
                )
            )
    except ImportError:
        pass

    entered = [manager.__enter__() for manager in patches]
    try:
        job_cls = entered[0]
        job_cls.return_value = MagicMock()
        result = runner.invoke(
            pyscf,
            cli_args,
            obj={"jobrunner": jobrunner},
            catch_exceptions=False,
        )
        settings = None
        if job_cls.call_args is not None:
            settings = job_cls.call_args.kwargs.get("settings")
        return result, settings
    finally:
        for manager in reversed(patches):
            manager.__exit__(None, None, None)


@pytest.mark.parametrize("jobtype", ["sp", "opt", "hess"])
def test_group_and_all_leaves_expose_help(single_molecule_xyz_file, jobtype):
    runner = CliRunner()

    group = runner.invoke(pyscf, ["--help"])
    leaf = runner.invoke(
        pyscf,
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            jobtype,
            "--help",
        ],
    )

    assert group.exit_code == 0, group.output
    assert leaf.exit_code == 0, leaf.output
    assert jobtype in group.output


@pytest.mark.parametrize(
    ("option", "field", "expected"),
    [
        (["-t", "receipt"], "title", "receipt"),
        (["-c", "1"], "charge", 1),
        (["-m", "2"], "multiplicity", 2),
        (["-x", "pbe0"], "functional", "pbe0"),
        (["-A", "hf"], "ab_initio", "hf"),
        (["-b", "def2-tzvp"], "basis", "def2-tzvp"),
        (
            ["--aux-basis", "def2-universal-jkfit"],
            "aux_basis",
            "def2-universal-jkfit",
        ),
        (["-d", "d3bj"], "dispersion", "d3bj"),
        (["--defgrid", "defgrid3"], "defgrid", "defgrid3"),
        (["--scf-tol", "1e-10"], "scf_tol", 1e-10),
        (["--scf-maxiter", "75"], "scf_maxiter", 75),
        (["--no-density-fit"], "density_fit", False),
        (["--opt-solver", "geometric"], "opt_solver", "geometric"),
        (["--opt-maxsteps", "40"], "opt_maxsteps", 40),
        (["-sm", "cpcm"], "solvent_model", "cpcm"),
        (["-si", "toluene"], "solvent_id", "toluene"),
    ],
)
def test_every_group_override_survives_the_merge(
    single_molecule_xyz_file, option, field, expected
):
    result, settings = _run_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            *option,
            "sp",
        ]
    )

    assert result.exit_code == 0, result.output
    assert settings is not None
    assert getattr(settings, field) == expected

    if field == "functional":
        assert settings.ab_initio is None
    if field == "ab_initio":
        assert settings.functional is None


@pytest.mark.parametrize(
    ("option", "expected"), [(["--gpu"], "gpu"), (["--no-gpu"], "cpu")]
)
def test_explicit_engine_override_wins(
    single_molecule_xyz_file, option, expected
):
    result, settings = _run_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            *option,
            "sp",
        ],
        num_gpus=1,
    )

    assert result.exit_code == 0, result.output
    assert settings.engine == expected


def test_explicit_gpu_rejects_zero_resolved_devices(single_molecule_xyz_file):
    with pytest.raises(ValueError, match="positive resolved NUM_GPUS"):
        _run_and_capture_settings(
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--gpu",
                "sp",
            ],
            num_gpus=0,
        )


@pytest.mark.parametrize("num_gpus", [0, 1, 4])
def test_omitted_gpu_flag_always_resolves_cpu(
    single_molecule_xyz_file, num_gpus
):
    result, settings = _run_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "sp",
        ],
        num_gpus=num_gpus,
    )

    assert result.exit_code == 0, result.output
    assert settings.engine == "cpu"


def test_remove_solvent_clears_project_single_point_solvent(
    single_molecule_xyz_file,
):
    result, settings = _run_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "--remove-solvent",
            "sp",
        ]
    )

    assert result.exit_code == 0, result.output
    assert (settings.solvent_model, settings.solvent_id) == (None, None)


@pytest.mark.parametrize("jobtype", ["sp", "opt", "hess"])
def test_leaf_selects_matching_project_settings(
    single_molecule_xyz_file, jobtype
):
    result, settings = _run_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            jobtype,
        ],
        jobtype=jobtype,
    )

    assert result.exit_code == 0, result.output
    assert settings.jobtype == jobtype
    assert settings.project_yaml_digest is not None


def test_sub_reconstruction_preserves_pyscf_run_settings(
    monkeypatch, single_molecule_xyz_file
):
    fake_server = Server(
        name="dummy", NUM_CORES=2, NUM_GPUS=0, MEM_GB=4
    )
    captured = {}

    def capture_submit(job, test=False, cli_args=None, **_kwargs):
        captured["submitted_job"] = job
        captured["test"] = test
        captured["cli_args"] = cli_args

    fake_server.submit = capture_submit
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )

    sub_result = CliRunner().invoke(
        sub,
        [
            "--test",
            "--server",
            "dummy",
            "--num-cores",
            "2",
            "--num-gpus",
            "0",
            "--mem-gb",
            "4",
            "pyscf",
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "--functional",
            "pbe0",
            "--basis",
            "def2-tzvp",
            "--scf-tol",
            "1e-10",
            "sp",
        ],
    )

    assert sub_result.exit_code == 0, sub_result.output
    assert captured["test"] is True
    assert "--test" not in captured["cli_args"]

    def capture_run_job(job, _jobrunner):
        captured["reconstructed_job"] = job

    monkeypatch.setattr(
        "chemsmart.cli.run._run_single_job", capture_run_job
    )
    run_result = CliRunner().invoke(run, captured["cli_args"])

    assert run_result.exit_code == 0, run_result.output
    submitted = captured["submitted_job"].settings
    reconstructed = captured["reconstructed_job"].settings
    for field in (
        "jobtype",
        "functional",
        "basis",
        "charge",
        "multiplicity",
        "scf_tol",
        "density_fit",
        "freq",
        "engine",
        "project_yaml_digest",
    ):
        assert getattr(reconstructed, field) == getattr(submitted, field)
