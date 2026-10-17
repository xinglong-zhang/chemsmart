"""Formatted checkpoint Hessians and conversion failure handling."""

import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.io.gaussian.fchk import GaussianFormattedCheckpoint
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.mecp import GaussianMECPJob
from chemsmart.jobs.gaussian.runner import GaussianJobRunner
from chemsmart.jobs.gaussian.settings import GaussianMECPJobSettings

FCHK = (
    "Test checkpoint\nSP RHF STO-3G\n"
    "Number of atoms                            I                1\n"
    "Cartesian Force Constants                  R   N=           6\n"
    " 1.123456789012D+00 2.0E-02 3.0E+00 -4.0E-02 5.0E-02\n"
    " 6.0E+00\n"
)


def test_fchk_preserves_precision_and_unpacks_lower_triangle(tmp_path):
    path = tmp_path / "state.fchk"
    path.write_text(FCHK)
    checkpoint = GaussianFormattedCheckpoint(path)
    assert checkpoint.num_atoms == 1
    assert checkpoint.filepath == str(path.resolve())
    matrix = checkpoint.read_fchk_hessian()
    np.testing.assert_array_equal(
        matrix,
        [[1.123456789012, 0.02, -0.04], [0.02, 3.0, 0.05], [-0.04, 0.05, 6.0]],
    )


@pytest.mark.parametrize(
    "contents",
    [
        FCHK.replace("I                1", "I                2"),
        FCHK.replace("N=           6", "N=           5"),
        FCHK.replace(" 6.0E+00\n", ""),
        FCHK.replace("6.0E+00", "nan"),
        FCHK.replace("6.0E+00", "Next field I 3"),
        "No Hessian\n",
        FCHK.replace("I                1", "I                0"),
    ],
)
def test_fchk_rejects_invalid_data(tmp_path, contents):
    path = tmp_path / "state.fchk"
    path.write_text(contents)
    with pytest.raises(RuntimeError):
        GaussianFormattedCheckpoint(path).read_fchk_hessian()


@pytest.mark.parametrize("failure", [None, "process", "invalid", "atom_count"])
def test_conversion_uses_environment_and_never_reuses_old_output(
    tmp_path, monkeypatch, failure
):
    import chemsmart.jobs.gaussian.runner as module

    runner = object.__new__(GaussianJobRunner)
    runner._get_executable = lambda: "/configured Gaussian/g16"
    env = {"PATH": "/configured Gaussian", "GAUSS_EXEDIR": "configured"}
    runner._update_os_environ = lambda job: env
    monkeypatch.setattr(module.os, "access", lambda path, mode: True)
    chk = tmp_path / "state.chk"
    chk.write_bytes(b"new checkpoint")
    destination = tmp_path / "state.fchk"
    destination.write_text("old result")
    job = SimpleNamespace(
        chkfile=str(chk),
        folder=str(tmp_path),
        label="state",
        molecule=SimpleNamespace(
            symbols=["H", "H"] if failure == "atom_count" else ["H"]
        ),
    )

    def run(command, **kwargs):
        assert command[:2] == ["/configured Gaussian/formchk", str(chk)]
        assert kwargs["env"] == env
        assert kwargs["check"] is True
        if failure == "process":
            raise subprocess.CalledProcessError(1, command)
        Path(command[2]).write_text(
            "invalid" if failure == "invalid" else FCHK
        )

    monkeypatch.setattr(module.subprocess, "run", run)
    if failure:
        with pytest.raises(RuntimeError):
            runner.format_hessian_checkpoint(job)
        assert destination.read_text() == "old result"
    else:
        result = runner.format_hessian_checkpoint(job)
        assert result[0, 0] == 1.123456789012
        assert destination.read_text() == FCHK
    assert not list(tmp_path.glob("formchk-*"))


def test_missing_reference_checkpoint_prevents_analytic_job(tmp_path):
    job = object.__new__(GaussianMECPJob)
    job.settings = GaussianMECPJobSettings(functional="b3lyp", basis="6-31g")
    job.molecule = Molecule(symbols=["H"], positions=np.zeros((1, 3)))
    job.label = "test"
    job._state_checkpoint_files = {
        ("check", "A"): str(tmp_path / "missing.chk")
    }
    with pytest.raises(
        RuntimeError, match="requires the reference checkpoint"
    ):
        job._run_state(
            np.zeros((1, 3)), "hessian_analytic", "A", "check", analytic=True
        )


def test_missing_analytic_checkpoint_prevents_conversion(tmp_path):
    runner = object.__new__(GaussianJobRunner)
    with pytest.raises(RuntimeError, match="Missing analytic checkpoint"):
        runner.format_hessian_checkpoint(
            SimpleNamespace(chkfile=str(tmp_path / "missing.chk"))
        )
