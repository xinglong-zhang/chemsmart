import shutil
import subprocess
from pathlib import Path

import pytest

from chemsmart.jobs.gromacs.job import GromacsEMJob
from chemsmart.jobs.gromacs.runner import GromacsJobRunner
from chemsmart.settings.gromacs import GromacsProjectSettings


pytestmark = pytest.mark.skipif(
    shutil.which("gmx") is None,
    reason="Real GROMACS executable 'gmx' is not available in PATH.",
)


def _write_demo_case(tmp_path):
    """
    Write a minimal prepared GROMACS EM demo case.

    This case is intentionally small so that it can be used for local
    integration testing. It uses one SPC water molecule and the built-in
    OPLS-AA force field files available in a standard GROMACS installation.
    """
    project_yaml = tmp_path / "project.yaml"
    mdp_file = tmp_path / "em.mdp"
    structure_file = tmp_path / "water.gro"
    top_file = tmp_path / "topol.top"

    project_yaml.write_text(
        """
project:
  name: water_prepared_em
  type: gromacs
  job_type: em
  mode: prepared

gromacs_settings:
  force_field: user_provided
  water_model: spc
  timestep: 0.001
  temperature: 300.0
  pressure: 1.0
  constraints: none
  grompp_maxwarn: 0
  mdrun_threads: 1

inputs:
  mdp_file: em.mdp
  structure_file: water.gro
  topology_file: topol.top
  tpr_file: em.tpr
""",
        encoding="utf-8",
    )

    mdp_file.write_text(
        """
integrator              = steep
nsteps                  = 50
emtol                   = 1000.0
emstep                  = 0.01

nstlog                  = 1
nstenergy               = 1

cutoff-scheme           = Verlet
nstlist                 = 1
coulombtype             = Cut-off
rcoulomb                = 0.8
vdwtype                 = Cut-off
rvdw                    = 0.8

pbc                     = xyz
constraints             = none
""",
        encoding="utf-8",
    )

    structure_file.write_text(
        """Single water molecule prepared EM demo
3
    1SOL     OW    1   0.000   0.000   0.000
    1SOL    HW1    2   0.096   0.000   0.000
    1SOL    HW2    3  -0.024   0.093   0.000
   2.00000   2.00000   2.00000
""",
        encoding="utf-8",
    )

    top_file.write_text(
        """
#include "oplsaa.ff/forcefield.itp"
#include "oplsaa.ff/spc.itp"

[ system ]
Single water molecule prepared EM demo

[ molecules ]
SOL 1
""",
        encoding="utf-8",
    )

    return project_yaml


def _make_real_runner(tmp_path):
    """
    Create a minimal real GROMACS runner for integration testing.

    __new__ is used to avoid depending on a real ChemSmart server object in
    this integration test. The goal here is to validate the GROMACS execution
    path, not the cluster submission layer.
    """
    runner = GromacsJobRunner.__new__(GromacsJobRunner)
    runner.gmx_executable = "gmx"
    runner.gmx_modules = []
    runner.gmx_env = {}
    runner.gmx_source_scripts = []
    runner.fake = False
    runner.running_directory = str(tmp_path)
    runner.job_outputfile = str(tmp_path / "water_prepared_em.out")
    runner.job_errfile = str(tmp_path / "water_prepared_em.err")
    return runner


def test_real_gromacs_prepared_em_from_project_yaml(tmp_path):
    """
    Test project.yaml -> settings -> job -> runner -> real grompp -> real mdrun.
    """
    project_yaml = _write_demo_case(tmp_path)

    settings = GromacsProjectSettings.from_yaml(project_yaml)
    settings.validate()

    job = GromacsEMJob.from_project_settings(
        settings=settings,
        molecule=None,
        jobrunner=None,
    )
    job.set_folder(str(tmp_path))

    runner = _make_real_runner(tmp_path)

    runner._assign_variables(job)
    runner._validate_gromacs_inputs(job)

    runner._assemble_tpr(job)
    assert (tmp_path / "em.tpr").exists()

    runner._run_command(
        runner._get_mdrun_command(job),
        cwd=tmp_path,
        check=True,
        stage="mdrun",
    )

    assert (tmp_path / "em.log").exists()
    assert (tmp_path / "em.edr").exists()
    assert (tmp_path / "em.gro").exists()


def test_real_gromacs_is_available():
    """
    Record the GROMACS version in the test log.
    """
    result = subprocess.run(
        ["gmx", "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )

    assert "GROMACS" in result.stdout