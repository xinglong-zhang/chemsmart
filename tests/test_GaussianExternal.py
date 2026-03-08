"""
Integration tests for Gaussian External ML potential calculator support.

Covers the full round-trip:
  GaussianExternalJobSettings  → correct route string
  ASEExternalCalculatorScript  → script written, executable, correct content
  GaussianExternalJob          → script written into job directory
  Generated script             → reads Gaussian External input, runs ASE
                                 calculator, writes correctly formatted output
"""

import os
import shutil
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from chemsmart.io.gaussian.external_script import ASEExternalCalculatorScript
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.external import GaussianExternalJob
from chemsmart.jobs.gaussian.settings import GaussianExternalJobSettings
from chemsmart.settings.gaussian import GaussianProjectSettings

# ── Stub calculator without todict() (used for the pickle-fallback test) ─
class _NoDictCalc:
    """Minimal pickle-able calculator stub that has no todict() method."""
    pass


# ── Unit conversions (shared with the generated script) ───────────────────
_BOHR_TO_ANG = 0.529177210903
_EV_PER_HA = 27.211386245988
_GRAD_CONV = _BOHR_TO_ANG / _EV_PER_HA  # (eV/Å) → (Ha/Bohr)


# ── Helpers ───────────────────────────────────────────────────────────────

def _write_gaussian_external_input(path, positions_ang, atomic_numbers, deriv=1):
    """Write a minimal Gaussian External input file (positions in Bohr)."""
    natoms = len(atomic_numbers)
    with open(path, "w") as fh:
        fh.write(f"{natoms:10d}{deriv:10d}{0:10d}{1:10d}\n")
        for z, (x, y, z_coord) in zip(atomic_numbers, positions_ang):
            xb = x / _BOHR_TO_ANG
            yb = y / _BOHR_TO_ANG
            zb = z_coord / _BOHR_TO_ANG
            fh.write(
                f"{z:10d}{xb:20.12f}{yb:20.12f}{zb:20.12f}{0.0:20.12f}\n"
            )


# ── Tests ─────────────────────────────────────────────────────────────────

class TestGaussianExternalIntegration:
    """End-to-end tests for the Gaussian External calculator interface."""

    # ── 1. Settings and route string ──────────────────────────────────────

    def test_route_string_opt(self):
        """opt job produces '# opt External=<script>' route."""
        s = GaussianExternalJobSettings(
            external_script="run_mace.py",
            jobtype="opt",
            charge=0,
            multiplicity=1,
        )
        route = s.route_string
        assert "opt" in route
        assert "External=run_mace.py" in route
        assert "functional" not in route.lower()
        assert "basis" not in route.lower()

    def test_route_string_sp(self):
        """sp job produces energy-only External route (no opt keyword)."""
        s = GaussianExternalJobSettings(
            external_script="run_mace.py",
            jobtype="sp",
            charge=0,
            multiplicity=1,
        )
        route = s.route_string
        assert "opt" not in route
        assert "External=run_mace.py" in route

    def test_route_string_extra_args_quoted(self):
        """Extra args produce a quoted External= value in Gaussian format."""
        s = GaussianExternalJobSettings(
            external_script="RunTink",
            extra_args="Amber",
            jobtype="opt",
            charge=0,
            multiplicity=1,
        )
        # Gaussian convention: External="<script> <extra_args>"
        assert 'External="RunTink Amber"' in s.route_string

    def test_project_settings_external(self):
        """GaussianProjectSettings.external_settings() returns correct type."""
        ps = GaussianProjectSettings()
        es = ps.external_settings()
        assert isinstance(es, GaussianExternalJobSettings)
        assert es.freq is False
        assert es.forces is False

    # ── 2. Script generation ──────────────────────────────────────────────

    def test_script_written_and_executable(self, tmp_path):
        """Code-generation path: script file is written and executable."""
        from ase.calculators.emt import EMT

        script = ASEExternalCalculatorScript(
            EMT, calc_kwargs={}, script_name="run_emt"
        )
        script_path = script.write(str(tmp_path))

        assert os.path.isfile(script_path)
        assert os.access(script_path, os.X_OK)
        content = open(script_path).read()
        assert "#!/usr/bin/env python" in content
        assert "from ase.calculators.emt import EMT" in content
        assert "_read_input" in content
        assert "_write_output" in content
        assert "GRAD_CONV" in content

    def test_pickle_fallback_writes_pkl(self, tmp_path):
        """Pickle path: .pkl file is written when calculator has no todict()."""
        calc = _NoDictCalc()
        script = ASEExternalCalculatorScript(calc, script_name="run_pkl")
        script.write(str(tmp_path))

        assert os.path.isfile(tmp_path / "run_pkl.py")
        assert os.path.isfile(tmp_path / "run_pkl.pkl")

    # ── 3. GaussianExternalJob.write_external_script() ───────────────────

    def test_job_write_external_script(
        self, tmp_path, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """GaussianExternalJob writes the script into its working directory."""
        from ase.calculators.emt import EMT

        settings = GaussianExternalJobSettings(
            external_script="run_emt.py",
            jobtype="opt",
            charge=0,
            multiplicity=1,
            title="EMT opt integration test",
        )
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        script_writer = ASEExternalCalculatorScript(
            EMT, calc_kwargs={}, script_name="run_emt"
        )
        job = GaussianExternalJob(
            molecule=mol,
            settings=settings,
            label="test_mol",
            jobrunner=gaussian_jobrunner_no_scratch,
            external_script=script_writer,
        )
        written = job.write_external_script(directory=str(tmp_path))

        assert written is not None
        assert os.path.isfile(written)
        assert os.access(written, os.X_OK)

    # ── 4. Full round-trip: script execution ─────────────────────────────

    def test_script_execution_energy_and_gradient(self, tmp_path):
        """
        End-to-end: generated EMT script reads Gaussian External input
        and writes correctly formatted energy + gradient output.

        Uses ASE's EMT (Effective Medium Theory) calculator — a bundled
        force field that requires no external model files.  Au2 is a
        valid EMT system.
        """
        from ase.calculators.emt import EMT

        # Write the script
        script = ASEExternalCalculatorScript(
            EMT, calc_kwargs={}, script_name="run_emt"
        )
        script_path = script.write(str(tmp_path))

        # Prepare file paths Gaussian would supply
        input_file = str(tmp_path / "ext.dat")
        output_file = str(tmp_path / "ext.out")
        msg_file = str(tmp_path / "ext.msg")
        fchk_file = str(tmp_path / "ext.fchk")
        matel_file = str(tmp_path / "ext.matel")

        # Au2 molecule: two gold atoms separated by 2.88 Å
        atomic_numbers = [79, 79]
        positions_ang = [[0.0, 0.0, 0.0], [0.0, 0.0, 2.88]]
        _write_gaussian_external_input(
            input_file, positions_ang, atomic_numbers, deriv=1
        )

        # Run the generated script exactly as Gaussian would
        result = subprocess.run(
            [sys.executable, script_path,
             "R", input_file, output_file, msg_file, fchk_file, matel_file],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"Script exited with code {result.returncode}.\n"
            f"stderr: {result.stderr}\n"
            f"msg: {open(msg_file).read() if os.path.exists(msg_file) else 'no msg file'}"
        )

        # Parse output: line 0 = energy+dipole (4×D20.12), lines 1-2 = gradient
        with open(output_file) as fh:
            lines = fh.readlines()

        assert len(lines) == 3, "Expected 1 energy line + 2 gradient lines"

        energy_line = lines[0].split()
        assert len(energy_line) == 4
        energy_ha = float(energy_line[0])
        assert np.isfinite(energy_ha), "Energy should be finite"

        for i, line in enumerate(lines[1:], start=1):
            grad = line.split()
            assert len(grad) == 3, f"Gradient line {i} should have 3 components"
            assert all(np.isfinite(float(g)) for g in grad)

        # Sanity-check units: EMT Au2 energy is ~few eV; in Ha that's ~0.1-1 Ha
        assert abs(energy_ha) < 10.0, "Energy unexpectedly large; check units"


# ── Real Gaussian run ─────────────────────────────────────────────────────

def _find_g16():
    """Return the g16 executable path, or None if Gaussian is not available."""
    on_path = shutil.which("g16")
    if on_path:
        return on_path
    gauss_exedir = os.environ.get("GAUSS_EXEDIR", "")
    candidate = os.path.join(gauss_exedir, "g16")
    return candidate if os.path.isfile(candidate) else None


_G16 = _find_g16()


@pytest.mark.skipif(
    _G16 is None,
    reason="Gaussian g16 not available (add to PATH or set GAUSS_EXEDIR)",
)
class TestGaussianExternalRealRun:
    """
    Integration tests that invoke the real g16 binary.

    Each test is skipped automatically when Gaussian is not installed.

    Note on calculator choice
    -------------------------
    ASE's ``EMT`` (Effective Medium Theory) is parametrised only for FCC
    metals (Al, Cu, Ag, Au, Ni, Pd, Pt) and cannot evaluate H.  For H2
    we use ASE's built-in ``LennardJones`` calculator, which is general-
    purpose and requires no external model files.
    """

    def test_h2_lj_external_singlepoint(self, tmp_path):
        """
        Full pipeline: LennardJones single-point on H2 via Gaussian External.

        Steps
        -----
        1. ``ASEExternalCalculatorScript`` writes ``run_lj.py`` (LJ params
           chosen so H2 at 0.74 Å sits in the well region).
        2. A minimal Gaussian ``.com`` file is written whose route line
           is ``# External="<python_exe> run_lj.py"``.  Using the explicit
           interpreter path avoids shebang / conda-environment mismatches.
        3. ``g16`` is invoked; it calls ``run_lj.py`` with the six standard
           Gaussian External arguments.
        4. The ``.log`` is checked for ``"Normal termination"`` and a
           reported energy value.
        """
        from ase.calculators.lj import LennardJones

        # ── 1. Write the ASE External script ─────────────────────────────
        # sigma=0.5 Å places H2 (r=0.74 Å) in the attractive LJ well.
        script_writer = ASEExternalCalculatorScript(
            LennardJones,
            calc_kwargs={"sigma": 0.5, "epsilon": 0.01},
            script_name="run_lj",
        )
        script_path = script_writer.write(str(tmp_path))

        # ── 2. Write the Gaussian .com input ─────────────────────────────
        # Embed the full interpreter path so Gaussian invokes the right Python,
        # regardless of what "python" resolves to in the shell that runs g16.
        # Gaussian calls: <python_exe> <script_path> layer In Out Msg Fchk MatEl
        external_value = f'"{sys.executable} {script_path}"'

        com_content = textwrap.dedent(f"""\
            %chk=h2_ext.chk
            %nprocshared=1
            %mem=1GB
            # External={external_value}

            H2 LJ external single-point

            0 1
            H   0.000000   0.000000   0.000000
            H   0.000000   0.000000   0.740000

        """)
        com_file = tmp_path / "h2_ext.com"
        com_file.write_text(com_content)

        # ── 3. Run g16 ───────────────────────────────────────────────────
        result = subprocess.run(
            [_G16, str(com_file)],
            capture_output=True,
            text=True,
            cwd=str(tmp_path),
            timeout=120,
        )

        log_file = tmp_path / "h2_ext.log"
        log_text = log_file.read_text() if log_file.exists() else "(no log)"

        # ── 4. Assert normal termination and reported energy ─────────────
        assert "Normal termination" in log_text, (
            "Gaussian did not terminate normally.\n"
            f"route: {external_value}\n"
            f"g16 stdout: {result.stdout[-400:]}\n"
            f"log tail:\n{log_text[-1500:]}"
        )

        # Gaussian echoes the energy it received from the External script.
        assert "Energy=" in log_text, (
            "No 'Energy=' line found in Gaussian log; "
            "External script may not have written output correctly.\n"
            f"log tail:\n{log_text[-1500:]}"
        )
