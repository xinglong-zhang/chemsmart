from types import SimpleNamespace

import numpy as np

from chemsmart.analysis.result_readers import reader_for
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.input import ORCAInput
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.writer import ORCAInputWriter


def test_orca_vpt2_yaml_input_and_result_round_trip(tmp_path):
    settings = ORCAJobSettings(
        ab_initio="rhf",
        basis="def2-SVP",
        scf_convergence="extremely tight",
        vpt2=True,
        vpt2_anharmonic_displacement=0.05,
        vpt2_hessian_cutoff=1e-12,
        jobtype="sp",
        charge=0,
        multiplicity=1,
    )
    molecule = Molecule(
        symbols=["O", "H", "H"],
        positions=np.array(
            [[0.0, 0.06, 0.06], [0.0, -0.06, 1.0], [0.0, 1.0, -0.06]]
        ),
    )
    job = SimpleNamespace(
        settings=settings,
        molecule=molecule,
        jobrunner=SimpleNamespace(num_cores=4, mem_gb=4),
        folder=str(tmp_path),
        label="water_vpt2",
    )
    ORCAInputWriter(job).write()

    input_path = tmp_path / "water_vpt2.inp"
    rendered = input_path.read_text()
    assert "!  VPT2 rhf def2-SVP" in rendered
    assert " Opt " not in rendered
    assert "  convergence extreme" in rendered
    assert "  AnharmDisp 0.05" in rendered
    assert "  HessianCutoff 1e-12" in rendered

    observed = ORCAInput(str(input_path)).read_settings()
    assert observed.jobtype == "sp"
    assert observed.vpt2 is True
    assert observed.vpt2_anharmonic_displacement == 0.05
    assert observed.vpt2_hessian_cutoff == 1e-12
    assert observed.scf_convergence == "extreme"

    output_path = tmp_path / "water_vpt2.out"
    output_path.write_text("""Fundamental transitions [1/cm]
-----------------------------------------
Mode     w(harm)      v(fund)      Diff
-----------------------------------------
  0     1750.523     1692.012     -58.510
  1     4148.947     3981.879    -167.068
  2     4245.104     4071.127    -173.977
-----------------------------------------

Zero-point ro-vibrational energy [1/cm]
---------------------------------------
Harmonic contribution:         5072.287
Anharmonic correction:          -73.729
Ro-vibrational correction:        3.448
---------------------------------------
Total:                         5002.005
""")
    output = ORCAOutput(str(output_path))
    reader = reader_for("orca")
    assert reader.read(output, "vpt2_harmonic_frequencies") == (
        [1750.523, 4148.947, 4245.104],
        "cm^-1",
    )
    assert reader.read(output, "vpt2_fundamental_frequencies") == (
        [1692.012, 3981.879, 4071.127],
        "cm^-1",
    )
    assert reader.read(output, "vpt2_zero_point_rovibrational_energy") == (
        5002.005,
        "cm^-1",
    )
