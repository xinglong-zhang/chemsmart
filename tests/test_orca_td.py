from types import SimpleNamespace

import numpy as np

from chemsmart.analysis.result_readers import reader_for
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.input import ORCAInput
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.writer import ORCAInputWriter


def test_orca_tda_input_and_spectrum_round_trip(tmp_path):
    settings = ORCAJobSettings(
        functional="BP86",
        basis="def2-SVP",
        aux_basis="def2/J",
        scf_tol="TightSCF",
        jobtype="td",
        charge=0,
        multiplicity=1,
        response_method="tda",
        nstates=3,
        state_manifold="singlet",
    )
    molecule = Molecule(
        symbols=["C", "O", "H", "H"],
        positions=np.array(
            [
                [0.0, 0.0, -0.523841],
                [0.0, 0.0, 0.676159],
                [0.0, 0.939547, -1.109428],
                [0.0, -0.939547, -1.109428],
            ]
        ),
    )
    job = SimpleNamespace(
        settings=settings,
        molecule=molecule,
        jobrunner=SimpleNamespace(num_cores=4, mem_gb=4),
        folder=str(tmp_path),
        label="formaldehyde_tda",
    )
    ORCAInputWriter(job).write()

    rendered = (tmp_path / "formaldehyde_tda.inp").read_text()
    assert "%tddft" in rendered
    assert "  NRoots 3" in rendered
    assert "  TDA true" in rendered
    assert "  Triplets false" in rendered
    assert " Opt" not in rendered

    observed = ORCAInput(
        str(tmp_path / "formaldehyde_tda.inp")
    ).read_settings()
    assert observed.jobtype == "td"
    assert observed.response_method == "tda"
    assert observed.nstates == 3
    assert observed.state_manifold == "singlet"

    output_path = tmp_path / "formaldehyde_tda.out"
    output_path.write_text("""|  1> ! BP86 def2-SVP def2/J TightSCF
|  2> %tddft
|  3> nroots 3
|  4> tda true
|  5> triplets false
|  6> end
----------------------------------------------------------------------------------------------------
                     ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
----------------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ
                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)
----------------------------------------------------------------------------------------------------
  0-1A  ->  1-1A    3.992771   32203.9   310.5   0.000000000   0.00000  0.0  0.0  0.0
  0-1A  ->  2-1A    7.677957   61926.9   161.5   0.144293775   0.76709  0.0  0.0  0.0
  0-1A  ->  3-1A    9.109781   73475.3   136.1   0.002190997   0.00982  0.0  0.0  0.0
----------------------------------------------------------------------------------------------------
                             ****ORCA TERMINATED NORMALLY****
""")
    output = ORCAOutput(str(output_path))
    reader = reader_for("orca")
    assert reader.read(output, "excitation_energies") == (
        [3.992771, 7.677957, 9.109781],
        "eV",
    )
    assert reader.read(output, "absorption_wavelengths") == (
        [310.5, 161.5, 136.1],
        "nm",
    )
    assert reader.read(output, "oscillator_strengths") == (
        [0.0, 0.144293775, 0.002190997],
        "",
    )
