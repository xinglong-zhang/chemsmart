import os.path

from chemsmart.io.gaussian.output import Gaussian16Output

class TestGaussian16Output:
    def test_normal_termination_with_forces_and_frequencies(self, td_outputfile):
        assert os.path.exists(td_outputfile)
        g16_output = Gaussian16Output(filename=td_outputfile)
        print(g16_output.tddft_transitions)
        print(g16_output.excitation_energies_eV)