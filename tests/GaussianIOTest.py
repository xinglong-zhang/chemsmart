import os.path

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.gaussian.cube import GaussianCubeFile

class TestGaussian16Output:
    def test_normal_termination_with_forces_and_frequencies(self, td_outputfile):
        assert os.path.exists(td_outputfile)
        g16_output = Gaussian16Output(filename=td_outputfile)
        print(g16_output.tddft_transitions)
        print(g16_output.excitation_energies_eV)


class TestGaussianCubeFile:
    def test_read_file_content(self, spin_cube_file):
        spin_cube = GaussianCubeFile(filename=spin_cube_file)

        print(spin_cube.filepath)

        assert spin_cube.cube_job_title == 'Gaussian job density'
        assert spin_cube.cube_job_description == 'Electron density from Total SCF Density'
        assert spin_cube.num_atoms == 2
        assert spin_cube.coordinate_origin == (-5.483229, -5.483229, -6.522947)
        assert type(spin_cube.coordinate_origin) is tuple
        assert spin_cube.grid_points == (9, 9, 11)
        assert type(spin_cube.grid_points) is tuple
        assert spin_cube.grid_increment_vector == ((1.2911, 0.0, 0.0), (0.0, 1.2911, 0.0), (0.0, 0.0, 1.2911))
        print(type(spin_cube.grid_increment_vector[0]))
        print(spin_cube.values_by_lines)

