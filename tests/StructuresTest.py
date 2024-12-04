import os
import numpy as np
from chemsmart.io.molecules.structure import CoordinateBlock
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.molecules.structure import XYZFile
from pyatoms.cli.traj import index
from tests.conftest import reference_genecp_txt_file_from_api


class TestCoordinateBlock:
    def test_read_coordinate_block(self):
        coordinates_string = """
C       -0.5448210000   -1.1694570000    0.0001270000
C        0.8378780000   -1.0476350000    0.0001900000
C        1.4329940000    0.2194290000    0.0001440000
C        0.6358350000    1.3657650000   -0.0000040000
C       -0.7521390000    1.2582750000    0.0000710000
C       -1.3283680000   -0.0113420000    0.0001560000
H       -1.0298620000   -2.1454490000    0.0001020000
H        1.4853190000   -1.9262430000    0.0002620000
H        1.1050940000    2.3527070000    0.0000530000
H       -1.3913950000    2.1406100000   -0.0000130000
C        2.9142600000    0.3363820000    0.0000140000
O        3.6625230000   -0.6037690000   -0.0002940000
H        3.3025560000    1.3842410000    0.0001510000
Cl      -3.0556310000   -0.1578960000   -0.0001400000
"""
        cb = CoordinateBlock(coordinate_block=coordinates_string)
        assert cb.symbols.get_chemical_formula() == "C7H5ClO"

    def test_read_gaussian_cb_with_tv(self):
        coordinates_string = """
C                  0.000000    0.000000    0.000000
C                  0.000000    1.429118    0.000000
TV                 2.475315    0.000000    0.000000
TV                -1.219952    2.133447    0.000000
"""
        cb = CoordinateBlock(coordinate_block=coordinates_string)
        assert cb.symbols.get_chemical_formula() == "C2"
        assert cb.translation_vectors == [
            [2.475315, 0.000000, 0.000000],
            [-1.219952, 2.133447, 0.000000],
        ]
        assert cb.molecule.pbc_conditions == [True, True, False]
        assert cb.molecule.translation_vectors == [
            [2.475315, 0.000000, 0.000000],
            [-1.219952, 2.133447, 0.000000],
        ]

    def test_read_gaussian_cb_frozen_atoms(self):
        coordinates_string = """
C        -1      -0.5448210000   -1.1694570000    0.0001270000
C        -1       0.8378780000   -1.0476350000    0.0001900000
C        -1       1.4329940000    0.2194290000    0.0001440000
C        -1       0.6358350000    1.3657650000   -0.0000040000
C        -1      -0.7521390000    1.2582750000    0.0000710000
C        -1      -1.3283680000   -0.0113420000    0.0001560000
H        -1      -1.0298620000   -2.1454490000    0.0001020000
H        -1       1.4853190000   -1.9262430000    0.0002620000
H        -1       1.1050940000    2.3527070000    0.0000530000
H        -1      -1.3913950000    2.1406100000   -0.0000130000
C        0       2.9142600000    0.3363820000    0.0000140000
O        0       3.6625230000   -0.6037690000   -0.0002940000
H        0       3.3025560000    1.3842410000    0.0001510000
Cl       0      -3.0556310000   -0.1578960000   -0.0001400000
"""
        cb = CoordinateBlock(coordinate_block=coordinates_string)
        assert cb.symbols.get_chemical_formula() == "C7H5ClO"
        assert cb.molecule.empirical_formula == "C7H5ClO"
        assert cb.molecule.frozen_atoms == [
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
        ]


class TestStructures:
    def test_read_molecule_from_single_molecule_xyz_file(
        self, single_molecule_xyz_file
    ):
        assert os.path.exists(single_molecule_xyz_file)
        assert os.path.isfile(single_molecule_xyz_file)

        xyz_file = XYZFile(filename=single_molecule_xyz_file)
        assert xyz_file.num_atoms == 71

        molecule = xyz_file.get_molecule(index="-1", return_list=False)
        assert isinstance(molecule, Molecule)

        molecule = xyz_file.get_molecule(index="-1", return_list=True)
        assert isinstance(molecule, list)
        assert len(molecule) == 1

        # molecule creation from path
        molecule = Molecule.from_filepath(
            single_molecule_xyz_file, return_list=False
        )
        assert isinstance(molecule, Molecule)

    def test_read_molecule_from_multiple_molecules_xyz_file(
        self, multiple_molecules_xyz_file
    ):
        assert os.path.exists(multiple_molecules_xyz_file)
        assert os.path.isfile(multiple_molecules_xyz_file)

        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        assert xyz_file.num_atoms == 71

        molecules = xyz_file.get_molecule(index=":", return_list=True)
        assert isinstance(molecules, list)
        assert len(molecules) == 18

        # obtain the last structure as molecule
        molecule = xyz_file.get_molecule(index="-1", return_list=False)
        assert isinstance(molecule, Molecule)
        assert molecule.empirical_formula == "C37H25Cl3N3O3"

        last_structure_positions = np.array(
            [
                [0.2264896660, -2.2726277143, -1.5005807047],
                [-1.1201033733, -1.8955122371, -1.3945085856],
                [0.9768480935, -2.8512670618, -0.3694870257],
                [-1.8258763303, -2.0106412689, -0.1062680111],
                [0.8617890058, -2.1432555829, -2.7332062296],
                [-1.7985851393, -1.4449652357, -2.5242164995],
                [0.4954580021, -3.9625367099, 0.3354063501],
                [2.2208677741, -2.3349087185, -0.0222807215],
                [0.1898630484, -1.6553784729, -3.8375759399],
                [-1.1484756939, -1.3091786134, -3.7357140585],
                [-2.9667193725, -2.7987082534, 0.0061646581],
                [-1.3141774390, -1.3966239933, 1.0467894815],
                [1.2339592377, -4.4830126857, 1.3905518197],
                [2.9482376828, -2.8511344152, 1.0317811065],
                [2.4461453107, -3.9241196042, 1.7474253947],
                [-3.5614407324, -3.0062761440, 1.2381549525],
                [-1.8827872905, -1.6603232771, 2.2933275850],
                [-3.0075412614, -2.4554483696, 2.3852692127],
                [1.8889420557, -2.4626415472, -2.8227646001],
                [-2.8491801264, -1.1996196224, -2.4455556214],
                [2.6144774694, -1.5031645558, -0.5857109595],
                [0.7001171594, -1.5677225836, -4.7842753217],
                [-1.6857700349, -0.9538644919, -4.6019587068],
                [-3.3572015271, -3.2797676010, -0.8772618514],
                [0.8506101699, -5.3402890161, 1.9275198344],
                [3.9019618146, -2.4188353570, 1.2894812560],
                [3.0020044917, -4.3366944551, 2.5745065021],
                [-4.4438861839, -3.6232549060, 1.3087192078],
                [-1.4460969280, -1.2127421123, 3.1725494581],
                [-3.4525784651, -2.6464262102, 3.3488105304],
                [-0.1865105937, -0.4538291540, 0.9887889655],
                [0.7688503822, -0.4765790832, 1.7264536163],
                [-0.6791027353, -4.5285401175, -0.0526063825],
                [-0.8344011095, -5.3335975978, 0.4570625557],
                [-0.2990824532, 0.6271820181, -0.0537026754],
                [0.7404894119, 1.2515803582, -0.6236564731],
                [-1.4005328482, 1.1173160931, -0.6640490756],
                [0.3377253026, 2.0635287872, -1.5990515584],
                [2.0966797201, 1.2303778522, -0.2473210677],
                [-0.9575730878, 1.9861507862, -1.6160079666],
                [-2.7967025857, 0.9928729668, -0.3155960225],
                [2.4790152573, 1.6497592409, 1.0328099240],
                [3.0833087680, 0.8900070398, -1.1793624563],
                [-1.9165764026, 2.7033420477, -2.4951424300],
                [-3.0596448982, -0.0671817865, -0.2414221278],
                [-3.7101642220, 1.6900284738, -1.3598720453],
                [-3.1584471615, 1.7036314935, 0.9716532507],
                [3.8095958751, 1.6183400087, 1.4096384082],
                [1.3108023677, 2.2858125753, 2.1092366311],
                [4.4161194014, 0.8616582054, -0.8065088124],
                [2.6310042273, 0.5028118970, -2.7888184432],
                [-2.0449496531, 2.1228196918, -3.4255189220],
                [-1.5334013669, 3.6929548225, -2.7515549748],
                [-3.1352385324, 2.8879140886, -1.8341783130],
                [-4.9525547499, 2.0893263022, -0.5508200660],
                [-3.9344004138, 1.0241590275, -2.2057714851],
                [-4.4076004213, 2.3020038695, 0.8310226638],
                [-2.4455006780, 1.7957561500, 2.1477016801],
                [4.7712350269, 1.2016499063, 0.4950718763],
                [4.0986677849, 1.9280440642, 2.4009935079],
                [5.1713381974, 0.5797256616, -1.5220545134],
                [-5.3941913179, 2.9926265804, -0.9679055458],
                [-5.6984020441, 1.2931485756, -0.5613787000],
                [-4.9682281649, 2.9891323043, 1.8887908748],
                [-3.0104573514, 2.4933743318, 3.2066378723],
                [-1.4664651823, 1.3511630586, 2.2521651249],
                [6.4120660120, 1.1438487037, 0.9633575644],
                [-4.2595413323, 3.0776346778, 3.0787905808],
                [-5.9348439371, 3.4565775747, 1.7874808818],
                [-2.4664057629, 2.5823290741, 4.1336857855],
                [-4.6820282248, 3.6172731919, 3.9117852027],
            ],
        )
        assert np.allclose(
            molecule.positions, last_structure_positions, rtol=1e-5
        )

        # obtain the last structure as a list
        molecule = xyz_file.get_molecule(index="-1", return_list=True)
        assert isinstance(molecule, list)
        assert len(molecule) == 1

        # obtain the last 10 structures
        molecules = xyz_file.get_molecule(index="-10:", return_list=True)
        assert isinstance(molecules, list)
        assert len(molecules) == 10

        # first coordinates of the 10th last structure
        assert np.allclose(
            molecules[0].positions[0],
            np.array([-0.5297471504, -3.4014322100, -1.3490458905]),
            rtol=1e-5,
        )

        # molecule creation from path
        molecules = Molecule.from_filepath(
            filepath=multiple_molecules_xyz_file, index=":", return_list=True
        )
        assert isinstance(molecules, list)
        assert len(molecules) == 18
        # obtain the last structure as molecule
        molecule = Molecule.from_filepath(
            filepath=multiple_molecules_xyz_file, index="-1", return_list=False
        )
        assert isinstance(molecule, Molecule)
        assert molecule.empirical_formula == "C37H25Cl3N3O3"
        assert np.allclose(
            molecule.positions, last_structure_positions, rtol=1e-5
        )
        # obtain the last structure as list
        molecule = Molecule.from_filepath(
            filepath=multiple_molecules_xyz_file, index="-1", return_list=True
        )
        assert isinstance(molecule, list)
        assert len(molecule) == 1
