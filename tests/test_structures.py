import os
from filecmp import cmp
from pathlib import Path

import networkx as nx
import numpy as np
import pytest
from ase import Atoms
from pymatgen.core.structure import Molecule as PMGMolecule
from rdkit import Chem
from rdkit.Chem.rdchem import Mol as RDKitMolecule

from chemsmart.io.file import CDXFile
from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.io.molecules.structure import (
    CoordinateBlock,
    Molecule,
    QMMMMolecule,
)
from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.utils.cluster import (
    is_pubchem_api_available,
    is_pubchem_network_available,
)
from chemsmart.utils.utils import cmp_with_ignore


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
        coordinates_string1 = """
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
        cb1 = CoordinateBlock(coordinate_block=coordinates_string1)
        assert cb1.symbols.get_chemical_formula() == "C7H5ClO"
        assert cb1.molecule.empirical_formula == "C7H5ClO"
        assert cb1.molecule.frozen_atoms == [
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
        # assert cb1.molecule.partition_level_strings is None

        coordinates_string2 = """ 
  F     -1.041506214819     0.000000000000    -2.126109488809 M
  F     -2.033681935634    -1.142892069126    -0.412218766901 M
  F     -2.033681935634     1.142892069126    -0.412218766901 M
  C     -1.299038105677     0.000000000000    -0.750000000000 M H 5
  C      0.000000000000     0.000000000000     0.000000000000 H
  H      0.000000000000     0.000000000000     1.100000000000 H
  O      1.125833024920     0.000000000000    -0.650000000000 H
 """
        cb2 = CoordinateBlock(coordinate_block=coordinates_string2)
        assert cb2.symbols.get_chemical_formula() == "C2HF3O"
        assert cb2.molecule.empirical_formula == "C2HF3O"
        assert cb2.molecule.partition_level_strings == [
            "M",
            "M",
            "M",
            "M",
            "H",
            "H",
            "H",
        ]
        assert cb2.molecule.high_level_atoms == [5, 6, 7]
        assert cb2.molecule.medium_level_atoms == [1, 2, 3, 4]

        # TODO:

    def test_coordinate_block_without_partitions_returns_molecule(self):
        """Non-ONIOM coordinate block should return Molecule, not QMMMMolecule."""
        normal_block = [
            "C   0.000  0.000  0.000",
            "H   1.089  0.000  0.000",
            "H  -0.363  1.027  0.000",
            "H  -0.363 -0.513  0.890",
            "H  -0.363 -0.513 -0.890",
        ]
        cb = CoordinateBlock(coordinate_block=normal_block)
        mol = cb.molecule
        assert isinstance(mol, Molecule)
        assert not isinstance(mol, QMMMMolecule)
        assert type(mol).__name__ == "Molecule"

    def test_coordinate_block_with_partitions_returns_qmmm_molecule(self):
        """ONIOM coordinate block should return QMMMMolecule."""
        oniom_block = [
            "C   0.000  0.000  0.000  H",
            "H   1.089  0.000  0.000  L",
            "H  -0.363  1.027  0.000  L",
            "H  -0.363 -0.513  0.890  L",
            "H  -0.363 -0.513 -0.890  L",
        ]
        cb = CoordinateBlock(coordinate_block=oniom_block)
        mol = cb.molecule
        assert isinstance(mol, QMMMMolecule)
        assert type(mol).__name__ == "QMMMMolecule"

    def test_read_molecule_from_single_molecule_xyz_file(
        self, single_molecule_xyz_file
    ):
        assert os.path.exists(single_molecule_xyz_file)
        assert os.path.isfile(single_molecule_xyz_file)

        xyz_file = XYZFile(filename=single_molecule_xyz_file)
        assert xyz_file.num_atoms == 71

        molecule = xyz_file.get_molecules(index="-1", return_list=False)
        assert isinstance(molecule, Molecule)
        assert len(molecule.chemical_symbols) == 71
        assert molecule.is_chiral
        assert molecule.is_aromatic

        # test conversion of molecule to RDKit molecule
        rdkit_mol = molecule.to_rdkit()
        assert isinstance(rdkit_mol, Chem.Mol)
        assert rdkit_mol.GetNumAtoms() == 71
        assert rdkit_mol.GetNumConformers() == 1
        assert rdkit_mol.GetConformer().GetPositions().shape == (71, 3)

        # molecule is chiral
        assert Chem.FindMolChiralCenters(rdkit_mol, force=True) != []

        # convert to smiles string
        smiles = molecule.to_smiles()
        assert isinstance(smiles, str)

        # test conversion of molecule to graph
        graph = molecule.to_graph()
        assert isinstance(graph, nx.Graph)
        assert len(graph.nodes) == 71
        assert len(graph.edges) == 78  # CH4 should have 4 bonds

        molecule = xyz_file.get_molecules(index="-1", return_list=True)
        assert isinstance(molecule, list)
        assert len(molecule) == 1

        # molecule creation from path
        molecule = Molecule.from_filepath(
            single_molecule_xyz_file, return_list=False
        )
        assert isinstance(molecule, Molecule)
        assert getattr(molecule, "partition_level_strings", None) is None
        assert len(molecule.chemical_symbols) == 71
        assert molecule.empirical_formula == "C37H25Cl3N3O3"
        assert np.isclose(molecule.mass, 665.982, atol=1e-2)
        assert molecule.energy == -126.2575508

        # test conversion to RDKit molecule
        rdkit_molecule = molecule.to_rdkit()
        assert isinstance(rdkit_molecule, RDKitMolecule)

    def test_read_molecule_energy_from_xyz_file(
        self,
        xtb_optimized_xyz_file,
        chemsmart_generated_xyz_file,
        extended_xyz_file,
    ):
        assert os.path.exists(xtb_optimized_xyz_file)
        assert os.path.isfile(xtb_optimized_xyz_file)
        xyz_file1 = XYZFile(filename=xtb_optimized_xyz_file)
        assert xyz_file1.num_atoms == 508
        molecule1 = xyz_file1.get_molecules(index="-1", return_list=False)
        assert isinstance(molecule1, Molecule)
        assert len(molecule1.chemical_symbols) == 508
        assert molecule1.energy == -978.449085030052

        assert os.path.exists(chemsmart_generated_xyz_file)
        assert os.path.isfile(chemsmart_generated_xyz_file)
        xyz_file2 = XYZFile(filename=chemsmart_generated_xyz_file)
        assert xyz_file2.num_atoms == 14
        molecule2 = xyz_file2.get_molecules(index="-1", return_list=False)
        assert isinstance(molecule2, Molecule)
        assert len(molecule2.chemical_symbols) == 14
        assert molecule2.energy == -804.614711

        assert os.path.exists(extended_xyz_file)
        assert os.path.isfile(extended_xyz_file)
        xyz_file3 = XYZFile(filename=extended_xyz_file)
        assert xyz_file3.num_atoms == 1
        molecule3 = xyz_file3.get_molecules(index="-1", return_list=False)
        assert isinstance(molecule3, Molecule)
        assert len(molecule3.chemical_symbols) == 1
        assert molecule3.energy == -157.72725320

    def test_read_molecule_from_multiple_molecules_xyz_file(
        self, multiple_molecules_xyz_file
    ):
        assert os.path.exists(multiple_molecules_xyz_file)
        assert os.path.isfile(multiple_molecules_xyz_file)

        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        assert xyz_file.num_atoms == 71

        all_molecules = xyz_file.get_molecules(index=":", return_list=True)

        # set correct charge and multiplicity for
        # molecules as needed by pymatgen checks
        molecules = []
        for molecule in all_molecules:
            molecule.charge = 1
            molecule.multiplicity = 1
            molecules.append(molecule)

        assert isinstance(molecules, list)
        assert len(molecules) == 18

        # test molecule bond orders
        first_mol = molecules[0]
        last_mol = molecules[-1]
        assert isinstance(first_mol, Molecule)
        assert isinstance(last_mol, Molecule)
        first_bond_orders = first_mol.get_bond_orders_from_rdkit_mol(
            bond_cutoff_buffer=0.0
        )
        last_bond_orders = last_mol.get_bond_orders_from_rdkit_mol(
            bond_cutoff_buffer=0.0
        )
        assert first_bond_orders == last_bond_orders

        # note that for conformers, due to the buffer values, the bond orders
        # may have changed, due to different distances in conformers,
        #  although, this is not supposed to be
        assert first_mol.bond_orders != last_mol.bond_orders

        first_rdkit_mol = first_mol.to_rdkit()
        last_rdkit_mol = last_mol.to_rdkit()
        assert isinstance(first_rdkit_mol, Chem.Mol)
        assert isinstance(last_rdkit_mol, Chem.Mol)

        # convert rdkit molecule back to Molecule object
        first_mol_conv = Molecule.from_rdkit_mol(first_rdkit_mol)
        last_mol_conv = Molecule.from_rdkit_mol(last_rdkit_mol)
        assert isinstance(first_mol_conv, Molecule)
        assert isinstance(last_mol_conv, Molecule)
        first_bond_orders = first_mol_conv.get_bond_orders_from_rdkit_mol(
            bond_cutoff_buffer=0.0
        )
        last_bond_orders = last_mol_conv.get_bond_orders_from_rdkit_mol(
            bond_cutoff_buffer=0.0
        )
        assert first_bond_orders == last_bond_orders
        assert np.all(first_mol.symbols == first_mol_conv.symbols)
        assert np.all(last_mol.symbols == last_mol_conv.symbols)
        assert np.all(first_mol.positions == first_mol_conv.positions)
        assert np.all(last_mol.positions == last_mol_conv.positions)
        assert first_mol.num_atoms == first_mol_conv.num_atoms == 71
        assert last_mol.num_atoms == last_mol_conv.num_atoms == 71

        # test conversion to ase Atoms
        first_ase_atoms = first_mol.to_ase()
        last_ase_atoms = last_mol.to_ase()
        assert isinstance(first_ase_atoms, Atoms)
        assert isinstance(last_ase_atoms, Atoms)

        # test conversion to pymatgen Structure
        assert first_mol.charge == 1
        assert first_mol.multiplicity == 1
        assert first_ase_atoms.charge == 1
        assert first_ase_atoms.multiplicity == 1
        first_py_structure = first_mol.to_pymatgen()
        last_py_structure = last_mol.to_pymatgen()
        assert isinstance(first_py_structure, PMGMolecule)
        assert isinstance(last_py_structure, PMGMolecule)
        assert first_py_structure.charge == 1
        assert first_py_structure.spin_multiplicity == 1

        # obtain the last structure as molecule
        molecule = xyz_file.get_molecules(index="-1", return_list=False)
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
        molecule = xyz_file.get_molecules(index="-1", return_list=True)
        assert isinstance(molecule, list)
        assert len(molecule) == 1

        # obtain the last 10 structures
        molecules = xyz_file.get_molecules(index="-10:", return_list=True)
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

        # test first molecule is 1-indexed
        molecule = Molecule.from_filepath(
            filepath=multiple_molecules_xyz_file, index="1", return_list=False
        )
        assert np.allclose(
            molecule.positions[0],
            np.array([-1.0440166707, -2.3921211654, -1.1765767093]),
            rtol=1e-5,
        )
        assert isinstance(molecule, Molecule)

    def test_molecular_geometry(self):
        """Test molecular geometry calculations."""
        mol = Molecule(
            symbols=["C", "O", "O"],
            positions=np.array([[-1.16, 0, 0], [0, 0, 0], [1.16, 0, 0]]),
        )
        assert np.isclose(mol.mass, 44.01, atol=1e-2)
        assert np.isclose(mol.get_distance(1, 2), 1.16)
        assert np.isclose(mol.get_distance(2, 3), 1.16)
        assert np.isclose(mol.get_angle(1, 2, 3), 180)
        assert np.isclose(mol.get_dihedral(0, 1, 2, 0), 0)
        assert mol.is_linear


class TestMoleculeAdvanced:
    def test_molecule_to_rdkit_conversion(self):
        """Test conversion of Molecule to RDKit Mol object."""
        mol = Molecule(
            symbols=["C", "O", "O"],
            positions=np.array([[0, 0, 0], [1.2, 0, 0], [0, 1.2, 0]]),
            charge=-1,
            multiplicity=2,
        )
        rdkit_mol = mol.to_rdkit()

        assert isinstance(rdkit_mol, Chem.Mol)
        assert rdkit_mol.GetNumAtoms() == 3
        assert rdkit_mol.GetNumConformers() == 1
        assert rdkit_mol.GetConformer().GetPositions().shape == (3, 3)

    def test_molecule_graph_generation(self):
        """Test molecular graph creation with bond detection."""
        mol = Molecule(
            symbols=["C", "H", "H", "H", "H"],
            positions=np.array(
                [
                    [0, 0, 0],
                    [1.09, 0, 0],
                    [-1.09, 0, 0],
                    [0, 1.09, 0],
                    [0, -1.09, 0],
                ]
            ),
        )

        assert not mol.is_chiral, "CH4 is not chiral"
        graph = mol.to_graph()

        assert isinstance(graph, nx.Graph)
        assert len(graph.nodes) == 5
        assert len(graph.edges) == 4  # CH4 should have 4 bonds
        assert all(
            graph.nodes[i]["element"] == ("C" if i == 0 else "H")
            for i in range(5)
        )

    def test_charge_and_multiplicity_handling(self):
        """Test preservation of charge and multiplicity states."""
        charged_mol = Molecule(
            symbols=["O"],
            positions=np.array([[0, 0, 0]]),
            charge=-1,
            multiplicity=2,
        )

        assert charged_mol.charge == -1
        assert charged_mol.multiplicity == 2

    def test_invalid_molecule_creation(self):
        """Test error handling for invalid molecule configurations."""
        with pytest.raises(ValueError):
            # Mismatched symbols and positions
            Molecule(symbols=["C", "H"], positions=np.array([[0, 0, 0]]))

    def test_empty_molecule_handling(self):
        """Test edge case of empty molecule."""
        with pytest.raises(ValueError):
            Molecule(symbols=[], positions=np.empty((0, 3)))

    def test_frozen_atoms_manipulation(self):
        """Test frozen atoms property handling."""
        mol = Molecule(
            symbols=["C", "O"],
            positions=np.array([[0, 0, 0], [1.2, 0, 0]]),
            frozen_atoms=[-1, 0],
        )

        assert mol.frozen_atoms == [-1, 0]
        assert not mol.is_chiral

    def test_convert_ase_atoms_with_constraints_to_molecule(
        self, constrained_atoms
    ):
        """Test conversion of ASE Atoms with constraints to Molecule."""
        from chemsmart.io.molecules.atoms import AtomsChargeMultiplicity

        mol = AtomsChargeMultiplicity.from_atoms(
            constrained_atoms
        ).to_molecule(charge=0, multiplicity=1)

        assert isinstance(mol, Molecule)
        assert np.all(mol.symbols == ["Ar", "Ar"])
        assert mol.energy == 0.0
        assert np.allclose(
            mol.forces, np.array([(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)])
        )
        assert np.allclose(
            mol.velocities, np.array([(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)])
        )
        assert mol.frozen_atoms == [
            -1,
            0,
        ]  # -1: frozen atom, 0: relaxed atom (Gaussian format)
        assert mol.charge == 0
        assert mol.multiplicity == 1

        # no charge and multiplicity specification
        mol_no_charge_mult = AtomsChargeMultiplicity.from_atoms(
            constrained_atoms
        ).to_molecule()
        assert isinstance(mol_no_charge_mult, Molecule)
        assert np.all(mol_no_charge_mult.symbols == ["Ar", "Ar"])
        assert mol_no_charge_mult.energy == 0.0
        assert np.allclose(
            mol_no_charge_mult.forces,
            np.array([(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]),
        )
        assert np.allclose(
            mol_no_charge_mult.velocities,
            np.array([(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)]),
        )
        assert mol_no_charge_mult.frozen_atoms == [
            -1,
            0,
        ]  # Frozen atoms should be 1-indexed
        assert mol_no_charge_mult.charge is None
        assert mol_no_charge_mult.multiplicity is None

    def test_molecule_from_db_with_pbc_and_constraints(
        self, constrained_pbc_db_file
    ):
        """Test creation of Molecule from database with PBC and constraints."""
        mol = Molecule.from_filepath(
            constrained_pbc_db_file,
            index="-1",
        )

        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 96
        assert mol.chemical_formula == "Co25Cu9Fe9Mo44Ni9"
        assert np.all(mol.pbc_conditions == [True, True, True])
        assert len(mol.translation_vectors) == 3
        assert mol.frozen_atoms == [
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
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
            0,
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
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
            0,
        ]
        my_list = mol.frozen_atoms
        indices = [i for i, x in enumerate(my_list) if x == -1]
        assert indices == [
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            34,
            35,
            36,
            37,
            43,
            44,
            45,
            46,
            47,
            48,
            49,
            50,
            51,
            52,
            53,
            54,
            55,
            56,
            57,
            58,
            87,
            88,
            89,
            90,
        ]

        mol2 = Molecule.from_filepath(
            constrained_pbc_db_file,
            index="-5",
        )
        assert isinstance(mol2, Molecule)
        assert mol2.num_atoms == 96
        assert mol2.chemical_formula == "Co25Cu9Fe9Mo44Ni9"
        assert np.all(mol2.pbc_conditions == [True, True, True])
        assert len(mol2.translation_vectors) == 3
        assert mol2.frozen_atoms == [
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
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
            0,
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
            -1,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        my_list = mol2.frozen_atoms
        indices = [i for i, x in enumerate(my_list) if x == -1]
        assert indices == [
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            25,
            26,
            27,
            28,
            34,
            35,
            36,
            37,
            43,
            44,
            45,
            46,
            47,
            48,
            49,
            50,
            51,
            52,
            53,
            54,
            55,
            56,
            87,
            88,
            89,
        ]

    def test_pbc_handling(self):
        """Test periodic boundary conditions handling."""
        mol = Molecule(
            symbols=["C"],
            positions=np.array([[0, 0, 0]]),
            pbc_conditions=[1, 1, 0],
            translation_vectors=[[2.5, 0, 0], [0, 2.5, 0]],
        )

        assert mol.pbc_conditions == [1, 1, 0]
        assert len(mol.translation_vectors) == 2

    def test_distance_calculation(self):
        """Test bond distance calculations."""
        mol = Molecule(
            symbols=["H", "H"], positions=np.array([[0, 0, 0], [1.0, 0, 0]])
        )
        distances = mol.get_all_distances()

        assert len(distances) == 1
        assert np.isclose(distances[0], 1.0)


class TestCoordinateBlockAdvanced:
    def test_mixed_coordinate_formats(self):
        """Test parsing of mixed coordinate formats."""
        mixed_block = """
        C       1.0 2.0 3.0
        TV      4.0 5.0 6.0
        H       7.0 8.0 9.0
        """
        cb = CoordinateBlock(mixed_block)
        mol = cb.molecule

        assert len(mol.symbols) == 2  # Should ignore TV line
        assert mol.translation_vectors == [[4.0, 5.0, 6.0]]

    def test_cube_file_format(self):
        """Test parsing of cube file format coordinates."""
        cube_block = """
        6    6.000000  -12.064399   -0.057172   -0.099010
        1    1.000000    1.234500    2.345600    3.456700
        """
        cb = CoordinateBlock(cube_block)

        assert all(cb.symbols == ["C", "H"])
        assert np.allclose(cb.positions[0], [-12.064399, -0.057172, -0.099010])


class TestFileHandlingAdvanced:
    def test_large_file_handling(self, tmpdir):
        """Test handling of large molecule files."""
        large_xyz = "\n".join(
            [
                "1000",
                "Large Molecule",
                *[f"H {i} {i} {i}" for i in range(1000)],
            ]
        )

        tmp_file = tmpdir.join("large.xyz")
        with open(tmp_file, "w") as f:
            f.write(large_xyz)

        mol = Molecule.from_filepath(tmp_file)
        assert mol.num_atoms == 1000

    def test_corrupted_file_handling(self, tmpdir):
        """Test error handling for corrupted files."""
        corrupted_content = "Invalid file content"
        tmp_file = tmpdir.join("corrupted.xyz")
        with open(tmp_file, "w") as f:
            f.write(corrupted_content)

        with pytest.raises(ValueError):
            Molecule.from_filepath(tmp_file)


class TestGraphFeatures:
    def test_graph_properties(self):
        """Test molecular graph properties."""
        mol = Molecule(
            symbols=["C", "C"], positions=np.array([[0, 0, 0], [1.5, 0, 0]])
        )
        graph = mol.to_graph()

        assert nx.number_connected_components(graph) == 1
        assert nx.diameter(graph) == 1

    def test_variable_bond_cutoffs(self):
        """Test bond detection with different cutoff buffers."""
        mol = Molecule(
            symbols=["H", "H"], positions=np.array([[0, 0, 0], [0.9, 0, 0]])
        )

        assert not mol.is_chiral

        # H has covalent radius of 0.31 Å from ase.data

        # With buffer too small
        tight_graph = mol.to_graph(bond_cutoff_buffer=0.0, adjust_H=False)
        assert len(tight_graph.edges) == 0

        # With reasonable buffer
        normal_graph = mol.to_graph(bond_cutoff_buffer=0.3, adjust_H=False)
        assert len(normal_graph.edges) == 1


class TestChemicalFeatures:
    def test_stereochemistry_handling(self, methyl3hexane_molecule):
        """Test preservation of stereochemical information."""
        methyl_3_hexane = methyl3hexane_molecule
        assert np.all(
            methyl_3_hexane.bond_orders
            == [
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
            ]
        )
        assert len(methyl_3_hexane.bond_orders) == 22
        assert all([i == 1 for i in methyl_3_hexane.bond_orders])
        assert methyl_3_hexane.is_chiral
        chiral_mol = Molecule(
            symbols=["C", "Cl", "F", "Br", "I"],
            positions=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.2, 0.0, -0.5],
                    [-0.6, 1.0, 0.5],
                    [-0.6, -1.0, 0.5],
                    [0.0, 0.0, 1.3],
                ]
            ),
        )

        # can't get the bond orders of this challenging molecule correctly
        # bond_orders = chiral_mol.get_bond_orders_from_rdkit_mol(bond_cutoff_buffer=-0.4)
        # print(bond_orders)
        assert chiral_mol.is_chiral

        rdkit_mol = chiral_mol.to_rdkit()
        assert Chem.FindMolChiralCenters(rdkit_mol) != []

    @pytest.mark.skipif(
        not is_pubchem_network_available() or not is_pubchem_api_available(),
        reason="Network to pubchem is unavailable",
    )
    def test_more_stereochemistry_handling(self):
        """Test preservation of stereochemical information with PubChem."""
        chiral_mol2 = Molecule.from_pubchem(
            "CC(C)(Oc1ccc(Cl)cc1)C(=O)N[C@H]1C2CCCC1C[C@@H](C(=O)O)C2"
        )
        assert chiral_mol2.is_chiral
        rdkit_mol2 = chiral_mol2.to_rdkit()
        assert Chem.FindMolChiralCenters(rdkit_mol2) != []

    def test_resonance_handling(
        self,
        gaussian_ozone_opt_outfile,
        gaussian_acetone_opt_outfile,
        gaussian_benzene_opt_outfile,
    ):
        """Test handling of resonance structures."""
        ozone = Molecule.from_filepath(gaussian_ozone_opt_outfile)
        assert ozone.get_chemical_formula() == "O3"
        assert ozone.chemical_formula == "O3"
        assert ozone.bond_orders == [
            1.5,
            1.5,
        ]  # correctly gets bond order of ozone as 1.5
        assert not ozone.is_chiral
        rdkit_mol = ozone.to_rdkit()
        assert Chem.FindMolChiralCenters(rdkit_mol) == []
        assert ozone.chemical_symbols == ["O", "O", "O"]
        assert ozone.atomic_radii_list == [0.66, 0.66, 0.66]
        assert ozone.vdw_radii_list == [1.52, 1.52, 1.52]

        graph = ozone.to_graph()
        assert any(
            bond["bond_order"] > 1 for bond in graph.edges.values()
        )  # Check for possible multiple bonds

        acetone = Molecule.from_filepath(gaussian_acetone_opt_outfile)
        assert acetone.bond_orders == [
            2.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
        ]
        # 1 double bond C2=O1, 2 C-C single bonds, 6 C-H single bonds

        graph = acetone.to_graph()
        assert any(bond["bond_order"] > 1 for bond in graph.edges.values())

        benzene = Molecule.from_filepath(gaussian_benzene_opt_outfile)
        assert benzene.bond_orders == [
            1.5,
            1.5,
            1.0,
            1.5,
            1.0,
            1.5,
            1.0,
            1.5,
            1.0,
            1.5,
            1.0,
            1.0,
        ]
        # check there are 6 aromatic C-C bonds
        # and 6 single C-H bonds in benzene
        assert len([bond for bond in benzene.bond_orders if bond == 1.5]) == 6
        assert len([bond for bond in benzene.bond_orders if bond == 1.0]) == 6

    def test_volume(
        self, gaussian_ozone_opt_outfile, gaussian_acetone_opt_outfile
    ):
        """Test volume calculation for molecules.

        Tests various volume calculation methods:
        - voronoi_dirichlet_occupied_volume
        - crude_volume_by_vdw_radii
        - crude_volume_by_atomic_radii
        - vdw_volume
        - vdw_volume_from_rdkit
        - voronoi_dirichlet_polyhedra_occupied_volume
        """
        ozone = Molecule.from_filepath(gaussian_ozone_opt_outfile)

        ozone_vd_vol = ozone.voronoi_dirichlet_occupied_volume
        assert ozone_vd_vol > 0
        assert np.isclose(ozone_vd_vol, 44.13068085447146, rtol=0.01)

        # Test other volume methods
        assert np.isclose(
            ozone.crude_volume_by_vdw_radii, 44.13068085447146, rtol=0.01
        )
        assert np.isclose(
            ozone.crude_volume_by_atomic_radii, 3.612781286145805, rtol=0.01
        )
        # vdw_volume uses exact lens-shaped
        # intersection formula for overlap correction
        assert np.isclose(ozone.vdw_volume, 29.65124427436735, rtol=0.01)
        # grid_vdw_volume uses grid-based integration (similar to RDKit)
        assert np.isclose(ozone.grid_vdw_volume, 31.464, rtol=0.05)
        assert np.isclose(
            ozone.vdw_volume_from_rdkit, 25.533471711063285, rtol=0.01
        )
        assert np.isclose(
            ozone.voronoi_dirichlet_polyhedra_occupied_volume,
            20.94252748967074,
            rtol=0.01,
        )

        acetone = Molecule.from_filepath(gaussian_acetone_opt_outfile)

        acetone_vd_vol = acetone.voronoi_dirichlet_occupied_volume
        assert acetone_vd_vol > 0
        assert np.isclose(acetone_vd_vol, 73.5820640781846, rtol=0.01)

        # Test other volume methods
        assert np.isclose(
            acetone.crude_volume_by_vdw_radii, 119.87818262306239, rtol=0.01
        )
        assert np.isclose(
            acetone.crude_volume_by_atomic_radii, 7.469325029468949, rtol=0.01
        )
        assert np.isclose(
            acetone.voronoi_dirichlet_polyhedra_occupied_volume,
            12.369068467588548,
            rtol=0.01,
        )
        # vdw_volume uses exact lens-shaped
        # intersection formula for overlap correction
        assert np.isclose(acetone.vdw_volume, 48.85540325089168, rtol=0.01)
        # grid_vdw_volume uses grid-based integration (similar to RDKit)
        assert np.isclose(acetone.grid_vdw_volume, 64.832, rtol=0.05)
        assert np.isclose(
            acetone.vdw_volume_from_rdkit, 61.98249809788294, rtol=0.01
        )


class TestStructuresFromGaussianInput:
    def test_read_molecule_from_gaussian_opt_input(
        self, tmpdir, gaussian_opt_inputfile
    ):
        assert os.path.exists(gaussian_opt_inputfile)
        assert os.path.isfile(gaussian_opt_inputfile)

        g16_file = Gaussian16Input(filename=gaussian_opt_inputfile)
        molecule = g16_file.molecule
        tmp_file = os.path.join(tmpdir, "tmp.txt")
        f = open(tmp_file, "w")
        molecule.write_coordinates(f, program="gaussian")
        f.close()
        coordinates = """C       -0.5448210000   -1.1694570000    0.0001270000
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
Cl      -3.0556310000   -0.1578960000   -0.0001400000"""
        with open(tmp_file, "r") as g:
            written_coordinates = g.read()
        assert os.path.exists(tmp_file)
        assert all([a == b for a, b in zip(coordinates, written_coordinates)])
        assert g16_file.num_atoms == 14
        assert isinstance(g16_file.molecule, Molecule)
        assert g16_file.molecule.charge == 0
        assert g16_file.molecule.multiplicity == 1
        assert g16_file.molecule.empirical_formula == "C7H5ClO"
        assert g16_file.molecule.frozen_atoms is None
        assert g16_file.molecule.pbc_conditions is None
        assert g16_file.molecule.translation_vectors is None

    def test_read_molecule_from_gaussian_modred(
        self, gaussian_modred_inputfile
    ):
        assert os.path.exists(gaussian_modred_inputfile)
        assert os.path.isfile(gaussian_modred_inputfile)

        g16_file = Gaussian16Input(filename=gaussian_modred_inputfile)
        assert g16_file.num_atoms == 14
        assert isinstance(g16_file.molecule, Molecule)
        assert g16_file.molecule.charge == 0
        assert g16_file.molecule.multiplicity == 1
        assert g16_file.molecule.empirical_formula == "C3H7NO3"
        assert g16_file.molecule.frozen_atoms is None
        assert g16_file.molecule.pbc_conditions is None
        assert g16_file.molecule.translation_vectors is None

    def test_read_molecule_from_gaussian_scan(self, gaussian_scan_inputfile):
        assert os.path.exists(gaussian_scan_inputfile)
        assert os.path.isfile(gaussian_scan_inputfile)

        g16_file = Gaussian16Input(filename=gaussian_scan_inputfile)
        assert g16_file.num_atoms == 14
        assert isinstance(g16_file.molecule, Molecule)
        assert g16_file.molecule.charge == 0
        assert g16_file.molecule.multiplicity == 1
        assert g16_file.molecule.empirical_formula == "C3H7NO3"
        assert g16_file.molecule.frozen_atoms is None
        assert g16_file.molecule.pbc_conditions is None
        assert g16_file.molecule.translation_vectors is None

    def test_read_molecule_from_opt_genecp(
        self, gaussian_opt_genecp_inputfile
    ):
        assert os.path.exists(gaussian_opt_genecp_inputfile)
        assert os.path.isfile(gaussian_opt_genecp_inputfile)

        g16_file = Gaussian16Input(filename=gaussian_opt_genecp_inputfile)
        assert g16_file.num_atoms == 15
        assert isinstance(g16_file.molecule, Molecule)
        assert g16_file.molecule.charge == 0
        assert g16_file.molecule.multiplicity == 1
        assert g16_file.molecule.empirical_formula == "C4H6O4Pd"
        assert g16_file.molecule.frozen_atoms is None
        assert g16_file.molecule.pbc_conditions is None
        assert g16_file.molecule.translation_vectors is None

    def test_read_molecule_from_modred_gen(self, modred_gen_inputfile):
        assert os.path.exists(modred_gen_inputfile)
        assert os.path.isfile(modred_gen_inputfile)

        g16_file = Gaussian16Input(filename=modred_gen_inputfile)
        assert g16_file.num_atoms == 15
        assert isinstance(g16_file.molecule, Molecule)
        assert g16_file.molecule.charge == 0
        assert g16_file.molecule.multiplicity == 1
        assert g16_file.molecule.empirical_formula == "C4H6BrO4"
        assert g16_file.molecule.frozen_atoms is None
        assert g16_file.molecule.pbc_conditions is None
        assert g16_file.molecule.translation_vectors is None

    def test_read_molecule_from_modred_genecp(self, modred_genecp_inputfile):
        assert os.path.exists(modred_genecp_inputfile)
        assert os.path.isfile(modred_genecp_inputfile)

        g16_file = Gaussian16Input(filename=modred_genecp_inputfile)
        assert g16_file.num_atoms == 15
        assert isinstance(g16_file.molecule, Molecule)
        assert g16_file.molecule.charge == 0
        assert g16_file.molecule.multiplicity == 1
        assert g16_file.molecule.empirical_formula == "C4H6O4Pd"
        assert g16_file.molecule.frozen_atoms is None
        assert g16_file.molecule.pbc_conditions is None
        assert g16_file.molecule.translation_vectors is None

    def test_read_molecule_from_modred_genecp_custom_solvent(
        self, modred_genecp_custom_solvent_inputfile
    ):
        assert os.path.exists(modred_genecp_custom_solvent_inputfile)
        assert os.path.isfile(modred_genecp_custom_solvent_inputfile)

        g16_file = Gaussian16Input(
            filename=modred_genecp_custom_solvent_inputfile
        )
        assert g16_file.num_atoms == 65
        assert isinstance(g16_file.molecule, Molecule)
        assert g16_file.molecule.charge == 0
        assert g16_file.molecule.multiplicity == 1
        assert g16_file.molecule.empirical_formula == "C28H26Cl2N2O4PdS2"
        assert g16_file.molecule.frozen_atoms is None
        assert g16_file.molecule.pbc_conditions is None
        assert g16_file.molecule.translation_vectors is None

    def test_read_molecule_from_gaussian_frozen_opt(
        self, tmpdir, gaussian_frozen_opt_inputfile
    ):
        assert os.path.exists(gaussian_frozen_opt_inputfile)
        assert os.path.isfile(gaussian_frozen_opt_inputfile)

        g16_file = Gaussian16Input(filename=gaussian_frozen_opt_inputfile)
        assert g16_file.num_atoms == 14
        assert isinstance(g16_file.molecule, Molecule)
        assert g16_file.molecule.charge == 0
        assert g16_file.molecule.multiplicity == 1
        assert g16_file.molecule.empirical_formula == "C7H5ClO"
        assert g16_file.molecule.frozen_atoms == [
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
        assert g16_file.molecule.pbc_conditions is None
        assert g16_file.molecule.translation_vectors is None
        molecule = g16_file.molecule
        tmp_file = os.path.join(tmpdir, "tmp.txt")
        f = open(tmp_file, "w")
        molecule.write_coordinates(f, program="gaussian")
        f.close()
        coordinates = """C         -1   -0.5448210000   -1.1694570000    0.0001270000
C         -1    0.8378780000   -1.0476350000    0.0001900000
C         -1    1.4329940000    0.2194290000    0.0001440000
C         -1    0.6358350000    1.3657650000   -0.0000040000
C         -1   -0.7521390000    1.2582750000    0.0000710000
C         -1   -1.3283680000   -0.0113420000    0.0001560000
H         -1   -1.0298620000   -2.1454490000    0.0001020000
H         -1    1.4853190000   -1.9262430000    0.0002620000
H         -1    1.1050940000    2.3527070000    0.0000530000
H         -1   -1.3913950000    2.1406100000   -0.0000130000
C          0    2.9142600000    0.3363820000    0.0000140000
O          0    3.6625230000   -0.6037690000   -0.0002940000
H          0    3.3025560000    1.3842410000    0.0001510000
Cl         0   -3.0556310000   -0.1578960000   -0.0001400000"""
        with open(tmp_file, "r") as g:
            written_coordinates = g.read()
        assert os.path.exists(tmp_file)
        assert all([a == b for a, b in zip(coordinates, written_coordinates)])

    def test_read_molecule_from_gaussian_pbc(
        self, tmpdir, gaussian_pbc_1d_inputfile
    ):
        assert os.path.exists(gaussian_pbc_1d_inputfile)
        assert os.path.isfile(gaussian_pbc_1d_inputfile)

        g16_file = Gaussian16Input(filename=gaussian_pbc_1d_inputfile)
        assert g16_file.num_atoms == 10
        assert isinstance(g16_file.molecule, Molecule)
        assert g16_file.molecule.charge == 0
        assert g16_file.molecule.multiplicity == 1
        assert g16_file.molecule.empirical_formula == "C4H5Cl"
        assert g16_file.molecule.frozen_atoms is None
        assert g16_file.molecule.pbc_conditions == [True, False, False]
        assert g16_file.molecule.translation_vectors == [
            [4.8477468928, 0.1714181332, 0.5112729831],
        ]
        molecule = g16_file.molecule
        tmp_file = os.path.join(tmpdir, "tmp.txt")
        f = open(tmp_file, "w")
        molecule.write_coordinates(f, program="gaussian")
        f.close()
        coordinates = """C       -1.9267226529    0.4060180273    0.0316702826
H       -2.3523143977    0.9206168644    0.9131400756
H       -1.8372739404    1.1548899113   -0.7707507970
C       -0.5737182157   -0.1434584477    0.3762843235
H       -0.5015912465   -0.7653394047    1.2791284293
C        0.5790889876    0.0220081655   -0.3005160849
C        1.9237098673   -0.5258773194    0.0966261209
H        1.7722344520   -1.2511397907    0.9159625120
H        2.3627869487   -1.0792380182   -0.7525115830
Cl       0.6209825739    0.9860944599   -1.7876398696
TV       4.8477468928    0.1714181332    0.5112729831"""
        with open(tmp_file, "r") as g:
            written_coordinates = g.read()
        assert os.path.exists(tmp_file)
        assert all([a == b for a, b in zip(coordinates, written_coordinates)])


class TestQMMMinMolecule:
    def test_atoms_in_levels_wrong_low_level(
        self, tmpdir, methyl3hexane_molecule
    ):
        methyl_3_hexane = QMMMMolecule(molecule=methyl3hexane_molecule)
        methyl_3_hexane.high_level_atoms = [1, 2, 3]
        methyl_3_hexane.medium_level_atoms = [4, 5, 6]
        methyl_3_hexane.low_level_atoms = [7, 8, 9]
        methyl_3_hexane.bonded_atoms = [(3, 4), (1, 7)]
        assert methyl_3_hexane.chemical_formula == "C7H16"
        assert methyl_3_hexane.num_atoms == 23
        with pytest.raises(ValueError):
            methyl_3_hexane.partition_level_strings
            # should raise error since high + medium +
            # low is not equal to total number of atoms

    def test_atoms_in_levels_default_low_level(
        self,
        tmpdir,
        qmmm_written_xyz_file,
        qmmm_written_xyz_only_file,
        methyl3hexane_molecule,
    ):
        methyl_3_hexane = QMMMMolecule(molecule=methyl3hexane_molecule)
        methyl_3_hexane.high_level_atoms = [1, 2, 3]
        methyl_3_hexane.medium_level_atoms = [4, 5, 6]
        methyl_3_hexane.bonded_atoms = [(3, 4), (1, 7)]
        assert methyl_3_hexane.chemical_formula == "C7H16"
        assert methyl_3_hexane.num_atoms == 23

        assert methyl_3_hexane.partition_level_strings == [
            "H",
            "H",
            "H",
            "M",
            "M",
            "M",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
            "L",
        ]

        written_input = os.path.join(tmpdir, "tmp.xyz")
        methyl_3_hexane.write(written_input, format="xyz", xyz_only=False)
        assert cmp_with_ignore(
            written_input, qmmm_written_xyz_only_file, ignore_string="tmp"
        )  # writes input file as expected

        written_input2 = os.path.join(tmpdir, "tmp_xyz_only.xyz")
        methyl_3_hexane.write(written_input2, format="xyz", xyz_only=True)
        assert cmp(
            written_input2, qmmm_written_xyz_only_file, shallow=False
        )  # writes input file as expected

    def test_qmmm_atoms_handling(self, tmpdir):
        """Test QM/MM atoms handling."""
        mol = QMMMMolecule(
            symbols=["O", "H", "H", "Cl"],
            positions=np.array(
                [
                    [-4.84098481, -0.56828899, 0.00000000],
                    [-3.88098484, -0.56804789, 0.00000000],
                    [-5.16121212, 0.33672729, 0.00000000],
                    [-1.93181817, -0.59090908, 0.00000000],
                ]
            ),
            high_level_atoms=[4],
            medium_level_atoms=[3],
            low_level_atoms=[1, 2],
            bonded_atoms=[(1, 3)],
            scale_factors={(1, 3): [0.9, 0.8, 0.7]},
        )

        assert mol.high_level_atoms == [4]
        assert mol.low_level_atoms == [1, 2]
        assert mol.medium_level_atoms == [3]
        assert mol.bonded_atoms == [(1, 3)]

        written_input = os.path.join(tmpdir, "tmp.xyz")
        with open(written_input, "w") as f:
            mol._write_gaussian_coordinates(f)
        with open(written_input, "r") as f:
            lines = [line.strip() for line in f.readlines()]

            expected_lines = [
                "O -4.8409848100 -0.5682889900 0.0000000000 L H 3  0.9 0.8 0.7",
                "H -3.8809848400 -0.5680478900 0.0000000000 L",
                "H -5.1612121200 0.3367272900 0.0000000000 M",
                "Cl -1.9318181700 -0.5909090800 0.0000000000 H",
            ]

            assert [" ".join(line.split()) for line in lines] == [
                " ".join(line.split()) for line in expected_lines
            ], f"Mismatch in written Gaussian coordinates:\nExpected: {expected_lines}\nGot: {lines}"
        if os.path.exists("tmp.xyz"):
            os.remove("tmp.xyz")


class TestSDFFile:
    def test_converts_sdf_string_to_molecule_object(self, tmpdir):
        sdf_string = """6999790
  -OEChem-03092302273D

 13 12  0     1  0  0  0  0  0999 V2000
   -1.1018    1.2385   -0.0496 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4517    0.0248    0.3184 C   0  0  1  0  0  0  0  0  0  0  0  0
   -1.3489   -1.1518   -0.0426 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8543   -0.0422   -0.4148 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0481   -0.0692    0.1885 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2932    0.0475    1.4026 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3037   -1.0804    0.4900 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5887   -1.1551   -1.1122 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8804   -2.1073    0.2141 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.8229   -0.0582   -1.5016 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9459    1.2770    0.4319 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.9591   -0.1121   -0.3983 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.1387   -0.0496    1.2689 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1 11  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2  6  1  0  0  0  0
  3  7  1  0  0  0  0
  3  8  1  0  0  0  0
  3  9  1  0  0  0  0
  4  5  2  0  0  0  0
  4 10  1  0  0  0  0
  5 12  1  0  0  0  0
  5 13  1  0  0  0  0
M  END
> <PUBCHEM_COMPOUND_CID>
6999790

> <PUBCHEM_CONFORMER_RMSD>
0.4

> <PUBCHEM_CONFORMER_DIVERSEORDER>
1
2
3

> <PUBCHEM_MMFF94_PARTIAL_CHARGES>
8
1 -0.68
10 0.15
11 0.4
12 0.15
13 0.15
2 0.42
4 -0.29
5 -0.3

> <PUBCHEM_EFFECTIVE_ROTOR_COUNT>
1

> <PUBCHEM_PHARMACOPHORE_FEATURES>
3
1 1 acceptor
1 1 donor
1 5 hydrophobe

> <PUBCHEM_HEAVY_ATOM_COUNT>
5

> <PUBCHEM_ATOM_DEF_STEREO_COUNT>
1

> <PUBCHEM_ATOM_UDEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_DEF_STEREO_COUNT>
0

> <PUBCHEM_BOND_UDEF_STEREO_COUNT>
0

> <PUBCHEM_ISOTOPIC_ATOM_COUNT>
0

> <PUBCHEM_COMPONENT_COUNT>
1

> <PUBCHEM_CACTVS_TAUTO_COUNT>
1

> <PUBCHEM_CONFORMER_ID>
006ACEEE00000001

> <PUBCHEM_MMFF94_ENERGY>
2.2534

> <PUBCHEM_FEATURE_SELFOVERLAP>
15.223

> <PUBCHEM_SHAPE_FINGERPRINT>
139733 1 9078577887532652523
16714656 1 18411707598977923407
20096714 4 18191872223362461888
21015797 1 9943532881608965411
21040471 1 18194406580373827108
29004967 10 18120096344745644827
5460574 1 9511472116931879378
5943 1 10672208095369030331

> <PUBCHEM_SHAPE_MULTIPOLES>
97.03
2.26
1.1
0.66
0.88
0.14
0
-0.36
0.09
-0.58
-0.02
0.05
-0.02
-0.01

> <PUBCHEM_SHAPE_SELFOVERLAP>
167.629

> <PUBCHEM_SHAPE_VOLUME>
65.3

> <PUBCHEM_COORDINATE_TYPE>
2
5
10

$$$$"""
        from chemsmart.io.file import SDFFile

        tmpfile = os.path.join(tmpdir, "test.sdf")
        with open(tmpfile, "w") as f:
            f.write(sdf_string)
        sdf_file = SDFFile(filename=tmpfile)
        sdf_molecule = sdf_file.molecule
        assert isinstance(sdf_molecule, Molecule)
        symbols = [
            "O",
            "C",
            "C",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
        ]
        assert all(
            [
                sdf_molecule.symbols[i] == symbols[i]
                for i in range(len(symbols))
            ]
        )

        assert sdf_molecule.empirical_formula == "C4H8O"
        structure_coords = np.array(
            [
                [-1.1018, 1.2385, -0.0496],
                [-0.4517, 0.0248, 0.3184],
                [-1.3489, -1.1518, -0.0426],
                [0.8543, -0.0422, -0.4148],
                [2.0481, -0.0692, 0.1885],
                [-0.2932, 0.0475, 1.4026],
                [-2.3037, -1.0804, 0.49],
                [-1.5887, -1.1551, -1.1122],
                [-0.8804, -2.1073, 0.2141],
                [0.8229, -0.0582, -1.5016],
                [-1.9459, 1.277, 0.4319],
                [2.9591, -0.1121, -0.3983],
                [2.1387, -0.0496, 1.2689],
            ]
        )

        assert np.all(sdf_molecule.positions == structure_coords)


class TestCDXFile:
    """Tests for ChemDraw file reading functionality."""

    def test_read_single_molecule_cdxml_file_benzene(
        self, single_molecule_cdxml_file_benzene
    ):
        """Test reading a single molecule from a CDXML file."""
        assert os.path.exists(single_molecule_cdxml_file_benzene)
        assert os.path.isfile(single_molecule_cdxml_file_benzene)

        cdx_file = CDXFile(filename=single_molecule_cdxml_file_benzene)
        molecules = cdx_file.molecules

        assert isinstance(molecules, list)
        assert len(molecules) == 1
        mol = molecules[0]
        assert isinstance(mol, Molecule)
        assert mol.chemical_formula == "C6H6"
        assert mol.num_atoms == 12  # benzene with hydrogens
        assert mol.is_aromatic

    def test_read_single_molecule_cdxml_file_methane(
        self, single_molecule_cdxml_file_methane
    ):
        """Test reading a single molecule from a CDXML file."""
        assert os.path.exists(single_molecule_cdxml_file_methane)
        assert os.path.isfile(single_molecule_cdxml_file_methane)

        cdx_file = CDXFile(filename=single_molecule_cdxml_file_methane)
        molecules = cdx_file.molecules

        assert isinstance(molecules, list)
        assert len(molecules) == 1
        mol = molecules[0]
        assert isinstance(mol, Molecule)
        assert mol.chemical_formula == "CH4"
        assert mol.num_atoms == 5  # benzene with hydrogens
        assert not mol.is_aromatic

    def test_read_multi_molecule_cdxml_file(self, multi_molecule_cdxml_file):
        """Test reading multiple molecules from a CDXML file."""
        assert os.path.exists(multi_molecule_cdxml_file)
        assert os.path.isfile(multi_molecule_cdxml_file)

        cdx_file = CDXFile(filename=multi_molecule_cdxml_file)
        molecules = cdx_file.molecules

        assert isinstance(molecules, list)
        assert len(molecules) == 2
        assert molecules[0].chemical_formula == "CH2O"  # formaldehyde
        assert molecules[1].chemical_formula == "N2"  # nitrogen

    def test_molecule_from_filepath_cdxml(
        self, single_molecule_cdxml_file_benzene
    ):
        """Test Molecule.from_filepath with CDXML file."""
        mol = Molecule.from_filepath(single_molecule_cdxml_file_benzene)

        assert isinstance(mol, Molecule)
        assert mol.chemical_formula == "C6H6"
        assert mol.num_atoms == 12
        assert mol.is_aromatic
        assert mol.positions is not None
        assert mol.positions.shape == (12, 3)

    def test_molecule_from_filepath_cdxml_pathlib(
        self, single_molecule_cdxml_file_benzene
    ):
        """Test Molecule.from_filepath with pathlib.Path."""
        path = Path(single_molecule_cdxml_file_benzene)
        mol = Molecule.from_filepath(path)

        assert isinstance(mol, Molecule)
        assert mol.chemical_formula == "C6H6"

    def test_molecule_from_filepath_cdxml_multi_molecules(
        self, multi_molecule_cdxml_file
    ):
        """Test reading multiple molecules from CDXML using from_filepath."""
        # Get all molecules
        molecules = Molecule.from_filepath(
            multi_molecule_cdxml_file, index=":", return_list=True
        )
        assert isinstance(molecules, list)
        assert len(molecules) == 2

        # Get first molecule (1-based indexing)
        mol1 = Molecule.from_filepath(multi_molecule_cdxml_file, index="1")
        assert mol1.chemical_formula == "CH2O"

        # Get last molecule
        mol_last = Molecule.from_filepath(
            multi_molecule_cdxml_file, index="-1"
        )
        assert mol_last.chemical_formula == "N2"

    def test_molecule_from_filepath_cdxml_return_list(
        self, single_molecule_cdxml_file_benzene
    ):
        """Test return_list parameter with single molecule CDXML file."""
        mol_single = Molecule.from_filepath(
            single_molecule_cdxml_file_benzene, return_list=False
        )
        assert isinstance(mol_single, Molecule)

        mol_list = Molecule.from_filepath(
            single_molecule_cdxml_file_benzene, return_list=True
        )
        assert isinstance(mol_list, list)
        assert len(mol_list) == 1
        assert isinstance(mol_list[0], Molecule)

    def test_cdxfile_get_molecules_index(self, multi_molecule_cdxml_file):
        """Test CDXFile.get_molecules with various index specifications."""
        cdx_file = CDXFile(filename=multi_molecule_cdxml_file)

        # Test index=":"
        all_mols = cdx_file.get_molecules(index=":")
        assert len(all_mols) == 2

        # Test index="-1"
        last_mol = cdx_file.get_molecules(index="-1")
        assert last_mol.chemical_formula == "N2"

        # Test 1-based indexing
        first_mol = cdx_file.get_molecules(index="1")
        assert first_mol.chemical_formula == "CH2O"

        second_mol = cdx_file.get_molecules(index="2")
        assert second_mol.chemical_formula == "N2"

    def test_cdx_molecule_to_rdkit_conversion(
        self, single_molecule_cdxml_file_benzene
    ):
        """Test that molecules from CDXML can be converted to RDKit."""
        mol = Molecule.from_filepath(single_molecule_cdxml_file_benzene)
        rdkit_mol = mol.to_rdkit()

        assert isinstance(rdkit_mol, Chem.Mol)
        assert rdkit_mol.GetNumAtoms() == 12
        assert rdkit_mol.GetNumConformers() == 1

    def test_cdx_molecule_to_graph(self, single_molecule_cdxml_file_benzene):
        """Test that molecules from CDXML can be converted to graph."""
        mol = Molecule.from_filepath(single_molecule_cdxml_file_benzene)
        graph = mol.to_graph()

        assert isinstance(graph, nx.Graph)
        assert len(graph.nodes) == 12  # benzene with hydrogens
        # Benzene should have 6 C-C bonds + 6 C-H bonds = 12 bonds
        assert len(graph.edges) == 12

    def test_read_single_molecule_cdx_file_imidazole(
        self, single_molecule_cdx_file_imidazole
    ):
        """Test reading a single molecule from a CDXML file."""
        assert os.path.exists(single_molecule_cdx_file_imidazole)
        assert os.path.isfile(single_molecule_cdx_file_imidazole)

        cdx_file = CDXFile(filename=single_molecule_cdx_file_imidazole)
        molecules = cdx_file.molecules

        assert isinstance(molecules, list)
        assert len(molecules) == 1
        mol = molecules[0]
        assert isinstance(mol, Molecule)
        assert mol.chemical_formula == "C8H10N2O"
        assert mol.num_atoms == 21
        assert mol.is_aromatic

    def test_read_complex_molecule_cdxml_file_(
        self, complex_molecule_cdxml_file
    ):
        """Test reading a single molecule from a CDXML file."""
        assert os.path.exists(complex_molecule_cdxml_file)
        assert os.path.isfile(complex_molecule_cdxml_file)

        cdx_file = CDXFile(filename=complex_molecule_cdxml_file)
        molecules = cdx_file.molecules

        assert isinstance(molecules, list)
        assert len(molecules) == 1
        mol = molecules[0]
        assert isinstance(mol, Molecule)
        assert mol.chemical_formula == "C32H31N5O5"
        assert mol.num_atoms == 73  # benzene with hydrogens
        assert mol.is_aromatic

    def test_read_metal_ligand_molecules_cdxml_file_(
        self, metal_ligand_molecules_cdxml_file
    ):
        """Test reading multiple organometallic molecules from a CDXML file with Cp and aromatic ligands."""
        assert os.path.exists(metal_ligand_molecules_cdxml_file)
        assert os.path.isfile(metal_ligand_molecules_cdxml_file)

        cdx_file = CDXFile(filename=metal_ligand_molecules_cdxml_file)
        molecules = cdx_file.molecules

        assert isinstance(molecules, list)
        assert len(molecules) == 7

        # Test molecule 0: Ti complex with Cp ligands
        mol = molecules[0]
        assert isinstance(mol, Molecule)
        assert mol.chemical_formula == "C14H30Ti"
        assert mol.num_atoms == 45
        assert mol.is_aromatic

        # Test molecule 5: Ir complex with aromatic benzene ligands
        mol = molecules[5]
        assert isinstance(mol, Molecule)
        assert mol.chemical_formula == "C35H31Cl2FeNO3P2"
        assert mol.num_atoms == 75
        assert mol.is_aromatic


def test_qmmm_partition_overlap_raises():
    """Creating a QMMMMolecule with overlapping
    partitions should raise a ValueError."""
    # Create a small dummy molecule
    symbols = ["C"] * 5
    positions = np.zeros((5, 3))
    m = Molecule(symbols=symbols, positions=positions)
    # High and medium overlap (atom index 2 appears in both)
    q = QMMMMolecule(
        molecule=m,
        high_level_atoms=[1, 2],
        medium_level_atoms=[2, 3],
        low_level_atoms=None,
    )
    with pytest.raises(ValueError) as exc:
        q._get_partition_levels()
    assert "Overlap" in str(exc.value)


def test_qmmm_partition_out_of_range_raises():
    """Specifying out-of-range atom indices should raise a ValueError."""
    symbols = ["C"] * 4
    positions = np.zeros((4, 3))
    m = Molecule(symbols=symbols, positions=positions)
    # index 10 out of range
    q = QMMMMolecule(
        molecule=m,
        high_level_atoms=[1],
        medium_level_atoms=[2],
        low_level_atoms=[10],
    )
    with pytest.raises(ValueError) as exc:
        q._get_partition_levels()
    assert "out of range" in str(exc.value)


class TestInChIKey:
    """Tests for Molecule.inchikey property (Open Babel backend)."""

    EXPECTED_NORMAL = "NNJYFTBCZFRDIO-UHFFFAOYSA-N"
    EXPECTED_R_ENANTIOMER = "YDCAVENCOFCEDV-HSZRJFAPSA-N"
    EXPECTED_S_ENANTIOMER = "YDCAVENCOFCEDV-QHCPKHFHSA-N"
    EXPECTED_LARGE_C3 = "WYLDIUSELJCHHK-MMELAICESA-M"
    EXPECTED_LARGE_C2 = "KRPJGRYSEYYRSW-YWQHEUOTSA-M"

    @staticmethod
    def _load_molecule(filepath):
        mol = Molecule.from_filepath(filepath)
        if isinstance(mol, list):
            mol = mol[-1]
        return mol

    def test_regression_inchikey(self, inchikey_normal_file):
        """InChIKey for a simple small molecule should be
        deterministic across repeated calls."""
        mol = self._load_molecule(inchikey_normal_file)
        for _ in range(3):
            assert mol.inchikey == self.EXPECTED_NORMAL

    def test_r_enantiomer_inchikey(self, inchikey_r_enantiomer_file):
        """InChIKey for the R-enantiomer should match the expected value."""
        mol = self._load_molecule(inchikey_r_enantiomer_file)
        assert mol.inchikey == self.EXPECTED_R_ENANTIOMER

    def test_s_enantiomer_inchikey(self, inchikey_s_enantiomer_file):
        """InChIKey for the S-enantiomer should match the expected value."""
        mol = self._load_molecule(inchikey_s_enantiomer_file)
        assert mol.inchikey == self.EXPECTED_S_ENANTIOMER

    def test_enantiomers_share_connectivity_layer(
        self, inchikey_r_enantiomer_file, inchikey_s_enantiomer_file
    ):
        """R and S enantiomers share the same first (connectivity) layer of
        the InChIKey (identical constitution) but differ in the stereo layer,
        confirming that Open Babel correctly resolves the axial chirality."""
        mol_r = self._load_molecule(inchikey_r_enantiomer_file)
        mol_s = self._load_molecule(inchikey_s_enantiomer_file)
        # First 14-character block: same connectivity
        assert mol_r.inchikey.split("-")[0] == mol_s.inchikey.split("-")[0]
        # Second block: stereo layer must differ for a chiral pair
        assert mol_r.inchikey.split("-")[1] != mol_s.inchikey.split("-")[1]
        # Overall InChIKeys are distinct
        assert mol_r.inchikey != mol_s.inchikey

    def test_large_molecule_c3_inchikey(self, inchikey_large_molecule_c3_file):
        """InChIKey for a large molecule (c3) should match the expected value."""
        mol = self._load_molecule(inchikey_large_molecule_c3_file)
        assert mol.inchikey == self.EXPECTED_LARGE_C3

    def test_large_molecule_c2_inchikey(self, inchikey_large_molecule_c2_file):
        """InChIKey for a large molecule (c2) should match the expected value."""
        mol = self._load_molecule(inchikey_large_molecule_c2_file)
        assert mol.inchikey == self.EXPECTED_LARGE_C2

    def test_large_molecules_differ(
        self, inchikey_large_molecule_c3_file, inchikey_large_molecule_c2_file
    ):
        """Two different large molecules should produce different InChIKeys."""
        mol_c3 = self._load_molecule(inchikey_large_molecule_c3_file)
        mol_c2 = self._load_molecule(inchikey_large_molecule_c2_file)
        assert mol_c3.inchikey != mol_c2.inchikey


class TestCXSMILES:
    """Tests for Molecule.cxsmiles property (RDKit backend)."""

    EXPECTED_NORMAL = (
        "[H]C([H])([H])C([H])([H])op(=O)oC([H])([H])C([H])([H])[H] "
        "|(3.6969,1.9448,0.2049;3.0842,1.2373,-0.3608;3.7313,0.4472,"
        "-0.7559;2.652,1.7556,-1.2232;1.9939,0.6488,0.5088;1.3583,"
        "1.4601,0.8855;2.4534,0.1313,1.3605;1.2255,-0.2622,-0.2653;"
        "0.0003,-0.9988,0.5013;0.0051,-2.2423,-0.2655;-1.2252,-0.2619,"
        "-0.2663;-1.9972,0.645,0.509;-2.4576,0.1238,1.3579;-1.364,"
        "1.4563,0.8899;-3.0867,1.2342,-0.3613;-3.702,1.9386,0.2054;"
        "-2.6535,1.756,-1.221;-3.7312,0.4439,-0.7604)|"
    )
    EXPECTED_R_ENANTIOMER = (
        "[H]c1c([H])c([H])c(P(=O)(c2c([H])c([H])c([H])c([H])c2[H])"
        "C([H])([H])[C@@]2(C([H])([H])[H])c(=O)n(C([H])([H])[H])"
        "c3c([H])c([H])c([H])c([H])c32)c([H])c1[H] "
        "|(0.657587,4.92454,-1.29892;0.276843,3.96138,-0.951759;"
        "0.485742,2.81258,-1.71847;1.03401,2.87447,-2.66093;"
        "0.006912,1.57925,-1.27484;0.204594,0.685275,-1.87284;"
        "-0.6776,1.49046,-0.055628;-1.25801,-0.0518,0.702519;"
        "-1.3091,0.044461,2.20006;-2.88635,-0.380829,-0.012479;"
        "-3.12941,-0.354113,-1.39219;-2.32549,-0.092551,-2.08777;"
        "-4.40354,-0.645056,-1.87792;-4.59613,-0.624486,-2.95272;"
        "-5.43533,-0.955747,-0.986524;-6.43361,-1.18057,-1.36882;"
        "-5.19575,-0.974224,0.388852;-6.00568,-1.21169,1.08208;"
        "-3.92009,-0.687184,0.879324;-3.70148,-0.688059,1.9505;"
        "-0.227699,-1.40463,0.052757;-0.119079,-1.24802,-1.03485;"
        "-0.807656,-2.33494,0.17008;1.15351,-1.5942,0.72979;"
        "0.99867,-2.27426,2.08799;0.308383,-1.69441,2.71697;"
        "1.97119,-2.35931,2.59369;0.59233,-3.28598,1.93715;"
        "1.93427,-2.48285,-0.25179;1.70453,-3.64544,-0.504979;"
        "2.90906,-1.70492,-0.846562;3.81634,-2.19747,-1.84937;"
        "3.73917,-1.61247,-2.77951;3.54494,-3.24141,-2.05378;"
        "4.8593,-2.15693,-1.49649;2.95902,-0.433075,-0.262488;"
        "3.8266,0.615031,-0.550651;4.58376,0.527564,-1.3328;"
        "3.69973,1.78628,0.210751;4.36465,2.62845,0.005447;"
        "2.75405,1.88809,1.23239;2.68234,2.80776,1.81661;"
        "1.88531,0.820237,1.50649;1.12177,0.893949,2.28612;"
        "1.97673,-0.328846,0.734248;-0.882212,2.64496,0.711945;"
        "-1.38448,2.55034,1.67858;-0.411788,3.878,0.260678;"
        "-0.570718,4.77518,0.863087),wU:23.24|"
    )
    EXPECTED_S_ENANTIOMER = (
        "[H]c1c([H])c([H])c(P(=O)(c2c([H])c([H])c([H])c([H])c2[H])"
        "C([H])([H])[C@]2(C([H])([H])[H])c(=O)n(C([H])([H])[H])"
        "c3c([H])c([H])c([H])c([H])c32)c([H])c1[H] "
        "|(-4.10545,4.56874,0.709448;-3.52977,3.64265,0.645345;"
        "-3.31587,2.87242,1.79018;-3.72201,3.19528,2.7513;"
        "-2.57592,1.69184,1.70894;-2.37818,1.07929,2.59263;"
        "-2.05266,1.27517,0.478011;-1.13514,-0.28848,0.507738;"
        "-0.75977,-0.677078,1.90906;-2.21055,-1.523,-0.263839;"
        "-2.718,-1.39784,-1.56435;-2.48778,-0.513904,-2.16641;"
        "-3.52676,-2.40254,-2.09389;-3.92258,-2.30619,-3.10717;"
        "-3.83217,-3.53062,-1.32596;-4.46746,-4.31554,-1.7426;"
        "-3.3289,-3.6555,-0.02996;-3.56957,-4.53716,0.56812;"
        "-2.51594,-2.65283,0.503124;-2.10519,-2.72167,1.51378;"
        "0.246332,-0.126442,-0.664672;0.507903,-1.14931,-0.980531;"
        "-0.124177,0.404384,-1.55997;1.49769,0.604701,-0.148795;"
        "1.1654,1.86783,0.660908;0.549974,2.55301,0.05834;"
        "2.09382,2.38433,0.946933;0.621442,1.59517,1.57793;"
        "2.29809,1.06731,-1.37867;1.86402,1.731,-2.29803;"
        "3.58921,0.612819,-1.24748;4.63005,0.881221,-2.20479;"
        "5.4733,1.41251,-1.73622;4.19794,1.50997,-2.994;"
        "5.00869,-0.052287,-2.65075;3.72297,-0.193357,-0.108558;"
        "4.85777,-0.860668,0.337892;5.8052,-0.793848,-0.201106;"
        "4.73765,-1.62678,1.50603;5.61099,-2.16479,1.88194;"
        "3.52468,-1.71469,2.1904;3.45722,-2.32196,3.09528;"
        "2.38859,-1.03029,1.72864;1.4267,-1.09209,2.24455;"
        "2.49776,-0.268364,0.573364;-2.26183,2.05492,-0.668208;"
        "-1.83522,1.76286,-1.63193;-3.00036,3.23568,-0.582173;"
        "-3.1564,3.84433,-1.47533),wU:23.24|"
    )
    EXPECTED_R_ROTAMER = (
        "[H]c1n=c(-c2c(-os(=O)(=O)C(F)(F)F)c([H])c([H])c3c([H])c([H])"
        "c([H])c([H])c23)c2c([H])c([H])c([H])c([H])c2c1[H] "
        "|(-0.329609,2.38122,-3.61273;-0.026158,2.27146,-2.56743;"
        "0.336125,1.02648,-2.1812;0.702396,0.825008,-0.934002;"
        "1.01703,-0.590455,-0.562362;-0.007714,-1.51067,-0.532668;"
        "-1.30548,-1.10616,-0.852577;-2.1204,-0.164611,0.168256;"
        "-2.55417,1.01923,-0.515766;-1.46291,-0.151398,1.44906;"
        "-3.57423,-1.28562,0.320271;-4.16542,-1.41565,-0.847695;"
        "-3.16121,-2.46571,0.747804;-4.41207,-0.76474,1.19438;"
        "0.181508,-2.86971,-0.206966;-0.687662,-3.52873,-0.187174;"
        "1.44931,-3.31282,0.079351;1.62504,-4.36094,0.333044;"
        "2.55451,-2.41778,0.057183;3.87361,-2.86277,0.348232;"
        "4.02499,-3.91674,0.595343;4.93563,-1.98941,0.318758;"
        "5.94394,-2.3428,0.544269;4.72459,-0.626471,-0.008915;"
        "5.57315,0.060292,-0.037079;3.46034,-0.164367,-0.295324;"
        "3.30174,0.885311,-0.552026;2.34077,-1.04278,-0.266984;"
        "0.743402,1.86566,0.050067;1.11738,1.6475,1.4052;"
        "1.36594,0.635744,1.73181;1.14034,2.69619,2.29419;"
        "1.42424,2.52279,3.3341;0.783222,4.00428,1.87401;"
        "0.805052,4.82508,2.59435;0.401105,4.23926,0.574492;"
        "0.114218,5.24117,0.245832;0.370036,3.17459,-0.368519;"
        "-0.025328,3.3532,-1.71845;-0.329238,4.34162,-2.07033)|"
    )
    EXPECTED_S_ROTAMER = (
        "[H]c1n=c(-c2c(-os(=O)(=O)C(F)(F)F)c([H])c([H])c3c([H])c([H])"
        "c([H])c([H])c23)c2c([H])c([H])c([H])c([H])c2c1[H] "
        "|(-0.329441,-2.38184,-3.61277;-0.026136,-2.27183,-2.56745;"
        "0.33626,-1.0268,-2.18149;0.702361,-0.825125,-0.934248;"
        "1.01698,0.590384,-0.562766;-0.007759,1.51064,-0.533148;"
        "-1.3055,1.10597,-0.853096;-2.12019,0.16467,0.168019;"
        "-1.46257,0.151709,1.44877;-2.55422,-1.01919,-0.515773;"
        "-3.57381,1.28579,0.320452;-3.16032,2.46582,0.747959;"
        "-4.16523,1.41612,-0.847356;-4.41149,0.765061,1.19469;"
        "0.181532,2.86961,-0.207313;-0.687458,3.52887,-0.187547;"
        "1.44935,3.31259,0.07922;1.62509,4.36068,0.33303;"
        "2.55452,2.41752,0.057133;3.87361,2.86235,0.348425;"
        "4.02509,3.91629,0.595642;4.93552,1.98887,0.319058;"
        "5.94384,2.34212,0.544769;4.7244,0.625956,-0.008734;"
        "5.5729,-0.060885,-0.036781;3.46015,0.164008,-0.295385;"
        "3.30143,-0.885637,-0.552158;2.3407,1.04257,-0.26718;"
        "0.743183,-1.86557,0.050012;1.11714,-1.64713,1.40512;"
        "1.36576,-0.635314,1.7315;1.1399,-2.69559,2.29437;"
        "1.42376,-2.522,3.33425;0.782587,-4.00373,1.87448;"
        "0.804249,-4.82437,2.595;0.400478,-4.23897,0.575013;"
        "0.11344,-5.24092,0.246622;0.369622,-3.17453,-0.368276;"
        "-0.025656,-3.35336,-1.71819;-0.329772,-4.34179,-2.06987)|"
    )

    @staticmethod
    def _load_molecule(filepath):
        mol = Molecule.from_filepath(filepath)
        if isinstance(mol, list):
            mol = mol[-1]
        return mol

    @staticmethod
    def _smiles_core(cxsmiles):
        """Extract the SMILES part before the CX coordinate extension."""
        return cxsmiles.split(" |")[0]

    @staticmethod
    def _load_expected(filepath):
        with open(filepath, "r") as f:
            return f.read().strip()

    def test_regression_cxsmiles(self, cxsmiles_normal_file):
        """CXSMILES for a simple molecule should be deterministic across
        repeated calls."""
        mol = self._load_molecule(cxsmiles_normal_file)
        for _ in range(3):
            assert mol.cxsmiles == self.EXPECTED_NORMAL

    def test_r_enantiomer_cxsmiles(self, cxsmiles_r_enantiomer_file):
        """CXSMILES for the R-enantiomer should match the expected value."""
        mol = self._load_molecule(cxsmiles_r_enantiomer_file)
        assert mol.cxsmiles == self.EXPECTED_R_ENANTIOMER

    def test_s_enantiomer_cxsmiles(self, cxsmiles_s_enantiomer_file):
        """CXSMILES for the S-enantiomer should match the expected value."""
        mol = self._load_molecule(cxsmiles_s_enantiomer_file)
        assert mol.cxsmiles == self.EXPECTED_S_ENANTIOMER

    def test_enantiomers_differ(
        self, cxsmiles_r_enantiomer_file, cxsmiles_s_enantiomer_file
    ):
        """R and S enantiomers must produce different CXSMILES.
        The SMILES core itself differs (@@/@ chirality annotation)."""
        mol_r = self._load_molecule(cxsmiles_r_enantiomer_file)
        mol_s = self._load_molecule(cxsmiles_s_enantiomer_file)
        core_r = self._smiles_core(mol_r.cxsmiles)
        core_s = self._smiles_core(mol_s.cxsmiles)
        # SMILES cores must differ (stereo annotation)
        assert core_r != core_s
        # Full CXSMILES must differ
        assert mol_r.cxsmiles != mol_s.cxsmiles

    def test_r_rotamer_cxsmiles(self, cxsmiles_r_rotamer_file):
        """CXSMILES for the R-rotamer should match the expected value."""
        mol = self._load_molecule(cxsmiles_r_rotamer_file)
        assert mol.cxsmiles == self.EXPECTED_R_ROTAMER

    def test_s_rotamer_cxsmiles(self, cxsmiles_s_rotamer_file):
        """CXSMILES for the S-rotamer should match the expected value."""
        mol = self._load_molecule(cxsmiles_s_rotamer_file)
        assert mol.cxsmiles == self.EXPECTED_S_ROTAMER

    def test_rotamers_differ(
        self, cxsmiles_r_rotamer_file, cxsmiles_s_rotamer_file
    ):
        """R and S rotamers must produce different CXSMILES.
        Rotamers share the same SMILES core (identical connectivity)
        but differ in the CX coordinate extension (3D geometry)."""
        mol_r = self._load_molecule(cxsmiles_r_rotamer_file)
        mol_s = self._load_molecule(cxsmiles_s_rotamer_file)
        core_r = self._smiles_core(mol_r.cxsmiles)
        core_s = self._smiles_core(mol_s.cxsmiles)
        # Rotamers share the same SMILES core
        assert core_r == core_s
        # Full CXSMILES must differ (different 3D coordinates)
        assert mol_r.cxsmiles != mol_s.cxsmiles

    def test_large_molecule_c2_cxsmiles(
        self, cxsmiles_large_molecule_c2_file, cxsmiles_expected_large_c2_file
    ):
        """CXSMILES for a large molecule (c2) should match the expected value."""
        expected = self._load_expected(cxsmiles_expected_large_c2_file)
        mol = self._load_molecule(cxsmiles_large_molecule_c2_file)
        assert mol.cxsmiles == expected

    def test_large_molecule_c3_cxsmiles(
        self, cxsmiles_large_molecule_c3_file, cxsmiles_expected_large_c3_file
    ):
        """CXSMILES for a large molecule (c3) should match the expected value."""
        expected = self._load_expected(cxsmiles_expected_large_c3_file)
        mol = self._load_molecule(cxsmiles_large_molecule_c3_file)
        assert mol.cxsmiles == expected

    def test_large_molecules_differ(
        self, cxsmiles_large_molecule_c2_file, cxsmiles_large_molecule_c3_file
    ):
        """Two different large molecules should produce different CXSMILES."""
        mol_c2 = self._load_molecule(cxsmiles_large_molecule_c2_file)
        mol_c3 = self._load_molecule(cxsmiles_large_molecule_c3_file)
        assert mol_c2.cxsmiles != mol_c3.cxsmiles
