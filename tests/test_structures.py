import os
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
)
from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.utils.cluster import is_pubchem_network_available


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

        # set correct charge and multiplicity for molecules as needed by pymatgen checks
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
    @pytest.mark.skipif(
        not is_pubchem_network_available(),
        reason="Network to pubchem is unavailable",
    )
    def test_stereochemistry_handling(self):
        """Test preservation of stereochemical information."""
        methyl_3_hexane = Molecule.from_pubchem("11507")
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
        # check there are 6 aromatic C-C bonds and 6 single C-H bonds in benzene
        assert len([bond for bond in benzene.bond_orders if bond == 1.5]) == 6
        assert len([bond for bond in benzene.bond_orders if bond == 1.0]) == 6

    def test_volume(
        self, gaussian_ozone_opt_outfile, gaussian_acetone_opt_outfile
    ):
        """Test volume calculation for molecules.

        Tests various volume calculation methods:
        - voronoi_dirichlet_occupied_volume (requires pyvoro, optional)
        - crude_volume_by_vdw_radii
        - crude_volume_by_atomic_radii
        - vdw_volume
        - vdw_volume_from_rdkit
        - voronoi_dirichlet_polyhedra_occupied_volume
        """
        ozone = Molecule.from_filepath(gaussian_ozone_opt_outfile)

        # Test pyvoro-based method (optional, may not be available in Python 3.12+)
        try:
            ozone_vd_vol = ozone.voronoi_dirichlet_occupied_volume
            assert ozone_vd_vol > 0
            assert np.isclose(ozone_vd_vol, 42.796979883456515, rtol=0.01)
        except ImportError:
            pass  # pyvoro not available, skip this test

        # Test other volume methods that don't require pyvoro
        assert np.isclose(
            ozone.crude_volume_by_vdw_radii, 44.13068085447146, rtol=0.01
        )
        assert np.isclose(
            ozone.crude_volume_by_atomic_radii, 3.612781286145805, rtol=0.01
        )
        # vdw_volume uses exact lens-shaped intersection formula for overlap correction
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

        # Test pyvoro-based method (optional)
        try:
            acetone_vd_vol = acetone.voronoi_dirichlet_occupied_volume
            assert acetone_vd_vol > 0
            assert np.isclose(acetone_vd_vol, 108.73483002110545, rtol=0.01)
        except ImportError:
            pass  # pyvoro not available, skip this test

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
        # vdw_volume uses exact lens-shaped intersection formula for overlap correction
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
        """Test reading a single molecule from a CDXML file."""
        assert os.path.exists(metal_ligand_molecules_cdxml_file)
        assert os.path.isfile(metal_ligand_molecules_cdxml_file)

        cdx_file = CDXFile(filename=metal_ligand_molecules_cdxml_file)
        molecules = cdx_file.molecules

        assert isinstance(molecules, list)
        assert len(molecules) == 7
        mol = molecules[0]
        assert isinstance(mol, Molecule)
        assert mol.chemical_formula == "C32H31N5O5"
        assert mol.num_atoms == 73  # benzene with hydrogens
        assert mol.is_aromatic
