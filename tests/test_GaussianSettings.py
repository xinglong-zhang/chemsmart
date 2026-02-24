import pytest

from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.io.molecules.structure import Molecule, QMMMMolecule
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.pka import GaussianpKaJob
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianpKaJobSettings,
    GaussianQMMMJobSettings,
)
from chemsmart.jobs.settings import read_molecular_job_yaml


class TestGaussianJobSettings:
    def test_merge_dict(self):
        settings = GaussianJobSettings.default()
        settings.charge = None
        settings.multiplicity = None

        merged_settings = settings.merge({"charge": -1, "multiplicity": 3})
        assert merged_settings.charge == -1
        assert merged_settings.multiplicity == 3

    def test_merge_other_settings(self):
        settings1 = GaussianJobSettings.default()
        settings1.charge = None
        settings1.multiplicity = None

        settings2 = GaussianJobSettings.default()
        settings2.charge = 0
        settings2.multiplicity = 1

        merged_settings = settings1.merge(
            settings2, keywords=("charge", "multiplicity")
        )
        assert merged_settings.charge == 0
        assert merged_settings.multiplicity == 1

    def test_get_settings_from_yaml_gas_solv(
        self, gaussian_yaml_settings_gas_solv
    ):
        all_project_settings = read_molecular_job_yaml(
            gaussian_yaml_settings_gas_solv, program="gaussian"
        )
        opt_settings_dict = all_project_settings["opt"]
        opt_settings = GaussianJobSettings.from_dict(opt_settings_dict)
        assert isinstance(opt_settings, GaussianJobSettings)
        assert opt_settings.chk is True  # from defaults file
        assert opt_settings.functional == "b3lyp empiricaldispersion=gd3bj"
        assert opt_settings.basis == "def2svp"
        assert opt_settings.solvent_model is None
        assert opt_settings.solvent_id is None

        modred_settings_dict = all_project_settings["modred"]
        modred_settings = GaussianJobSettings.from_dict(modred_settings_dict)
        assert isinstance(modred_settings, GaussianJobSettings)
        assert modred_settings.chk is True  # from defaults file
        assert modred_settings.functional == "b3lyp empiricaldispersion=gd3bj"
        assert modred_settings.basis == "def2svp"
        assert modred_settings.solvent_model is None
        assert modred_settings.solvent_id is None

        sp_settings_dict = all_project_settings["sp"]
        sp_settings = GaussianJobSettings.from_dict(sp_settings_dict)
        assert isinstance(sp_settings, GaussianJobSettings)
        assert sp_settings.chk is True  # from defaults file
        assert sp_settings.functional == "b3lyp empiricaldispersion=gd3bj"
        assert sp_settings.basis == "def2tzvp"
        assert sp_settings.solvent_model == "smd"
        assert sp_settings.solvent_id == "toluene"

    def test_get_settings_from_yaml_solv(self, gaussian_yaml_settings_solv):
        all_project_settings = read_molecular_job_yaml(
            gaussian_yaml_settings_solv, program="gaussian"
        )
        opt_settings_dict = all_project_settings["opt"]
        opt_settings = GaussianJobSettings.from_dict(opt_settings_dict)
        assert isinstance(opt_settings, GaussianJobSettings)
        assert opt_settings.chk is True  # from defaults file
        assert opt_settings.functional == "m062x"
        assert opt_settings.basis == "def2tzvp"
        assert opt_settings.solvent_model == "smd"
        assert opt_settings.solvent_id == "toluene"

        modred_settings_dict = all_project_settings["modred"]
        modred_settings = GaussianJobSettings.from_dict(modred_settings_dict)
        assert isinstance(modred_settings, GaussianJobSettings)
        assert modred_settings.chk is True  # from defaults file
        assert modred_settings.functional == "m062x"
        assert modred_settings.basis == "def2tzvp"
        assert modred_settings.solvent_model == "smd"
        assert modred_settings.solvent_id == "toluene"

        sp_settings_dict = all_project_settings["sp"]
        sp_settings = GaussianJobSettings.from_dict(sp_settings_dict)
        assert isinstance(sp_settings, GaussianJobSettings)
        assert sp_settings.chk is True  # from defaults file
        assert sp_settings.functional == "m062x"
        assert sp_settings.basis == "def2tzvp"
        assert sp_settings.solvent_model == "smd"
        assert sp_settings.solvent_id == "toluene"

    def test_read_gaussian_hf_comfile(self, hf_com_filepath):
        settings = GaussianJobSettings.from_comfile(hf_com_filepath)
        assert isinstance(settings, GaussianJobSettings)
        assert settings.ab_initio == "hf"
        assert settings.functional is None
        assert settings.basis == "6-31g"
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_read_gaussian_settings_from_orca_inp(self, water_sp_input_path):
        settings = GaussianJobSettings.from_inpfile(water_sp_input_path)
        assert isinstance(settings, GaussianJobSettings)
        assert settings.ab_initio == "hf"
        assert settings.functional is None
        assert settings.basis == "def2-svp"
        assert settings.solvent_model is None
        assert settings.solvent_id is None


class TestGaussianQMMMJobSettings:
    def test_qmmm_settings(self):
        settings1 = GaussianQMMMJobSettings(
            high_level_functional="b3lyp",
            high_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="opt",
            freq=True,
        )
        assert settings1.route_string == "# opt freq oniom(b3lyp/6-31g(d):uff)"

        settings2 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="sp",
        )
        assert (
            settings2.route_string
            == "# oniom(mn15/def2svp:b3lyp/6-31g(d):uff)"
        )

        settings3 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            high_level_force_field="uff",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="sp",
        )
        # assert settings3.route_string == "#
        # oniom(mn15/def2svp:uff:b3lyp/6-31g(d):uff)"
        # ValueError: For high level of theory, one should
        # specify only functional/basis or force field!
        with pytest.raises(ValueError):
            settings3.route_string

        # settings with solvent specification
        settings4 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            solvent_model="smd",
            solvent_id="toluene",
            parent_jobtype="sp",
        )
        assert (
            settings4.route_string
            == "# oniom(mn15/def2svp:b3lyp/6-31g(d):uff) scrf=(smd,solvent=toluene)"
        )

        # settings with solvent specification for opt job
        settings5 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="opt",
            freq=True,
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings5.route_string
            == "# opt freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) scrf=(smd,solvent=toluene)"
        )

        # settings with solvent specification for ts job
        settings5 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="ts",
            freq=True,
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings5.route_string
            == "# opt=(ts,calcfc,noeigentest) freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) scrf=(smd,solvent=toluene)"
        )

        # settings with solvent specification for ts job
        settings6 = GaussianQMMMJobSettings(
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            high_level_atoms=[1, 2, 3],
            parent_jobtype="ts",
            freq=False,
            numfreq=True,
            solvent_model="smd",
            solvent_id="toluene",
        )
        assert (
            settings6.route_string
            == "# opt=(ts,calcfc,noeigentest) freq oniom(mn15/def2svp:b3lyp/6-31g(d):uff) scrf=(smd,solvent=toluene)"
        )

    def test_qmmm_settings_for_atoms(
        self,
        gaussian_inputs_test_directory,
        gaussian_semiempirical_pm6_output_file,
    ):
        mol1 = QMMMMolecule(
            molecule=Molecule.from_filepath(
                gaussian_semiempirical_pm6_output_file
            )
        )

        settings1 = QMMMMolecule(
            symbols=mol1.symbols,
            positions=mol1.positions,
            high_level_atoms=[1, 2, 3, 8, 9, 10],
            bonded_atoms=[[1, 2], [2, 3], [8, 9], [9, 10]],
        )
        assert settings1.high_level_atoms == [1, 2, 3, 8, 9, 10]
        assert settings1.partition_level_strings == [
            "H",
            "H",
            "H",
            "L",
            "L",
            "L",
            "L",
            "H",
            "H",
            "H",
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

        # Test for 3-layer ONIOM calculation with example
        # input dppeFeCl2_phenyldioxazolone_qmmm.com
        # mol2 = QMMM(molecule=Molecule._read_gaussian_inputfile(
        #     os.path.join(
        #         gaussian_inputs_test_directory,
        #         "qmmm/dppeFeCl2_phenyldioxazolone_qmmm.com",
        #     ),
        # ))
        # settings2 = QMMM(
        #     symbols=mol2.symbols,
        #     positions=mol2.positions,
        #     high_level_atoms="[18-28,29-39,40-50,51-61,62-72]",
        #     medium_level_atoms=[1, 2, 3, 16],
        #     bonded_atoms=[[2, 18], [2, 29], [1, 40], [1, 51], [16, 62]],
        # )
        # assert settings2.high_level_atoms == list(range(18, 29)) + list(
        #     range(29, 40)
        # ) + list(range(40, 51)) + list(range(51, 62)) + list(range(62, 73))
        # assert settings2.medium_level_atoms == [1, 2, 3, 16]
        # assert settings2.low_level_atoms == list(range(4, 16)) + [17]

    def test_qmmm_settings_for_charge_and_multiplicity(self):
        # test cases for 3-layer ONIOM model
        settings1 = GaussianQMMMJobSettings(
            jobtype="opt",
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
        )
        assert (
            settings1.charge_and_multiplicity_string
            == "0 1 0 1 0 1 0 1 0 1 0 1"
        )

        settings2 = GaussianQMMMJobSettings(
            jobtype="sp",
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            int_charge=1,
            int_multiplicity=3,
        )
        assert (
            settings2.charge_and_multiplicity_string
            == "0 1 1 3 1 3 1 3 1 3 1 3"
        )

        settings3 = GaussianQMMMJobSettings(
            jobtype="ts",
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            low_level_force_field="uff",
            real_charge=0,
            real_multiplicity=1,
            int_charge=-1,
            int_multiplicity=2,
            model_charge=0,
            model_multiplicity=1,
        )
        assert (
            settings3.charge_and_multiplicity_string
            == "0 1 -1 2 -1 2 0 1 0 1 0 1"
        )

        # test cases for 2-layer ONIOM model

        settings4 = GaussianQMMMJobSettings(
            jobtype="opt",
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            real_charge=0,
            real_multiplicity=1,
        )
        assert settings4.charge_and_multiplicity_string == "0 1 0 1 0 1"

        settings5 = GaussianQMMMJobSettings(
            jobtype="sp",
            high_level_functional="mn15",
            high_level_basis="def2svp",
            medium_level_functional="b3lyp",
            medium_level_basis="6-31g(d)",
            real_charge=0,
            real_multiplicity=1,
            model_charge=-1,
            model_multiplicity=2,
        )
        assert settings5.charge_and_multiplicity_string == "0 1 -1 2 -1 2"


class TestGaussianRoute:
    gas_route = "opt=(ts,calcfc,noeigentest) freq m062x def2svp"
    gas_route2 = "opt b3lyp 6-31G(d) empiricaldispersion=gd3bj"
    solv_route = "mn15 def2svp scrf=(smd,solvent=DiChloroMethane)"
    solv_route2 = "# opt=(calcfc,ts,noeigentest,maxstep=5) freq\nscrf=(smd,solvent=chlorobenzene) nosymm scf=qc def2svp m06"

    def test_settings_from_route(self):
        route_object1 = GaussianRoute(route_string=self.gas_route)
        assert isinstance(route_object1, object)
        assert route_object1.jobtype == "ts"
        assert route_object1.freq is True
        assert route_object1.functional == "m062x"
        assert route_object1.basis == "def2svp"
        assert route_object1.solvent_model is None
        assert route_object1.solvent_id is None

    def test_settings_from_route2(self):
        route_object2 = GaussianRoute(route_string=self.gas_route2)
        assert isinstance(route_object2, object)
        assert route_object2.jobtype == "opt"
        assert route_object2.freq is False
        assert route_object2.functional == "b3lyp empiricaldispersion=gd3bj"
        assert route_object2.basis == "6-31G(d)".lower()
        assert route_object2.solvent_model is None
        assert route_object2.solvent_id is None

    def test_settings_from_route3(self):
        route_object3 = GaussianRoute(route_string=self.solv_route)
        assert isinstance(route_object3, object)
        assert route_object3.jobtype == "sp"
        assert route_object3.freq is False
        assert route_object3.functional == "mn15"
        assert route_object3.basis == "def2svp"
        assert route_object3.solvent_model == "smd"
        assert route_object3.solvent_id == "DiChloroMethane".lower()

    def test_settings_from_route4(self):
        route_object4 = GaussianRoute(route_string=self.solv_route2)
        assert isinstance(route_object4, object)
        assert route_object4.jobtype == "ts"
        assert route_object4.freq is True
        assert route_object4.functional == "m06"
        assert route_object4.basis == "def2svp"
        assert route_object4.solvent_model == "smd"
        assert route_object4.solvent_id == "chlorobenzene"


class TestGaussianJobFromComFile:
    def test_reads_com_file(self, gaussian_opt_inputfile):
        com_settings = GaussianJobSettings.from_comfile(gaussian_opt_inputfile)
        assert com_settings.chk is True
        assert com_settings.jobtype == "opt"
        assert com_settings.freq is True
        assert com_settings.functional == "m062x"
        assert com_settings.basis == "def2svp"
        assert com_settings.charge == 0
        assert com_settings.multiplicity == 1
        assert com_settings.solvent_model is None
        assert com_settings.solvent_id is None

    def test_update_solvent(self, gaussian_opt_inputfile, tmpdir):
        com_settings = GaussianJobSettings.from_comfile(gaussian_opt_inputfile)
        assert com_settings.solvent_model is None
        assert com_settings.solvent_id is None
        com_settings.update_solvent(solvent_model="smd", solvent_id="toluene")
        assert com_settings.solvent_model == "smd"
        assert com_settings.solvent_id == "toluene"

    def test_include_solvent(self, gaussian_opt_inputfile, tmpdir):
        com_settings = GaussianJobSettings.from_filepath(
            gaussian_opt_inputfile
        )
        assert com_settings.solvent_model is None
        assert com_settings.solvent_id is None

        com_settings.modify_solvent(
            remove_solvent=False, solvent_model="cpcm", solvent_id="water"
        )
        assert com_settings.solvent_model == "cpcm"
        assert com_settings.solvent_id == "water"

        com_settings.modify_solvent(remove_solvent=True)
        assert com_settings.solvent_model is None
        assert com_settings.solvent_id is None


class TestGaussianJobFromLogFile:
    def test_accumulates_settings(self, tmpdir, gaussian_ts_genecp_outfile):
        settings = GaussianJobSettings.from_logfile(gaussian_ts_genecp_outfile)
        assert settings.functional == "mn15"
        assert settings.basis == "genecp"
        assert settings.solvent_model is None
        assert settings.solvent_id is None

        settings.functional = "b3lyp"
        assert settings.solvent_model is None
        assert settings.functional == "b3lyp"

        settings.basis = "def2tzvp"
        assert settings.functional == "b3lyp"
        assert settings.basis == "def2tzvp"

        settings.update_solvent(solvent_model="smd", solvent_id="toluene")
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_reads_gaussian_ts_genecp_outfile(
        self, tmpdir, gaussian_ts_genecp_outfile
    ):
        settings = GaussianJobSettings.from_logfile(gaussian_ts_genecp_outfile)
        assert settings.jobtype == "ts"
        assert settings.functional == "mn15"
        assert settings.basis == "genecp"
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_reads_gaussian_pm6_outfile(
        self, tmpdir, gaussian_semiempirical_pm6_output_file
    ):
        settings = GaussianJobSettings.from_logfile(
            gaussian_semiempirical_pm6_output_file
        )
        assert settings.jobtype == "opt"
        assert settings.ab_initio is None
        assert settings.functional is None
        assert settings.basis is None
        assert settings.semiempirical == "PM6"
        assert settings.solvent_model is None
        assert settings.solvent_id is None


class TestGaussianPBCJob:
    def test_writes_gaussian_input_from_pbc_comfile(
        self, tmpdir, gaussian_pbc_3d_outputfile
    ):
        settings = GaussianJobSettings.from_filepath(
            filepath=gaussian_pbc_3d_outputfile
        )
        assert settings.jobtype == "sp"
        assert settings.functional.lower() == "pbepbe"
        assert settings.basis.lower() == "6-31g(d,p)/auto"
        assert settings.additional_route_parameters.lower() == "scf=tight"


class TestGaussianpKaJobSettings:
    """Tests for GaussianpKaJobSettings class."""

    def test_init_custom_values(self):
        """Test initialization with custom values."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            reference="acetic_acid",
            solvent_model="PCM",
            solvent_id="water",
            thermodynamic_cycle="isodesmic",
            charge=0,  # Protonated form charge (inherited from parent)
            multiplicity=1,  # Protonated form multiplicity (inherited from parent)
            conjugate_base_charge=-1,
            conjugate_base_multiplicity=1,
            functional="B3LYP",
            basis="6-311+G(d,p)",
        )
        assert settings.proton_index == 10
        assert settings.reference == "acetic_acid"
        assert settings.solvent_model == "PCM"
        assert settings.solvent_id == "water"
        assert settings.thermodynamic_cycle == "isodesmic"
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.protonated_charge == 0
        assert settings.protonated_multiplicity == 1
        assert settings.conjugate_base_charge == -1
        assert settings.conjugate_base_multiplicity == 1
        assert settings.functional == "B3LYP"
        assert settings.basis == "6-311+G(d,p)"

    def test_gas_phase_optimization_settings(self, single_molecule_xyz_file):
        """Test that gas phase optimization has no solvent."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
            solvent_model="SMD",
            solvent_id="water",
        )

        prot_settings, conj_base_settings = (
            settings._create_gas_phase_job_settings(mol)
        )

        # Gas phase should have no solvent
        assert prot_settings.solvent_model is None
        assert prot_settings.solvent_id is None
        assert conj_base_settings.solvent_model is None
        assert conj_base_settings.solvent_id is None
        # Should use same functional/basis
        assert prot_settings.functional == "B3LYP"
        assert prot_settings.basis == "6-31G*"

    def test_solution_phase_sp_settings(self, single_molecule_xyz_file):
        """Test that solution phase SP uses same level of theory with solvent."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
            solvent_model="SMD",
            solvent_id="water",
        )

        prot_sp_settings, conj_base_sp_settings = (
            settings._create_solution_phase_sp_settings(mol)
        )

        # Solution phase should have solvent
        assert prot_sp_settings.solvent_model == "SMD"
        assert prot_sp_settings.solvent_id == "water"
        assert conj_base_sp_settings.solvent_model == "SMD"
        assert conj_base_sp_settings.solvent_id == "water"
        # Should use SAME functional/basis as gas phase for error cancellation
        assert prot_sp_settings.functional == "B3LYP"
        assert prot_sp_settings.basis == "6-31G*"
        assert conj_base_sp_settings.functional == "B3LYP"
        assert conj_base_sp_settings.basis == "6-31G*"

    def test_protonated_charge_multiplicity_properties(self):
        """Test that protonated_charge/multiplicity are aliases for charge/multiplicity."""
        settings = GaussianpKaJobSettings(
            proton_index=10,
            charge=2,
            multiplicity=3,
        )
        # Properties should return the same values
        assert settings.protonated_charge == settings.charge
        assert settings.protonated_multiplicity == settings.multiplicity

        # Setting via property should update the underlying attribute
        settings.protonated_charge = 5
        assert settings.charge == 5
        settings.protonated_multiplicity = 4
        assert settings.multiplicity == 4

    def test_create_conjugate_base_molecule(self, single_molecule_xyz_file):
        """Test creating conjugate base molecule by removing a proton."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        original_num_atoms = len(mol)

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        assert (
            len(h_indices) > 0
        ), "Test molecule must have at least one hydrogen"

        proton_index = h_indices[0]
        settings = GaussianpKaJobSettings(proton_index=proton_index)

        conjugate_base = settings._create_conjugate_base_molecule(mol)

        # Check that one atom was removed
        assert len(conjugate_base) == original_num_atoms - 1
        # Check that charge decreased by 1
        assert conjugate_base.charge == -1
        # Check that multiplicity is preserved
        assert conjugate_base.multiplicity == 1

    def test_create_conjugate_base_molecule_custom_charge(
        self, single_molecule_xyz_file
    ):
        """Test creating conjugate base with custom charge/multiplicity."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 1
        mol.multiplicity = 2

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            conjugate_base_charge=0,
            conjugate_base_multiplicity=1,
        )

        conjugate_base = settings._create_conjugate_base_molecule(mol)

        # Custom values should override defaults
        assert conjugate_base.charge == 0
        assert conjugate_base.multiplicity == 1

    def test_create_conjugate_base_molecule_no_proton_index(
        self, single_molecule_xyz_file
    ):
        """Test that error is raised when proton_index is not specified."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        settings = GaussianpKaJobSettings()

        with pytest.raises(ValueError, match="proton_index must be specified"):
            settings._create_conjugate_base_molecule(mol)

    def test_create_conjugate_base_molecule_invalid_index(
        self, single_molecule_xyz_file
    ):
        """Test that error is raised for out-of-range proton index."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        settings = GaussianpKaJobSettings(proton_index=999)

        with pytest.raises(ValueError, match="out of range"):
            settings._create_conjugate_base_molecule(mol)

    def test_create_conjugate_base_molecule_not_hydrogen(
        self, single_molecule_xyz_file
    ):
        """Test that error is raised when index is not a hydrogen."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)

        # Find a non-hydrogen atom index
        non_h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s != "H"]
        assert len(non_h_indices) > 0

        settings = GaussianpKaJobSettings(proton_index=non_h_indices[0])

        with pytest.raises(ValueError, match="not hydrogen"):
            settings._create_conjugate_base_molecule(mol)

    def test_create_job_settings(self, single_molecule_xyz_file):
        """Test creating gas phase job settings for both forms."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
            solvent_model="SMD",
            solvent_id="water",
        )

        prot_settings, conj_base_settings = settings._create_job_settings(mol)

        # Check protonated settings - GAS PHASE (no solvent)
        assert isinstance(prot_settings, GaussianJobSettings)
        assert prot_settings.charge == 0
        assert prot_settings.multiplicity == 1
        assert prot_settings.functional == "B3LYP"
        assert prot_settings.basis == "6-31G*"
        assert prot_settings.jobtype == "opt"
        assert prot_settings.freq is True
        assert prot_settings.solvent_model is None  # Gas phase
        assert prot_settings.solvent_id is None

        # Check conjugate base settings - GAS PHASE (no solvent)
        assert isinstance(conj_base_settings, GaussianJobSettings)
        assert conj_base_settings.charge == -1
        assert conj_base_settings.multiplicity == 1
        assert conj_base_settings.functional == "B3LYP"
        assert conj_base_settings.basis == "6-31G*"
        assert conj_base_settings.jobtype == "opt"
        assert conj_base_settings.freq is True
        assert conj_base_settings.solvent_model is None  # Gas phase
        assert conj_base_settings.solvent_id is None

    def test_create_molecules(self, single_molecule_xyz_file):
        """Test creating both protonated and conjugate base molecules."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        original_num_atoms = len(mol)

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        prot_mol, conj_base_mol = settings._create_molecules(mol)

        # Check protonated molecule
        assert len(prot_mol) == original_num_atoms
        assert prot_mol.charge == 0
        assert prot_mol.multiplicity == 1

        # Check conjugate base molecule
        assert len(conj_base_mol) == original_num_atoms - 1
        assert conj_base_mol.charge == -1
        assert conj_base_mol.multiplicity == 1

    def test_conjugate_base_molecule_method(self, single_molecule_xyz_file):
        """Test the public conjugate_base_molecule method."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        conjugate_base = settings.conjugate_base_molecule(mol)

        assert len(conjugate_base) == len(mol) - 1
        assert conjugate_base.charge == -1

    def test_conjugate_pair_molecules_method(self, single_molecule_xyz_file):
        """Test the public conjugate_pair_molecules method."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        prot_mol, conj_base_mol = settings.conjugate_pair_molecules(mol)

        assert len(prot_mol) == len(mol)
        assert len(conj_base_mol) == len(mol) - 1
        assert prot_mol.charge == 0
        assert conj_base_mol.charge == -1

    def test_conjugate_pair_job_settings_method(
        self, single_molecule_xyz_file
    ):
        """Test the public conjugate_pair_job_settings method returns gas phase settings."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
            solvent_model="SMD",
            solvent_id="water",
        )

        prot_settings, conj_base_settings = (
            settings.conjugate_pair_job_settings(mol)
        )

        assert isinstance(prot_settings, GaussianJobSettings)
        assert isinstance(conj_base_settings, GaussianJobSettings)
        assert prot_settings.charge == 0
        assert conj_base_settings.charge == -1
        # Should be gas phase (no solvent for optimization)
        assert prot_settings.solvent_model is None
        assert conj_base_settings.solvent_model is None


class TestGaussianpKaJob:
    """Tests for GaussianpKaJob class."""

    def test_init_valid_settings(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test initialization with valid pKa settings."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
        )

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="test_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianpKaJob)
        assert job.TYPE == "g16pka"
        assert job.label == "test_pka"

    def test_init_invalid_settings_type(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test that error is raised for non-pKa settings."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)

        settings = GaussianJobSettings(functional="B3LYP", basis="6-31G*")

        with pytest.raises(
            ValueError, match="must be instance of GaussianpKaJobSettings"
        ):
            GaussianpKaJob(
                molecule=mol,
                settings=settings,
                label="test_pka",
                jobrunner=gaussian_jobrunner_no_scratch,
            )

    def test_init_no_proton_index(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test that error is raised when proton_index is not specified."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)

        settings = GaussianpKaJobSettings(functional="B3LYP", basis="6-31G*")

        with pytest.raises(ValueError, match="proton_index must be specified"):
            GaussianpKaJob(
                molecule=mol,
                settings=settings,
                label="test_pka",
                jobrunner=gaussian_jobrunner_no_scratch,
            )

    def test_pka_jobs_property(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test that pka_jobs returns both jobs."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
        )

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="test_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        pka_jobs = job.pka_jobs
        assert len(pka_jobs) == 2
        assert isinstance(pka_jobs[0], GaussianOptJob)
        assert isinstance(pka_jobs[1], GaussianOptJob)

    def test_protonated_job_property(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test protonated_job property."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
        )

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="test_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        protonated_job = job.protonated_job
        assert isinstance(protonated_job, GaussianOptJob)
        assert protonated_job.label == "test_pka_HA"
        assert protonated_job.settings.charge == 0

    def test_conjugate_base_job_property(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test conjugate_base_job property."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            functional="B3LYP",
            basis="6-31G*",
        )

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="test_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        conjugate_base_job = job.conjugate_base_job
        assert isinstance(conjugate_base_job, GaussianOptJob)
        assert conjugate_base_job.label == "test_pka_A"
        assert conjugate_base_job.settings.charge == -1

    def test_protonated_molecule_property(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test protonated_molecule property."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        original_num_atoms = len(mol)

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="test_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        protonated_mol = job.protonated_molecule
        assert len(protonated_mol) == original_num_atoms
        assert protonated_mol.charge == 0

    def test_conjugate_base_molecule_property(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test conjugate_base_molecule property."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        original_num_atoms = len(mol)

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="test_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        conjugate_base_mol = job.conjugate_base_molecule
        assert len(conjugate_base_mol) == original_num_atoms - 1
        assert conjugate_base_mol.charge == -1

    def test_job_labels(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        """Test that job labels are correctly generated."""
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1

        h_indices = [i + 1 for i, s in enumerate(mol.symbols) if s == "H"]
        proton_index = h_indices[0]

        settings = GaussianpKaJobSettings(proton_index=proton_index)

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="acetic_acid_pka",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        protonated_job, conjugate_base_job = job.pka_jobs

        assert protonated_job.label == "acetic_acid_pka_HA"
        assert conjugate_base_job.label == "acetic_acid_pka_A"

    def test_settings_class(self):
        """Test settings_class class method."""
        assert GaussianpKaJob.settings_class() == GaussianpKaJobSettings
