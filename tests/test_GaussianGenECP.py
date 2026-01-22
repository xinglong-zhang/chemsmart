from ase import Atoms

from chemsmart.io.gaussian.gengenecp import GenGenECPSection
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.utils import two_files_have_similar_contents


class TestGaussianGenGenECP:
    def test_genecp_from_base_api(
        self, tmpdir, reference_genecp_txt_file_from_api
    ):
        genecp_section = GenGenECPSection.from_bse_api(
            light_elements=["C", "H", "O"],
            light_elements_basis="def2-SVP",
            heavy_elements=["Pd"],
            heavy_elements_basis="def2-TZVPPD",
        )

        assert genecp_section.heavy_elements == ["Pd"]
        assert genecp_section.heavy_elements_basis == "def2-TZVPPD".lower()
        assert genecp_section.light_elements == ["H", "C", "O"]
        assert (
            genecp_section.light_elements_basis == "def2svp"
        )  # no hyphen in this one, since this is the required format for Gausian input

        # temp file for writing genecp section
        genecp_txt_file_from_api_strings = tmpdir.join("genecp_from_api.txt")
        with open(genecp_txt_file_from_api_strings, "w") as f:
            for line in genecp_section.string_list:
                f.write(f"{line}\n")

        assert two_files_have_similar_contents(
            reference_file=reference_genecp_txt_file_from_api,
            output_file=genecp_txt_file_from_api_strings,
            ignored_string="Version",
        )

    def test_genecp_from_comfile(
        self, tmpdir, gaussian_opt_genecp_inputfile, genecp_txt_file_from_web
    ):
        genecp_section = GenGenECPSection.from_comfile(
            gaussian_opt_genecp_inputfile
        )
        assert genecp_section.heavy_elements == ["Pd"]
        assert genecp_section.heavy_elements_basis == "Def2-TZVPPD".lower()
        assert set(genecp_section.light_elements) == {"C", "H", "O"}
        assert genecp_section.light_elements_basis == "def2svp"

        # temp file for writing genecp section
        genecp_txt_file_from_web_strings = tmpdir.join("genecp_from_api.txt")
        with open(genecp_txt_file_from_web_strings, "w") as f:
            for line in genecp_section.string_list:
                f.write(f"{line}\n")

        assert two_files_have_similar_contents(
            reference_file=genecp_txt_file_from_web,
            output_file=genecp_txt_file_from_web_strings,
            ignored_string="Version",
        )

    def test_genecp_from_modred_gen_comfile(
        self, tmpdir, modred_gen_inputfile, gen_txt_file_from_web
    ):
        genecp_section = GenGenECPSection.from_comfile(modred_gen_inputfile)
        assert genecp_section.heavy_elements == ["Br"]
        assert genecp_section.heavy_elements_basis == "Def2-TZVPPD".lower()
        assert set(genecp_section.light_elements) == {"C", "H", "O"}
        assert genecp_section.light_elements_basis == "def2svp"

        # temp file for writing genecp section
        genecp_txt_file_from_web_strings = tmpdir.join("genecp_from_api.txt")
        with open(genecp_txt_file_from_web_strings, "w") as f:
            for line in genecp_section.string_list:
                f.write(f"{line}\n")

        assert two_files_have_similar_contents(
            reference_file=gen_txt_file_from_web,
            output_file=genecp_txt_file_from_web_strings,
            ignored_string="Version",
        )

    def test_genecp_from_modred_genecp_comfile(self, modred_genecp_inputfile):
        genecp_section = GenGenECPSection.from_comfile(modred_genecp_inputfile)
        assert genecp_section.heavy_elements == ["Pd"]
        assert genecp_section.heavy_elements_basis == "Def2-TZVPPD".lower()
        assert set(genecp_section.light_elements) == {"C", "H", "O"}
        assert genecp_section.light_elements_basis == "def2svp"

    def test_genecp_from_modred_genecp_solvent_comfile(
        self, modred_genecp_custom_solvent_inputfile
    ):
        genecp_section = GenGenECPSection.from_comfile(
            modred_genecp_custom_solvent_inputfile
        )
        assert genecp_section.heavy_elements == ["Pd"]
        assert genecp_section.heavy_elements_basis == "Def2-TZVPPD".lower()
        assert set(genecp_section.light_elements) == {
            "C",
            "H",
            "Cl",
            "N",
            "O",
            "S",
        }
        assert genecp_section.light_elements_basis == "def2svp"


class TestGenGenECPBasisDetermination:
    """Test automatic determination of gen vs genecp basis keywords
    based on elements present."""

    def test_periodic_table_requires_ecp(self):
        """Test that requires_ecp correctly identifies elements needing ECPs."""
        pt = PeriodicTable()

        # Elements with Z <= 36 should not require ECP
        assert not pt.requires_ecp("H")  # Z=1
        assert not pt.requires_ecp("C")  # Z=6
        assert not pt.requires_ecp("Br")  # Z=35
        assert not pt.requires_ecp("Kr")  # Z=36 (last element in period 4)

        # Elements with Z > 36 should require ECP
        assert pt.requires_ecp("Rb")  # Z=37 (first element in period 5)
        assert pt.requires_ecp("I")  # Z=53
        assert pt.requires_ecp("Pd")  # Z=46
        assert pt.requires_ecp("Ir")  # Z=77

    def test_determine_basis_keyword_with_br(self):
        """Test that molecules with only Br use 'gen' instead of 'genecp'."""
        # Create molecule with Br
        br_molecule = Atoms(
            "BrC2H4O2",
            positions=[
                [0, 0, 0],
                [1, 0, 0],
                [2, 0, 0],
                [3, 0, 0],
                [4, 0, 0],
                [5, 0, 0],
                [6, 0, 0],
                [7, 0, 0],
                [8, 0, 0],
            ],
        )
        symbols = br_molecule.get_chemical_symbols()
        positions = br_molecule.get_positions()
        br_mol = Molecule(
            symbols=symbols, positions=positions, charge=0, multiplicity=1
        )

        # Create settings with genecp basis
        settings = GaussianJobSettings(
            functional="mn15",
            basis="genecp",
            heavy_elements=["Ir", "Br"],
            heavy_elements_basis="def2-SVPD",
            light_elements_basis="def2SVP",
            charge=0,
            multiplicity=1,
        )

        # Determine basis should return 'gen' since Br (Z=35) <= 36
        determined_basis = settings.determine_basis_keyword(br_mol)
        assert determined_basis == "gen"

    def test_determine_basis_keyword_with_i(self):
        """Test that molecules with I use 'genecp'."""
        # Create molecule with I
        i_molecule = Atoms(
            "IC2H4O2",
            positions=[
                [0, 0, 0],
                [1, 0, 0],
                [2, 0, 0],
                [3, 0, 0],
                [4, 0, 0],
                [5, 0, 0],
                [6, 0, 0],
                [7, 0, 0],
                [8, 0, 0],
            ],
        )
        symbols = i_molecule.get_chemical_symbols()
        positions = i_molecule.get_positions()
        i_mol = Molecule(
            symbols=symbols, positions=positions, charge=0, multiplicity=1
        )

        # Create settings with genecp basis
        settings = GaussianJobSettings(
            functional="mn15",
            basis="genecp",
            heavy_elements=["Ir", "I"],
            heavy_elements_basis="def2-SVPD",
            light_elements_basis="def2SVP",
            charge=0,
            multiplicity=1,
        )

        # Determine basis should return 'genecp' since I (Z=53) > 36
        determined_basis = settings.determine_basis_keyword(i_mol)
        assert determined_basis == "genecp"

    def test_determine_basis_keyword_with_pd(self):
        """Test that molecules with Pd use 'genecp'."""
        # Create molecule with Pd
        pd_molecule = Atoms(
            "PdC2H4O2",
            positions=[
                [0, 0, 0],
                [1, 0, 0],
                [2, 0, 0],
                [3, 0, 0],
                [4, 0, 0],
                [5, 0, 0],
                [6, 0, 0],
                [7, 0, 0],
                [8, 0, 0],
            ],
        )
        symbols = pd_molecule.get_chemical_symbols()
        positions = pd_molecule.get_positions()
        pd_mol = Molecule(
            symbols=symbols, positions=positions, charge=0, multiplicity=1
        )

        # Create settings with genecp basis
        settings = GaussianJobSettings(
            functional="mn15",
            basis="genecp",
            heavy_elements=["Pd"],
            heavy_elements_basis="def2-TZVPPD",
            light_elements_basis="def2SVP",
            charge=0,
            multiplicity=1,
        )

        # Determine basis should return 'genecp' since Pd (Z=46) > 36
        determined_basis = settings.determine_basis_keyword(pd_mol)
        assert determined_basis == "genecp"

    def test_determine_basis_keyword_no_heavy_elements(self):
        """Test that molecules with no heavy elements use light elements basis."""
        # Create molecule with only light elements
        h2o_molecule = Atoms(
            "H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        )
        symbols = h2o_molecule.get_chemical_symbols()
        positions = h2o_molecule.get_positions()
        h2o_mol = Molecule(
            symbols=symbols, positions=positions, charge=0, multiplicity=1
        )

        # Create settings with genecp basis
        settings = GaussianJobSettings(
            functional="mn15",
            basis="genecp",
            heavy_elements=["Ir", "Br"],
            heavy_elements_basis="def2-SVPD",
            light_elements_basis="def2SVP",
            charge=0,
            multiplicity=1,
        )

        # Determine basis should use light elements since no heavy elements present
        determined_basis = settings.determine_basis_keyword(h2o_mol)
        assert determined_basis == "def2svp"

    def test_determine_basis_keyword_non_gen_basis(self):
        """Test that non-gen/genecp basis keywords are returned unchanged."""
        # Create molecule
        br_molecule = Atoms(
            "BrC2H4O2",
            positions=[
                [0, 0, 0],
                [1, 0, 0],
                [2, 0, 0],
                [3, 0, 0],
                [4, 0, 0],
                [5, 0, 0],
                [6, 0, 0],
                [7, 0, 0],
                [8, 0, 0],
            ],
        )
        symbols = br_molecule.get_chemical_symbols()
        positions = br_molecule.get_positions()
        br_mol = Molecule(
            symbols=symbols, positions=positions, charge=0, multiplicity=1
        )

        # Create settings with regular basis
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="6-31g(d)",
            charge=0,
            multiplicity=1,
        )

        # Determine basis should return original basis since it's not gen/genecp
        determined_basis = settings.determine_basis_keyword(br_mol)
        assert determined_basis == "6-31g(d)"

    def test_determine_basis_keyword_mixed_elements(self):
        """Test molecule with both gen and genecp elements uses genecp."""
        # Create molecule with both Br (gen) and I (genecp)
        mixed_molecule = Atoms(
            "BrIC2H4O2",
            positions=[
                [0, 0, 0],
                [1, 0, 0],
                [2, 0, 0],
                [3, 0, 0],
                [4, 0, 0],
                [5, 0, 0],
                [6, 0, 0],
                [7, 0, 0],
                [8, 0, 0],
                [9, 0, 0],
            ],
        )
        symbols = mixed_molecule.get_chemical_symbols()
        positions = mixed_molecule.get_positions()
        mixed_mol = Molecule(
            symbols=symbols, positions=positions, charge=0, multiplicity=1
        )

        # Create settings with both Br and I as heavy elements
        settings = GaussianJobSettings(
            functional="mn15",
            basis="genecp",
            heavy_elements=["Br", "I"],
            heavy_elements_basis="def2-SVPD",
            light_elements_basis="def2SVP",
            charge=0,
            multiplicity=1,
        )

        # Should use 'genecp' since I requires ECP (even though Br doesn't)
        determined_basis = settings.determine_basis_keyword(mixed_mol)
        assert determined_basis == "genecp"
