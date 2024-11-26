from chemsmart.io.gaussian.gengenecp import GenGenECPSection
from chemsmart.utils.utils import two_files_have_similar_contents

class TestGaussianGenGenECP:
    def test_genecp_from_base_api(self, tmpdir, reference_genecp_txt_file_from_api):
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

    def test_genecp_from_comfile(self, tmpdir, gaussian_opt_genecp_inputfile, genecp_txt_file_from_web):
        genecp_section = GenGenECPSection.from_comfile(gaussian_opt_genecp_inputfile)
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
        self, tmpdir, modred_gen_comfile, gen_txt_file_from_web
    ):
        genecp_section = GenGenECPSection.from_comfile(modred_gen_comfile)
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

    def test_genecp_from_modred_genecp_comfile(self, modred_genecp_comfile):
        genecp_section = GenGenECPSection.from_comfile(modred_genecp_comfile)
        assert genecp_section.heavy_elements == ["Pd"]
        assert genecp_section.heavy_elements_basis == "Def2-TZVPPD".lower()
        assert set(genecp_section.light_elements) == {"C", "H", "O"}
        assert genecp_section.light_elements_basis == "def2svp"

    def test_genecp_from_modred_genecp_solvent_comfile(
        self, modred_genecp_custom_solvent_comfile
    ):
        genecp_section = GenGenECPSection.from_comfile(
            modred_genecp_custom_solvent_comfile
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
