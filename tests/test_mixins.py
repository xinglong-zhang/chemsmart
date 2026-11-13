import os
from types import SimpleNamespace

import pytest

from chemsmart.utils.mixins import (
    GaussianFileMixin,
    ORCAFileMixin,
    XTBFileMixin,
)
from tests.helpers import (
    BaseRegistry,
    SubRegistry1,
    SubRegistry2,
    make_dummy_file,
    make_dummy_folder,
    make_dummy_folder_output,
    make_dummy_gaussian_file,
    make_mixin_instance,
    write_text_file,
)


class TestFileMixin:
    def test_file_properties(self, temp_text_file):
        dummy = make_dummy_file(temp_text_file)
        assert os.path.abspath(temp_text_file) == dummy.filepath
        assert (
            os.path.basename(temp_text_file)
            == dummy.base_filename_with_extension
        )
        assert (
            dummy.basename
            == os.path.splitext(os.path.basename(temp_text_file))[0]
        )
        assert dummy.contents == ["Line1", "Line2"]
        assert dummy.content_lines_string == "Line1\nLine2\n"
        assert dummy.forces_in_eV_per_angstrom == [None, None]
        assert dummy.input_translation_vectors == []
        assert dummy.num_energies == 2

    def test_input_translation_vectors_from_block(self, temp_text_file):
        dummy = make_dummy_file(
            temp_text_file,
            input_coordinates_block=SimpleNamespace(
                translation_vectors=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
            ),
        )
        assert dummy.input_translation_vectors == [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]

    def test_input_translation_vectors_block_without_vectors(
        self, temp_text_file
    ):
        dummy = make_dummy_file(
            temp_text_file,
            input_coordinates_block=SimpleNamespace(translation_vectors=None),
        )
        assert dummy.input_translation_vectors == []


class TestFrontierOrbitalProperties:
    """Direct tests for FileMixin's frontier-orbital cached_properties,
    isolating the multiplicity==1 vs. open-shell, and populated vs.
    empty eigenvalue-list branches independently of any real output
    parser."""

    def test_somo_energies_none_for_closed_shell(
        self, make_frontier_orbital_file
    ):
        dummy = make_frontier_orbital_file(multiplicity=1)
        assert dummy.somo_energies is None

    def test_somo_energies_none_when_no_alpha_occ(
        self, make_frontier_orbital_file
    ):
        dummy = make_frontier_orbital_file(
            multiplicity=3, alpha_occ_eigenvalues=None
        )
        assert dummy.somo_energies is None

    def test_somo_energies_open_shell(self, make_frontier_orbital_file):
        dummy = make_frontier_orbital_file(
            multiplicity=3,
            alpha_occ_eigenvalues=[-10.0, -5.0, -1.0],
        )
        assert dummy.somo_energies == [-5.0, -1.0]

    def test_lowest_and_highest_somo_energy(self, make_frontier_orbital_file):
        dummy = make_frontier_orbital_file(
            multiplicity=3,
            alpha_occ_eigenvalues=[-10.0, -5.0, -1.0],
        )
        assert dummy.lowest_somo_energy == -5.0
        assert dummy.highest_somo_energy == -1.0

    def test_lowest_and_highest_somo_energy_none_for_closed_shell(
        self, make_frontier_orbital_file
    ):
        dummy = make_frontier_orbital_file(multiplicity=1)
        assert dummy.lowest_somo_energy is None
        assert dummy.highest_somo_energy is None

    def test_alpha_beta_homo_energy_none_without_eigenvalues(
        self, make_frontier_orbital_file
    ):
        dummy = make_frontier_orbital_file()
        assert dummy.alpha_homo_energy is None
        assert dummy.beta_homo_energy is None

    def test_alpha_beta_homo_energy_populated(
        self, make_frontier_orbital_file
    ):
        dummy = make_frontier_orbital_file(
            alpha_occ_eigenvalues=[-10.0, -5.0],
            beta_occ_eigenvalues=[-9.0, -4.0],
        )
        assert dummy.alpha_homo_energy == -5.0
        assert dummy.beta_homo_energy == -4.0

    def test_alpha_beta_lumo_energy_none_without_eigenvalues(
        self, make_frontier_orbital_file
    ):
        dummy = make_frontier_orbital_file()
        assert dummy.alpha_lumo_energy is None
        assert dummy.beta_lumo_energy is None

    def test_alpha_beta_lumo_energy_populated(
        self, make_frontier_orbital_file
    ):
        dummy = make_frontier_orbital_file(
            alpha_virtual_eigenvalues=[1.0, 2.0],
            beta_virtual_eigenvalues=[1.5, 2.5],
        )
        assert dummy.alpha_lumo_energy == 1.0
        assert dummy.beta_lumo_energy == 1.5

    def test_homo_lumo_energy_none_when_no_eigenvalues(
        self, make_frontier_orbital_file
    ):
        dummy = make_frontier_orbital_file(multiplicity=1)
        assert dummy.homo_energy is None
        assert dummy.lumo_energy is None

    def test_homo_lumo_energy_none_for_open_shell(
        self, make_frontier_orbital_file
    ):
        dummy = make_frontier_orbital_file(
            multiplicity=3,
            alpha_occ_eigenvalues=[-5.0],
            alpha_virtual_eigenvalues=[1.0],
        )
        assert dummy.homo_energy is None
        assert dummy.lumo_energy is None

    def test_homo_lumo_energy_closed_shell(self, make_frontier_orbital_file):
        dummy = make_frontier_orbital_file(
            multiplicity=1,
            alpha_occ_eigenvalues=[-10.0, -5.0],
            alpha_virtual_eigenvalues=[1.0, 2.0],
        )
        assert dummy.homo_energy == -5.0
        assert dummy.lumo_energy == 1.0

    def test_fmo_gap_none_when_missing_data(self, make_frontier_orbital_file):
        dummy = make_frontier_orbital_file(multiplicity=1)
        assert dummy.fmo_gap is None

        dummy2 = make_frontier_orbital_file(multiplicity=3)
        assert dummy2.fmo_gap is None

    def test_fmo_gap_closed_shell(self, make_frontier_orbital_file):
        dummy = make_frontier_orbital_file(
            multiplicity=1,
            alpha_occ_eigenvalues=[-10.0, -5.0],
            alpha_virtual_eigenvalues=[1.0, 2.0],
        )
        assert dummy.fmo_gap == 6.0

    def test_fmo_gap_open_shell(self, make_frontier_orbital_file):
        dummy = make_frontier_orbital_file(
            multiplicity=3,
            alpha_occ_eigenvalues=[-10.0, -5.0, -1.0],
            alpha_virtual_eigenvalues=[2.0],
            beta_virtual_eigenvalues=[3.0],
        )
        assert dummy.fmo_gap == 3.0  # min(2.0, 3.0) - (-1.0)

    def test_alpha_beta_fmo_gap_none_without_data(
        self, make_frontier_orbital_file
    ):
        dummy = make_frontier_orbital_file()
        assert dummy.alpha_fmo_gap is None
        assert dummy.beta_fmo_gap is None

    def test_alpha_beta_fmo_gap_populated(self, make_frontier_orbital_file):
        dummy = make_frontier_orbital_file(
            alpha_occ_eigenvalues=[-10.0, -5.0],
            alpha_virtual_eigenvalues=[1.0],
            beta_occ_eigenvalues=[-9.0, -4.0],
            beta_virtual_eigenvalues=[1.5],
        )
        assert dummy.alpha_fmo_gap == 6.0
        assert dummy.beta_fmo_gap == 5.5


class TestFrequencyValidation:
    def test_opt_job_no_imaginary_freqs_is_valid(self):
        dummy = make_dummy_file("test.log")
        dummy.vibrational_frequencies = [100.0, 200.0]
        dummy.jobtype = "opt"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "OPT"
        assert report["total_imaginary_frequencies"] == 0
        assert report["is_valid_minimum"] is True
        assert report["is_valid_ts"] is False

    def test_opt_job_with_imaginary_freqs_is_invalid(self):
        dummy = make_dummy_file("test.log")
        dummy.vibrational_frequencies = [-50.0, 200.0]
        dummy.jobtype = "opt"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "OPT"
        assert report["total_imaginary_frequencies"] == 1
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is False

    def test_ts_job_one_imaginary_freq_is_valid(self):
        dummy = make_dummy_file("test.log")
        dummy.vibrational_frequencies = [-150.0, 200.0, 300.0]
        dummy.jobtype = "ts"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "TS"
        assert report["total_imaginary_frequencies"] == 1
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is True

    def test_ts_job_no_imaginary_freqs_is_invalid(self):
        dummy = make_dummy_file("test.log")
        dummy.vibrational_frequencies = [100.0, 200.0]
        dummy.jobtype = "ts"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "TS"
        assert report["total_imaginary_frequencies"] == 0
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is False

    def test_ts_job_multiple_imaginary_freqs_is_invalid(self):
        dummy = make_dummy_file("test.log")
        dummy.vibrational_frequencies = [-100.0, -200.0]
        dummy.jobtype = "ts"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "TS"
        assert report["total_imaginary_frequencies"] == 2
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is False

    def test_no_frequencies_returns_valid_for_opt(self):
        dummy = make_dummy_file("test.log")
        dummy.vibrational_frequencies = None
        dummy.jobtype = "opt"
        report = dummy.validate_frequencies()
        assert report["is_valid_minimum"] is True
        assert report["total_imaginary_frequencies"] == 0

    def test_ignore_threshold(self):
        dummy = make_dummy_file("test.log")
        dummy.vibrational_frequencies = [-10.0, -20.0]
        dummy.jobtype = "ts"
        report = dummy.validate_frequencies(ignore_threshold=-15.0)
        assert report["total_imaginary_frequencies"] == 1
        assert report["is_valid_ts"] is True

    def test_unknown_job_type(self):
        dummy = make_dummy_file("test.log")
        dummy.vibrational_frequencies = [-50.0]
        dummy.jobtype = "sp"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "SP"
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is False

    def test_none_job_type(self):
        dummy = make_dummy_file("test.log")
        dummy.vibrational_frequencies = [-50.0]
        dummy.jobtype = None
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "UNKNOWN"
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is False


class TestGaussianFileMixin:
    def test_gaussian_file_properties(self):
        dummy = make_dummy_gaussian_file()
        assert dummy.chk is True
        assert dummy._get_mem() == 32
        assert dummy._get_nproc() == 8

    def test_route_forwarding_properties(self):
        dummy = make_dummy_gaussian_file(
            route_string=(
                "#p opt freq=numer force b3lyp/6-31g(d) "
                "scrf=(smd,solvent=water)"
            )
        )
        assert dummy.dieze_tag == "#p"
        assert dummy.numfreq is True
        assert dummy.method == "b3lyp"
        assert dummy.force is True
        assert dummy.solvent_on is True
        assert dummy.solvent_model == "smd"
        assert dummy.solvent_id == "water"
        assert dummy.additional_solvent_options is None

    def test_route_object_returns_none_on_type_error(self, capsys):
        dummy = make_dummy_gaussian_file(route_string=None)
        assert dummy.route_object is None
        assert "TypeError" not in capsys.readouterr().out

    def test_get_route_raises_not_implemented(self):
        bare = make_mixin_instance(GaussianFileMixin, filename="bare.gjf")
        with pytest.raises(NotImplementedError):
            bare._get_route()

    def test_get_chk_false_when_no_directive(self):
        no_chk = make_mixin_instance(
            GaussianFileMixin,
            contents=["%mem=8GB", "#p opt"],
            filename="nochk.gjf",
        )
        assert no_chk._get_chk() is False

    def test_get_version_found(self):
        versioned = make_mixin_instance(
            GaussianFileMixin,
            contents=[
                "some header",
                "******************************************",
                "Gaussian 16:  ES64L-G16RevB.01 20-Dec-2017",
                "more text",
            ],
            filename="v.log",
        )
        assert versioned._get_version() == "G16RevB.01"

    def test_get_version_not_found(self):
        no_version = make_mixin_instance(
            GaussianFileMixin,
            contents=["some header", "no marker here"],
            filename="v.log",
        )
        assert no_version._get_version() is None

    def test_get_version_marker_without_gaussian_line_keeps_looking(self):
        """Covers the `if "Gaussian" in next_line:` False arm: a
        "****" marker line followed by an unrelated line must not
        return, and should keep scanning for a later real match."""
        multi_marker = make_mixin_instance(
            GaussianFileMixin,
            contents=[
                "******************************************",
                "not a version line",
                "******************************************",
                "Gaussian 16:  ES64L-G16RevB.01 20-Dec-2017",
            ],
            filename="v.log",
        )
        assert multi_marker._get_version() == "G16RevB.01"

    def test_file_date_empty_contents_returns_none(self):
        empty = make_mixin_instance(
            GaussianFileMixin, contents=[], filename="empty.log"
        )
        assert empty.file_date is None

    def test_file_date_malformed_date_returns_none(self):
        malformed_date = make_mixin_instance(
            GaussianFileMixin,
            contents=[
                " Job cpu time: 0 days  0 hours  0 minutes  0.0 seconds."
            ],
            filename="bad.log",
        )
        assert malformed_date.file_date is None

    def test_file_date_matched_but_strptime_fails(self):
        bad_date = make_mixin_instance(
            GaussianFileMixin,
            contents=[" Normal termination of Gaussian at nonsense-date."],
            filename="bd.log",
        )
        assert bad_date.file_date is None

    def test_read_settings_builds_gaussian_job_settings(self):
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        dummy = make_dummy_gaussian_file(
            route_string="#p opt freq b3lyp/6-31g(d)",
            filename="mytest.gjf",
            charge=0,
            multiplicity=1,
            heavy_elements=None,
            heavy_elements_basis=None,
            light_elements_basis=None,
            custom_solvent=None,
        )

        settings = dummy.read_settings()

        assert isinstance(settings, GaussianJobSettings)
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.functional == "b3lyp"
        assert "mytest.gjf" in settings.title

    def test_modred_scan_coords_skips_non_scan_lines(self):
        """The dummy Gaussian file's modredundant_group mixes a frozen
        ("F") line with a scan ("S") line; _get_modred_scan_coords must
        skip the non-scan line (continue) rather than trying to parse it
        as a scan spec.

        Note: the `self.jobtype = "scan"` assignment inside
        _get_modredundant_conditions doesn't actually persist, since
        `jobtype`'s setter writes through the non-cached `route_object`
        property (a fresh GaussianRoute is built on every access) --
        see BUGS_FOUND.md #61. This test only asserts the returned
        modred dict, not the jobtype mutation."""
        dummy = make_dummy_gaussian_file()
        modred = dummy.modred
        assert modred["coords"] == [[9, 2]]
        assert modred["num_steps"] == 10
        assert modred["step_size"] == 0.05


class TestORCAFileMixin:
    def test_orca_file_properties(self, dummy_orca_file):
        assert dummy_orca_file.mdci_cutoff == "1e-5"
        assert dummy_orca_file.mdci_density == "1e-6"
        assert dummy_orca_file.solvent_model == "smd"
        assert dummy_orca_file.solvent_id == "water"

    def test_solvent_on_true_when_model_and_id_present(self, dummy_orca_file):
        assert dummy_orca_file.solvent_on is True

    def test_solvent_on_false_without_solvent_block(self):
        no_solvent = make_mixin_instance(
            ORCAFileMixin,
            contents=["! B3LYP def2-SVP"],
            route_string="! B3LYP def2-SVP",
            filename="nosolv.inp",
        )
        assert no_solvent.solvent_on is False
        assert no_solvent.solvent_model is None
        assert no_solvent.solvent_id is None

    def test_solvent_model_cpcm_without_smd(self):
        cpcm_only = make_mixin_instance(
            ORCAFileMixin,
            contents=["%cpcm", "  epsilon 80.4", "%end"],
            filename="cpcm.inp",
        )
        assert cpcm_only.solvent_model == "cpcm"

    def test_contents_string_joins_lines(self, dummy_orca_file):
        assert dummy_orca_file.contents_string == "\n".join(
            dummy_orca_file.contents
        )

    def test_mdci_cutoff_and_density_none_without_mdci_block(self):
        no_mdci = make_mixin_instance(
            ORCAFileMixin,
            contents=["! B3LYP def2-SVP"],
            filename="no_mdci.inp",
        )
        assert no_mdci.mdci_cutoff is None
        assert no_mdci.mdci_density is None

    def test_mdci_cutoff_absent_within_block_returns_none(self):
        no_cutoff = make_mixin_instance(
            ORCAFileMixin,
            contents=["%mdci", "  density 1e-6", "%end"],
            filename="nc.inp",
        )
        assert no_cutoff.mdci_cutoff is None

    def test_mdci_cutoff_skips_non_cutoff_lines_before_match(self):
        delayed_cutoff = make_mixin_instance(
            ORCAFileMixin,
            contents=["%mdci", "  density 1e-6", "  cutoff 1e-5", "%end"],
            filename="dc.inp",
        )
        assert delayed_cutoff.mdci_cutoff == "1e-5"

    def test_mdci_density_absent_within_block_returns_none(self):
        no_density = make_mixin_instance(
            ORCAFileMixin,
            contents=["%mdci", "  cutoff 1e-5", "%end"],
            filename="nd2.inp",
        )
        assert no_density.mdci_density is None

    def test_get_version_not_found(self):
        no_version = make_mixin_instance(
            ORCAFileMixin,
            contents=["no version marker here"],
            filename="v.out",
        )
        assert no_version._get_version() is None

    def test_file_date_parses_starting_time(self):
        dated = make_mixin_instance(
            ORCAFileMixin,
            contents=["* Starting time: Mon Jan  5 12:34:56 2026"],
            filename="d.out",
        )
        assert dated.file_date == "2026-01-05 12:34:56"

    def test_file_date_regex_matches_but_strptime_fails(self):
        """A line matching orca_date_pattern's shape but with an
        invalid weekday/month combination that datetime.strptime
        still rejects should be treated as "continue looking" rather
        than raising."""
        bad_date = make_mixin_instance(
            ORCAFileMixin,
            contents=["* Starting time: Xxx Jan 99 99:99:99 2026"],
            filename="bd.out",
        )
        assert bad_date.file_date is None

    def test_file_date_no_starting_time_returns_none(self):
        no_date = make_mixin_instance(
            ORCAFileMixin,
            contents=["nothing relevant here"],
            filename="nd.out",
        )
        assert no_date.file_date is None

    def test_file_date_starting_time_substring_without_full_pattern(self):
        """ "Starting time:" is present, but without the leading "* "
        marker orca_date_pattern requires, so the regex search itself
        must fail (as opposed to strptime failing on a match)."""
        unmatched_marker = make_mixin_instance(
            ORCAFileMixin,
            contents=["Starting time: not actually parseable"],
            filename="um.out",
        )
        assert unmatched_marker.file_date is None

    @pytest.mark.parametrize(
        "keyword, expected",
        [
            ("COSMORS(water)", "cosmors"),
            ("CPCMC", "cpcmc"),
            ("SMD(water)", "smd"),
            ("CPCM", "cpcm"),
        ],
    )
    def test_solvent_model_route_fallback_variants(self, keyword, expected):
        route = f"! B3LYP def2-SVP {keyword}"
        route_file = make_mixin_instance(
            ORCAFileMixin,
            contents=[route],
            route_string=route,
            filename="r.inp",
        )
        assert route_file.solvent_model == expected

    def test_solvent_model_route_string_not_implemented_is_swallowed(self):
        no_route_string = make_mixin_instance(
            ORCAFileMixin,
            contents=["! B3LYP def2-SVP"],
            filename="nr.inp",
        )
        assert no_route_string.solvent_model is None

    def test_solvent_id_from_solvent_name_line(self):
        solvent_name = make_mixin_instance(
            ORCAFileMixin,
            contents=["Solvent name: water"],
            filename="sn.out",
        )
        assert solvent_name.solvent_id == "water"

    def test_solvent_id_raises_when_not_quoted(self):
        unquoted_solvent = make_mixin_instance(
            ORCAFileMixin,
            contents=["  solvent water"],
            filename="us.inp",
        )
        with pytest.raises(Exception, match="not in quotes"):
            unquoted_solvent.solvent_id

    def test_solvent_id_route_fallback_and_not_implemented(self):
        route_solvent = make_mixin_instance(
            ORCAFileMixin,
            contents=["! B3LYP def2-SVP SMD(cyclohexane)"],
            route_string="! B3LYP def2-SVP SMD(cyclohexane)",
            filename="rs.inp",
        )
        assert route_solvent.solvent_id == "cyclohexane"

        no_route_string = make_mixin_instance(
            ORCAFileMixin,
            contents=["! B3LYP def2-SVP"],
            filename="nrs.inp",
        )
        assert no_route_string.solvent_id is None

    def test_route_string_not_overridden_raises(self):
        bare = make_mixin_instance(ORCAFileMixin, filename="bare.inp")
        with pytest.raises(NotImplementedError):
            bare.route_string

    def test_method_dispersion_scf_algorithm_forwarding(self):
        route = "! B3LYP D3BJ def2-SVP DIIS"
        dummy = make_mixin_instance(
            ORCAFileMixin,
            contents=[route],
            route_string=route,
            filename="rich.inp",
        )
        assert dummy.method == "b3lyp"
        assert dummy.dispersion == "d3bj"
        assert dummy.scf_algorithm == "diis"

    def test_read_settings_builds_orca_job_settings(self):
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        route = "! B3LYP def2-SVP"
        settings_file = make_mixin_instance(
            ORCAFileMixin,
            contents=[route],
            route_string=route,
            filename="settings.inp",
            charge=0,
            multiplicity=1,
            scf_maxiter=None,
            scf_convergence=None,
            dipole=False,
            quadrupole=False,
        )

        settings = settings_file.read_settings()
        assert isinstance(settings, ORCAJobSettings)
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.functional == "b3lyp"


class TestXTBFileMixin:
    def test_xtb_file_properties(self, dummy_xtb_file):
        dummy = dummy_xtb_file
        assert dummy.jobtype == "opt"
        assert dummy.optimization_level == "loose"
        assert dummy.gfn_version == "gfn2"
        assert dummy.solvent_model == "alpb"
        assert dummy.solvent_id == "toluene"
        assert dummy.charge == 0
        assert dummy.uhf == 0
        assert dummy.freq is False
        assert dummy.grad is True

    def test_get_route_raises_not_implemented(self):
        bare = make_mixin_instance(XTBFileMixin, filename="bare.out")
        with pytest.raises(NotImplementedError):
            bare._get_route()

    def test_custom_solvent_always_none(self, dummy_xtb_file):
        assert dummy_xtb_file.custom_solvent is None

    @pytest.mark.parametrize(
        "lines, expected",
        [
            (["  * xtb version 6.5.1"], "6.5.1"),
            (["  * xtb version"], None),
            (["no marker"], None),
            # "xtb version" substring matches, but whitespace-splitting
            # doesn't produce an exact "version" token (it's part of a
            # larger word), so the inner `if "version" in parts:` is False.
            (["xtb versioning 1.0"], None),
        ],
    )
    def test_get_version_found_and_not_found(self, lines, expected):
        versioned = make_mixin_instance(
            XTBFileMixin, contents=lines, filename="v.out"
        )
        assert versioned._get_version() == expected

    def test_file_date_parses_finished_run_line(self):
        dated = make_mixin_instance(
            XTBFileMixin,
            contents=["* finished run on 2026/01/05 at 12:34:56"],
            filename="d.out",
        )
        assert dated.file_date == "2026-01-05 12:34:56"

    def test_file_date_no_marker_returns_none(self):
        no_date = make_mixin_instance(
            XTBFileMixin, contents=["nothing relevant"], filename="nd.out"
        )
        assert no_date.file_date is None

    def test_file_date_marker_present_but_unmatched_pattern(self):
        unmatched_date = make_mixin_instance(
            XTBFileMixin,
            contents=["* finished run on today at noon"],
            filename="ud.out",
        )
        assert unmatched_date.file_date is None

    def test_file_date_matched_but_strptime_fails(self):
        bad_date = make_mixin_instance(
            XTBFileMixin,
            contents=["* finished run on 2026/99/99 at 99:99:99"],
            filename="bd.out",
        )
        assert bad_date.file_date is None

    def test_method_returns_gfn_directly_when_present_in_route(self):
        gfn_route = make_mixin_instance(
            XTBFileMixin,
            route_string="xtb structure.xyz --opt --gfn 2",
            filename="g.out",
        )
        assert gfn_route.method == "gfn2"

    def test_method_falls_back_to_hamiltonian(self):
        hamiltonian = make_mixin_instance(
            XTBFileMixin,
            route_string="xtb structure.xyz --opt",
            filename="h.out",
            hamiltonian="GFN2-xTB",
        )
        assert hamiltonian.method == "gfn2"

    def test_method_none_without_gfn_or_hamiltonian(self):
        no_method = make_mixin_instance(
            XTBFileMixin,
            route_string="xtb structure.xyz --opt",
            filename="nm.out",
            hamiltonian=None,
        )
        assert no_method.method is None


class TestCRESTFileMixin:
    def test_crest_file_properties(self, dummy_crest_file):
        dummy = dummy_crest_file
        assert dummy.jobtype == "conformers"
        assert dummy.optimization_level == "tight"
        assert dummy.gfn_version == "gfn2"
        assert dummy.solvent_model is None
        assert dummy.solvent_id is None
        assert dummy.charge == 0
        assert dummy.uhf == 0
        assert dummy.constrained is True


class TestYAMLFileMixin:
    def test_yaml_file_properties(self, dummy_yaml_file):
        dummy = dummy_yaml_file
        assert dummy.yaml_contents_dict == {"key1": "value1", "key2": "value2"}
        assert "key1" in dummy.yaml_contents_keys
        assert "value1" in dummy.yaml_contents_values
        assert dummy.yaml_contents_by_key("key1") == "value1"


class TestRegistryMixin:
    def test_subclasses(self):
        subclasses = BaseRegistry.subclasses()
        assert SubRegistry1 in subclasses
        assert SubRegistry2 in subclasses

    def test_non_registerable_subclass_is_not_registered(
        self,
    ):
        # Defining the subclass is the act under test: RegistryMixin's
        # __init_subclass__ must skip it.
        class NonRegisterableSub(BaseRegistry):
            REGISTERABLE = False

        assert NonRegisterableSub not in BaseRegistry._REGISTRY


class TestFolderMixin:
    def test_get_all_files_by_suffix(self, temp_folder_with_files):
        folder, file1, file2 = temp_folder_with_files
        dummy = make_dummy_folder(folder)
        txt_files = dummy.get_all_files_in_current_folder_by_suffix(".txt")
        assert file1 in txt_files

    def test_suffix_without_dot_excludes_compound_extensions(self, tmp_path):
        """filetype='xyz' must match .xyz but not .extxyz."""
        content = "3\n\nO 0 0 0\nH 1 0 0\nH 0 1 0\n"
        xyz = str(write_text_file(tmp_path / "water.xyz", content))
        extxyz = str(write_text_file(tmp_path / "crystal.extxyz", content))

        dummy = make_dummy_folder(str(tmp_path))
        files = dummy.get_all_files_in_current_folder_by_suffix("xyz")
        assert xyz in files
        assert extxyz not in files

        # Leading-dot form remains supported.
        files_dotted = dummy.get_all_files_in_current_folder_by_suffix(".xyz")
        assert xyz in files_dotted
        assert extxyz not in files_dotted

    def test_get_all_files_by_regex(self, temp_folder_with_files):
        folder, file1, file2 = temp_folder_with_files
        dummy = make_dummy_folder(folder)
        log_files = dummy.get_all_files_in_current_folder_matching_regex(
            r".*\.log"
        )
        assert file2 in log_files
        assert file1 not in log_files

    def test_folderpath_is_absolute(self, tmp_path):
        dummy = make_dummy_folder(str(tmp_path))
        assert dummy.folderpath == os.path.abspath(str(tmp_path))

    def test_get_all_output_files_unsupported_program_raises(self, tmp_path):
        dummy = make_dummy_folder(str(tmp_path))
        with pytest.raises(ValueError, match="Unsupported program"):
            dummy.get_all_output_files_in_current_folder_by_program(
                "bogus_program"
            )

    def test_get_all_output_files_by_program_non_recursive(self, tmp_path):
        gaussian_log = str(
            write_text_file(tmp_path / "job.log", "Entering Gaussian System\n")
        )
        write_text_file(tmp_path / "empty.log", "")
        write_text_file(
            tmp_path / "subdir" / "nested.log", "Entering Gaussian System\n"
        )
        write_text_file(tmp_path / "notes.txt", "just some text\n")
        write_text_file(tmp_path / "unknown.log", "nothing recognizable\n")

        dummy = make_dummy_folder(str(tmp_path))
        gaussian_files = (
            dummy.get_all_output_files_in_current_folder_by_program("gaussian")
        )
        assert gaussian_files == [gaussian_log]

        all_files = dummy.get_all_output_files_in_current_folder_by_program()
        assert gaussian_log in all_files
        assert all(f.endswith(".log") or f.endswith(".out") for f in all_files)

    def test_get_all_output_files_by_program_recursive(self, tmp_path):
        nested_log = str(
            write_text_file(
                tmp_path / "sub" / "nested.log", "Entering Gaussian System\n"
            )
        )
        top_log = str(
            write_text_file(tmp_path / "top.log", "Entering Gaussian System\n")
        )
        # An empty file and a non-matching-suffix file, both within the
        # recursive walk, must be skipped without raising.
        write_text_file(tmp_path / "sub" / "empty.log", "")
        write_text_file(tmp_path / "sub" / "notes.txt", "irrelevant\n")
        # A broken symlink appears in os.walk's "files" list but fails
        # os.path.isfile, covering that skip branch too.
        os.symlink(
            tmp_path / "sub" / "does_not_exist.log",
            tmp_path / "sub" / "broken_link.log",
        )

        dummy = make_dummy_folder(str(tmp_path))
        found = dummy.get_all_output_files_in_current_folder_and_subfolders_by_program(
            "gaussian"
        )
        assert set(found) == {nested_log, top_log}

    def test_get_all_output_files_recursive_program_none(self, tmp_path):
        nested_log = str(
            write_text_file(
                tmp_path / "sub" / "nested.log", "Entering Gaussian System\n"
            )
        )
        write_text_file(tmp_path / "sub" / "notes.dat", "irrelevant\n")

        dummy = make_dummy_folder(str(tmp_path))
        found = (
            dummy.get_all_output_files_in_current_folder_and_subfolders_by_program()
        )
        assert found == [nested_log]

    def test_is_program_calculation_directory(self, tmp_path):
        write_text_file(tmp_path / "job.log", "Entering Gaussian System\n")
        dummy = make_dummy_folder(str(tmp_path))
        assert dummy.is_program_calculation_directory("gaussian") is True
        assert dummy.is_program_calculation_directory("orca") is False

    def test_is_program_calculation_directory_missing_folder(self, tmp_path):
        dummy = make_dummy_folder(str(tmp_path / "does_not_exist"))
        assert dummy.is_program_calculation_directory("gaussian") is False

    def test_get_program_type_from_folder_unknown(self, tmp_path):
        write_text_file(tmp_path / "notes.txt", "nothing relevant\n")
        dummy = make_dummy_folder(str(tmp_path))
        assert dummy.get_program_type_from_folder() == "unknown"

    def test_get_program_type_from_folder_single_program(self, tmp_path):
        write_text_file(tmp_path / "job.out", "x T B\nsome xtb output\n")
        dummy = make_dummy_folder(str(tmp_path))
        assert dummy.get_program_type_from_folder() == "xtb"

    def test_get_program_type_from_folder_mixed(self, tmp_path):
        write_text_file(tmp_path / "xtb_job.out", "x T B\nsome xtb output\n")
        write_text_file(
            tmp_path / "crest_job.out",
            "C R E S T\nsome crest output\n",
        )
        dummy = make_dummy_folder(str(tmp_path))
        assert dummy.get_program_type_from_folder() == "mixed"

    def test_get_all_files_by_suffix_skips_empty_files(self, tmp_path):
        kept = str(write_text_file(tmp_path / "kept.txt", "content"))
        write_text_file(tmp_path / "empty.txt", "")
        dummy = make_dummy_folder(str(tmp_path))
        assert dummy.get_all_files_in_current_folder_by_suffix(".txt") == [
            kept
        ]

    def test_get_all_files_by_regex_skips_empty_files(self, tmp_path):
        kept = str(write_text_file(tmp_path / "kept.log", "content"))
        write_text_file(tmp_path / "empty.log", "")
        dummy = make_dummy_folder(str(tmp_path))
        matched = dummy.get_all_files_in_current_folder_matching_regex(
            r".*\.log"
        )
        assert matched == [kept]

    def test_get_all_files_and_subfolders_by_suffix(self, tmp_path):
        top = str(write_text_file(tmp_path / "top.log", "content"))
        nested = str(
            write_text_file(tmp_path / "sub" / "nested.log", "content")
        )
        write_text_file(tmp_path / "other.txt", "content")

        dummy = make_dummy_folder(str(tmp_path))
        found = dummy.get_all_files_in_current_folder_and_subfolders_by_suffix(
            ".log"
        )
        assert set(found) == {top, nested}

    def test_get_all_files_by_program_and_suffix_non_recursive(self, tmp_path):
        gaussian_log = str(
            write_text_file(tmp_path / "job.log", "Entering Gaussian System\n")
        )
        write_text_file(tmp_path / "empty.log", "")
        write_text_file(tmp_path / "unknown.log", "nothing recognizable\n")
        write_text_file(
            tmp_path / "job.dat", "Entering Gaussian System\n"
        )  # wrong suffix
        (tmp_path / "subdir").mkdir()

        dummy = make_dummy_folder(str(tmp_path))
        found = dummy.get_all_files_in_current_folder_by_program_and_suffix(
            "gaussian", ".log"
        )
        assert found == [gaussian_log]

    def test_get_all_files_by_program_and_suffix_recursive(self, tmp_path):
        nested = str(
            write_text_file(
                tmp_path / "sub" / "nested.log", "Entering Gaussian System\n"
            )
        )
        top = str(
            write_text_file(tmp_path / "top.log", "Entering Gaussian System\n")
        )
        write_text_file(
            tmp_path / "sub" / "wrong_suffix.dat",
            "Entering Gaussian System\n",
        )
        write_text_file(
            tmp_path / "sub" / "wrong_program.log", "nothing recognizable\n"
        )

        dummy = make_dummy_folder(str(tmp_path))
        found = dummy.get_all_files_in_current_folder_and_subfolders_by_program_and_suffix(
            "gaussian", ".log"
        )
        assert set(found) == {nested, top}

    def test_get_all_files_and_subfolders_matching_regex(self, tmp_path):
        top = str(write_text_file(tmp_path / "top.log", "content"))
        nested = str(
            write_text_file(tmp_path / "sub" / "nested.log", "content")
        )
        write_text_file(tmp_path / "other.txt", "content")

        dummy = make_dummy_folder(str(tmp_path))
        found = dummy.get_all_files_in_current_folder_and_subfolders_matching_regex(
            r".*\.log"
        )
        assert set(found) == {top, nested}


class TestFolderOutputMixin:
    def test_delegates_to_first_parser(self):
        output = make_dummy_folder_output(
            main_out=SimpleNamespace(total_energy=-123.45),
            secondary_out=SimpleNamespace(total_energy=-123),
        )
        assert output.total_energy == -123.45

    def test_falls_back_to_next_parser(self):
        output = make_dummy_folder_output(
            main_out=SimpleNamespace(),
            secondary_out=SimpleNamespace(charge=0),
        )
        assert output.charge == 0

    def test_skips_none_parser(self):
        output = make_dummy_folder_output(
            main_out=None,
            secondary_out=SimpleNamespace(multiplicity=1),
        )
        assert output.multiplicity == 1

    def test_raises_attribute_error_when_not_found(self):
        output = make_dummy_folder_output()
        with pytest.raises(AttributeError, match="no attribute 'missing'"):
            _ = output.missing
