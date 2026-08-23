import os

import pytest

from chemsmart.utils.mixins import (
    CRESTFileMixin,
    FileMixin,
    FolderMixin,
    FolderOutputMixin,
    GaussianFileMixin,
    ORCAFileMixin,
    RegistryMixin,
    XTBFileMixin,
)


class DummyFile(FileMixin):
    def __init__(self, filename):
        self.filename = filename
        self.forces = None
        self.energies = [1.0, 2.0]
        self.input_coordinates_block = None


class TestFileMixin:
    def test_file_properties(self, temp_text_file):
        dummy = DummyFile(temp_text_file)
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
        class DummyBlock:
            translation_vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]

        dummy = DummyFile(temp_text_file)
        dummy.input_coordinates_block = DummyBlock()
        assert dummy.input_translation_vectors == [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]

    def test_input_translation_vectors_block_without_vectors(
        self, temp_text_file
    ):
        class DummyBlockNoVectors:
            translation_vectors = None

        dummy = DummyFile(temp_text_file)
        dummy.input_coordinates_block = DummyBlockNoVectors()
        assert dummy.input_translation_vectors == []


class TestFrontierOrbitalProperties:
    """Direct tests for FileMixin's frontier-orbital cached_properties,
    isolating the multiplicity==1 vs. open-shell, and populated vs.
    empty eigenvalue-list branches independently of any real output
    parser."""

    def _make(self, temp_text_file, **overrides):
        dummy = DummyFile(temp_text_file)
        dummy.multiplicity = overrides.pop("multiplicity", 1)
        dummy.alpha_occ_eigenvalues = overrides.pop(
            "alpha_occ_eigenvalues", None
        )
        dummy.beta_occ_eigenvalues = overrides.pop(
            "beta_occ_eigenvalues", None
        )
        dummy.alpha_virtual_eigenvalues = overrides.pop(
            "alpha_virtual_eigenvalues", None
        )
        dummy.beta_virtual_eigenvalues = overrides.pop(
            "beta_virtual_eigenvalues", None
        )
        return dummy

    def test_somo_energies_none_for_closed_shell(self, temp_text_file):
        dummy = self._make(temp_text_file, multiplicity=1)
        assert dummy.somo_energies is None

    def test_somo_energies_none_when_no_alpha_occ(self, temp_text_file):
        dummy = self._make(
            temp_text_file, multiplicity=3, alpha_occ_eigenvalues=None
        )
        assert dummy.somo_energies is None

    def test_somo_energies_open_shell(self, temp_text_file):
        dummy = self._make(
            temp_text_file,
            multiplicity=3,
            alpha_occ_eigenvalues=[-10.0, -5.0, -1.0],
        )
        assert dummy.somo_energies == [-5.0, -1.0]

    def test_lowest_and_highest_somo_energy(self, temp_text_file):
        dummy = self._make(
            temp_text_file,
            multiplicity=3,
            alpha_occ_eigenvalues=[-10.0, -5.0, -1.0],
        )
        assert dummy.lowest_somo_energy == -5.0
        assert dummy.highest_somo_energy == -1.0

    def test_lowest_and_highest_somo_energy_none_for_closed_shell(
        self, temp_text_file
    ):
        dummy = self._make(temp_text_file, multiplicity=1)
        assert dummy.lowest_somo_energy is None
        assert dummy.highest_somo_energy is None

    def test_alpha_beta_homo_energy_none_without_eigenvalues(
        self, temp_text_file
    ):
        dummy = self._make(temp_text_file)
        assert dummy.alpha_homo_energy is None
        assert dummy.beta_homo_energy is None

    def test_alpha_beta_homo_energy_populated(self, temp_text_file):
        dummy = self._make(
            temp_text_file,
            alpha_occ_eigenvalues=[-10.0, -5.0],
            beta_occ_eigenvalues=[-9.0, -4.0],
        )
        assert dummy.alpha_homo_energy == -5.0
        assert dummy.beta_homo_energy == -4.0

    def test_alpha_beta_lumo_energy_none_without_eigenvalues(
        self, temp_text_file
    ):
        dummy = self._make(temp_text_file)
        assert dummy.alpha_lumo_energy is None
        assert dummy.beta_lumo_energy is None

    def test_alpha_beta_lumo_energy_populated(self, temp_text_file):
        dummy = self._make(
            temp_text_file,
            alpha_virtual_eigenvalues=[1.0, 2.0],
            beta_virtual_eigenvalues=[1.5, 2.5],
        )
        assert dummy.alpha_lumo_energy == 1.0
        assert dummy.beta_lumo_energy == 1.5

    def test_homo_lumo_energy_none_when_no_eigenvalues(self, temp_text_file):
        dummy = self._make(temp_text_file, multiplicity=1)
        assert dummy.homo_energy is None
        assert dummy.lumo_energy is None

    def test_homo_lumo_energy_none_for_open_shell(self, temp_text_file):
        dummy = self._make(
            temp_text_file,
            multiplicity=3,
            alpha_occ_eigenvalues=[-5.0],
            alpha_virtual_eigenvalues=[1.0],
        )
        assert dummy.homo_energy is None
        assert dummy.lumo_energy is None

    def test_homo_lumo_energy_closed_shell(self, temp_text_file):
        dummy = self._make(
            temp_text_file,
            multiplicity=1,
            alpha_occ_eigenvalues=[-10.0, -5.0],
            alpha_virtual_eigenvalues=[1.0, 2.0],
        )
        assert dummy.homo_energy == -5.0
        assert dummy.lumo_energy == 1.0

    def test_fmo_gap_none_when_missing_data(self, temp_text_file):
        dummy = self._make(temp_text_file, multiplicity=1)
        assert dummy.fmo_gap is None

        dummy2 = self._make(temp_text_file, multiplicity=3)
        assert dummy2.fmo_gap is None

    def test_fmo_gap_closed_shell(self, temp_text_file):
        dummy = self._make(
            temp_text_file,
            multiplicity=1,
            alpha_occ_eigenvalues=[-10.0, -5.0],
            alpha_virtual_eigenvalues=[1.0, 2.0],
        )
        assert dummy.fmo_gap == 6.0

    def test_fmo_gap_open_shell(self, temp_text_file):
        dummy = self._make(
            temp_text_file,
            multiplicity=3,
            alpha_occ_eigenvalues=[-10.0, -5.0, -1.0],
            alpha_virtual_eigenvalues=[2.0],
            beta_virtual_eigenvalues=[3.0],
        )
        assert dummy.fmo_gap == 3.0  # min(2.0, 3.0) - (-1.0)

    def test_alpha_beta_fmo_gap_none_without_data(self, temp_text_file):
        dummy = self._make(temp_text_file)
        assert dummy.alpha_fmo_gap is None
        assert dummy.beta_fmo_gap is None

    def test_alpha_beta_fmo_gap_populated(self, temp_text_file):
        dummy = self._make(
            temp_text_file,
            alpha_occ_eigenvalues=[-10.0, -5.0],
            alpha_virtual_eigenvalues=[1.0],
            beta_occ_eigenvalues=[-9.0, -4.0],
            beta_virtual_eigenvalues=[1.5],
        )
        assert dummy.alpha_fmo_gap == 6.0
        assert dummy.beta_fmo_gap == 5.5


class TestFrequencyValidation:
    def test_opt_job_no_imaginary_freqs_is_valid(self):
        dummy = DummyFile("test.log")
        dummy.vibrational_frequencies = [100.0, 200.0]
        dummy.jobtype = "opt"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "OPT"
        assert report["total_imaginary_frequencies"] == 0
        assert report["is_valid_minimum"] is True
        assert report["is_valid_ts"] is False

    def test_opt_job_with_imaginary_freqs_is_invalid(self):
        dummy = DummyFile("test.log")
        dummy.vibrational_frequencies = [-50.0, 200.0]
        dummy.jobtype = "opt"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "OPT"
        assert report["total_imaginary_frequencies"] == 1
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is False

    def test_ts_job_one_imaginary_freq_is_valid(self):
        dummy = DummyFile("test.log")
        dummy.vibrational_frequencies = [-150.0, 200.0, 300.0]
        dummy.jobtype = "ts"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "TS"
        assert report["total_imaginary_frequencies"] == 1
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is True

    def test_ts_job_no_imaginary_freqs_is_invalid(self):
        dummy = DummyFile("test.log")
        dummy.vibrational_frequencies = [100.0, 200.0]
        dummy.jobtype = "ts"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "TS"
        assert report["total_imaginary_frequencies"] == 0
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is False

    def test_ts_job_multiple_imaginary_freqs_is_invalid(self):
        dummy = DummyFile("test.log")
        dummy.vibrational_frequencies = [-100.0, -200.0]
        dummy.jobtype = "ts"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "TS"
        assert report["total_imaginary_frequencies"] == 2
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is False

    def test_no_frequencies_returns_valid_for_opt(self):
        dummy = DummyFile("test.log")
        dummy.vibrational_frequencies = None
        dummy.jobtype = "opt"
        report = dummy.validate_frequencies()
        assert report["is_valid_minimum"] is True
        assert report["total_imaginary_frequencies"] == 0

    def test_ignore_threshold(self):
        dummy = DummyFile("test.log")
        dummy.vibrational_frequencies = [-10.0, -20.0]
        dummy.jobtype = "ts"
        report = dummy.validate_frequencies(ignore_threshold=-15.0)
        assert report["total_imaginary_frequencies"] == 1
        assert report["is_valid_ts"] is True

    def test_unknown_job_type(self):
        dummy = DummyFile("test.log")
        dummy.vibrational_frequencies = [-50.0]
        dummy.jobtype = "sp"
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "SP"
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is False

    def test_none_job_type(self):
        dummy = DummyFile("test.log")
        dummy.vibrational_frequencies = [-50.0]
        dummy.jobtype = None
        report = dummy.validate_frequencies()
        assert report["detected_job_type"] == "UNKNOWN"
        assert report["is_valid_minimum"] is False
        assert report["is_valid_ts"] is False


class DummyGaussianFile(GaussianFileMixin):
    def __init__(self, filename):
        self.filename = filename
        self._route_string = "modred"

    @property
    def contents(self):
        return [
            "%chk=test.chk",
            "%mem=32GB",
            "%nproc=8",
            "#p opt freq",
            "modred",
        ]

    @property
    def route_string(self):
        return self._route_string

    @property
    def modredundant_group(self):
        return ["B 2 12 F", "B 9 2 S 10 0.05"]


class TestGaussianFileMixin:
    def test_gaussian_file_properties(self):
        dummy = DummyGaussianFile("test.gjf")
        assert dummy.chk is True
        assert dummy._get_mem() == 32
        assert dummy._get_nproc() == 8

    def test_route_forwarding_properties(self):
        dummy = DummyGaussianFile("test.gjf")
        dummy._route_string = (
            "#p opt freq=numer force b3lyp/6-31g(d) "
            "scrf=(smd,solvent=water)"
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
        dummy = DummyGaussianFile("test.gjf")
        dummy._route_string = None
        assert dummy.route_object is None
        assert "TypeError" not in capsys.readouterr().out

    def test_get_route_raises_not_implemented(self):
        import pytest

        class BareGaussianFile(GaussianFileMixin):
            def __init__(self):
                self.filename = "bare.gjf"

        bare = BareGaussianFile()
        with pytest.raises(NotImplementedError):
            bare._get_route()

    def test_get_chk_false_when_no_directive(self):
        class NoChkFile(GaussianFileMixin):
            def __init__(self):
                self.filename = "nochk.gjf"

            @property
            def contents(self):
                return ["%mem=8GB", "#p opt"]

        assert NoChkFile()._get_chk() is False

    def test_get_version_found(self):
        class VersionedFile(GaussianFileMixin):
            def __init__(self):
                self.filename = "v.log"

            @property
            def contents(self):
                return [
                    "some header",
                    "******************************************",
                    "Gaussian 16:  ES64L-G16RevB.01 20-Dec-2017",
                    "more text",
                ]

        assert VersionedFile()._get_version() == "G16RevB.01"

    def test_get_version_not_found(self):
        class NoVersionFile(GaussianFileMixin):
            def __init__(self):
                self.filename = "v.log"

            @property
            def contents(self):
                return ["some header", "no marker here"]

        assert NoVersionFile()._get_version() is None

    def test_get_version_marker_without_gaussian_line_keeps_looking(self):
        """Covers the `if "Gaussian" in next_line:` False arm: a
        "****" marker line followed by an unrelated line must not
        return, and should keep scanning for a later real match."""

        class MultiMarkerFile(GaussianFileMixin):
            def __init__(self):
                self.filename = "v.log"

            @property
            def contents(self):
                return [
                    "******************************************",
                    "not a version line",
                    "******************************************",
                    "Gaussian 16:  ES64L-G16RevB.01 20-Dec-2017",
                ]

        assert MultiMarkerFile()._get_version() == "G16RevB.01"

    def test_file_date_empty_contents_returns_none(self):
        class EmptyFile(GaussianFileMixin):
            def __init__(self):
                self.filename = "empty.log"

            @property
            def contents(self):
                return []

        assert EmptyFile().file_date is None

    def test_file_date_malformed_date_returns_none(self):
        class MalformedDateFile(GaussianFileMixin):
            def __init__(self):
                self.filename = "bad.log"

            @property
            def contents(self):
                return [
                    " Job cpu time: 0 days  0 hours  0 minutes  0.0 seconds."
                ]

        assert MalformedDateFile().file_date is None

    def test_file_date_matched_but_strptime_fails(self):
        class BadDateFile(GaussianFileMixin):
            def __init__(self):
                self.filename = "bd.log"

            @property
            def contents(self):
                return [" Normal termination of Gaussian at nonsense-date."]

        assert BadDateFile().file_date is None

    def test_read_settings_builds_gaussian_job_settings(self):
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        dummy = DummyGaussianFile("mytest.gjf")
        dummy._route_string = "#p opt freq b3lyp/6-31g(d)"
        dummy.charge = 0
        dummy.multiplicity = 1
        dummy.heavy_elements = None
        dummy.heavy_elements_basis = None
        dummy.light_elements_basis = None
        dummy.custom_solvent = None

        settings = dummy.read_settings()

        assert isinstance(settings, GaussianJobSettings)
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.functional == "b3lyp"
        assert "mytest.gjf" in settings.title

    def test_modred_scan_coords_skips_non_scan_lines(self):
        """DummyGaussianFile's modredundant_group mixes a frozen ("F")
        line with a scan ("S") line; _get_modred_scan_coords must skip
        the non-scan line (continue) rather than trying to parse it as
        a scan spec.

        Note: the `self.jobtype = "scan"` assignment inside
        _get_modredundant_conditions doesn't actually persist, since
        `jobtype`'s setter writes through the non-cached `route_object`
        property (a fresh GaussianRoute is built on every access) --
        see BUGS_FOUND.md #61. This test only asserts the returned
        modred dict, not the jobtype mutation."""
        dummy = DummyGaussianFile("test.gjf")
        modred = dummy.modred
        assert modred["coords"] == [[9, 2]]
        assert modred["num_steps"] == 10
        assert modred["step_size"] == 0.05


class DummyORCAFile(ORCAFileMixin):
    def __init__(self, filename):
        self.filename = filename

    @property
    def contents(self):
        return [
            "%mdci",
            "  cutoff 1e-5",
            "  density 1e-6",
            "%cpcm",
            "  smd true",
            '  solvent "water"',
        ]

    @property
    def route_string(self):
        return "! B3LYP def2-SVP"


class TestORCAFileMixin:
    def test_orca_file_properties(self):
        dummy = DummyORCAFile("test.inp")
        assert dummy.mdci_cutoff == "1e-5"
        assert dummy.mdci_density == "1e-6"
        assert dummy.solvent_model == "smd"
        assert dummy.solvent_id == "water"

    def test_solvent_on_true_when_model_and_id_present(self):
        dummy = DummyORCAFile("test.inp")
        assert dummy.solvent_on is True

    def test_solvent_on_false_without_solvent_block(self):
        class NoSolventFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "nosolv.inp"

            @property
            def contents(self):
                return ["! B3LYP def2-SVP"]

            @property
            def route_string(self):
                return "! B3LYP def2-SVP"

        assert NoSolventFile().solvent_on is False
        assert NoSolventFile().solvent_model is None
        assert NoSolventFile().solvent_id is None

    def test_solvent_model_cpcm_without_smd(self):
        class CpcmOnlyFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "cpcm.inp"

            @property
            def contents(self):
                return ["%cpcm", "  epsilon 80.4", "%end"]

        assert CpcmOnlyFile().solvent_model == "cpcm"

    def test_contents_string_joins_lines(self):
        dummy = DummyORCAFile("test.inp")
        assert dummy.contents_string == "\n".join(dummy.contents)

    def test_mdci_cutoff_and_density_none_without_mdci_block(self):
        class NoMdciFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "no_mdci.inp"

            @property
            def contents(self):
                return ["! B3LYP def2-SVP"]

        assert NoMdciFile().mdci_cutoff is None
        assert NoMdciFile().mdci_density is None

    def test_mdci_cutoff_absent_within_block_returns_none(self):
        class NoCutoffFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "nc.inp"

            @property
            def contents(self):
                return ["%mdci", "  density 1e-6", "%end"]

        assert NoCutoffFile().mdci_cutoff is None

    def test_mdci_cutoff_skips_non_cutoff_lines_before_match(self):
        class DelayedCutoffFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "dc.inp"

            @property
            def contents(self):
                return ["%mdci", "  density 1e-6", "  cutoff 1e-5", "%end"]

        assert DelayedCutoffFile().mdci_cutoff == "1e-5"

    def test_mdci_density_absent_within_block_returns_none(self):
        class NoDensityFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "nd2.inp"

            @property
            def contents(self):
                return ["%mdci", "  cutoff 1e-5", "%end"]

        assert NoDensityFile().mdci_density is None

    def test_get_version_not_found(self):
        class NoVersionFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "v.out"

            @property
            def contents(self):
                return ["no version marker here"]

        assert NoVersionFile()._get_version() is None

    def test_file_date_parses_starting_time(self):
        class DatedFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "d.out"

            @property
            def contents(self):
                return [
                    "* Starting time: Mon Jan  5 12:34:56 2026",
                ]

        assert DatedFile().file_date == "2026-01-05 12:34:56"

    def test_file_date_regex_matches_but_strptime_fails(self):
        """A line matching orca_date_pattern's shape but with an
        invalid weekday/month combination that datetime.strptime
        still rejects should be treated as "continue looking" rather
        than raising."""

        class BadDateFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "bd.out"

            @property
            def contents(self):
                return [
                    "* Starting time: Xxx Jan 99 99:99:99 2026",
                ]

        assert BadDateFile().file_date is None

    def test_file_date_no_starting_time_returns_none(self):
        class NoDateFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "nd.out"

            @property
            def contents(self):
                return ["nothing relevant here"]

        assert NoDateFile().file_date is None

    def test_file_date_starting_time_substring_without_full_pattern(self):
        """ "Starting time:" is present, but without the leading "* "
        marker orca_date_pattern requires, so the regex search itself
        must fail (as opposed to strptime failing on a match)."""

        class UnmatchedMarkerFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "um.out"

            @property
            def contents(self):
                return ["Starting time: not actually parseable"]

        assert UnmatchedMarkerFile().file_date is None

    def test_solvent_model_route_fallback_variants(self):
        for keyword, expected in [
            ("COSMORS(water)", "cosmors"),
            ("CPCMC", "cpcmc"),
            ("SMD(water)", "smd"),
            ("CPCM", "cpcm"),
        ]:

            class RouteFile(ORCAFileMixin):
                def __init__(self, route):
                    self.filename = "r.inp"
                    self._route = route

                @property
                def contents(self):
                    return [f"! B3LYP def2-SVP {self._route}"]

                @property
                def route_string(self):
                    return f"! B3LYP def2-SVP {self._route}"

            assert RouteFile(keyword).solvent_model == expected

    def test_solvent_model_route_string_not_implemented_is_swallowed(self):
        class NoRouteStringFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "nr.inp"

            @property
            def contents(self):
                return ["! B3LYP def2-SVP"]

        assert NoRouteStringFile().solvent_model is None

    def test_solvent_id_from_solvent_name_line(self):
        class SolventNameFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "sn.out"

            @property
            def contents(self):
                return ["Solvent name: water"]

        assert SolventNameFile().solvent_id == "water"

    def test_solvent_id_raises_when_not_quoted(self):
        import pytest

        class UnquotedSolventFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "us.inp"

            @property
            def contents(self):
                return ["  solvent water"]

        with pytest.raises(Exception, match="not in quotes"):
            UnquotedSolventFile().solvent_id

    def test_solvent_id_route_fallback_and_not_implemented(self):
        class RouteSolventFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "rs.inp"

            @property
            def contents(self):
                return ["! B3LYP def2-SVP SMD(cyclohexane)"]

            @property
            def route_string(self):
                return "! B3LYP def2-SVP SMD(cyclohexane)"

        assert RouteSolventFile().solvent_id == "cyclohexane"

        class NoRouteStringFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "nrs.inp"

            @property
            def contents(self):
                return ["! B3LYP def2-SVP"]

        assert NoRouteStringFile().solvent_id is None

    def test_route_string_not_overridden_raises(self):
        import pytest

        class BareOrcaFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "bare.inp"

        with pytest.raises(NotImplementedError):
            BareOrcaFile().route_string

    def test_method_dispersion_scf_algorithm_forwarding(self):
        class RichRouteFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "rich.inp"

            @property
            def contents(self):
                return [self.route_string]

            @property
            def route_string(self):
                return "! B3LYP D3BJ def2-SVP DIIS"

        dummy = RichRouteFile()
        assert dummy.method == "b3lyp"
        assert dummy.dispersion == "d3bj"
        assert dummy.scf_algorithm == "diis"

    def test_read_settings_builds_orca_job_settings(self):
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        class SettingsFile(ORCAFileMixin):
            def __init__(self):
                self.filename = "settings.inp"
                self.charge = 0
                self.multiplicity = 1
                self.scf_maxiter = None
                self.scf_convergence = None
                self.dipole = False
                self.quadrupole = False

            @property
            def contents(self):
                return [self.route_string]

            @property
            def route_string(self):
                return "! B3LYP def2-SVP"

        settings = SettingsFile().read_settings()
        assert isinstance(settings, ORCAJobSettings)
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert settings.functional == "b3lyp"


class DummyXTBFile(XTBFileMixin):
    def __init__(self, filename):
        self.filename = filename

    @property
    def route_string(self):
        return "xtb p_benzyne.xyz --opt loose --gfn 2 --alpb toluene --chrg 0 --uhf 0 --grad"


class TestXTBFileMixin:
    def test_xtb_file_properties(self):
        dummy = DummyXTBFile("test.out")
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
        import pytest

        class BareXTBFile(XTBFileMixin):
            def __init__(self):
                self.filename = "bare.out"

        with pytest.raises(NotImplementedError):
            BareXTBFile()._get_route()

    def test_custom_solvent_always_none(self):
        assert DummyXTBFile("test.out").custom_solvent is None

    def test_get_version_found_and_not_found(self):
        class VersionedXTBFile(XTBFileMixin):
            def __init__(self, lines):
                self.filename = "v.out"
                self._lines = lines

            @property
            def contents(self):
                return self._lines

        assert (
            VersionedXTBFile(["  * xtb version 6.5.1"])._get_version()
            == "6.5.1"
        )
        assert VersionedXTBFile(["  * xtb version"])._get_version() is None
        assert VersionedXTBFile(["no marker"])._get_version() is None
        # "xtb version" substring matches, but whitespace-splitting
        # doesn't produce an exact "version" token (it's part of a
        # larger word), so the inner `if "version" in parts:` is False.
        assert VersionedXTBFile(["xtb versioning 1.0"])._get_version() is None

    def test_file_date_parses_finished_run_line(self):
        class DatedXTBFile(XTBFileMixin):
            def __init__(self):
                self.filename = "d.out"

            @property
            def contents(self):
                return ["* finished run on 2026/01/05 at 12:34:56"]

        assert DatedXTBFile().file_date == "2026-01-05 12:34:56"

    def test_file_date_no_marker_returns_none(self):
        class NoDateXTBFile(XTBFileMixin):
            def __init__(self):
                self.filename = "nd.out"

            @property
            def contents(self):
                return ["nothing relevant"]

        assert NoDateXTBFile().file_date is None

    def test_file_date_marker_present_but_unmatched_pattern(self):
        class UnmatchedDateXTBFile(XTBFileMixin):
            def __init__(self):
                self.filename = "ud.out"

            @property
            def contents(self):
                return ["* finished run on today at noon"]

        assert UnmatchedDateXTBFile().file_date is None

    def test_file_date_matched_but_strptime_fails(self):
        class BadDateXTBFile(XTBFileMixin):
            def __init__(self):
                self.filename = "bd.out"

            @property
            def contents(self):
                return ["* finished run on 2026/99/99 at 99:99:99"]

        assert BadDateXTBFile().file_date is None

    def test_method_returns_gfn_directly_when_present_in_route(self):
        class GfnRouteXTBFile(XTBFileMixin):
            def __init__(self):
                self.filename = "g.out"

            @property
            def route_string(self):
                return "xtb structure.xyz --opt --gfn 2"

        assert GfnRouteXTBFile().method == "gfn2"

    def test_method_falls_back_to_hamiltonian(self):
        class HamiltonianXTBFile(XTBFileMixin):
            def __init__(self):
                self.filename = "h.out"
                self.hamiltonian = "GFN2-xTB"

            @property
            def route_string(self):
                return "xtb structure.xyz --opt"

        assert HamiltonianXTBFile().method == "gfn2"

    def test_method_none_without_gfn_or_hamiltonian(self):
        class NoMethodXTBFile(XTBFileMixin):
            def __init__(self):
                self.filename = "nm.out"
                self.hamiltonian = None

            @property
            def route_string(self):
                return "xtb structure.xyz --opt"

        assert NoMethodXTBFile().method is None


class DummyCRESTFile(CRESTFileMixin):
    def __init__(self, filename):
        self.filename = filename

    @property
    def route_string(self):
        return "crest 1a.xyz --cinp constraints.inp --gfn2 --chrg 0 --uhf 0 --optlev tight"


class TestCRESTFileMixin:
    def test_crest_file_properties(self):
        dummy = DummyCRESTFile("test.out")
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


class BaseRegistry(RegistryMixin):
    pass


class SubRegistry1(BaseRegistry):
    pass


class SubRegistry2(BaseRegistry):
    pass


class TestRegistryMixin:
    def test_subclasses(self):
        subclasses = BaseRegistry.subclasses()
        assert SubRegistry1 in subclasses
        assert SubRegistry2 in subclasses

    def test_non_registerable_subclass_is_not_registered(self):
        class NonRegisterableSub(BaseRegistry):
            REGISTERABLE = False

        assert NonRegisterableSub not in BaseRegistry._REGISTRY


class DummyFolder(FolderMixin):
    def __init__(self, folder):
        self.folder = folder


class TestFolderMixin:
    def test_get_all_files_by_suffix(self, temp_folder_with_files):
        folder, file1, file2 = temp_folder_with_files
        dummy = DummyFolder(folder)
        txt_files = dummy.get_all_files_in_current_folder_by_suffix(".txt")
        assert file1 in txt_files

    def test_suffix_without_dot_excludes_compound_extensions(self, tmpdir):
        """filetype='xyz' must match .xyz but not .extxyz."""
        xyz = os.path.join(str(tmpdir), "water.xyz")
        extxyz = os.path.join(str(tmpdir), "crystal.extxyz")
        for path, content in (
            (xyz, "3\n\nO 0 0 0\nH 1 0 0\nH 0 1 0\n"),
            (extxyz, "3\n\nO 0 0 0\nH 1 0 0\nH 0 1 0\n"),
        ):
            with open(path, "w") as f:
                f.write(content)

        dummy = DummyFolder(str(tmpdir))
        files = dummy.get_all_files_in_current_folder_by_suffix("xyz")
        assert xyz in files
        assert extxyz not in files

        # Leading-dot form remains supported.
        files_dotted = dummy.get_all_files_in_current_folder_by_suffix(".xyz")
        assert xyz in files_dotted
        assert extxyz not in files_dotted

    def test_get_all_files_by_regex(self, temp_folder_with_files):
        folder, file1, file2 = temp_folder_with_files
        dummy = DummyFolder(folder)
        log_files = dummy.get_all_files_in_current_folder_matching_regex(
            r".*\.log"
        )
        assert file2 in log_files
        assert file1 not in log_files

    def test_folderpath_is_absolute(self, tmp_path):
        dummy = DummyFolder(str(tmp_path))
        assert dummy.folderpath == os.path.abspath(str(tmp_path))

    def test_get_all_output_files_unsupported_program_raises(self, tmp_path):
        import pytest

        dummy = DummyFolder(str(tmp_path))
        with pytest.raises(ValueError, match="Unsupported program"):
            dummy.get_all_output_files_in_current_folder_by_program(
                "bogus_program"
            )

    def _write(self, path, content):
        path.write_text(content)
        return str(path)

    def test_get_all_output_files_by_program_non_recursive(self, tmp_path):
        gaussian_log = self._write(
            tmp_path / "job.log", "Entering Gaussian System\n"
        )
        self._write(tmp_path / "empty.log", "")
        (tmp_path / "subdir").mkdir()
        self._write(
            tmp_path / "subdir" / "nested.log", "Entering Gaussian System\n"
        )
        self._write(tmp_path / "notes.txt", "just some text\n")
        self._write(tmp_path / "unknown.log", "nothing recognizable\n")

        dummy = DummyFolder(str(tmp_path))
        gaussian_files = (
            dummy.get_all_output_files_in_current_folder_by_program("gaussian")
        )
        assert gaussian_files == [gaussian_log]

        all_files = dummy.get_all_output_files_in_current_folder_by_program()
        assert gaussian_log in all_files
        assert all(f.endswith(".log") or f.endswith(".out") for f in all_files)

    def test_get_all_output_files_by_program_recursive(self, tmp_path):
        (tmp_path / "sub").mkdir()
        nested_log = self._write(
            tmp_path / "sub" / "nested.log", "Entering Gaussian System\n"
        )
        top_log = self._write(
            tmp_path / "top.log", "Entering Gaussian System\n"
        )
        # An empty file and a non-matching-suffix file, both within the
        # recursive walk, must be skipped without raising.
        self._write(tmp_path / "sub" / "empty.log", "")
        self._write(tmp_path / "sub" / "notes.txt", "irrelevant\n")
        # A broken symlink appears in os.walk's "files" list but fails
        # os.path.isfile, covering that skip branch too.
        os.symlink(
            tmp_path / "sub" / "does_not_exist.log",
            tmp_path / "sub" / "broken_link.log",
        )

        dummy = DummyFolder(str(tmp_path))
        found = dummy.get_all_output_files_in_current_folder_and_subfolders_by_program(
            "gaussian"
        )
        assert set(found) == {nested_log, top_log}

    def test_get_all_output_files_recursive_program_none(self, tmp_path):
        (tmp_path / "sub").mkdir()
        nested_log = self._write(
            tmp_path / "sub" / "nested.log", "Entering Gaussian System\n"
        )
        self._write(tmp_path / "sub" / "notes.dat", "irrelevant\n")

        dummy = DummyFolder(str(tmp_path))
        found = (
            dummy.get_all_output_files_in_current_folder_and_subfolders_by_program()
        )
        assert found == [nested_log]

    def test_is_program_calculation_directory(self, tmp_path):
        self._write(tmp_path / "job.log", "Entering Gaussian System\n")
        dummy = DummyFolder(str(tmp_path))
        assert dummy.is_program_calculation_directory("gaussian") is True
        assert dummy.is_program_calculation_directory("orca") is False

    def test_is_program_calculation_directory_missing_folder(self, tmp_path):
        dummy = DummyFolder(str(tmp_path / "does_not_exist"))
        assert dummy.is_program_calculation_directory("gaussian") is False

    def test_get_program_type_from_folder_unknown(self, tmp_path):
        self._write(tmp_path / "notes.txt", "nothing relevant\n")
        dummy = DummyFolder(str(tmp_path))
        assert dummy.get_program_type_from_folder() == "unknown"

    def test_get_program_type_from_folder_single_program(self, tmp_path):
        self._write(tmp_path / "job.out", "x T B\nsome xtb output\n")
        dummy = DummyFolder(str(tmp_path))
        assert dummy.get_program_type_from_folder() == "xtb"

    def test_get_program_type_from_folder_mixed(self, tmp_path):
        self._write(tmp_path / "xtb_job.out", "x T B\nsome xtb output\n")
        self._write(
            tmp_path / "crest_job.out",
            "C R E S T\nsome crest output\n",
        )
        dummy = DummyFolder(str(tmp_path))
        assert dummy.get_program_type_from_folder() == "mixed"

    def test_get_all_files_by_suffix_skips_empty_files(self, tmp_path):
        kept = self._write(tmp_path / "kept.txt", "content")
        self._write(tmp_path / "empty.txt", "")
        dummy = DummyFolder(str(tmp_path))
        assert dummy.get_all_files_in_current_folder_by_suffix(".txt") == [
            kept
        ]

    def test_get_all_files_by_regex_skips_empty_files(self, tmp_path):
        kept = self._write(tmp_path / "kept.log", "content")
        self._write(tmp_path / "empty.log", "")
        dummy = DummyFolder(str(tmp_path))
        matched = dummy.get_all_files_in_current_folder_matching_regex(
            r".*\.log"
        )
        assert matched == [kept]

    def test_get_all_files_and_subfolders_by_suffix(self, tmp_path):
        top = self._write(tmp_path / "top.log", "content")
        (tmp_path / "sub").mkdir()
        nested = self._write(tmp_path / "sub" / "nested.log", "content")
        self._write(tmp_path / "other.txt", "content")

        dummy = DummyFolder(str(tmp_path))
        found = dummy.get_all_files_in_current_folder_and_subfolders_by_suffix(
            ".log"
        )
        assert set(found) == {top, nested}

    def test_get_all_files_by_program_and_suffix_non_recursive(self, tmp_path):
        gaussian_log = self._write(
            tmp_path / "job.log", "Entering Gaussian System\n"
        )
        self._write(tmp_path / "empty.log", "")
        self._write(tmp_path / "unknown.log", "nothing recognizable\n")
        self._write(
            tmp_path / "job.dat", "Entering Gaussian System\n"
        )  # wrong suffix
        (tmp_path / "subdir").mkdir()

        dummy = DummyFolder(str(tmp_path))
        found = dummy.get_all_files_in_current_folder_by_program_and_suffix(
            "gaussian", ".log"
        )
        assert found == [gaussian_log]

    def test_get_all_files_by_program_and_suffix_recursive(self, tmp_path):
        (tmp_path / "sub").mkdir()
        nested = self._write(
            tmp_path / "sub" / "nested.log", "Entering Gaussian System\n"
        )
        top = self._write(tmp_path / "top.log", "Entering Gaussian System\n")
        self._write(
            tmp_path / "sub" / "wrong_suffix.dat",
            "Entering Gaussian System\n",
        )
        self._write(
            tmp_path / "sub" / "wrong_program.log", "nothing recognizable\n"
        )

        dummy = DummyFolder(str(tmp_path))
        found = dummy.get_all_files_in_current_folder_and_subfolders_by_program_and_suffix(
            "gaussian", ".log"
        )
        assert set(found) == {nested, top}

    def test_get_all_files_and_subfolders_matching_regex(self, tmp_path):
        top = self._write(tmp_path / "top.log", "content")
        (tmp_path / "sub").mkdir()
        nested = self._write(tmp_path / "sub" / "nested.log", "content")
        self._write(tmp_path / "other.txt", "content")

        dummy = DummyFolder(str(tmp_path))
        found = dummy.get_all_files_in_current_folder_and_subfolders_matching_regex(
            r".*\.log"
        )
        assert set(found) == {top, nested}


class DummyParser:
    def __init__(self, **attrs):
        for name, value in attrs.items():
            setattr(self, name, value)


class DummyFolderOutput(FolderOutputMixin):
    FILE_PARSERS = ("main_out", "secondary_out")

    def __init__(self, main_out=None, secondary_out=None):
        self.main_out = main_out
        self.secondary_out = secondary_out


class TestFolderOutputMixin:
    def test_delegates_to_first_parser(self):
        output = DummyFolderOutput(
            main_out=DummyParser(total_energy=-123.45),
            secondary_out=DummyParser(total_energy=-123),
        )
        assert output.total_energy == -123.45

    def test_falls_back_to_next_parser(self):
        output = DummyFolderOutput(
            main_out=DummyParser(),
            secondary_out=DummyParser(charge=0),
        )
        assert output.charge == 0

    def test_skips_none_parser(self):
        output = DummyFolderOutput(
            main_out=None,
            secondary_out=DummyParser(multiplicity=1),
        )
        assert output.multiplicity == 1

    def test_raises_attribute_error_when_not_found(self):
        output = DummyFolderOutput()
        with pytest.raises(AttributeError, match="no attribute 'missing'"):
            _ = output.missing
