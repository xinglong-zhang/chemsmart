import os

from chemsmart.utils.mixins import (
    FileMixin,
    FolderMixin,
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
        return ["F 1 2 3", "S 1 2 10 0.05"]


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


class DummyFolder(FolderMixin):
    def __init__(self, folder):
        self.folder = folder


class TestFolderMixin:
    def test_get_all_files_by_suffix(self, temp_folder_with_files):
        folder, file1, file2 = temp_folder_with_files
        dummy = DummyFolder(folder)
        txt_files = dummy.get_all_files_in_current_folder_by_suffix(".txt")
        assert file1 in txt_files

    def test_get_all_files_by_regex(self, temp_folder_with_files):
        folder, file1, file2 = temp_folder_with_files
        dummy = DummyFolder(folder)
        log_files = dummy.get_all_files_in_current_folder_matching_regex(
            r".*\.log"
        )
        assert file2 in log_files
        assert file1 not in log_files
