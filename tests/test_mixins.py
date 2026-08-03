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
