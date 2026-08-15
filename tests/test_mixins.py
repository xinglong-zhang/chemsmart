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
