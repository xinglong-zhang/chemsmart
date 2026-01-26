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
