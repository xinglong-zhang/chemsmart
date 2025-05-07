import inspect
import os
import re
from functools import cached_property

from ase import units

from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.io.orca.route import ORCARoute


class FileMixin:
    """Mixin class for files that can be opened and read"""

    @property
    def filepath(self):
        return os.path.abspath(self.filename)

    @property
    def filepath_directory(self):
        return os.path.split(self.filepath)[0]

    @property
    def base_filename_with_extension(self):
        return os.path.split(self.filepath)[1]

    @property
    def basename(self):
        return self.base_filename_with_extension.split(".")[0]

    @cached_property
    def contents(self):
        with open(self.filepath, "r") as f:
            return [line.strip() for line in f.readlines()]

    @cached_property
    def content_lines_string(self):
        with open(self.filepath, "r") as f:
            return f.read()

    @cached_property
    def forces_in_eV_per_angstrom(self):
        """Convert forces from Hartrees/Bohr to eV/Angstrom."""
        if self.forces is None:
            return [None] * len(self.energies)
        forces_in_eV_per_A = []
        for forces in self.forces:
            forces_in_eV_per_A.append(forces * units.Hartree / units.Bohr)
        return forces_in_eV_per_A

    @cached_property
    def input_translation_vectors(self):
        """Obtain the translation vectors from the input that is printed in the outputfile."""
        tvs = []
        if self.input_coordinates_block is not None:
            cb = self.input_coordinates_block
            if cb.translation_vectors is not None:
                tvs = cb.translation_vectors
        return tvs

    @cached_property
    def symbols(self):
        return self.input_coordinates_block.chemical_symbols

    @cached_property
    def list_of_pbc_conditions(self):
        return self.input_coordinates_block.pbc_conditions

    @cached_property
    def energies_in_eV(self):
        """Convert energies from Hartree to eV."""
        return [energy * units.Hartree for energy in self.energies]

    @property
    def num_energies(self):
        return len(self.energies)


class GaussianFileMixin(FileMixin):
    """Mixin class for Gaussian files."""

    def _get_chk(self):
        for line in self.contents:
            if line.startswith("%chk"):
                return True
        return False

    def _get_mem(self):
        mem = 20  # default value: 20 GB
        for line in self.contents:
            if line.startswith("%mem"):
                mem = int(line.split("=")[-1].split("GB")[0])
        return mem

    def _get_nproc(self):
        nproc = 16  # default value
        for line in self.contents:
            if line.startswith("%nproc"):
                nproc = int(line.split("=")[-1])
        return nproc

    @property
    def route_string(self):
        """Route string for Gaussian file.
        Returned by individual subclasses."""
        return self._get_route()

    def _get_route(self):
        """Default implementation. Subclasses must override this method."""
        raise NotImplementedError("Subclasses must implement `_get_route`.")

    @property
    def modred(self):
        return self._get_modredundant_conditions()

    def _get_modredundant_conditions(self):
        modred = None
        if (
            "modred" in self.route_string
            and self.modredundant_group is not None
        ):
            for line in self.modredundant_group:
                if "F" in line or "f" in line:
                    modred = self._get_modred_frozen_coords(
                        self.modredundant_group
                    )
                    self.job_type = "modred"
                elif "S" in line or "s" in line:
                    modred = self._get_modred_scan_coords(
                        self.modredundant_group
                    )
                    self.job_type = "scan"
                return modred

    @staticmethod
    def _get_modred_frozen_coords(modred_list_of_string):
        modred = []
        for raw_line in modred_list_of_string:
            line = raw_line[2:-2]
            line_elems = line.split()
            assert all(
                line_elem.isdigit() for line_elem in line_elems
            ), f"modred coordinates should be integers, but is {line_elems} instead."
            each_modred_list = [int(line_elem) for line_elem in line_elems]
            modred.append(each_modred_list)
        return modred

    @staticmethod
    def _get_modred_scan_coords(modred_list_of_string):
        modred = {}
        coords = []
        # modred = {'num_steps': 10, 'step_size': 0.05, 'coords': [[1, 2], [3, 4]]}
        for raw_line in modred_list_of_string:
            line = raw_line.strip()[2:]
            line_elems = line.split("S")

            # obtain coords
            coords_string = line_elems[0]
            each_coords_list = coords_string.split()
            assert all(
                line_elem.isdigit() for line_elem in each_coords_list
            ), f"modred coordinates should be integers, but is {line_elems[0]} instead."
            each_modred_list = [
                int(line_elem) for line_elem in each_coords_list
            ]
            coords.append(each_modred_list)
            modred["coords"] = coords

            # obtain num_steps and step_size (assumed the same for each scan coordinate)
            steps_string = line_elems[-1]
            steps_list = steps_string.split()
            num_steps = int(steps_list[0])
            step_size = float(steps_list[1])
            modred["num_steps"] = num_steps
            modred["step_size"] = step_size

        return modred

    @property
    def route_object(self):
        try:
            route_object = GaussianRoute(route_string=self.route_string)
            return route_object
        except TypeError as err:
            print(err)

    @property
    def dieze_tag(self):
        return self.route_object.dieze_tag

    @property
    def job_type(self):
        return self.route_object.job_type

    @job_type.setter
    def job_type(self, value):
        self.route_object.job_type = value

    @property
    def chk(self):
        return self._get_chk()

    @property
    def freq(self):
        return self.route_object.freq

    @property
    def numfreq(self):
        return self.route_object.numfreq

    @property
    def ab_initio(self):
        return self.route_object.ab_initio

    @property
    def functional(self):
        return self.route_object.functional

    @property
    def basis(self):
        return self.route_object.basis

    @property
    def force(self):
        return self.route_object.force

    @property
    def solvent_on(self):
        return self.route_object.solv

    @property
    def solvent_model(self):
        return self.route_object.solvent_model

    @property
    def solvent_id(self):
        return self.route_object.solvent_id

    @property
    def additional_solvent_options(self):
        return self.route_object.additional_solvent_options

    @property
    def additional_opt_options_in_route(self):
        return self.route_object.additional_opt_options_in_route

    @property
    def additional_route_parameters(self):
        return self.route_object.additional_route_parameters

    def read_settings(self):
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        filename = os.path.basename(self.filename)

        title = f"Job prepared from Gaussian file {filename}"

        return GaussianJobSettings(
            ab_initio=self.ab_initio,
            functional=self.functional,
            basis=self.basis,
            charge=self.charge,
            multiplicity=self.multiplicity,
            chk=self.chk,
            job_type=self.job_type,
            title=title,
            freq=self.freq,
            numfreq=self.numfreq,
            dieze_tag=self.dieze_tag,
            solvent_model=self.solvent_model,
            solvent_id=self.solvent_id,
            additional_solvent_options=self.additional_solvent_options,
            additional_opt_options_in_route=self.additional_opt_options_in_route,
            additional_route_parameters=self.additional_route_parameters,
            route_to_be_written=None,
            modred=self.modred,
            gen_genecp_file=None,
            heavy_elements=self.heavy_elements,
            heavy_elements_basis=self.heavy_elements_basis,
            light_elements_basis=self.light_elements_basis,
            custom_solvent=self.custom_solvent,
            append_additional_info=None,
            forces=False,
        )


class ORCAFileMixin(FileMixin):
    """Mixin class for ORCA files."""

    @cached_property
    def contents_string(self):
        return "\n".join(self.contents)

    @property
    def mdci_cutoff(self):
        for i, line in enumerate(self.contents):
            if "%mdci" in line.lower():
                # mdci cutoff in %mdci block
                next_lines = self.contents[i + 1 :]
                for next_line in next_lines:
                    if "cutoff" in next_line.lower():
                        # Find the string prior to it. This is assuming the input is written by
                        # this program where the comment on the cutoff is also written
                        l_elem = next_line.split()
                        c_idx = l_elem.index("cutoff")
                        return l_elem[c_idx - 1]
        return None

    @property
    def mdci_density(self):
        for i, line in enumerate(self.contents):
            if "%mdci" in line.lower():
                # mdci density in %mdci block
                next_lines = self.contents[i + 1 :]
                for next_line in next_lines:
                    next_line_lower = next_line.lower()
                    if "density" in next_line_lower.lower():
                        # Find the string after it
                        l_elem = next_line_lower.split()
                        c_idx = l_elem.index("density")
                        return l_elem[c_idx + 1]
        return None

    @property
    def solvent_model(self):
        cpcm = False
        smd = False

        for i, line in enumerate(self.contents):
            # not even needed to test the route string for solvent
            # if '!' in line and 'cpcm' in line:
            #    cpcm = True

            # solvent specification not in the route string but in the %cpcm block (ORCA_Test_0829.inp)
            if "%cpcm" in line.lower():
                cpcm = True
                next_lines = self.contents[i + 1 :]
                for _, next_line in enumerate(next_lines):
                    if re.search(r"\bsmd\s+true\b", next_line, re.IGNORECASE):
                        smd = True

        if cpcm:
            if smd:
                return "smd"
            return "cpcm"
        return None

    @property
    def solvent_id(self):
        pattern = re.compile(
            r'"([^"]*)"'
        )  # pattern to find text between double quotes
        for line in self.contents:
            line_lower = line.lower()
            if "solvent" in line_lower:
                if not pattern.search(line_lower):
                    raise Exception(
                        "Your input file specifies solvent but solvent is not in quotes, "
                        "thus, your input file is not valid to run for ORCA!"
                    )

                # Find all matches of the pattern in the line
                matches = pattern.findall(line_lower)
                if len(matches) == 1:
                    return matches[0]

                raise Exception(
                    f"{len(matches)} solvents found! Only can specify 1 solvent!"
                )
        return None

    # properties from orca route string
    @property
    def route_string(self):
        """Route string for ORCA file."""
        raise NotImplementedError

    ## properties derived from route string
    @property
    def route_object(self):
        return ORCARoute(route_string=self.route_string)

    @property
    def route_keywords(self):
        return self.route_object.route_keywords

    @property
    def functional(self):
        return self.route_object.functional

    @property
    def ab_initio(self):
        return self.route_object.ab_initio

    @property
    def dispersion(self):
        return self.route_object.dispersion

    @property
    def basis(self):
        return self.route_object.basis

    @property
    def aux_basis(self):
        return self.route_object.auxiliary_basis

    @property
    def extrapolation_basis(self):
        return self.route_object.extrapolation_basis

    @property
    def defgrid(self):
        return self.route_object.defgrid

    @property
    def scf_tol(self):
        return self.route_object.scf_tol

    @property
    def scf_algorithm(self):
        return self.route_object.scf_algorithm

    @property
    def job_type(self):
        return self.route_object.job_type

    @property
    def freq(self):
        return self.route_object.freq

    @property
    def numfreq(self):
        return self.route_object.numfreq

    def read_settings(self):
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        dv = ORCAJobSettings.default()
        return ORCAJobSettings(
            ab_initio=self.ab_initio,
            functional=self.functional,
            dispersion=self.dispersion,
            basis=self.basis,
            aux_basis=self.aux_basis,
            extrapolation_basis=self.extrapolation_basis,
            defgrid=self.defgrid,
            scf_tol=self.scf_tol,
            scf_algorithm=self.scf_algorithm,
            scf_maxiter=self.scf_maxiter,
            scf_convergence=self.scf_convergence,
            charge=self.charge,
            multiplicity=self.multiplicity,
            gbw=dv.gbw,
            freq=self.freq,
            numfreq=self.numfreq,
            dipole=self.dipole,
            quadrupole=self.quadrupole,
            mdci_cutoff=self.mdci_cutoff,
            mdci_density=self.mdci_density,
            job_type=self.job_type,
            solvent_model=self.solvent_model,
            solvent_id=self.solvent_id,
            additional_route_parameters=dv.additional_route_parameters,
            route_to_be_written=dv.route_to_be_written,
            modred=dv.modred,
            gen_genecp_file=dv.gen_genecp_file,
            heavy_elements=dv.heavy_elements,
            heavy_elements_basis=dv.heavy_elements_basis,
            light_elements_basis=dv.light_elements_basis,
            custom_solvent=dv.custom_solvent,
            forces=dv.forces,
        )


class YAMLFileMixin(FileMixin):
    """Mixin class for YAML files."""

    @cached_property
    def yaml_contents_dict(self):
        import yaml

        return yaml.safe_load(self.content_lines_string)

    @property
    def yaml_contents_keys(self):
        return self.yaml_contents_dict.keys()

    @property
    def yaml_contents_values(self):
        return self.yaml_contents_dict.values()

    def yaml_contents_by_key(self, key):
        return self.yaml_contents_dict[key]


class RegistryMeta(type):
    """Metaclass to ensure all subclasses are registered in the root class's _REGISTRY."""

    def __init__(cls, name, bases, dct):
        super().__init__(name, bases, dct)
        # Only initialize _REGISTRY in the root parent class
        if not hasattr(cls, "_REGISTRY"):
            cls._REGISTRY = []


class RegistryMixin(metaclass=RegistryMeta):
    """Mixin to register subclasses in a shared registry."""

    REGISTERABLE = True

    @classmethod
    def subclasses(cls, allow_abstract=False):
        return cls._subclasses(cls, cls._REGISTRY, allow_abstract)

    @staticmethod
    def _subclasses(parent_cls, registry, allow_abstract):
        return [
            c
            for c in registry
            if issubclass(c, parent_cls)
            and c != parent_cls
            and (not inspect.isabstract(c) or allow_abstract)
        ]

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.REGISTERABLE:
            # Append the subclass to the root _REGISTRY
            cls._REGISTRY.append(cls)


# class BlockMixin:
#     """Mixin class for files that can be opened and read"""
#     def write(self, f):
#         for line in self.contents:
#             f.write(line)


class FolderMixin:
    """Mixin class for folders."""

    def get_all_files_in_current_folder_by_suffix(self, filetype):
        """Obtain a list of files of specified type in the folder."""

        all_files = []
        for file in os.listdir(self.folder):
            # check that the file is not empty:
            if os.stat(os.path.join(self.folder, file)).st_size == 0:
                continue
            # collect files of specified type
            if file.endswith(filetype):
                all_files.append(os.path.join(self.folder, file))
        return all_files

    def get_all_files_in_current_folder_and_subfolders_by_suffix(
        self, filetype
    ):
        """Obtain a list of files of specified type in the folder and subfolders."""
        all_files = []
        for subdir, _dirs, files in os.walk(self.folder):
            # subdir is the full path to the subdirectory
            for file in files:
                if file.endswith(filetype):
                    all_files.append(os.path.join(subdir, file))
        return all_files

    def get_all_files_in_current_folder_and_subfolders_matching_regex(
        self, regex
    ):
        """Obtain a list of files matching the regex in the folder and subfolders."""
        all_files = []
        for subdir, _dirs, files in os.walk(self.folder):
            # subdir is the full path to the subdirectory
            for file in files:
                if re.match(regex, file):
                    all_files.append(os.path.join(subdir, file))
        return all_files

    def get_all_files_in_current_folder_matching_regex(self, regex):
        """Obtain a list of files matching the regex in the folder."""
        all_files = []
        for file in os.listdir(self.folder):
            # check that the file is not empty:
            if os.stat(os.path.join(self.folder, file)).st_size == 0:
                continue
            # collect files of specified type
            if re.match(regex, file):
                all_files.append(os.path.join(self.folder, file))
        return all_files


class BaseFolder(FolderMixin):
    def __init__(self, folder):
        self.folder = folder
