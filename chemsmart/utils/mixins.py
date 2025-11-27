"""
Mixin classes for file handling and computational chemistry operations.

This module provides mixin classes that add common functionality to
computational chemistry file readers and data processors. Includes
mixins for file operations, Gaussian/ORCA specific functionality,
YAML file handling, and folder operations.

Key mixin classes:
- FileMixin: Basic file reading and property extraction
- GaussianFileMixin: Gaussian-specific file parsing
- ORCAFileMixin: ORCA-specific file parsing
- YAMLFileMixin: YAML file handling
- RegistryMixin: Automatic subclass registration
- FolderMixin: Directory and file search operations
"""

import inspect
import os
import re
from functools import cached_property

from ase import units

from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.io.orca.route import ORCARoute


class FileMixin:
    """
    Mixin class for files that can be opened and read.

    Provides common file operations and property extraction methods
    for computational chemistry files. Includes path handling,
    content reading, and unit conversion utilities.
    """

    @property
    def filepath(self):
        """
        Get the absolute path of the file.

        Returns:
            str: Absolute file path.
        """
        return os.path.abspath(self.filename)

    @property
    def filepath_directory(self):
        """
        Get the directory containing the file.

        Returns:
            str: Directory path containing the file.
        """
        return os.path.split(self.filepath)[0]

    @property
    def base_filename_with_extension(self):
        """
        Get the filename with extension (no directory path).

        Returns:
            str: Filename with extension.
        """
        return os.path.split(self.filepath)[1]

    @property
    def basename(self):
        """
        Get the filename without extension.

        Returns:
            str: Base filename without extension.
        """
        return self.base_filename_with_extension.split(".")[0]

    @cached_property
    def contents(self):
        """
        Read and cache file contents as list of stripped lines.

        Returns:
            list: List of strings, each representing a line from the file.
        """
        with open(self.filepath, "r") as f:
            return [line.strip() for line in f.readlines()]

    @cached_property
    def content_lines_string(self):
        """
        Read and cache file contents as a single string.

        Returns:
            str: Complete file contents as a single string.
        """
        with open(self.filepath, "r") as f:
            return f.read()

    @cached_property
    def forces_in_eV_per_angstrom(self):
        """
        Convert forces from Hartrees/Bohr to eV/Angstrom.

        Returns:
            list: List of force arrays converted to eV/Angstrom units.
        """
        if self.forces is None:
            return [None] * len(self.energies)
        forces_in_eV_per_A = []
        for forces in self.forces:
            forces_in_eV_per_A.append(forces * units.Hartree / units.Bohr)
        return forces_in_eV_per_A

    @cached_property
    def input_translation_vectors(self):
        """
        Obtain translation vectors from input printed in the output file.

        Returns:
            list: Translation vectors from input coordinates block.
        """
        tvs = []
        if self.input_coordinates_block is not None:
            cb = self.input_coordinates_block
            if cb.translation_vectors is not None:
                tvs = cb.translation_vectors
        return tvs

    @cached_property
    def symbols(self):
        """
        Get chemical symbols from input coordinates block.

        Returns:
            list: List of chemical symbols.
        """
        return self.input_coordinates_block.chemical_symbols

    @cached_property
    def list_of_pbc_conditions(self):
        """
        Get periodic boundary conditions from input coordinates.

        Returns:
            list: Periodic boundary conditions.
        """
        return self.input_coordinates_block.pbc_conditions

    @cached_property
    def energies_in_eV(self):
        """
        Convert energies from Hartree to eV.

        Returns:
            list: List of energies converted to eV units.
        """
        return [energy * units.Hartree for energy in self.energies]

    @property
    def num_energies(self):
        """
        Get the number of energy values.

        Returns:
            int: Number of energy values.
        """
        return len(self.energies)


class GaussianFileMixin(FileMixin):
    """
    Mixin class for Gaussian computational chemistry files.

    Extends FileMixin with Gaussian-specific functionality including
    route string parsing, job type detection, and settings extraction.
    Handles Gaussian input/output file formats and job parameters.
    """

    def _get_chk(self):
        """
        Check if checkpoint file directive is present.

        Returns:
            bool: True if %chk directive found, False otherwise.
        """
        for line in self.contents:
            if line.startswith("%chk"):
                return True
        return False

    def _get_mem(self):
        """
        Extract memory allocation from file.

        Returns:
            int: Memory allocation in GB (default 20 if not found).
        """
        mem = 20  # Default value: 20 GB
        for line in self.contents:
            if line.startswith("%mem"):
                mem = int(line.split("=")[-1].split("GB")[0])
        return mem

    def _get_nproc(self):
        """
        Extract number of processors from file.

        Returns:
            int: Number of processors (default 16 if not found).
        """
        nproc = 16  # Default value
        for line in self.contents:
            if line.startswith("%nproc"):
                nproc = int(line.split("=")[-1])
        return nproc

    @property
    def freq(self):
        """
        Get frequency calculation setting from route.

        Returns:
            bool: True if frequency calculations are requested.
        """
        return self.route_object.freq

    @property
    def route_string(self):
        """
        Get the route string for Gaussian file.

        Returns the computational route string as defined in the
        Gaussian input file. Implementation is provided by subclasses.

        Returns:
            str: Route string for Gaussian calculations.
        """
        return self._get_route()

    def _get_route(self):
        """
        Get route string from file contents.

        Default implementation that must be overridden by subclasses
        to provide specific route string extraction logic.

        Raises:
            NotImplementedError: Must be implemented by subclasses.
        """
        raise NotImplementedError("Subclasses must implement `_get_route`.")

    @property
    def is_link(self):
        """
        Check if file contains a Gaussian link job.

        Searches for link job indicators such as 'stable=opt' in
        route lines to determine if this is a multi-step job.

        Returns:
            bool: True if link job detected, False otherwise.
        """
        for line in self.contents:
            if line.startswith("#") and "stable=opt" in line:
                return True
        return False

    @property
    def modred(self):
        """
        Get modified redundant coordinate specifications.

        Extracts modredundant coordinate definitions from the file
        for constrained optimizations or coordinate scans.

        Returns:
            dict or list or None: Modredundant coordinate specifications.
        """
        return self._get_modredundant_conditions()

    def _get_modredundant_conditions(self):
        """
        Extract modredundant coordinate conditions from route.

        Parses modredundant coordinate specifications for frozen
        coordinates (F) or scan coordinates (S) from the file.

        Returns:
            dict or list or None: Modredundant conditions or None.
        """
        modred = None
        if (
            "modred" in self.route_string
            and self.modredundant_group is not None
        ):
            # Check if any line is a scan coordinate
            is_scan = False
            for line in self.modredundant_group:
                if "S" in line or "s" in line:
                    is_scan = True
                    break

            if is_scan:
                modred = self._get_modred_scan_coords(self.modredundant_group)
                self.job_type = "scan"
            else:
                modred = self._get_modred_frozen_coords(
                    self.modredundant_group
                )
                self.job_type = "modred"
            return modred

    @staticmethod
    def _get_modred_frozen_coords(modred_list_of_string):
        """
        Parse frozen coordinate specifications from modredundant block.

        Extracts integer coordinate indices for frozen internal
        coordinates from modredundant input lines.

        Args:
            modred_list_of_string (list): Raw modredundant input lines.

        Returns:
            list: List of frozen coordinate specifications.
        """
        modred = []
        for raw_line in modred_list_of_string:
            if "F" not in raw_line and "f" not in raw_line:
                continue
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
        """
        Parse scan coordinate specifications from modredundant block.

        Extracts coordinate scan parameters including atom indices,
        number of steps, and step size from modredundant input.

        Args:
            modred_list_of_string (list): Raw modredundant input lines.

        Returns:
            dict: Scan coordinate specifications with keys: 'coords',
                  'num_steps', 'step_size'.
        """
        modred = {}
        coords = []
        # modred = {'num_steps': 10, 'step_size': 0.05, 'coords': [[1, 2], [3, 4]]}
        for raw_line in modred_list_of_string:
            if "S" not in raw_line and "s" not in raw_line:
                continue
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
        """
        Get parsed route object from route string.

        Creates a GaussianRoute object from the route string to
        provide structured access to calculation parameters.

        Returns:
            GaussianRoute or None: Parsed route object or None on error.
        """
        try:
            route_object = GaussianRoute(route_string=self.route_string)
            return route_object
        except TypeError as err:
            print(err)

    @property
    def dieze_tag(self):
        """
        Get dieze tag from route object.

        Returns the dieze tag (calculation level indicator) from
        the parsed route string.

        Returns:
            str: Dieze tag indicating calculation level.
        """
        return self.route_object.dieze_tag

    @property
    def job_type(self):
        """
        Get job type from route object.

        Returns the type of calculation (e.g., 'opt', 'freq', 'sp')
        as determined from the route string.

        Returns:
            str: Job type specification.
        """
        return self.route_object.job_type

    @job_type.setter
    def job_type(self, value):
        """
        Set job type in route object.

        Updates the job type in the route object for dynamic
        job type modification during parsing.

        Args:
            value (str): New job type to set.
        """
        self.route_object.job_type = value

    @property
    def chk(self):
        """
        Get checkpoint file usage status.

        Returns whether checkpoint file directive is specified
        in the input file.

        Returns:
            bool: True if checkpoint file is used, False otherwise.
        """
        return self._get_chk()

    @property
    def numfreq(self):
        """
        Get numerical frequency calculation setting.

        Returns whether numerical frequency calculations are
        requested in the route string.

        Returns:
            bool: True if numerical frequencies requested.
        """
        return self.route_object.numfreq

    @property
    def ab_initio(self):
        """
        Get ab initio method from route string.

        Returns the ab initio quantum chemistry method specified
        in the calculation route.

        Returns:
            str or None: Ab initio method name or None if not specified.
        """
        return self.route_object.ab_initio

    @property
    def functional(self):
        """
        Get DFT functional from route string.

        Returns the density functional theory functional specified
        in the calculation route.

        Returns:
            str or None: DFT functional name or None if not specified.
        """
        return self.route_object.functional

    @property
    def basis(self):
        """
        Get basis set from route string.

        Returns the basis set specification from the calculation
        route string.

        Returns:
            str or None: Basis set name or None if not specified.
        """
        return self.route_object.basis

    @property
    def semiempirical(self):
        """
        Get semiempirical method from route string.

        Returns the semiempirical quantum chemistry method
        specified in the calculation route.

        Returns:
            str or None: Semiempirical method or None if not specified.
        """
        return self.route_object.semiempirical

    @property
    def force(self):
        """
        Get force calculation setting from route string.

        Returns whether force calculations are requested in
        the calculation route.

        Returns:
            bool: True if force calculations requested.
        """
        return self.route_object.force

    @property
    def solvent_on(self):
        """
        Get solvent model activation status.

        Returns whether solvent model calculations are enabled
        in the calculation route.

        Returns:
            bool: True if solvent model is active.
        """
        return self.route_object.solv

    @property
    def solvent_model(self):
        """
        Get solvent model type from route string.

        Returns the type of implicit solvent model specified
        in the calculation route.

        Returns:
            str or None: Solvent model name or None if not specified.
        """
        return self.route_object.solvent_model

    @property
    def solvent_id(self):
        """
        Get solvent identifier from route string.

        Returns the specific solvent identifier used in implicit
        solvent calculations.

        Returns:
            str or None: Solvent identifier or None if not specified.
        """
        return self.route_object.solvent_id

    @property
    def additional_solvent_options(self):
        """
        Get additional solvent options from route string.

        Returns any additional solvent-related parameters
        specified in the calculation route.

        Returns:
            str or None: Additional solvent options or None.
        """
        return self.route_object.additional_solvent_options

    @property
    def additional_opt_options_in_route(self):
        """
        Get additional optimization options from route string.

        Returns any additional optimization parameters specified
        in the calculation route.

        Returns:
            str or None: Additional optimization options or None.
        """
        return self.route_object.additional_opt_options_in_route

    @property
    def additional_route_parameters(self):
        """
        Get additional route parameters from route string.

        Returns any additional calculation parameters not covered
        by standard route parsing.

        Returns:
            str or None: Additional route parameters or None.
        """
        return self.route_object.additional_route_parameters

    def read_settings(self):
        """
        Create GaussianJobSettings from file parameters.

        Extracts all relevant calculation parameters from the file
        and creates a comprehensive settings object for job submission.

        Returns:
            GaussianJobSettings: Complete job settings configuration.
        """
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        filename = os.path.basename(self.filename)

        title = f"Job prepared from Gaussian file {filename}"

        return GaussianJobSettings(
            ab_initio=self.ab_initio,
            functional=self.functional,
            basis=self.basis,
            semiempirical=self.semiempirical,
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
    """
    Mixin class for ORCA computational chemistry files.

    Extends FileMixin with ORCA-specific functionality including
    route string parsing, solvent model detection, and MDCI parameter
    extraction. Handles ORCA input/output file formats and job settings.
    """

    @cached_property
    def contents_string(self):
        """
        Get file contents as a single joined string.

        Returns:
            str: Complete file contents with newlines preserved.
        """
        return "\n".join(self.contents)

    @property
    def mdci_cutoff(self):
        """
        Extract MDCI cutoff parameter from input file.

        Searches for MDCI cutoff specification in the %mdci block
        and returns the cutoff value for multireference calculations.

        Returns:
            str or None: MDCI cutoff value or None if not found.
        """
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
        """
        Extract MDCI density parameter from input file.

        Searches for MDCI density specification in the %mdci block
        and returns the density value for multireference calculations.

        Returns:
            str or None: MDCI density value or None if not found.
        """
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
        """
        Determine solvent model type from ORCA input file.

        Scans the "%cpcm" block to detect CPCM/SMD usage. Returns "smd" if a
        "%cpcm" block is present and contains "smd true"; returns "cpcm" if a
        "%cpcm" block is present without the SMD flag; returns None otherwise.
        Note: the route line is not parsed for solvent specification here.

        Returns:
            str or None: "cpcm", "smd", or None if no solvent model found.
        """
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
        """
        Extract solvent identifier from ORCA input file.

        Searches for solvent specification in quoted strings and
        validates that exactly one solvent is specified.

        Returns:
            str or None: Solvent identifier or None if not found.

        Raises:
            Exception: If solvent not in quotes or multiple solvents found.
        """
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
        """
        Route string for ORCA file.
        """
        raise NotImplementedError

    ## properties derived from route string
    @property
    def route_object(self):
        """
        Get parsed ORCA route object from route string.

        Creates an ORCARoute object from the route string to
        provide structured access to calculation parameters.

        Returns:
            ORCARoute: Parsed ORCA route object.
        """
        return ORCARoute(route_string=self.route_string)

    @property
    def route_keywords(self):
        """
        Get route keywords from ORCA route object.

        Returns the parsed route keywords from the ORCA
        calculation specification.

        Returns:
            list: List of route keywords.
        """
        return self.route_object.route_keywords

    @property
    def functional(self):
        """
        Get DFT functional from ORCA route string.

        Returns the density functional theory functional specified
        in the ORCA calculation route.

        Returns:
            str or None: DFT functional name or None if not specified.
        """
        return self.route_object.functional

    @property
    def ab_initio(self):
        """
        Get ab initio method from ORCA route string.

        Returns the ab initio quantum chemistry method specified
        in the ORCA calculation route.

        Returns:
            str or None: Ab initio method name or None if not specified.
        """
        return self.route_object.ab_initio

    @property
    def dispersion(self):
        """
        Get dispersion correction from ORCA route string.

        Returns the dispersion correction method specified
        in the ORCA calculation route.

        Returns:
            str or None: Dispersion correction or None if not specified.
        """
        return self.route_object.dispersion

    @property
    def basis(self):
        """
        Get basis set from ORCA route string.

        Returns the basis set specification from the ORCA
        calculation route string.

        Returns:
            str or None: Basis set name or None if not specified.
        """
        return self.route_object.basis

    @property
    def aux_basis(self):
        """
        Get auxiliary basis set from ORCA route string.

        Returns the auxiliary basis set for density fitting
        calculations from the ORCA route.

        Returns:
            str or None: Auxiliary basis set or None if not specified.
        """
        return self.route_object.auxiliary_basis

    @property
    def extrapolation_basis(self):
        """
        Get extrapolation basis from ORCA route string.

        Returns the basis set used for extrapolation schemes
        in the ORCA calculation route.

        Returns:
            str or None: Extrapolation basis or None if not specified.
        """
        return self.route_object.extrapolation_basis

    @property
    def defgrid(self):
        """
        Get integration grid specification from ORCA route.

        Returns the integration grid definition specified
        in the ORCA calculation route.

        Returns:
            str or None: Grid specification or None if not specified.
        """
        return self.route_object.defgrid

    @property
    def scf_tol(self):
        """
        Get SCF convergence tolerance from ORCA route.

        Returns the self-consistent field convergence tolerance
        specified in the ORCA calculation route.

        Returns:
            str or None: SCF tolerance or None if not specified.
        """
        return self.route_object.scf_tol

    @property
    def scf_algorithm(self):
        """
        Get SCF algorithm from ORCA route string.

        Returns the self-consistent field algorithm specified
        in the ORCA calculation route.

        Returns:
            str or None: SCF algorithm or None if not specified.
        """
        return self.route_object.scf_algorithm

    @property
    def job_type(self):
        """
        Get job type from ORCA route string.

        Returns the type of calculation (e.g., 'opt', 'freq', 'sp')
        as determined from the ORCA route string.

        Returns:
            str: Job type specification.
        """
        return self.route_object.job_type

    @property
    def freq(self):
        """
        Get frequency calculation setting from ORCA route.

        Returns whether frequency calculations are requested
        in the ORCA calculation route.

        Returns:
            bool: True if frequency calculations requested.
        """
        return self.route_object.freq

    @property
    def numfreq(self):
        """
        Get numerical frequency calculation setting from ORCA route.

        Returns whether numerical frequency calculations are
        requested in the ORCA route string.

        Returns:
            bool: True if numerical frequencies requested.
        """
        return self.route_object.numfreq

    def read_settings(self):
        """
        Create ORCAJobSettings from file parameters.

        Extracts all relevant ORCA calculation parameters from the file
        and creates a comprehensive settings object for job submission.

        Returns:
            ORCAJobSettings: Complete ORCA job settings configuration.
        """
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
    """
    Mixin class for YAML file handling and parsing.

    Extends FileMixin with YAML-specific functionality including
    content parsing, key/value extraction, and dictionary access.
    Provides convenient methods for working with YAML configuration files.
    """

    @cached_property
    def yaml_contents_dict(self):
        """
        Parse YAML file contents into a Python object.

        Uses `yaml.safe_load` to read the YAML root (typically a mapping).
        The return type depends on the file contents: dict (common), list,
        scalar, or None for empty files.

        Returns:
            Any: Parsed YAML root object.
        """
        import yaml

        return yaml.safe_load(self.content_lines_string)

    @property
    def yaml_contents_keys(self):
        """
        Get all keys from the parsed YAML dictionary.

        Returns the top-level keys from the YAML file contents
        for iteration and inspection.

        Returns:
            dict_keys: Dictionary keys from YAML contents.
        """
        return self.yaml_contents_dict.keys()

    @property
    def yaml_contents_values(self):
        """
        Get all values from the parsed YAML dictionary.

        Returns the top-level values from the YAML file contents
        for iteration and inspection.

        Returns:
            dict_values: Dictionary values from YAML contents.
        """
        return self.yaml_contents_dict.values()

    def yaml_contents_by_key(self, key):
        """
        Get YAML content value by specific key.

        Retrieves the value associated with a specific key from
        the parsed YAML dictionary.

        Args:
            key: The key to look up in the YAML dictionary.

        Returns:
            Any: Value associated with the key.
        """
        return self.yaml_contents_dict[key]


class RegistryMeta(type):
    """
    Metaclass that seeds a shared subclass registry on the root class.

    Initializes a `_REGISTRY` list on the first (root) class in a hierarchy.
    Actual automatic registration of subclasses happens in
    `RegistryMixin.__init_subclass__`, which appends new subclasses to
    the root's `_REGISTRY` when `REGISTERABLE` is True.
    """

    def __init__(cls, name, bases, dct):
        """
        Initialize class and set up registry if needed.

        Args:
            name (str): Class name.
            bases (tuple): Base classes.
            dct (dict): Class dictionary.
        """
        super().__init__(name, bases, dct)
        # Only initialize _REGISTRY in the root parent class
        if not hasattr(cls, "_REGISTRY"):
            cls._REGISTRY = []


class RegistryMixin(metaclass=RegistryMeta):
    """
    Mixin to automatically register subclasses in a shared registry.

    Provides automatic subclass registration and discovery functionality.
    Useful for implementing plugin systems and factory patterns where
    subclasses need to be dynamically discovered and instantiated.
    """

    # Flag to control whether this class should be registered in the registry
    REGISTERABLE = True

    @classmethod
    def subclasses(cls, allow_abstract=False):
        """
        Get all registered subclasses of this class.

        Args:
            allow_abstract (bool): Whether to include abstract classes.
                Defaults to False.

        Returns:
            list: List of subclass types.
        """
        return cls._subclasses(cls, cls._REGISTRY, allow_abstract)

    @staticmethod
    def _subclasses(parent_cls, registry, allow_abstract):
        """
        Filter registry for subclasses of the parent class.

        Args:
            parent_cls (type): Parent class to filter by.
            registry (list): Registry of all classes.
            allow_abstract (bool): Whether to include abstract classes.

        Returns:
            list: Filtered list of subclasses.
        """
        return [
            c
            for c in registry
            if issubclass(c, parent_cls)
            and c != parent_cls
            and (not inspect.isabstract(c) or allow_abstract)
        ]

    def __init_subclass__(cls, **kwargs):
        """
        Automatically register subclass in the registry.

        Called when a subclass is created to automatically add it
        to the shared registry if REGISTERABLE is True.

        Args:
            **kwargs: Additional keyword arguments passed to parent.
        """
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
    """
    Mixin class for folder operations and file discovery.

    Provides methods for searching and filtering files within directories
    based on file extensions, regular expressions, and other criteria.
    Supports both current directory and recursive subdirectory searches.

    Attributes:
        folder (str): Path to the base directory. Consumers are expected to
            set this attribute (e.g., via BaseFolder or a subclass).
    """

    def get_all_files_in_current_folder_by_suffix(self, filetype):
        """
        Obtain a list of files of specified type in the current folder.

        Non-recursively lists files in `self.folder` whose names end with
        `filetype`. Empty files are excluded.

        Args:
            filetype (str): File name suffix to match (e.g., '.log', '.out').

        Returns:
            list[str]: Full file paths matching the suffix.
        """
        all_files = []
        for file in os.listdir(self.folder):
            # Check that the file is not empty:
            if os.stat(os.path.join(self.folder, file)).st_size == 0:
                continue
            # Collect files of specified type
            if file.endswith(filetype):
                all_files.append(os.path.join(self.folder, file))
        return all_files

    def get_all_files_in_current_folder_and_subfolders_by_suffix(
        self, filetype
    ):
        """
        Obtain files of specified type in folder and all subfolders.

        Recursively searches `self.folder` for files whose names end with
        `filetype`. Unlike the non-recursive variant, empty files are not
        filtered out here.

        Args:
            filetype (str): File name suffix to match.

        Returns:
            list[str]: Full file paths matching the suffix.
        """
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
        """
        Obtain files matching regex pattern in folder and subfolders.

        Recursively searches `self.folder` for files whose names match
        the given regular expression using `re.match` (anchored at the
        beginning of the filename). Empty files are not filtered out here.

        Args:
            regex (str | re.Pattern): Regular expression pattern to match
                against filenames (not full paths).

        Returns:
            list[str]: Full file paths whose basenames match the pattern.
        """
        all_files = []
        for subdir, _dirs, files in os.walk(self.folder):
            # subdir is the full path to the subdirectory
            for file in files:
                if re.match(regex, file):
                    all_files.append(os.path.join(subdir, file))
        return all_files

    def get_all_files_in_current_folder_matching_regex(self, regex):
        """
        Obtain files matching regex pattern in the current folder only.

        Searches only `self.folder` for files whose names match the given
        regular expression using `re.match` (anchored at the beginning of
        the filename). Empty files are excluded.

        Args:
            regex (str | re.Pattern): Regular expression pattern to match
                against filenames (not full paths).

        Returns:
            list[str]: Full file paths whose basenames match the pattern.
        """
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
    """
    Base class for folder management with path operations.

    Provides basic folder initialization and validation functionality.
    Combines folder path management with file discovery capabilities
    from the FolderMixin class.
    """

    def __init__(self, folder):
        """
        Initialize the BaseFolder with a specified path.

        Args:
            folder (str): Path to the folder to be managed.
        """
        self.folder = folder
