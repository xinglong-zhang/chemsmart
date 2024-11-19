import os
import re
from functools import cached_property
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


class ORCAFileMixin:
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
                        # Find the string prior to it. This is assuming the input is written by pyatoms
                        # where the comment on the cutoff is also written
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


# class BlockMixin:
#     """Mixin class for files that can be opened and read"""
#     def write(self, f):
#         for line in self.contents:
#             f.write(line)
