from functools import cached_property

from chemsmart.io.gaussian.gengenecp import GenGenECPSection
from chemsmart.utils.mixins import GaussianFileMixin
from chemsmart.utils.utils import content_blocks_by_paragraph


class Gaussian16Input(GaussianFileMixin):
    def __init__(self, filename):
        self.filename = filename

    @property
    def num_content_blocks(self):
        return len(self.content_groups)

    @cached_property
    def content_groups(self):
        """Gaussian input file with content groups.
        # content_groups[0] gives the header information and the route string
        # content_groups[1] gives the title
        # content_groups[2] gives the charge/multiplicity and the xyz coordinates
        # content_groups[3:] gives everything else that are appended at the end of the coordinates:
        # modred, followed by gen/genecp, then custom solvent definitions - the details vary as it
        depends on the actual calculation
        # may need to be updated if there are job specific sections appended at the end.
        """
        return content_blocks_by_paragraph(string_list=self.contents)

    @property
    def coordinate_block(self):
        from chemsmart.io.molecules.structure import CoordinateBlock

        cb = CoordinateBlock(coordinate_block=self.content_groups[2])
        return cb

    @property
    def modredundant_group(self):
        if (
            "modred" in self.route_string or "modred" in self.route_string
        ) and self.num_content_blocks > 3:
            # in case the input .com file has opt=modred
            # in route but no modred section at the end
            return self.content_groups[3]

    @property
    def is_pbc(self):
        for line in self.contents:
            line_elements = line.split()
            if len(line_elements) == 4 and line_elements[0].upper() == "TV":
                return True
        return False

    @property
    def translation_vectors(self):
        return self.coordinate_block.translation_vectors

    @property
    def mem(self):
        return self._get_mem()

    @property
    def nproc(self):
        return self._get_nproc()

    @property
    def charge(self):
        charge, _ = self._get_charge_and_multiplicity()
        return charge

    @property
    def multiplicity(self):
        _, multiplicity = self._get_charge_and_multiplicity()
        return multiplicity

    @property
    def route_string(self):
        return self._get_route()

    @property
    def gen_genecp_group(self):
        """Block of strings in the input file specifying the gen/genecp group."""
        return self._get_gen_genecp_group()

    @property
    def genecp_section(self):
        if self.gen_genecp_group:
            return GenGenECPSection.from_genecp_group(
                genecp_group=self.gen_genecp_group
            )
        return None

    @property
    def light_elements(self):
        if self.genecp_section:
            return self.genecp_section.light_elements
        return None

    @property
    def light_elements_basis(self):
        if self.genecp_section:
            return self.genecp_section.light_elements_basis
        return None

    @property
    def heavy_elements(self):
        if self.genecp_section:
            return self.genecp_section.heavy_elements
        return None

    @property
    def heavy_elements_basis(self):
        if self.genecp_section:
            return self.genecp_section.heavy_elements_basis
        return None

    @property
    def custom_solvent(self):
        """Get the custom solvent string."""
        return self._get_custom_solvent_string()

    @property
    def custom_solvent_group(self):
        return self._get_custom_solvent_group()

    @property
    def molecule(self):
        molecule = self.coordinate_block.molecule
        # update the charge and multiplicity of the molecule
        molecule.charge = self.charge
        molecule.multiplicity = self.multiplicity
        return molecule

    @property
    def num_atoms(self):
        return self.molecule.num_atoms

    @property
    def constrained_atoms(self):
        return self.coordinate_block.constrained_atoms

    @constrained_atoms.setter
    def constrained_atoms(self, value):
        self.constrained_atoms = value

    def _get_route(self):
        """Obtain route string that may span over multiple lines
        and convert route to lower case."""
        concatenated_string = ""
        for line in self.content_groups[0]:
            found_hash = False
            if line.startswith("#"):
                concatenated_string += (
                    line.strip()
                )  # Remove the '#' and any leading/trailing whitespace
                found_hash = True
            elif found_hash:
                concatenated_string += (
                    " " + line.strip()
                )  # Concatenate with a space separator
            else:
                continue
            return concatenated_string.lower()

    def _get_charge_and_multiplicity(self):
        for line in self.contents:
            line_elements = line.split()
            if (
                len(line_elements) == 2
                and line_elements[0].replace("-", "").isdigit()
                and line_elements[1].isdigit()
            ):
                charge = int(line_elements[0])
                multiplicity = int(line_elements[1])
                return charge, multiplicity

    def _get_gen_genecp_group(self):
        if self.basis is None:  # this happens for semi-empirical calculations
            return None
        if "gen" not in self.basis:
            return None
        if (
            "modred" in self.route_string
            and "solvent=generic" in self.route_string
        ):
            return self.content_groups[4:-1]
        if (
            "modred" in self.route_string
            and "solvent=generic" not in self.route_string
        ):
            return self.content_groups[
                4:
            ]  # need to change if there are additional append info after these
        if (
            "modred" not in self.route_string
            and "solvent=generic" in self.route_string
        ):
            return self.content_groups[3:-1]
        if (
            "modred" not in self.route_string
            and "solvent=generic" not in self.route_string
        ):
            return self.content_groups[3:]
        return None

    def _get_custom_solvent_group(self):
        """Get the custom solvent group from the content groups.
        Custom solvent is always the last content group in the input file.
        """
        if "solvent=generic" in self.route_string:
            return self.content_groups[-1]
        return None

    def _get_custom_solvent_string(self):
        if self.custom_solvent_group is not None:
            return "\n".join(self.custom_solvent_group)
        return None
