from functools import cached_property

from chemsmart.io.gaussian.gengenecp import GenGenECPSection
from chemsmart.utils.mixins import GaussianFileMixin
from chemsmart.utils.utils import (
    content_blocks_by_paragraph,
    get_range_from_list,
)


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

    @cached_property
    def num_content_groups(self):
        return len(self.content_groups)

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
    def oniom_charge(self):
        oniom_charge, _ = self._get_oniom_charge_and_multiplicity()
        return oniom_charge

    @property
    def oniom_multiplicity(self):
        _, oniom_multiplicity = self._get_oniom_charge_and_multiplicity()
        return oniom_multiplicity

    @property
    def route_string(self):
        return self._get_route()

    @property
    def has_frozen_coordinates(self):
        return self.coordinate_block.constrained_atoms

    @property
    def frozen_coordinate_indices(self):
        frozen_coordinate_indices = []
        if self.has_frozen_coordinates:
            for i, frozen_mask in enumerate(list(self.molecule.frozen_atoms)):
                if frozen_mask == -1:
                    # use 1-index throughout
                    frozen_coordinate_indices.append(i + 1)
        if len(frozen_coordinate_indices) == 0:
            return None
        return frozen_coordinate_indices

    @property
    def free_coordinate_indices(self):
        """Obtain list of free coordinate indices from the input format."""
        if self.frozen_coordinate_indices is None:
            return None
        return [
            i
            for i in range(self.num_atoms)
            if i not in self.frozen_coordinate_indices
        ]

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
                (
                    len(line_elements) == 2
                    or len(line_elements) == 12
                    or len(line_elements) == 6
                )
                and line_elements[0].replace("-", "").isdigit()
                and line_elements[1].isdigit()
            ):
                charge = int(line_elements[0])
                multiplicity = int(line_elements[1])
                return charge, multiplicity

    def _get_oniom_charge_and_multiplicity(self):
        # line = self.contents[5]
        line_elements = []
        for line in self.contents:
            line_elements = line.split()
            if (
                all(element.isdigit() for element in line_elements)
                and len(line_elements) > 0
            ):
                break
        charge_multiplicity_list = [
            "real_charge",
            "real_multiplicity",
            "int_charge",
            "int_multiplicity",
            "model_charge",
            "model_multiplicity",
        ]
        oniom_charge = {}
        oniom_multiplicity = {}
        full_line = 12
        if len(self.partition) == 2:
            charge_multiplicity_list = charge_multiplicity_list[0:1, 4:5]
            full_line = 6
        for j in range(0, int(full_line) - len(line_elements)):
            line_elements.append("Not specified, will use default value.")
        for charge in range(0, len(charge_multiplicity_list), 2):
            oniom_charge[charge_multiplicity_list[charge]] = line_elements[
                charge
            ]
        for multiplicity in range(1, len(charge_multiplicity_list), 2):
            oniom_multiplicity[charge_multiplicity_list[multiplicity]] = (
                line_elements[multiplicity]
            )
        return oniom_charge, oniom_multiplicity

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


class Gaussian16QMMMInput(Gaussian16Input):
    """This class has all the properties of Gaussian16Input but also additional ones
    for Gaussian16QMMMInput."""

    @property
    def oniom_charge(self):
        oniom_charge, _ = self._get_oniom_charge_and_multiplicity()
        return oniom_charge

    @property
    def real_charge(self):
        oniom_charge, _ = self._get_oniom_charge_and_multiplicity()
        return int(oniom_charge["real_charge"])

    @property
    def int_charge(self):
        oniom_charge, _ = self._get_oniom_charge_and_multiplicity()
        return int(oniom_charge["int_charge"])

    @property
    def model_charge(self):
        oniom_charge, _ = self._get_oniom_charge_and_multiplicity()
        return int(oniom_charge["model_charge"])

    @property
    def oniom_multiplicity(self):
        _, oniom_multiplicity = self._get_oniom_charge_and_multiplicity()
        return oniom_multiplicity

    @property
    def real_multiplicity(self):
        _, oniom_multiplicity = self._get_oniom_charge_and_multiplicity()
        return int(oniom_multiplicity["real_multiplicity"])

    @property
    def int_multiplicity(self):
        _, oniom_multiplicity = self._get_oniom_charge_and_multiplicity()
        return int(oniom_multiplicity["int_multiplicity"])

    @property
    def model_multiplicity(self):
        _, oniom_multiplicity = self._get_oniom_charge_and_multiplicity()
        return int(oniom_multiplicity["model_multiplicity"])

    @property
    def partition(self):
        """Get the partition string."""
        partition = {}
        for key, val in [
            ("high level atoms", self.molecule.high_level_atoms),
            ("medium level atoms", self.molecule.medium_level_atoms),
            ("low level atoms", self.molecule.low_level_atoms),
        ]:
            if val is not None:
                partition[key] = get_range_from_list(val)
        return partition

    def _get_oniom_charge_and_multiplicity(self):
        line_elements = []
        for line in self.contents:
            line_elements = line.split()
            if (
                all(element.isdigit() for element in line_elements)
                and len(line_elements) > 0
            ):
                break
        charge_multiplicity_list = [
            "real_charge",
            "real_multiplicity",
            "int_charge",
            "int_multiplicity",
            "model_charge",
            "model_multiplicity",
        ]
        oniom_charge = {}
        oniom_multiplicity = {}
        full_line = 12
        if len(self.partition) == 2:
            charge_multiplicity_list = charge_multiplicity_list[0:1, 4:5]
            full_line = 6
        for j in range(0, int(full_line) - len(line_elements)):
            line_elements.append("Not specified, will use default value.")
        for charge in range(0, len(charge_multiplicity_list), 2):
            oniom_charge[charge_multiplicity_list[charge]] = line_elements[
                charge
            ]
        for multiplicity in range(1, len(charge_multiplicity_list), 2):
            oniom_multiplicity[charge_multiplicity_list[multiplicity]] = (
                line_elements[multiplicity]
            )
        return oniom_charge, oniom_multiplicity
