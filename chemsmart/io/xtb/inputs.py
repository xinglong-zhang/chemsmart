from functools import cached_property
from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import content_blocks_by_paragraph
from chemsmart.io.gaussian.route import GaussianRoute


class XTBInput(FileMixin):
    """XTB input file."""
    def __init__(self, filename):
        self.filename = filename

        from chemsmart.io.molecules.structure import CoordinateBlock

        cb = CoordinateBlock(coordinate_block=self.content_groups[2])
        self.cb = cb



    @property
    def num_content_blocks(self):
        return len(content_blocks_by_paragraph(string_list=self.contents))

    @cached_property
    def content_groups(self):
        """Gaussian input file with content groups.
        # content_groups[0] gives the header information and the route string
        # content_groups[1] gives the title
        # content_groups[2] gives the charge/multiplicity and the xyz coordinates
        # content_groups[3:] gives everything else that are appended at the end of the coordinates:
        # modred, followed by gen/genecp, then custom solvent definitions - the details vary as it
        depends on the actual calculation
        """
        return content_blocks_by_paragraph(string_list=self.contents)

    @property
    def modredundant_group(self):
        if (
            "modred" in self.route_string and self.num_content_blocks > 3
        ):  # in case the input .com file has opt=modredundant
            # in route but no modredundant section at the end
            return self.content_groups[3]

    @property
    def modredundant(self):
        return self._get_modredundant_conditions()

    @property
    def route_object(self):
        try:
            route_object = GaussianRoute(route_string=self.route_string)
            return route_object
        except TypeError as err:
            print(err)

    @property
    def is_pbc(self):
        for line in self.contents:
            line_elements = line.split()
            if len(line_elements) == 4 and line_elements[0].upper() == "TV":
                return True
        return False

    @property
    def translation_vectors(self):
        return self.cb.translation_vectors

    @property
    def chk(self):
        return self._get_chk()

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
    def dieze_tag(self):
        return self.route_object.dieze_tag

    @property
    def job_type(self):
        return self.route_object.job_type

    @job_type.setter
    def job_type(self, value):
        self.route_object.job_type = value

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
    def solv_on(self):
        return self.route_object.solv

    @property
    def solvent_model(self):
        return self.route_object.solvent_model

    @property
    def solvent_id(self):
        return self.route_object.solvent_id

    @property
    def additional_opt_options_in_route(self):
        return self.route_object.additional_opt_options_in_route

    @property
    def additional_route_parameters(self):
        return self.route_object.additional_route_parameters

    @property
    def molecule(self):
        return self.cb.molecule

    @property
    def constrained_atoms(self):
        self.cb.constrained_atoms

    @constrained_atoms.setter
    def constrained_atoms(self, value):
        self.cb.constrained_atoms = value

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

    def _get_route(self):
        """Obtain route string that may span over multiple lines."""
        concatenated_string = ""
        found_hash = False
        for line in self.content_groups[0]:
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

    def _get_modred_frozen_coords(self, modred_list_of_string):
        modred = []
        for raw_line in modred_list_of_string:
            line = raw_line[2:-2]
            line_elems = line.split()
            assert all(
                line_elem.isdigit() for line_elem in line_elems
            ), f"modredundant coordinates should be integers, but is {line_elems} instead."
            each_modred_list = [int(line_elem) for line_elem in line_elems]
            modred.append(each_modred_list)
        return modred

    def _get_modred_scan_coords(self, modred_list_of_string):
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
            ), f"modredundant coordinates should be integers, but is {line_elems[0]} instead."
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
