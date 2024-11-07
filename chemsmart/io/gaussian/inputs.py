from functools import cached_property
from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import content_blocks_by_paragraph
from chemsmart.io.gaussian.route import GaussianRoute


class Gaussian16Input(FileMixin):
    def __init__(self, filename):
        self.filename = filename
        self.route_object = GaussianRoute(route_string=self.route_string)

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
    def is_pbc(self):
        for line in self.contents:
            line_elements = line.split()
            if len(line_elements) == 4 and line_elements[0].upper() == "TV":
                return True
        return False

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

        return None

    @property
    def dieze_tag(self):
        return self.route_object.dieze_tag

    @property
    def job_type(self):
        return self.route_object.job_type

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
        from chemsmart.io.molecules.structure import CoordinateBlock

        cb = CoordinateBlock(coordinate_block=self.content_groups[2])
        return cb.molecule

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
        pass

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
