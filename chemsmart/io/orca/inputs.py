import logging

import numpy as np

from pyatoms.io.ase.atoms import AtomsWrapper
from pyatoms.io.orca.utils import ORCAFile
from pyatoms.utils.periodictable import to_element
from pyatoms.utils.utils import is_float

logger = logging.getLogger(__name__)


class ORCAInput(ORCAFile):
    def __init__(self, inpfile):
        super().__init__(filename=inpfile)

    @property
    def route_string(self):
        """Route string for ORCA file, convert to lower case."""
        route_string_lines = []
        for line in self.contents:
            if line.startswith("!"):
                new_line = line.replace("!", "")
                route_string_lines.append(new_line)

        route_string = " ".join(route_string_lines)
        return f"! {route_string}".lower()

    @property
    def coordinate_type(self):
        for line in self.contents:
            if line.startswith("*") and len(line) > 1:
                line_elem = line.split()
                return line_elem[1]
        return None

    @property
    def charge(self):
        for line in self.contents:
            if line.startswith("*") and len(line) > 1:
                line_elem = line.split()
                return int(line_elem[2])
        return None

    @property
    def multiplicity(self):
        for line in self.contents:
            if line.startswith("*") and len(line) > 1:
                line_elem = line.split()
                return int(line_elem[3])
        return None

    @property
    def coordinate_lines(self):
        coordinate_lines = []
        for line in self.contents:
            line_elem = line.split()
            if (
                len(line_elem) == 4
                and line_elem[0].isalpha()
                and is_float(line_elem[1])
                and is_float(line_elem[2])
                and is_float(line_elem[3])
            ):
                coordinate_lines.append(line)
        return coordinate_lines

    @property
    def natoms(self):
        return len(self.coordinate_lines)

    @property
    def list_of_symbols(self):
        """List of elements to follow periodic table with proper caps/cases."""
        elements = []
        for line in self.coordinate_lines:
            element = line.split()[0]
            element = to_element(element)
            elements.append(element)
        return elements

    @property
    def atoms(self):
        positions = []
        for line in self.coordinate_lines:
            line_elements = line.split()
            each_coord = [
                float(line_elements[1]),
                float(line_elements[2]),
                float(line_elements[3]),
            ]
            positions.append(each_coord)
        positions = np.array(positions)
        return AtomsWrapper.from_positions(
            symbols=self.list_of_symbols, positions=positions
        )

    @property
    def empirical_formula(self):
        return self.atoms.get_chemical_formula(empirical=True)

    @property
    def scf_maxiter(self):
        for i, raw_line in enumerate(self.contents):
            line = raw_line.lower()
            if "maxiter" in line:
                try:
                    # Find the index of "maxiter" in the same line
                    index = line.index("maxiter")
                    # Find the substring immediately following "maxiter"
                    num_maxiter = line[index + len("maxiter") :].split()[0]
                    return int(num_maxiter)
                except ValueError:
                    # If "maxiter" is not found, search the next line for it
                    next_lines = self.contents[i + 1 :]
                    for line in next_lines:
                        if "maxiter" in line:
                            # Find the index of "maxiter" in the same line
                            index = line.index("maxiter")
                            # Find the substring immediately following "maxiter"
                            num_maxiter = line[index + len("maxiter") :].split()[0]
                            return int(num_maxiter)
        return None

    @property
    def scf_convergence(self):
        """Within %scf block."""
        for i, line in enumerate(self.contents):
            if "%scf" not in line:
                continue

            if "convergence" in line:
                # Find the index of "convergence" in the same line
                index = line.index("convergence")
                # Find the substring immediately following "convergence"
                return line[index + len("convergence") :].split()[0]

            # If "convergence" is not found, search the next line for it
            next_lines = self.contents[i + 1 :]
            for next_line in next_lines:
                if "convergence" in next_line:
                    # Find the index of "convergence" in the same line
                    index = next_line.index("convergence")
                    # Find the substring immediately following "convergence"
                    return next_line[index + len("convergence") :].split()[0]
        return None

    @property
    def dipole(self):
        # in elprop block
        for line in self.contents:
            if "dipole" in line:
                return line.split()[-1]
        return None

    @property
    def quadrupole(self):
        # in elprop block
        for line in self.contents:
            if "quadrupole" in line:
                return line.split()[-1]
        return None

    def read_settings(self):
        from pyatoms.jobs.orca.settings import ORCAJobSettings

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
            gen_genecp=dv.gen_genecp,
            heavy_elements=dv.heavy_elements,
            heavy_elements_basis=dv.heavy_elements_basis,
            light_elements_basis=dv.light_elements_basis,
            custom_solvent=dv.custom_solvent,
            forces=dv.forces,
        )
