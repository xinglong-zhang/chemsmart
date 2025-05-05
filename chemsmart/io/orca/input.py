import logging
import re

from chemsmart.io.molecules.structure import CoordinateBlock
from chemsmart.io.orca.route import ORCARoute
from chemsmart.utils.mixins import ORCAFileMixin
from chemsmart.utils.repattern import (
    orca_qm_h_bond_length_pattern,
    standard_coord_pattern,
)

logger = logging.getLogger(__name__)


class ORCAInput(ORCAFileMixin):
    def __init__(self, filename):
        self.filename = filename

        cb = CoordinateBlock(coordinate_block=self.coordinate_lines)
        self.cb = cb

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
    def route_object(self):
        try:
            route_object = ORCARoute(route_string=self.route_string)
            return route_object
        except TypeError as err:
            print(err)

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
        pattern = re.compile(standard_coord_pattern)
        coordinate_lines = []
        for line in self.contents:
            if pattern.match(line):
                coordinate_lines.append(line)
        return coordinate_lines

    @property
    def molecule(self):
        molecule = self.cb.molecule
        # update charge and spin multiplicity
        molecule.charge = self.charge
        molecule.spin_multiplicity = self.multiplicity
        return molecule

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
                            num_maxiter = line[
                                index + len("maxiter") :
                            ].split()[0]
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

class ORCAQMMMInput(ORCAInput):

    @property
    def qm_atoms(self):
        """Get QM atoms from the QMMM block."""
        return self._get_qmmm_block()["qmatoms"]

    @property
    def qm_force_field(self):
        """Get QM atoms from the QMMM block."""
        return self._get_qmmm_block()["force field"]

    @property
    def qm_functional(self):
        """Get QM functional from the QMMM block."""
        # todo: need to double check
        method_block = []
        advanced_method = False
        for line in self.contents:
            if "%method" in line.lower():
                advanced_method = True
                continue
            if advanced_method:
                method_block.append(line)
                if "end" in line.lower():
                    break
        if len(method_block) > 0:
            return method_block
        else:
            return self.functional

    @property
    def qm_active_atoms(self):
        return self._get_active_atoms()

    @property
    def qm2_atoms(self):
        return self._get_qmmm_block()["qm2atoms"]


    @property
    def qm2_level_of_theory(self):
        qm2_level_of_theory, _, _ = self._get_qm2_level_of_theory()
        if qm2_level_of_theory is not None:
            return qm2_level_of_theory

    @property
    def qm2_functional(self):
        _, qm2_functional, _ = self._get_qm2_level_of_theory()
        if qm2_functional is not None:
            return qm2_functional

    @property
    def qm2_basis(self):
        _, _, qm2_basis = self._get_qm2_level_of_theory()
        if qm2_basis is not None:
            return qm2_basis

    @property
    def qm_opt_region_fixed_atoms(self):
        if self._get_qmmm_block().get("opt region fixed atoms") is not None:
            return self._get_qmmm_block()["opt region fixed atoms"]

    @property
    def qm_h_bond_length(self):
        return self._get_h_bond_length()

    @property
    def qm_boundary_interaction(self):
        boundary_interaction, _, _ = self._get_qm_boundary_interaction()
        return boundary_interaction

    @property
    def qm_embedding_type(self):
        _, embedding_type,_ = self._get_qm_boundary_interaction()
        return embedding_type

    @property
    def qm_qm2_boundary_treatment(self):
        _, _, boundary_treatment = self._get_qm_boundary_interaction()
        return boundary_treatment

    @property
    def qm_charge(self):
        return self.charge

    @property
    def qm_multiplicity(self):
        return self.multiplicity

    @property
    def qm2_charge(self):
        qm2_charge, _, qm2_layer = self._get_qmmm_charge_and_multiplicity()
        assert qm2_layer == True, f"QM2 layer has not been specified!"
        return qm2_charge

    @property
    def qm2_multiplicity(self):
        _, qm2_multiplicity, qm2_layer = (
            self._get_qmmm_charge_and_multiplicity()
        )
        assert qm2_layer == True, f"QM2 region has not been specified!"
        return qm2_multiplicity

    @property
    def qm_total_charge(self):
        total_charge, _, qm2_layer = self._get_qmmm_charge_and_multiplicity()
        assert (
            qm2_layer == False
        ), f"Only charge of QM and medium region will be specified!"
        return total_charge

    @property
    def qm_total_multiplicity(self):
        _, total_multiplicity, qm2_layer = (
            self._get_qmmm_charge_and_multiplicity()
        )
        assert (
            qm2_layer == False
        ), f"Only multiplicity of QM and medium region will be specified!"
        return total_multiplicity

    def _get_active_atoms(self):
        active_atoms = None
        for line in self.contents:
            line = line.lower()
            if "use_qm_infofrompdb" in line and "true" in line:
                active_atoms = (
                    "Will use active atoms information from PDB file."
                )
            elif "activeatoms" in line:
                active_atoms = re.search(r"\{(.+?)\}", line).group(1).split()
        return active_atoms

    def _get_h_bond_length(self):
        h_bond_pattern = re.compile(orca_qm_h_bond_length_pattern.lower())
        h_bond = []
        for line in self.contents:
            line = line.lower()
            if "h_dist_filename" in line:
                h_bond = line.split()[-1]
                break
            elif h_bond_pattern.match(line):
                line = re.split(r"[_\s]+", line)
                atom1 = line[1]
                atom2 = line[2]
                bond_length = line[3]
                h_bond.append((atom1, atom2, bond_length))
        return h_bond

    def _get_qm_boundary_interaction(self):
        boundary_interaction = ""
        boundary_treatment = "xtb"
        embedding = "electrostatic"  # default
        for line in self.contents:
            line = line.lower()
            if "deleteladoublecounting" in line:
                if "true" in line:
                    boundary_interaction += "Will neglect bends at QM2-QM1-MM1 and torsions at QM3-QM2-QM1-MM1 boundary.\n"
                else:
                    boundary_interaction += "Will include bends at QM2-QM1-MM1 and torsions at QM3-QM2-QM1-MM1 boundary.\n"
            if "deletelabonddoublecounting" in line:
                if "true" in line:
                    boundary_interaction += (
                        "Will neglect bonds at QM1-MM1 boundary.\n"
                    )
                else:
                    boundary_interaction += (
                        "Will include bonds at QM1-MM1 boundary.\n"
                    )
            if "embedding" in line:
                embedding = line.split()[-1]
            if 'autoff qm2 method' in line:
                boundary_treatment = line.replace('autoff qm2 method', "").strip()
        return boundary_interaction, embedding, boundary_treatment

    def _get_qm2_level_of_theory(self):
        """Get QM2 level of theory from the QMMM block."""
        qm2_level_of_theory = None
        qm2_custom_functional = None
        qm2_custom_basis = None
        for line in self.contents:
            line = line.lower()
            if "qm2customfile" in line:
                qm2_level_of_theory = line.replace("qm2customfile", "").strip()
            elif "qm2custommethod" in line:
                qm2_custom_functional = line.replace(
                    "qm2custommethod", ""
                ).strip()
            elif "qm2custombasis" in line:
                qm2_custom_basis = line.replace("qm2custombasis", "").strip()
        return qm2_level_of_theory, qm2_custom_functional, qm2_custom_basis

    def _get_qmmm_charge_and_multiplicity(self):
        charge = None
        multiplicity = None
        qm2_layer = False
        for line in self.contents:
            line = line.lower()
            if "charge_total" in line:
                charge = line.split()[1]
            if "mult_total" in line:
                multiplicity = line.split()[1]
            if "charge_medium" in line:
                charge = line.split()[1]
                qm2_layer = True
            if "mult_medium" in line:
                multiplicity = line.split()[1]
        return int(charge), int(multiplicity), qm2_layer

    def _get_qmmm_block(self):
        # todo: need to refactor this
        "return QMMM block as a dictionary"
        block = {}
        for line in self.contents:
            line = line.lower()
            if "qmatoms" in line:
                qm_atoms = re.search(r"\{(.+?)\}", line).group(1).split()
                block["qmatoms"] = qm_atoms
            if "qm2atoms" in line:
                qm2_atoms = re.search(r"\{(.+?)\}", line).group(1).split()
                block["qm2atoms"] = qm2_atoms
            if "optregion_fixedatoms" in line:
                opt_region_fixed_atoms = (
                    re.search(r"\{(.+?)\}", line).group(1).split()
                )
                block["opt region fixed atoms"] = opt_region_fixed_atoms
            if "ORCAFFFilename" in line:
                force_field = line.split()[-1]
                block["force field"] = force_field

        return block
