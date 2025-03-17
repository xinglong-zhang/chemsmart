import os
from functools import cached_property

from chemsmart.utils.mixins import FileMixin


class XTBInput(FileMixin):
    def __init__(self, filename):
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File {filename} not found.")
        self.filename = filename

    @cached_property
    def content_groups(self):
        groups = {}
        current_group = None
        for line in self.contents:
            if line.startswith("$"):
                parts = line.split(maxsplit=1)
                current_group = parts[0].strip("$")
                groups[current_group] = []
                if len(parts) > 1:
                    groups[current_group].append(parts[1].strip())
            elif current_group:
                groups[current_group].append(line)
        return groups

    def _get_key(self, group, key, default):
        if group not in self.content_groups:
            return default
        for line in self.content_groups[group]:
            parts = line.split("=", maxsplit=1)
            if len(parts) > 1 and key == parts[0].strip():
                return parts[1].strip()
        return default

    # Content group: chrg

    @property
    def charge(self):
        """Set the charge of the molecule, default is 0"""
        if "chrg" in self.content_groups:
            try:
                return int(self.content_groups["chrg"][0])
            except ValueError:
                return 0
        return 0

    # Content group: spin

    @property
    def spin(self):
        """Set the Nalpha-Nbeta of the molecule, default is 0"""
        if "spin" in self.content_groups:
            try:
                return int(self.content_groups["spin"][0])
            except ValueError:
                return 0
        return 0

    # Content group: gfn

    @property
    def method(self):
        """Set the version of the GFN Hamiltonian, default is GFN2-xTB"""
        key = self._get_key("gfn", "method", "2")
        return {"0": "GFN0-xTB", "1": "GFN1-xTB", "2": "GFN2-xTB"}.get(
            key, "GFN2-xTB"
        )

    @property
    def scc(self):
        """Switch self-consistent charge (SCC) calculation, default is True"""
        return self._get_key("gfn", "scc", "true").lower() == "true"

    @property
    def periodic(self):
        """Use periodic boundary condition, default is False"""
        return self._get_key("gfn", "periodic", "false").lower() == "true"

    @property
    def dispersion_energy_scale(self):
        """Set the scale dispersion energy of GFN-FF, default is 1.0"""
        try:
            return float(self._get_key("gfn", "dispscale", "1.0"))
        except ValueError:
            return 1.0

    # Content group: scc

    @property
    def max_iterations(self):
        """Set the maximum number of SCC iterations, default is 250"""
        try:
            return int(self._get_key("scc", "maxiterations", "250"))
        except ValueError:
            return 250

    @property
    def electronic_temperature(self):
        """Set the electronic temperature for the Fermi smearing, default is 300.0 K"""
        try:
            return float(self._get_key("scc", "temp", "300.0"))
        except ValueError:
            return 300.0

    @property
    def broyden_damping(self):
        """Set the damping for the Broyden convergence accelerator, default is 0.4"""
        try:
            return float(self._get_key("scc", "broydamp", "0.4"))
        except ValueError:
            return 0.4

    @property
    def guess_charges(self):
        """Specify the guess charges for GFN2-xTB SCC calculation"""
        return self._get_key("scc", "guess", None)

    # Content group: opt

    @property
    def engine(self):
        """Set the optimization engines, include rf (ANCopt), lbfgs (L-ANCopt) and inertial (FIRE), default is rf"""
        return self._get_key("opt", "engine", "rf")

    @property
    def optimization_level(self):
        """Specify the convergence thresholds for optimization, default is normal"""
        optlevel = self._get_key("opt", "optlevel", "normal")
        if optlevel.lower() in {
            "crude",
            "sloppy",
            "loose",
            "normal",
            "tight",
            "verytight",
            "extreme",
        }:
            return optlevel.lower()
        try:
            int_optlevel = int(optlevel)
            return {
                -3: "crude",
                -2: "sloppy",
                -1: "loose",
                0: "normal",
                1: "tight",
                2: "verytight",
                3: "extreme",
            }.get(int_optlevel, "normal")
        except ValueError:
            return "normal"

    @property
    def anc_microcycles(self):
        """Set the number of optimization cycles before new ANC are made, default is 20"""
        try:
            return int(self._get_key("opt", "microcycle", "20"))
        except ValueError:
            return 20

    @property
    def max_optcycles(self):
        """Set the total number of optimization cycles, default is 0, which means automatically determined"""
        try:
            return int(self._get_key("opt", "maxcycle", "0"))
        except ValueError:
            return 0

    @property
    def max_displacement(self):
        """Set the maximum coordinate displacement, default is 1.0"""
        try:
            return float(self._get_key("opt", "maxdispl", "1.0"))
        except ValueError:
            return 1.0

    @property
    def low_frequency_cutoff(self):
        """Set the lowest force constant in ANC generation, default is 0.01"""
        try:
            return float(self._get_key("opt", "hlow", "0.01"))
        except ValueError:
            return 0.01

    @property
    def hessian_model(self):
        """Set the model Hessian for generation of ANC used in optimization"""
        return self._get_key("opt", "hessian", "old")

    @property
    def s6_in_model_hessian(self):
        """Specify the dispersion scaling in ANC generation, default is 20.0"""
        try:
            return float(self._get_key("opt", "s6", "20.0"))
        except ValueError:
            return 20.0

    @property
    def stretch_force_constant(self):
        """Set the stretch force constant in model Hessian, default is 0.4"""
        try:
            return float(self._get_key("opt", "kstretch", "0.4"))
        except ValueError:
            return 0.4

    @property
    def bend_force_constant(self):
        """Set the bend force constant in model Hessian, default is 0.13"""
        try:
            return float(self._get_key("opt", "kbend", "0.13"))
        except ValueError:
            return 0.13

    @property
    def torsion_force_constant(self):
        """Set the torsion force constant in model Hessian, default is 0.0075"""
        try:
            return float(self._get_key("opt", "ktorsion", "0.0075"))
        except ValueError:
            return 0.0075

    @property
    def out_of_plane_force_constant(self):
        """Set the out-of-plane force constant in model Hessian, default is 0.0"""
        try:
            return float(self._get_key("opt", "koutofp", "0.0"))
        except ValueError:
            return 0.0

    @property
    def additional_vdw_contribution(self):
        """Set additional vdW-contribution, default is 0.0"""
        try:
            return float(self._get_key("opt", "kvdw", "0.0"))
        except ValueError:
            return 0.0

    @property
    def electrostatic_contribution(self):
        """Set electrostatic contribution to model Hessian by EEQ model, default is 0.0"""
        try:
            return float(self._get_key("opt", "kes", "0.0"))
        except ValueError:
            return 0.0

    @property
    def distance_cutoff(self):
        """Set the distance cutoff for bonds in model Hessian, default is 8.366600265340756"""
        try:
            return float(self._get_key("opt", "rcut", "8.366600265340756"))
        except ValueError:
            return 8.366600265340756

    @property
    def exact_rational_function(self):
        """Use better solver during the rational function optimization, default is False"""
        return self._get_key("opt", "exact rf", "false").lower() == "true"

    @property
    def average_convergence(self):
        """average the energy and gradient before checking for convergence, default is False"""
        return self._get_key("opt", "average conv", "false").lower() == "true"

    # Content group: thermo

    @property
    def temperature(self):
        """Set the temperature for thermostatistical calculation, default is 298.15 K"""
        try:
            return float(self._get_key("thermo", "temp", "298.15"))
        except ValueError:
            return 298.15

    @property
    def rotor_cutoff(self):
        """Set the rotor cut-off, default is 50.0 cm-1"""
        try:
            return float(self._get_key("thermo", "sthr", "50.0"))
        except ValueError:
            return 50.0

    @property
    def imaginary_frequency_cutoff(self):
        """Set the threshold for inverting imaginary frequencies, default is -20.0 cm-1"""
        try:
            return float(self._get_key("thermo", "imagthr", "-20.0"))
        except ValueError:
            return -20.0

    @property
    def scaling_factor(self):
        """Set the scaling factor for frequencies in vibrational partition function, default is 1.0"""
        try:
            return float(self._get_key("thermo", "scale", "1.0"))
        except ValueError:
            return 1.0
