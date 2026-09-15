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

    def _get_key(self, group, key, default=None):
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

    # Content group: cube

    @property
    def cube_step(self):
        """Specify the grid spacing for cube file, default is 0.4"""
        try:
            return float(self._get_key("cube", "step", "0.4"))
        except ValueError:
            return 0.4

    @property
    def density_matrix_threshold(self):
        """Set the density matrix neglect threshold, default is 0.05"""
        try:
            return float(self._get_key("cube", "pthr", "0.05"))
        except ValueError:
            return 0.05

    @property
    def boundary_offset(self):
        """Set the grid boundary offset, default is 3.0"""
        try:
            return float(self._get_key("cube", "boff", "3.0"))
        except ValueError:
            return 3.0

    @property
    def cube_output(self):
        """Set cube output, 0 = writing molden file, 1 = writing cube file (default), 2 = no file output"""
        try:
            return int(self._get_key("cube", "cal", "1"))
        except ValueError:
            return 1

    # Content group: embedding

    @property
    def atom_type(self):
        """Set the default atom type for point charges"""
        at = self._get_key("embedding", "at")
        if at is not None:
            try:
                return int(at)
            except ValueError:
                return None
        return None

    @property
    def isotropic_electrostatic(self):
        """Determine whether to use isotropic electrostatic with point charges, default is True"""
        return self._get_key("embedding", "es", "true").lower() == "true"

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

    # Content group: hess

    @property
    def hess_scc_accuracy(self):
        """Set SCC accuracy level in Hessian runs, default is 0.3"""
        try:
            return float(self._get_key("hess", "sccacc", "0.3"))
        except ValueError:
            return 0.3

    @property
    def hess_step(self):
        """Set the Cartesian displacement increment for numerical Hessian, default is 0.005"""
        try:
            return float(self._get_key("hess", "step", "0.005"))
        except ValueError:
            return 0.005

    @property
    def hess_scale(self):
        """Set the scaling factor for the hessian elements, default is 1.0"""
        try:
            return float(self._get_key("hess", "scale", "1.0"))
        except ValueError:
            return 1.0

    # Content group: md

    @property
    def md_temperature(self):
        """Set the temperature for MD thermostat or GBSA solvation in K, default is 298.15 K"""
        try:
            return float(self._get_key("md", "temp", "298.15"))
        except ValueError:
            return 298.15

    @property
    def md_time(self):
        """Set the MD run time in ps, default is 50.0 ps"""
        try:
            return float(self._get_key("md", "time", "50.0"))
        except ValueError:
            return 50.0

    @property
    def dump_structure(self):
        """Set the interval for trajectory printout in fs, default is 50 fs"""
        try:
            return float(self._get_key("md", "dump", "50.0"))
        except ValueError:
            return 50.0

    @property
    def velocity_in_trj(self):
        """Determine whether velocities should be included in the trajectory (trj) file, default is False"""
        return int(self._get_key("md", "velo", "0")) == 1

    @property
    def nvt_ensemble(self):
        """Specify whether to use a thermostat (perform simulation in NVT ensemble), or the NVE ensemble is used, default is True"""
        return int(self._get_key("md", "nvt", "1")) == 1

    @property
    def skip_interval(self):
        """Define the skip interval for -mdav and -mdopt calculations, default is 500"""
        try:
            return int(self._get_key("md", "skip", "500"))
        except ValueError:
            return 500

    @property
    def md_step(self):
        """Set MD time step in fs, default is 4.0 fs"""
        try:
            return float(self._get_key("md", "step", "4.0"))
        except ValueError:
            return 4.0

    @property
    def hydrogen_mass(self):
        """Set hydrogen mass to this value in amu, default is 4 amu"""
        try:
            return int(self._get_key("md", "hmass", "4"))
        except ValueError:
            return 4

    @property
    def shake_algorithm(self):
        """Control SHAKE algorithm to constrain bonds, 0 = off, 1 = X-H only, 2 = all bonds (default)"""
        try:
            return int(self._get_key("md", "shake", "2"))
        except ValueError:
            return 2

    @property
    def md_scc_accuracy(self):
        """Set SCC accuracy level in MD, default is 2.0"""
        try:
            return float(self._get_key("md", "sccacc", "2.0"))
        except ValueError:
            return 2.0

    @property
    def force_writing_restart(self):
        """Determine whether to force the writing of a restart file at each dump step, default is False"""
        return int(self._get_key("md", "forcewrrestart", "0")) == 1

    # Content group: modef

    @property
    def modef_n(self):
        """Specify the number of points along the normal mode path scan"""
        n = self._get_key("modef", "n")
        if n is not None:
            try:
                return int(n)
            except ValueError:
                return None
        return None

    @property
    def modef_step(self):
        """Set the step lengths for scan, default is 1.0"""
        try:
            return float(self._get_key("modef", "step", "1.0"))
        except ValueError:
            return 1.0

    @property
    def modef_update(self):
        """Update the search mode by a fraction of the displacement at every step, 0.0 means no update, default is 0.2"""
        try:
            return float(self._get_key("modef", "updat", "0.2"))
        except ValueError:
            return 0.2

    @property
    def modef_local(self):
        """Specify the type of normal modes used, 0 = canonical normal modes (default), 1 = Pipek-Mezey localized modes"""
        try:
            return int(self._get_key("modef", "local", "0"))
        except ValueError:
            return 0

    @property
    def modef_threshold(self):
        """Set the frequency threshold up to which frequency modes are used for mode based conformer search, default is 0.0"""
        try:
            return float(self._get_key("modef", "vthr", "0.0"))
        except ValueError:
            return 0.0

    @property
    def projected_mode(self):
        """Set the number of second mode which should be projected out in mode following, default is 0"""
        try:
            return int(self._get_key("modef", "prj", "0"))
        except ValueError:
            return 0

    @property
    def mode_following(self):
        """Specify the mode number to follow"""
        mode = self._get_key("modef", "mode")
        if mode is not None:
            try:
                return int(mode)
            except ValueError:
                return None
        return None

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

    # Content group: path

    @property
    def pathfinder_runs(self):
        """Set the number of runs for pathfinder, default is 3"""
        try:
            return int(self._get_key("path", "nrun", "3"))
        except ValueError:
            return 3

    @property
    def path_points(self):
        """Set the number of points on the path to optimize, default is 50"""
        try:
            return int(self._get_key("path", "npoint", "50"))
        except ValueError:
            return 50

    @property
    def path_optimization_steps(self):
        """Set the number of steps to optimize the points on the path, default is 3"""
        try:
            return int(self._get_key("path", "anopt", "3"))
        except ValueError:
            return 3

    @property
    def rmsd_push_factor(self):
        """Set the factor for RMSD criterium pushing away from the reactant structure, default is 0.05"""
        try:
            return float(self._get_key("path", "kpush", "0.05"))
        except ValueError:
            return 0.05

    @property
    def rmsd_pull_factor(self):
        """Set the factor for RMSD criterium pulling towards the product structure, default is -0.04"""
        try:
            return float(self._get_key("path", "kpull", "-0.04"))
        except ValueError:
            return -0.04

    @property
    def rmsd_width(self):
        """Set the width of the RMSD criterium, default is 0.7"""
        try:
            return float(self._get_key("path", "alp", "0.7"))
        except ValueError:
            return 0.7

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
        """Set the electronic temperature for the Fermi smearing in K, default is 300.0 K"""
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
        return self._get_key("scc", "guess")

    # Content group: symmetry

    @property
    def symmetry_threshold(self):
        """Set the point group symmetrization threshold, default is 0.1"""
        try:
            return float(self._get_key("symmetry", "desy", "0.1"))
        except ValueError:
            return 0.1

    @property
    def symmetry_max_atoms(self):
        """Set the maximum number of atoms for point group determination.
        If the number of atoms exceeds this value, symmetry determination is skipped.
        Default is 200, 0 switches it off"""
        try:
            return int(self._get_key("symmetry", "maxat", "200"))
        except ValueError:
            return 200

    # Content group: thermo

    @property
    def thermo_temperature(self):
        """Set the temperature for thermostatistical calculation in K, default is 298.15 K"""
        try:
            return float(self._get_key("thermo", "temp", "298.15"))
        except ValueError:
            return 298.15

    @property
    def rotor_cutoff(self):
        """Set the rotor cut-off in cm-1, default is 50.0 cm-1"""
        try:
            return float(self._get_key("thermo", "sthr", "50.0"))
        except ValueError:
            return 50.0

    @property
    def imaginary_frequency_cutoff(self):
        """Set the threshold for inverting imaginary frequencies in cm-1, default is -20.0 cm-1"""
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

    # Content group: wall

    @property
    def wall_potential(self):
        """Set kind of wall potential used, default is polynomial"""
        return self._get_key("wall", "potential", "polynomial")

    @property
    def wall_potential_exponent(self):
        """Specify the exponent of polynomial wall potential, default is 30"""
        try:
            return int(self._get_key("wall", "alpha", "30"))
        except ValueError:
            return 30

    @property
    def logfermi_bias_exponent(self):
        """Specify the exponent of logfermi bias potential, default is 6.0"""
        try:
            return float(self._get_key("wall", "beta", "6.0"))
        except ValueError:
            return 6.0

    @property
    def wall_temperature(self):
        """Set the temperature of the logfermi wall, default is 300.0 K"""
        try:
            return float(self._get_key("wall", "temp", "300.0"))
        except ValueError:
            return 300.0

    @property
    def auto_scale(self):
        """Set the scales axis of automatic determined wall potentials by real, default is 1.0"""
        try:
            return float(self._get_key("wall", "autoscale", "1.0"))
        except ValueError:
            return 1.0

    @property
    def axis_shift(self):
        """Set the constant offset used in automatic dermined wall potential axis, default is 3.5"""
        try:
            return float(self._get_key("wall", "axisshift", "3.5"))
        except ValueError:
            return 3.5
