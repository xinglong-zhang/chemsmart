from functools import cached_property

from chemsmart.utils.constants import kcal_per_mol_to_hartree
from chemsmart.utils.mixins import CRESTFileMixin, FileMixin


class CRESTMainOut(CRESTFileMixin):
    """
    Parser for the main CREST output log file.

    Extracts run results from the CREST stdout log including:
    - Termination status and quality indicators
    - Number of conformers / rotamers found
    - Final ensemble thermodynamics
    - CREGEN filtering parameters
    - Energy summary and wall/CPU time

    Args:
        filename (str): Path to the CREST output file (e.g., crest.out, water_conformers.out)
    """

    def __init__(self, filename):
        self.filename = filename

    @property
    def normal_termination(self):
        """Whether the output contains 'CREST terminated normally.'"""
        for line in reversed(self.contents):
            if "CREST terminated normally." in line:
                return True
        return False

    def _get_route(self):
        """Extract the CREST command line from the log.

        Looks for the '$ crest ...' or '> crest ...' line in the
        "Command line input:" section.
        """
        for line in self.contents:
            stripped = line.strip()
            if stripped.startswith("$ crest "):
                return stripped[2:]
            if stripped.startswith("> crest "):
                return stripped[2:]
        return None

    @cached_property
    def topology_mismatch(self):
        """Whether a topology change was detected during the run.

        CREST prints '*WARNING* Change in topology detected!' when
        the optimized geometry has a different connectivity than the input.
        """
        for line in self.contents:
            if "Change in topology detected" in line:
                return True
        return False

    @cached_property
    def num_conformers(self):
        """Number of unique conformers found."""
        for line in reversed(self.contents):
            if "number of unique conformers for further calc" in line:
                parts = line.strip().split()
                try:
                    return int(parts[-1])
                except (ValueError, IndexError):
                    continue
        return None

    @property
    def single_conformer(self):
        """Whether only one unique conformer was found."""
        return self.normal_termination and self.num_conformers == 1

    @cached_property
    def num_rotamers(self):
        """Total number of unique rotamers (pre-CREGEN unique points)."""
        for line in reversed(self.contents):
            if "total number unique points considered further" in line:
                parts = line.strip().split()
                try:
                    return int(parts[-1])
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def num_atoms(self):
        """Number of atoms."""
        for line in reversed(self.contents):
            if "number of atoms" in line and ":" in line:
                try:
                    return int(line.split(":")[-1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def rmsd_threshold(self):
        """CREGEN RMSD threshold in Angstrom."""
        for line in reversed(self.contents):
            if "RMSD threshold" in line and ":" in line and "(Ang" not in line:
                try:
                    return float(line.split(":")[-1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def bconst_threshold(self):
        """CREGEN rotational-constant threshold."""
        for line in reversed(self.contents):
            if "Bconst threshold" in line and ":" in line:
                try:
                    return float(line.split(":")[-1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def population_threshold(self):
        """CREGEN Boltzmann population threshold."""
        for line in reversed(self.contents):
            if "population threshold" in line and ":" in line:
                try:
                    return float(line.split(":")[-1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def energy_window(self):
        """CREGEN sorting energy window (EWIN) in kcal/mol."""
        for line in reversed(self.contents):
            if "sorting energy window" in line and ":" in line:
                try:
                    return float(line.split(":")[-1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def reference_state_energy(self):
        """Reference-state total energy (Etot) in Hartree."""
        for line in reversed(self.contents):
            if "reference state Etot" in line and ":" in line:
                try:
                    return float(line.split(":")[-1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def temperature(self):
        """CREGEN temperature in Kelvin."""
        for line in reversed(self.contents):
            if line.strip().startswith("T /K") and ":" in line:
                try:
                    return float(line.split(":")[-1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def lowest_energy(self):
        """Lowest energy in Hartree."""
        for line in reversed(self.contents):
            if line.strip().startswith("E lowest") and ":" in line:
                value_str = line.split(":")[-1].strip()
                try:
                    return float(value_str.split()[0])
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def average_energy_ensemble(self):
        """Absolute Boltzmann-averaged ensemble energy in Hartree.

        CREST prints a relative ensemble average energy (kcal/mol) with
        respect to the CREGEN reference state. This returns
        E_ensemble = E_reference + <E_rel>
        """
        if self.reference_state_energy is None:
            return None
        for line in reversed(self.contents):
            if "ensemble average energy" in line and ":" in line:
                try:
                    relative_hartree = (
                        float(line.split(":")[-1].strip().split()[0])
                        * kcal_per_mol_to_hartree
                    )
                    return self.reference_state_energy + relative_hartree
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def entropy_ensemble(self):
        """Ensemble entropy, S_ensemble, in Hartree/K."""
        for line in reversed(self.contents):
            if "ensemble entropy" in line and ":" in line:
                try:
                    entropy_hartree = (
                        float(line.split(":")[-1].strip().split()[-1])
                        * 1e-3
                        * kcal_per_mol_to_hartree
                    )
                    return entropy_hartree
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def free_energy_ensemble(self):
        """Absolute ensemble free energy in Hartree.

        CREST reports the ensemble free-energy contribution (-T*S)
        (kcal/mol). This returns
        G_ensemble = E_reference + <E_rel> - T*S_ensemble.
        """
        if self.average_energy_ensemble is None:
            return None
        for line in reversed(self.contents):
            if "ensemble free energy" in line and ":" in line:
                try:
                    relative_free_hartree = (
                        float(line.split(":")[-1].strip().split()[0])
                        * kcal_per_mol_to_hartree
                    )
                    return self.average_energy_ensemble + relative_free_hartree
                except (ValueError, IndexError):
                    continue
        return None

    @cached_property
    def population_of_lowest(self):
        """Boltzmann population of the lowest conformer in percent."""
        for line in reversed(self.contents):
            if "population of lowest" in line and ":" in line:
                try:
                    return float(line.split(":")[-1].strip().split()[0])
                except (ValueError, IndexError):
                    continue
        return None

    @staticmethod
    def sum_time_hours(line):
        """Parse a time string from CREST output and return the total time in hours."""
        n_days = float(line.split(" d,")[0].split()[-1])
        n_hours = float(line.split(" h,")[0].split()[-1])
        n_minutes = float(line.split(" min,")[0].split()[-1])
        n_seconds = float(line.split(" sec")[0].split()[-1])
        total_seconds = (
            n_days * 24 * 60 * 60
            + n_hours * 60 * 60
            + n_minutes * 60
            + n_seconds
        )
        total_hours = total_seconds / 3600
        return total_hours

    @cached_property
    def wall_time(self):
        """Total wall time in hours."""
        for line in reversed(self.contents):
            if "wall-time:" in line:
                return self.sum_time_hours(line)
        return None

    @cached_property
    def cpu_time(self):
        """Total CPU time in hours."""
        for line in reversed(self.contents):
            if "cpu-time:" in line:
                return self.sum_time_hours(line)
        return None


class CRESTEnergiesFile(FileMixin):
    """
    Parser for the crest.energies file.

    This file contains relative energies (in kcal/mol) for all unique
    conformers, one per line in the format:
        INDEX    RELATIVE_ENERGY_KCAL
    """

    def __init__(self, filename):
        self.filename = filename

    @cached_property
    def relative_energies(self):
        """Relative conformer energies in kcal/mol.

        Returns:
            list[float]: Relative energies for each conformer.
        """
        energies = []
        for line in self.contents:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    energies.append(float(parts[1]))
                except ValueError:
                    continue
        return energies

    @property
    def num_conformers(self):
        """Number of conformers in the energies file."""
        if self.relative_energies:
            return len(self.relative_energies)
        return None

    @property
    def energy_range(self):
        """Energy spread (max - min) in kcal/mol."""
        if not self.relative_energies:
            return None
        return max(self.relative_energies) - min(self.relative_energies)
