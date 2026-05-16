"""
Gaussian DI-AS (Distortion interaction-Activation strain)
calculation job implementation.

This module provides the GaussianDIASJob class for performing
DI-AS analysis calculations using Gaussian.
These calculations analyze molecular fragmentation along reaction
coordinates by computing energies of fragments and whole molecules.
"""

from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.utils.utils import get_list_from_string_range


class GaussianDIASJob(GaussianJob):
    """
    Gaussian job class for DI-AS fragmentation analysis calculations.

    Performs distortion interaction-Activation strain (DI-AS) analysis
    by calculating energies of molecular fragments and complete molecules
    along a reaction coordinate. Supports both IRC and transition state
    modes with configurable fragment definitions and sampling.

    The calculation involves three sets of jobs:
    1. Complete molecule calculations at sampled points
    2. Fragment 1 calculations at the same points
    3. Fragment 2 calculations at the same points

    Attributes:
        TYPE (str): Job type identifier ('g16dias').
        all_molecules (list): Molecule objects along reaction coordinate.
        fragment_indices (str): String range defining fragment 1 atoms.
        every_n_points (int): Sampling frequency along coordinate.
        mode (str): Calculation mode ('irc' or 'ts').
        charge_of_fragment1 (int, optional): Charge for fragment 1.
        multiplicity_of_fragment1 (int, optional): Multiplicity for fragment 1.
        charge_of_fragment2 (int, optional): Charge for fragment 2.
        multiplicity_of_fragment2 (int, optional): Multiplicity for fragment 2.
        fragment1_settings (GaussianJobSettings): Settings used for fragment 1
            calculations (derived from base settings).
        fragment2_settings (GaussianJobSettings): Settings used for fragment 2
            calculations (derived from base settings).
    """

    TYPE = "g16dias"

    def __init__(
        self,
        molecules,
        settings,
        label,
        jobrunner,
        fragment_indices,
        every_n_points,
        mode,
        charge_of_fragment1=None,
        multiplicity_of_fragment1=None,
        charge_of_fragment2=None,
        multiplicity_of_fragment2=None,
        **kwargs,
    ):
        """
        Initialize a Gaussian DI-AS fragmentation analysis calculation.

        Sets up DI-AS analysis with molecular structures along a reaction
        coordinate, fragment definitions, and separate settings for each
        fragment. Configures sampling frequency and calculation mode.

        Args:
            molecules (list): Molecule objects along reaction coordinate.
            settings (GaussianJobSettings): Base calculation settings.
            label (str): Base label for job identification.
            jobrunner: Job execution handler.
            fragment_indices (str): String range defining fragment 1 atoms
                (e.g., "1-10", "1,3,5-8").
            every_n_points (int): Sample every nth point along coordinate.
            mode (str): Calculation mode ('irc' or 'ts').
            charge_of_fragment1 (int, optional): Charge for fragment 1.
            multiplicity_of_fragment1 (int, optional): Multiplicity for
                fragment 1.
            charge_of_fragment2 (int, optional): Charge for fragment 2.
            multiplicity_of_fragment2 (int, optional): Multiplicity for
                fragment 2.
            **kwargs: Additional keyword arguments for parent class.
        """
        super().__init__(
            molecule=molecules[0],  # Use the first molecule as a placeholder
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
        self.all_molecules = molecules  # alone IRC coordinate
        self.fragment_indices = fragment_indices
        self.every_n_points = every_n_points
        self.settings.freq = False  # turn off freq calc for DI-AS
        self.mode = mode
        self.charge_of_fragment1 = charge_of_fragment1
        self.multiplicity_of_fragment1 = multiplicity_of_fragment1
        self.charge_of_fragment2 = charge_of_fragment2
        self.multiplicity_of_fragment2 = multiplicity_of_fragment2

        fragment1_settings = self.settings.copy()
        fragment2_settings = self.settings.copy()

        if self.charge_of_fragment1 is not None:
            fragment1_settings.charge = self.charge_of_fragment1
        if self.multiplicity_of_fragment1 is not None:
            fragment1_settings.multiplicity = self.multiplicity_of_fragment1

        if self.charge_of_fragment2 is not None:
            fragment2_settings.charge = self.charge_of_fragment2
        if self.multiplicity_of_fragment2 is not None:
            fragment2_settings.multiplicity = self.multiplicity_of_fragment2

        self.fragment1_settings = fragment1_settings
        self.fragment2_settings = fragment2_settings

    @property
    def num_molecules(self):
        """
        Get the total number of molecules along the reaction coordinate.

        Returns:
            int: Number of molecular structures in the coordinate series.
        """
        return len(self.all_molecules)

    def _fragment_structure(self, molecule):
        """
        Split a molecule into two fragments based on atom indices.

        Divides the molecular structure into fragment 1 (specified by
        fragment_indices) and fragment 2 (remaining atoms). This
        fragmentation is used for DI-AS energy analysis.

        Args:
            molecule: Molecule object to fragment.

        Returns:
            tuple: (fragment1_atoms, fragment2_atoms) as separate
                molecular structures.
        """
        fragment1_indices = get_list_from_string_range(self.fragment_indices)
        fragment2_indices = [
            i + 1
            for i in range(len(molecule))
            if i + 1 not in fragment1_indices
        ]

        fragment1_atoms = molecule[fragment1_indices]
        fragment2_atoms = molecule[fragment2_indices]
        return fragment1_atoms, fragment2_atoms

    @property
    def fragment1_atoms(self):
        """
        Get fragment 1 structures for all molecules along coordinate.

        Returns:
            list: Fragment 1 molecular structures for each point
                along the reaction coordinate.
        """
        return [
            self._fragment_structure(molecule=molecule)[0]
            for molecule in self.all_molecules
        ]

    @property
    def fragment2_atoms(self):
        """
        Get fragment 2 structures for all molecules along coordinate.

        Returns:
            list: Fragment 2 molecular structures for each point
                along the reaction coordinate.
        """
        return [
            self._fragment_structure(molecule=molecule)[1]
            for molecule in self.all_molecules
        ]

    def _sample_molecules(self, molecules):
        """
        Sample molecules at regular intervals along the coordinate.

        Selects every nth molecule according to every_n_points setting,
        ensuring the last molecule is always included for proper
        endpoint coverage.

        Args:
            molecules (list): Complete list of molecular structures.

        Returns:
            list: Sampled subset of molecular structures.
        """
        filtered_molecules = molecules[0 :: self.every_n_points]
        if (self.num_molecules - 1) / self.every_n_points != 0:
            filtered_molecules.append(molecules[-1])
        return filtered_molecules

    @property
    def fragment1_jobs(self):
        """
        Generate calculation jobs for fragment 1 structures.

        Creates GaussianGeneralJob objects for fragment 1 calculations
        based on the specified mode. For IRC mode, samples multiple
        points along the coordinate. For TS mode, uses only the final
        (transition state) structure.

        Returns:
            list: GaussianGeneralJob objects for fragment 1 calculations.

        Raises:
            ValueError: If mode is not 'irc' or 'ts'.
        """
        if self.mode.lower() == "irc":
            # using IRC log file
            molecules = self._sample_molecules(self.fragment1_atoms)
            jobs = []
            for i, molecule in enumerate(molecules):
                label = f"{self.label}_p{i}_f1"
                jobs += [
                    GaussianGeneralJob(
                        molecule=molecule,
                        settings=self.fragment1_settings,
                        label=label,
                        jobrunner=self.jobrunner,
                        skip_completed=self.skip_completed,
                    )
                ]
            return jobs
        elif self.mode.lower() == "ts":
            # using TS log file
            molecule = self.fragment1_atoms[-1]
            label = f"{self.label}_p1_f1"
            return [
                GaussianGeneralJob(
                    molecule=molecule,
                    settings=self.fragment1_settings,
                    label=label,
                    jobrunner=self.jobrunner,
                    skip_completed=self.skip_completed,
                )
            ]
        else:
            raise ValueError(
                f"Invalid mode: {self.mode}. Must be 'irc' or 'ts'."
            )

    @property
    def fragment2_jobs(self):
        """
        Generate calculation jobs for fragment 2 structures.

        Creates GaussianGeneralJob objects for fragment 2 calculations
        using the same sampling strategy as fragment 1. The mode
        determines whether to use IRC sampling or TS endpoint.

        Returns:
            list: GaussianGeneralJob objects for fragment 2 calculations.

        Raises:
            ValueError: If mode is not 'irc' or 'ts'.
        """
        if self.mode.lower() == "irc":
            molecules = self._sample_molecules(self.fragment2_atoms)
            jobs = []
            for i, molecule in enumerate(molecules):
                label = f"{self.label}_p{i}_f2"
                jobs += [
                    GaussianGeneralJob(
                        molecule=molecule,
                        settings=self.fragment2_settings,
                        label=label,
                        jobrunner=self.jobrunner,
                        skip_completed=self.skip_completed,
                    )
                ]
            return jobs
        elif self.mode.lower() == "ts":
            molecule = self.fragment2_atoms[-1]
            label = f"{self.label}_p1_f2"
            return [
                GaussianGeneralJob(
                    molecule=molecule,
                    settings=self.fragment2_settings,
                    label=label,
                    jobrunner=self.jobrunner,
                    skip_completed=self.skip_completed,
                )
            ]
        else:
            raise ValueError(
                f"Invalid mode: {self.mode}. Must be 'irc' or 'ts'."
            )

    @property
    def all_molecules_jobs(self):
        """
        Generate calculation jobs for complete molecular structures.

        Creates GaussianGeneralJob objects for the complete molecules
        (before fragmentation) at sampled points along the reaction
        coordinate. These provide reference energies for DI-AS analysis.

        Returns:
            list: GaussianGeneralJob objects for complete molecule
                calculations.
        """
        if self.mode.lower() == "irc":
            molecules = self._sample_molecules(self.all_molecules)
            jobs = []
            for i, molecule in enumerate(molecules):
                label = f"{self.label}_p{i}"
                jobs += [
                    GaussianGeneralJob(
                        molecule=molecule,
                        settings=self.settings,
                        label=label,
                        jobrunner=self.jobrunner,
                        skip_completed=self.skip_completed,
                    )
                ]
            return jobs
        elif self.mode.lower() == "ts":
            molecule = self.all_molecules[-1]
            label = f"{self.label}_p1"
            return [
                GaussianGeneralJob(
                    molecule=molecule,
                    settings=self.settings,
                    label=label,
                    jobrunner=self.jobrunner,
                    skip_completed=self.skip_completed,
                )
            ]

    def _run_all_molecules_jobs(self):
        """
        Execute all complete molecule calculation jobs.

        Runs all jobs for complete molecular structures at sampled
        points along the reaction coordinate. These calculations
        provide reference energies for DI-AS analysis.
        """
        for job in self.all_molecules_jobs:
            job.run()

    def _run_fragment1_jobs(self):
        """
        Execute all fragment 1 calculation jobs.

        Runs all jobs for fragment 1 structures at sampled points.
        Fragment 1 energies are used with fragment 2 and complete
        molecule energies to compute dissociation energies.
        """
        for job in self.fragment1_jobs:
            job.run()

    def _run_fragment2_jobs(self):
        """
        Execute all fragment 2 calculation jobs.

        Runs all jobs for fragment 2 structures at sampled points.
        These calculations complete the energy data needed for
        DI-AS fragmentation analysis.
        """
        for job in self.fragment2_jobs:
            job.run()

    def _run(self, **kwargs):
        """
        Execute the complete DI-AS calculation workflow.

        Runs all three sets of calculations in sequence:
        1. Complete molecules at sampled points
        2. Fragment 1 structures at the same points
        3. Fragment 2 structures at the same points

        Args:
            **kwargs: Additional keyword arguments for job execution.
        """
        self._run_all_molecules_jobs()
        self._run_fragment1_jobs()
        self._run_fragment2_jobs()

    def is_complete(self):
        """
        Check if all DI-AS calculation jobs are complete.

        Verifies that all three sets of calculations (complete molecules,
        fragment 1, and fragment 2) have finished successfully.

        Returns:
            bool: True if all jobs are complete, False otherwise.
        """
        return (
            self._run_all_molecules_jobs_are_complete()
            and self._run_fragment1_jobs_are_complete()
            and self._run_fragment2_jobs_are_complete()
        )

    def _run_all_molecules_jobs_are_complete(self):
        """
        Check completion status of all complete molecule jobs.

        Returns:
            bool: True if all complete molecule calculations are
                finished, False otherwise.
        """
        return all(job.is_complete() for job in self.all_molecules_jobs)

    def _run_fragment1_jobs_are_complete(self):
        """
        Check completion status of all fragment 1 jobs.

        Returns:
            bool: True if all fragment 1 calculations are finished,
                False otherwise.
        """
        return all(job.is_complete() for job in self.fragment1_jobs)

    def _run_fragment2_jobs_are_complete(self):
        """
        Check completion status of all fragment 2 jobs.

        Returns:
            bool: True if all fragment 2 calculations are finished,
                False otherwise.
        """
        return all(job.is_complete() for job in self.fragment2_jobs)
