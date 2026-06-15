"""
Gaussian DI-AS (Distortion interaction-Activation strain)
calculation job implementation.

This module provides the GaussianDIASJob class for performing
DI-AS analysis calculations using Gaussian.
These calculations analyze molecular fragmentation along reaction
coordinates by computing energies of fragments and whole molecules.
"""

import logging
import re

from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.qrc import GaussianQRCJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.utils.utils import get_list_from_string_range

logger = logging.getLogger(__name__)


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
        reactant_opt_settings=None,
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
        if reactant_opt_settings is None:
            reactant_opt_settings = self.settings
        # Normalize to GaussianJobSettings to avoid carrying IRC-specific
        # route generation logic into fragment optimization jobs.
        reactant_opt_settings = GaussianJobSettings(
            **reactant_opt_settings.__dict__
        )

        self.fragment1_reactant_opt_settings = reactant_opt_settings.copy()
        self.fragment1_reactant_opt_settings.charge = self.fragment1_settings.charge
        self.fragment1_reactant_opt_settings.multiplicity = (
            self.fragment1_settings.multiplicity
        )
        self.fragment1_reactant_opt_settings.jobtype = "opt"
        self.fragment1_reactant_opt_settings.freq = True
        self.fragment1_reactant_sp_settings = self.fragment1_settings.copy()
        self.fragment1_reactant_sp_settings.jobtype = "sp"
        self.fragment1_reactant_sp_settings.freq = False
        self.fragment2_reactant_opt_settings = reactant_opt_settings.copy()
        self.fragment2_reactant_opt_settings.charge = self.fragment2_settings.charge
        self.fragment2_reactant_opt_settings.multiplicity = (
            self.fragment2_settings.multiplicity
        )
        self.fragment2_reactant_opt_settings.jobtype = "opt"
        self.fragment2_reactant_opt_settings.freq = True
        self.fragment2_reactant_sp_settings = self.fragment2_settings.copy()
        self.fragment2_reactant_sp_settings.jobtype = "sp"
        self.fragment2_reactant_sp_settings.freq = False

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

    @property
    def fragment1_reactant_opt_job(self):
        label = f"{self.label}_fragment1_opt"
        return GaussianOptJob(
            molecule=self.fragment1_atoms[-1],
            settings=self.fragment1_reactant_opt_settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    @property
    def fragment2_reactant_opt_job(self):
        label = f"{self.label}_fragment2_opt"
        return GaussianOptJob(
            molecule=self.fragment2_atoms[-1],
            settings=self.fragment2_reactant_opt_settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    @property
    def fragment1_reactant_sp_job(self):
        molecule = self.fragment1_reactant_opt_job.optimized_structure()
        if molecule is None:
            molecule = self.fragment1_atoms[-1]
        label = f"{self.label}_fragment1_r1"
        return GaussianGeneralJob(
            molecule=molecule,
            settings=self.fragment1_reactant_sp_settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    @property
    def fragment2_reactant_sp_job(self):
        molecule = self.fragment2_reactant_opt_job.optimized_structure()
        if molecule is None:
            molecule = self.fragment2_atoms[-1]
        label = f"{self.label}_fragment2_i2"
        return GaussianGeneralJob(
            molecule=molecule,
            settings=self.fragment2_reactant_sp_settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    @staticmethod
    def _contains_imaginary_frequencies(opt_job):
        output = opt_job._output()
        if output is None or not output.normal_termination:
            return True
        frequencies = output.vibrational_frequencies
        if not frequencies:
            return True
        return any(freq < 0.0 for freq in frequencies)

    def _retry_reactant_opt_job(self, opt_job):
        retry_settings = opt_job.settings.copy()
        retry_route_string = retry_settings.route_string

        opt_options_match = re.search(
            r"\bopt\s*=\s*\(([^)]*)\)", retry_route_string, flags=re.IGNORECASE
        )
        if opt_options_match is not None:
            opt_options = [
                option.strip()
                for option in opt_options_match.group(1).split(",")
                if option.strip()
            ]
            if not any(option.lower() == "maxstep=5" for option in opt_options):
                opt_options.append("maxstep=5")
                retry_route_string = (
                    f"{retry_route_string[:opt_options_match.start()]}"
                    f"opt=({','.join(opt_options)})"
                    f"{retry_route_string[opt_options_match.end():]}"
                )
        else:
            opt_keyword_match = re.search(
                r"\bopt(?:\s*=\s*(?!\()([^\s]+))?\b",
                retry_route_string,
                flags=re.IGNORECASE,
            )
            if opt_keyword_match is not None:
                existing_opt = opt_keyword_match.group(1)
                if existing_opt is None:
                    replacement = "opt=(maxstep=5)"
                elif existing_opt.lower() == "maxstep=5":
                    replacement = "opt=(maxstep=5)"
                else:
                    replacement = f"opt=({existing_opt},maxstep=5)"
                retry_route_string = (
                    f"{retry_route_string[:opt_keyword_match.start()]}"
                    f"{replacement}"
                    f"{retry_route_string[opt_keyword_match.end():]}"
                )

        if re.search(r"\bscf\s*=\s*qc\b", retry_route_string, flags=re.IGNORECASE):
            retry_settings.route_to_be_written = retry_route_string
        else:
            retry_settings.route_to_be_written = (
                f"{retry_route_string} scf=qc"
            )

        retry_settings.additional_opt_options_in_route = None
        retry_settings.additional_route_parameters = None
        return GaussianOptJob(
            molecule=opt_job.molecule,
            settings=retry_settings,
            label=opt_job.label,
            jobrunner=opt_job.jobrunner,
            skip_completed=opt_job.skip_completed,
        )

    @staticmethod
    def _first_imaginary_mode_index(job):
        output = job._output()
        if output is None or not output.normal_termination:
            return 1
        frequencies = output.vibrational_frequencies or []
        for idx, freq in enumerate(frequencies, start=1):
            if freq < 0.0:
                return idx
        return 1

    def _run_qrc_for_reactant_opt_job(self, opt_job):
        mode_idx = self._first_imaginary_mode_index(opt_job)
        qrc_label = f"{opt_job.label}_qrc"
        qrc_job = GaussianQRCJob(
            molecule=opt_job.optimized_structure() or opt_job.molecule,
            settings=opt_job.settings.copy(),
            label=qrc_label,
            jobrunner=opt_job.jobrunner,
            mode_idx=mode_idx,
            skip_completed=opt_job.skip_completed,
        )
        logger.info(
            f"Running QRC fallback for {opt_job.label} with mode index {mode_idx}."
        )
        for qrc_subjob in qrc_job.both_qrc_jobs:
            qrc_subjob.run()
            if not self._contains_imaginary_frequencies(qrc_subjob):
                logger.info(
                    f"QRC fallback produced frequency-validated structure: "
                    f"{qrc_subjob.label}"
                )
                return qrc_subjob
        logger.error(
            f"Imaginary frequencies remained after QRC fallback for {opt_job.label}. "
            "Fragment SP will not be run."
        )
        raise ValueError(
            f"Imaginary frequencies remained after QRC fallback for {opt_job.label}. "
            "Fragment SP was not started."
        )

    def _run_reactant_opt_and_sp_job(
        self, opt_job, sp_settings, sp_label, fallback_molecule
    ):
        logger.info(f"Running reactant fragment optimization: {opt_job.label}")
        opt_job.run()
        if self._contains_imaginary_frequencies(opt_job):
            logger.warning(
                f"Imaginary frequencies detected in {opt_job.label}. "
                "Retrying optimization with maxstep=5 and scf=qc."
            )
            opt_job = self._retry_reactant_opt_job(opt_job)
            logger.info(
                f"Running reactant fragment optimization retry: {opt_job.label}"
            )
            opt_job.run()
        else:
            logger.info(
                f"Reactant fragment optimization passed frequency check: "
                f"{opt_job.label}"
            )

        if self._contains_imaginary_frequencies(opt_job):
            logger.warning(
                f"Imaginary frequencies remained after retry for {opt_job.label}. "
                "Running QRC fallback."
            )
            opt_job = self._run_qrc_for_reactant_opt_job(opt_job)

        molecule = opt_job.optimized_structure()
        if molecule is None:
            logger.warning(
                f"Optimized structure unavailable for {opt_job.label}; using "
                f"fallback geometry for {sp_label}."
            )
            molecule = fallback_molecule
        logger.info(f"Running reactant fragment SP job: {sp_label}")
        sp_job = GaussianGeneralJob(
            molecule=molecule,
            settings=sp_settings,
            label=sp_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        sp_job.run()
        logger.info(f"Completed reactant fragment SP job: {sp_label}")

    def _run_fragment_reactant_jobs(self):
        logger.info(
            "Running optimized reactant fragment jobs with frequency checks."
        )
        self._run_reactant_opt_and_sp_job(
            opt_job=self.fragment1_reactant_opt_job,
            sp_settings=self.fragment1_reactant_sp_settings,
            sp_label=f"{self.label}_fragment1_r1",
            fallback_molecule=self.fragment1_atoms[-1],
        )
        self._run_reactant_opt_and_sp_job(
            opt_job=self.fragment2_reactant_opt_job,
            sp_settings=self.fragment2_reactant_sp_settings,
            sp_label=f"{self.label}_fragment2_i2",
            fallback_molecule=self.fragment2_atoms[-1],
        )

    def _run(self, **kwargs):
        """
        Execute the complete DI-AS calculation workflow.

        Runs all three sets of calculations in sequence:
        1. Complete molecules at sampled points
        2. Fragment 1 structures at the same points
        3. Fragment 2 structures at the same points
        4. Optimized isolated fragments and SP energies

        Args:
            **kwargs: Additional keyword arguments for job execution.
        """
        self._run_all_molecules_jobs()
        self._run_fragment1_jobs()
        self._run_fragment2_jobs()
        self._run_fragment_reactant_jobs()

    def is_complete(self):
        """
        Check if all DI-AS calculation jobs are complete.

        Verifies that all three sets of calculations (complete molecules,
        fragment 1, fragment 2, and reactant fragment jobs) have finished
        successfully.

        Returns:
            bool: True if all jobs are complete, False otherwise.
        """
        return (
            self._run_all_molecules_jobs_are_complete()
            and self._run_fragment1_jobs_are_complete()
            and self._run_fragment2_jobs_are_complete()
            and self._run_fragment_reactant_jobs_are_complete()
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

    def _run_fragment_reactant_jobs_are_complete(self):
        return (
            self.fragment1_reactant_opt_job.is_complete()
            and self.fragment1_reactant_sp_job.is_complete()
            and self.fragment2_reactant_opt_job.is_complete()
            and self.fragment2_reactant_sp_job.is_complete()
        )
