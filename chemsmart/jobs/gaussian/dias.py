from chemsmart.jobs.gaussian.job import GaussianGeneralJob, GaussianJob
from chemsmart.utils.utils import get_list_from_string_range


class GaussianDIASJob(GaussianJob):
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
        return len(self.all_molecules)

    def _fragment_structure(self, molecule):
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
        return [
            self._fragment_structure(molecule=molecule)[0]
            for molecule in self.all_molecules
        ]

    @property
    def fragment2_atoms(self):
        return [
            self._fragment_structure(molecule=molecule)[1]
            for molecule in self.all_molecules
        ]

    def _sample_molecules(self, molecules):
        """Samples provided molecules every self.every_n_points."""
        filtered_molecules = molecules[0 :: self.every_n_points]
        if (self.num_molecules - 1) / self.every_n_points != 0:
            filtered_molecules.append(molecules[-1])
        return filtered_molecules

    @property
    def fragment1_jobs(self):
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
        for job in self.all_molecules_jobs:
            job.run()

    def _run_fragment1_jobs(self):
        for job in self.fragment1_jobs:
            job.run()

    def _run_fragment2_jobs(self):
        for job in self.fragment2_jobs:
            job.run()

    def _run(self, jobrunner, **kwargs):
        self._run_all_molecules_jobs()
        self._run_fragment1_jobs()
        self._run_fragment2_jobs()

    def is_complete(self):
        return (
            self._run_all_molecules_jobs_are_complete()
            and self._run_fragment1_jobs_are_complete()
            and self._run_fragment2_jobs_are_complete()
        )

    def _run_all_molecules_jobs_are_complete(self):
        return all(job.is_complete() for job in self.all_molecules_jobs)

    def _run_fragment1_jobs_are_complete(self):
        return all(job.is_complete() for job in self.fragment1_jobs)

    def _run_fragment2_jobs_are_complete(self):
        return all(job.is_complete() for job in self.fragment2_jobs)
