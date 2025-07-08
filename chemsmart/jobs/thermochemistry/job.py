import logging
import os
from typing import Type

from chemsmart.analysis.thermochemistry import Thermochemistry
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.thermochemistry.settings import ThermochemistryJobSettings

logger = logging.getLogger(__name__)


class ThermochemistryJob(Job):
    PROGRAM = "Thermochemistry"
    TYPE = "thermochemistry"

    def __init__(
        self,
        molecule=None,
        settings=None,
        label=None,
        jobrunner=None,
        filename=None,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        if settings is not None and not isinstance(
            settings, ThermochemistryJobSettings
        ):
            raise ValueError(
                f"Settings must be instance of {ThermochemistryJobSettings} for {self}, but is {settings} instead!"
            )

        if molecule is not None and not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is {molecule} instead!"
            )

        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = (
            settings.copy()
            if settings is not None
            else ThermochemistryJobSettings()
        )
        self.filename = filename

        if label is None:
            if filename is not None:
                label = os.path.splitext(os.path.basename(filename))[0]
            elif molecule is not None:
                label = molecule.get_chemical_formula(empirical=True)
            else:
                label = "thermochemistry_job"
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[ThermochemistryJobSettings]:
        return ThermochemistryJobSettings

    @property
    def inputfile(self):
        if self.filename:
            return os.path.abspath(self.filename)
        return None

    @property
    def outputfile(self):
        outputfile = self.label + ".dat"
        return os.path.join(self.folder, outputfile)

    @property
    def errfile(self):
        errfile = self.label + ".err"
        return os.path.join(self.folder, errfile)

    def _backup_files(self, **kwargs):
        folder = self._create_backup_folder_name()
        self.backup_file(self.outputfile, folder=folder, **kwargs)

    def _output(self):
        if not os.path.exists(self.outputfile):
            return None
        return os.path.abspath(self.outputfile)

    def _job_is_complete(self):
        return os.path.exists(self.outputfile)

    def _run(self, **kwargs):
        """Run the thermochemistry analysis job."""
        logger.info(
            f"Running ThermochemistryJob {self} with jobrunner {self.jobrunner}"
        )
        self.jobrunner.run(self, **kwargs)

    @classmethod
    def from_filename(
        cls,
        filename,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """Create a ThermochemistryJob from a Gaussian or ORCA output file."""
        if not filename.endswith((".log", ".out")):
            raise ValueError(
                f"Unsupported file extension for '{filename}'. Only .log or .out files are accepted."
            )

        logger.info(f"Reading molecule from file: {filename}.")
        molecule = Molecule.from_filepath(filename)

        if settings is None:
            settings = ThermochemistryJobSettings()

        if label is None:
            label = os.path.splitext(os.path.basename(filename))[0]

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecule,
                    settings=settings,
                    label=label,
                    filename=filename,
                    jobrunner=None,
                    **kwargs,
                ),
                server=kwargs.get("server"),
                scratch=kwargs.get("scratch"),
                fake=kwargs.get("fake", False),
                **kwargs,
            )

        return cls(
            molecule=molecule,
            settings=settings,
            label=label,
            filename=filename,
            jobrunner=jobrunner,
            **kwargs,
        )

    def compute_thermochemistry(self):
        """Perform the thermochemistry calculation and save results."""
        if not self.filename:
            raise ValueError(
                "No input file provided for thermochemistry calculation."
            )

        try:
            thermochemistry = Thermochemistry(
                filename=self.filename,
                temperature=self.settings.temperature,
                concentration=self.settings.concentration,
                pressure=self.settings.pressure,
                use_weighted_mass=self.settings.use_weighted_mass,
                alpha=self.settings.alpha,
                s_freq_cutoff=self.settings.s_freq_cutoff,
                h_freq_cutoff=self.settings.h_freq_cutoff,
            )
            thermochemistry.compute_thermochemistry()

            # # Energy conversion based on settings.units
            # if self.settings.energy_units == "ev":
            #     unit_conversion = joule_per_mol_to_eV  # J/mol to eV
            #     energy_unit = "eV"
            # elif self.settings.energy_units == "kcal/mol":
            #     unit_conversion = (
            #         joule_per_mol_to_kcal_per_mol  # J/mol to kcal/mol
            #     )
            #     energy_unit = "kcal/mol"
            # elif self.settings.energy_units == "kj/mol":
            #     unit_conversion = (
            #         joule_per_mol_to_kJ_per_mol  # J/mol to kJ/mol
            #     )
            #     energy_unit = "kJ/mol"
            # else:
            #     unit_conversion = joule_per_mol_to_hartree  # J/mol to Hartree
            #     energy_unit = "Hartree"
            #
            # # Calculate thermochemical properties
            # structure = os.path.splitext(os.path.basename(self.filename))[0]
            # energy = thermochemistry.electronic_energy * unit_conversion
            # zero_point_energy = (
            #     thermochemistry.zero_point_energy * unit_conversion
            # )
            # enthalpy = thermochemistry.enthalpy * unit_conversion
            # entropy_times_temperature = (
            #     thermochemistry.entropy_times_temperature * unit_conversion
            # )
            # gibbs_free_energy = (
            #     thermochemistry.gibbs_free_energy * unit_conversion
            # )
            #
            # qrrho_enthalpy = (
            #     thermochemistry.qrrho_enthalpy * unit_conversion
            #     if self.settings.h_freq_cutoff
            #     else None
            # )
            # qrrho_entropy_times_temperature = (
            #     thermochemistry.qrrho_entropy_times_temperature
            #     * unit_conversion
            #     if self.settings.s_freq_cutoff
            #     else None
            # )
            # if self.settings.s_freq_cutoff and self.settings.h_freq_cutoff:
            #     qrrho_gibbs_free_energy = (
            #         thermochemistry.qrrho_gibbs_free_energy * unit_conversion
            #     )
            # elif self.settings.h_freq_cutoff:
            #     qrrho_gibbs_free_energy = (
            #         thermochemistry.qrrho_gibbs_free_energy_qh
            #         * unit_conversion
            #     )
            # elif self.settings.s_freq_cutoff:
            #     qrrho_gibbs_free_energy = (
            #         thermochemistry.qrrho_gibbs_free_energy_qs
            #         * unit_conversion
            #     )
            # else:
            #     qrrho_gibbs_free_energy = None
            #
            # # Check for imaginary frequencies
            # if thermochemistry.imaginary_frequencies:
            #     if (
            #         thermochemistry.job_type == "ts"
            #         and len(thermochemistry.imaginary_frequencies) == 1
            #     ):
            #         logger.info(
            #             f"Correct Transition State detected: only 1 imaginary frequency\n"
            #             f"Imaginary frequency excluded for thermochemistry calculation in {self.filename}."
            #         )
            #     else:
            #         raise ValueError(
            #             f"Invalid number of imaginary frequencies for {self.filename}. "
            #             f"Expected 0 for optimization or 1 for TS, but found "
            #             f"{len(thermochemistry.imaginary_frequencies)} for job: {thermochemistry.job_type}!"
            #         )
            #
            # # Write results to output file
            # with open(self.outputfile, "w") as out:
            #     out.write(f"Thermochemistry Results for {structure}:\n\n")
            #     out.write(f"Temperature: {self.settings.temperature:.2f} K\n")
            #     if self.settings.concentration is not None:
            #         out.write(
            #             f"Concentration: {self.settings.concentration:.1f} mol/L\n"
            #         )
            #     else:
            #         out.write(f"Pressure: {self.settings.pressure:.1f} atm\n")
            #     if self.settings.s_freq_cutoff:
            #         out.write(
            #             f"Entropy Frequency Cut-off: {self.settings.s_freq_cutoff:.1f} cm^-1\n"
            #         )
            #     if self.settings.h_freq_cutoff:
            #         out.write(
            #             f"Enthalpy Frequency Cut-off: {self.settings.h_freq_cutoff:.1f} cm^-1\n"
            #         )
            #     if self.settings.s_freq_cutoff or self.settings.h_freq_cutoff:
            #         out.write(
            #             f"Damping Function Exponent: {self.settings.alpha}\n"
            #         )
            #     out.write(
            #         f"Mass Weighted: {'Most Abundant Masses' if not self.settings.use_weighted_mass else 'Natural Abundance Weighted Masses'}\n"
            #     )
            #     out.write(f"Energy Unit: {energy_unit}\n\n")
            #
            #     if self.settings.h_freq_cutoff or self.settings.s_freq_cutoff:
            #         out.write(qrrho_header)
            #         out.write(head_gordon_damping_function_ref)
            #         if self.settings.s_freq_cutoff:
            #             out.write(grimme_quasi_rrho_entropy_ref)
            #         if self.settings.h_freq_cutoff:
            #             out.write(head_gordon_quasi_rrho_enthalpy_ref)
            #         out.write("\n")
            #
            #     # Write table header
            #     if self.settings.h_freq_cutoff and self.settings.s_freq_cutoff:
            #         out.write(
            #             "{:<39} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} {:>13}\n".format(
            #                 "Structure",
            #                 "E",
            #                 "ZPE",
            #                 "H",
            #                 "qh-H",
            #                 "T.S",
            #                 "T.qh-S",
            #                 "G(T)",
            #                 "qh-G(T)",
            #             )
            #         )
            #         out.write("=" * 142 + "\n")
            #         out.write(
            #             "{:39} {:13.6f} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
            #                 structure,
            #                 energy,
            #                 zero_point_energy,
            #                 enthalpy,
            #                 qrrho_enthalpy,
            #                 entropy_times_temperature,
            #                 qrrho_entropy_times_temperature,
            #                 gibbs_free_energy,
            #                 qrrho_gibbs_free_energy,
            #             )
            #         )
            #     elif self.settings.s_freq_cutoff:
            #         out.write(
            #             "{:<39} {:>13} {:>10} {:>13} {:>10} {:>10} {:>13} {:>13}\n".format(
            #                 "Structure",
            #                 "E",
            #                 "ZPE",
            #                 "H",
            #                 "T.S",
            #                 "T.qh-S",
            #                 "G(T)",
            #                 "qh-G(T)",
            #             )
            #         )
            #         out.write("=" * 128 + "\n")
            #         out.write(
            #             "{:39} {:13.6f} {:10.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
            #                 structure,
            #                 energy,
            #                 zero_point_energy,
            #                 enthalpy,
            #                 entropy_times_temperature,
            #                 qrrho_entropy_times_temperature,
            #                 gibbs_free_energy,
            #                 qrrho_gibbs_free_energy,
            #             )
            #         )
            #     elif self.settings.h_freq_cutoff:
            #         out.write(
            #             "{:<39} {:>13} {:>10} {:>13} {:>13} {:>10} {:>13} {:>13}\n".format(
            #                 "Structure",
            #                 "E",
            #                 "ZPE",
            #                 "H",
            #                 "qh-H",
            #                 "T.S",
            #                 "G(T)",
            #                 "qh-G(T)",
            #             )
            #         )
            #         out.write("=" * 131 + "\n")
            #         out.write(
            #             "{:39} {:13.6f} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
            #                 structure,
            #                 energy,
            #                 zero_point_energy,
            #                 enthalpy,
            #                 qrrho_enthalpy,
            #                 entropy_times_temperature,
            #                 gibbs_free_energy,
            #                 qrrho_gibbs_free_energy,
            #             )
            #         )
            #     else:
            #         out.write(
            #             "{:<39} {:>13} {:>10} {:>13} {:>10} {:>13}\n".format(
            #                 "Structure", "E", "ZPE", "H", "T.S", "G(T)"
            #             )
            #         )
            #         out.write("=" * 103 + "\n")
            #         out.write(
            #             "{:39} {:13.6f} {:10.6f} {:13.6f} {:10.6f} {:13.6f}\n".format(
            #                 structure,
            #                 energy,
            #                 zero_point_energy,
            #                 enthalpy,
            #                 entropy_times_temperature,
            #                 gibbs_free_energy,
            #             )
            #         )
            #
            # logger.info(f"Thermochemistry results saved to {self.outputfile}")

        except Exception as e:
            logger.error(f"Error processing {self.filename}: {e}")
            raise
