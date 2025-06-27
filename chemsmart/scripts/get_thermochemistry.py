#!/usr/bin/env python
import glob
import logging
import os

import click

from chemsmart.analysis.thermochemistry import Thermochemistry

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)-7s - [%(name)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-f",
    "--filenames",
    type=str,
    multiple=True,
    default=None,
    help="Gaussian or ORCA output files for parsing thermochemistry.",
)
@click.option(
    "-ft",
    "--filetype",
    default=None,
    help="Type of file to be converted, if directory is specified.",
)
@click.option(
    "-c",
    "--cutoff",
    default=100.0,
    type=float,
    show_default=True,
    help="Cutoff frequency for both entropy and enthalpy in wavenumbers",
)
@click.option(
    "-cs",
    "--entropy-cutoff",
    default=100.0,
    type=float,
    show_default=True,
    help="Cutoff frequency for entropy in wavenumbers",
)
@click.option(
    "-ch",
    "--enthalpy-cutoff",
    default=100.0,
    type=float,
    show_default=True,
    help="Cutoff frequency for enthalpy in wavenumbers",
)
@click.option(
    "-c",
    "--concentration",
    default=1.0,
    type=float,
    show_default=True,
    help="Concentration in mol/L",
)
@click.option(
    "-t",
    "--temperature",
    required=True,
    default=None,
    type=float,
    help="Temperature in Kelvin.",
)
@click.option(
    "-a",
    "--alpha",
    default=4,
    type=int,
    show_default=True,
    help="Interpolator exponent used in the quasi-RRHO approximation.",
)
@click.option(
    "-w",
    "--weighted",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use natural abundance weighted masses (True) or use most abundant masses (False).\n"
    "Default to False, i.e., use single isotopic mass.",
)
@click.option(
    "-q",
    "--quasi-rrho",
    is_flag=True,
    default=False,
    show_default=True,
    help="Quasi-RRHO approximation for both entropy and enthalpy.",
)
@click.option(
    "-qs",
    "--quasi-rrho-entropy",
    is_flag=True,
    default=False,
    show_default=True,
    help="Apply quasi-RRHO approximation for entropy.",
)
@click.option(
    "-qh",
    "--quasi-rrho-enthalpy",
    is_flag=True,
    default=False,
    show_default=True,
    help="Apply quasi-RRHO approximation for enthalpy.",
)
@click.option(
    "-u",
    "--units",
    default="hartree",
    show_default=True,
    type=click.Choice(
        ["hartree", "eV", "kcal/mol", "kJ/mol"], case_sensitive=False
    ),
    help="Units of energetic values.",
)
@click.option(
    "--scf",
    is_flag=True,
    default=False,
    show_default=True,
    help="Only output SCF energy without thermochemical corrections.",
)
def entry_point(
    filenames,
    filetype,
    cutoff,
    entropy_cutoff,
    enthalpy_cutoff,
    concentration,
    temperature,
    alpha,
    weighted,
    quasi_rrho,
    quasi_rrho_entropy,
    quasi_rrho_enthalpy,
    units,
    scf,
):
    """Thermochemistry calculation script using quasi-RRHO approximation."""

    def log(message, output="thermochemistry.dat"):
        # log is a math function for logarithm
        logger.info(message)
        with open(output, "a") as out:
            out.write(message)

    if cutoff != 100.0:
        entropy_cutoff = cutoff
        enthalpy_cutoff = cutoff

    # Energy Conversion
    if units.lower() == "ev":
        energy_unit = "eV"
        unit_conversion = (
            1.0364269574711572e-05  # 1 J/mol = 1.0364269574711572e-05 eV
        )
    elif units.lower() == "kcal/mol":
        energy_unit = "kcal/mol"
        unit_conversion = (
            0.0002390057361376673  # 1 J/mol = 0.0002390057361376673 kcal/mol
        )
    elif units.lower() == "kj/mol":
        energy_unit = "kJ/mol"
        unit_conversion = 0.001  # 1 J/mol = 0.001 kJ/mol
    else:
        energy_unit = "Hartree"
        unit_conversion = (
            3.8087991196914175e-07  # 1 J/mol = 3.8087991196914175e-07 Eh
        )

    if not filenames:
        logger.info("No filenames provided.")
        logger.info(
            "Will use all files in the folder with suffixes .log or .out."
        )
        assert filetype, "Type of file to be converted must be specified."
        if filetype == "log":
            filenames = glob.glob("*.log")
        elif filetype == "out":
            filenames = glob.glob("*.out")
        else:
            logger.error(
                "Please provide calculation output files on the command line."
            )
            logger.error("Try 'get_thermochemistry.py --help' for help.")
            return

    # Error Handling
    for file in filenames:
        if not file.endswith((".log", ".out")):
            logger.error(
                f"Unsupported file extension for '{file}'. Only .log or .out files are accepted."
            )
            logger.error("Try 'get_thermochemistry.py --help' for help.")
            return

    # ASCII Arts for CHEMSMART
    log("\n")
    log(
        "   "
        + " " * 25
        + "  ____ _   _ _____ __  __ ____  __  __    _    ____ _____ "
    )
    log(
        "   "
        + " " * 25
        + " / ___| | | | ____|  \/  / ___||  \/  |  / \  |  _ \_   _|"
    )
    log(
        "   "
        + " " * 25
        + "| |   | |_| |  _| | |\/| \___ \| |\/| | / _ \ | |_) || | "
    )
    log(
        "   "
        + " " * 25
        + "| |___|  _  | |___| |  | |___) | |  | |/ ___ \|  _ < | | "
    )
    log(
        "   "
        + " " * 25
        + " \____|_| |_|_____|_|  |_|____/|_|  |_/_/   \_\_| \_\|_| \n"
    )

    # Only output SCF energy without thermochemical corrections.
    if scf:
        log(f" * SCF Energies (in {energy_unit}): \n\n")
        log(
            "   "
            + "{:<60} {:>22} {:>22}\n".format("Structure", "E", "Job Type")
        )
        log("   " + "=" * 106 + "\n")
        index = 1
        for file in filenames:
            thermochemistry = Thermochemistry(file, temperature=temperature)
            structure = os.path.splitext(os.path.basename(file))[0]
            energy = thermochemistry.energies * unit_conversion
            job_type = thermochemistry.job_type
            log(
                "{:2} {:60} {:22.6f} {:>22}\n".format(
                    index,
                    structure,
                    energy,
                    job_type,
                )
            )
            index += 1
        return

    # Output Thermochemistry analysis.
    log("   " + "┌" + "─" * 106 + "┐" + "\n")
    log(
        "   "
        + "├"
        + " " * 41
        + "Thermochemistry Settings"
        + " " * 41
        + "┤"
        + "\n"
    )
    log("   " + "└" + "─" * 106 + "┘" + "\n")
    log("   " + f"Temperature                : {temperature:.2f} K" + "\n")
    log(
        "   "
        + f"Concentration              : {concentration:.1f} mol/L"
        + "\n"
    )
    if quasi_rrho or quasi_rrho_entropy:
        log(
            "   "
            + f"Entropy Frequency Cut-off  : {entropy_cutoff:.1f} cm-1"
            + "\n"
        )
    if quasi_rrho or quasi_rrho_enthalpy:
        log(
            "   "
            + f"Enthalpy Frequency Cut-off : {enthalpy_cutoff:.1f} cm-1"
            + "\n"
        )
    if quasi_rrho or quasi_rrho_entropy or quasi_rrho_enthalpy:
        log("   " + f"Damping Function Exponent  : {alpha}" + "\n")
    log(
        "   "
        + f"Mass Weighted              : {'Most Abundant Masses' if not weighted else 'Natural Abundance Weighted Masses'}"
        + "\n"
    )
    log("   " + f"Energy Unit                : {energy_unit}" + "\n\n")
    if quasi_rrho or quasi_rrho_entropy or quasi_rrho_enthalpy:
        log("   " + "-" * 108 + "\n")
        log(
            "   "
            + " " * 32
            + "Quasi-Rigid-Rotor-Harmonic-Oscillator Scheme"
            + "\n"
        )
        log("   " + "-" * 108 + "\n")
        log("   - Damping function: Chai and Head-Gordon\n")
        log(
            "     REF: Chai, J.-D.; Head-Gordon, M. Phys. Chem. Chem. Phys. 2008, 10, 6615–6620\n\n"
        )
    if quasi_rrho or quasi_rrho_entropy:
        log("   - Entropic quasi-harmonic treatment: Grimme\n")
        log("     REF: Grimme, S. Chem. Eur. J. 2012, 18, 9955-9964\n\n")
    if quasi_rrho or quasi_rrho_enthalpy:
        log("   - Enthalpy quasi-harmonic treatment: Head-Gordon\n")
        log(
            "     REF: Li, Y.; Gomes, J.; Sharada, S. M.; Bell, A. T.; Head-Gordon, M. J. Phys. Chem. C 2015, 119, 1840-1850\n\n"
        )
    log("\n")

    log(f" * Thermochemistry Results (in {energy_unit}): \n\n")
    if quasi_rrho:
        log(
            "   {:<39} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} {:>13}\n".format(
                "Structure",
                "E",
                "ZPE",
                "H",
                "qh-H",
                "T.S",
                "T.qh-S",
                "G(T)",
                "qh-G(T)",
            )
        )
        log("   " + "=" * 142 + "\n")
    elif quasi_rrho_entropy:
        log(
            "   {:<39} {:>13} {:>10} {:>13} {:>10} {:>10} {:>13} {:>13}\n".format(
                "Structure",
                "E",
                "ZPE",
                "H",
                "T.S",
                "T.qh-S",
                "G(T)",
                "qh-G(T)",
            )
        )
        log("   " + "=" * 128 + "\n")
    elif quasi_rrho_enthalpy:
        log(
            "   {:<39} {:>13} {:>10} {:>13} {:>13} {:>10} {:>13} {:>13}\n".format(
                "Structure", "E", "ZPE", "H", "qh-H", "T.S", "G(T)", "qh-G(T)"
            )
        )
        log("   " + "=" * 131 + "\n")
    else:
        log(
            "   {:<39} {:>13} {:>10} {:>13} {:>10} {:>13}\n".format(
                "Structure", "E", "ZPE", "H", "T.S", "G(T)"
            )
        )
        log("   " + "=" * 103 + "\n")

    index = 1
    for file in filenames:
        try:
            thermochemistry = Thermochemistry(
                file,
                temperature=temperature,
                concentration=concentration,
                use_weighted_mass=weighted,
                alpha=alpha,
                s_freq_cutoff=entropy_cutoff,
                h_freq_cutoff=enthalpy_cutoff,
            )
            structure = os.path.splitext(os.path.basename(file))[0]
            energy = thermochemistry.energies * unit_conversion
            zero_point_energy = (
                thermochemistry.zero_point_energy * unit_conversion
            )
            enthalpy = thermochemistry.enthalpy * unit_conversion
            qrrho_enthalpy = thermochemistry.qrrho_enthalpy * unit_conversion
            entropy_times_temperature = (
                    thermochemistry.entropy_times_temperature_concentration * unit_conversion
            )
            qrrho_entropy_times_temperature = (
                thermochemistry.qrrho_entropy_times_temperature_concentration
                * unit_conversion
            )
            gibbs_free_energy = (
                    thermochemistry.gibbs_free_energy_concentration * unit_conversion
            )
            if quasi_rrho:
                qrrho_gibbs_free_energy = (
                        thermochemistry.qrrho_gibbs_free_energy_concentration * unit_conversion
                )
            elif quasi_rrho_enthalpy:
                qrrho_gibbs_free_energy = (
                    thermochemistry.qrrho_gibbs_free_energy_concentration_qh
                    * unit_conversion
                )
            elif quasi_rrho_entropy:
                qrrho_gibbs_free_energy = (
                    thermochemistry.qrrho_gibbs_free_energy_concentration_qs
                    * unit_conversion
                )

            # Warning for imaginary frequency.
            if thermochemistry.imaginary_frequencies:
                if (
                    thermochemistry.job_type == "ts"
                    and thermochemistry.vibrational_frequencies[0] < 0.0
                ):
                    if thermochemistry.num_replaced_frequencies > 0:
                        logger.error(
                            f"!! Detected multiple imaginary frequencies in transition state for {file} — aborting.\n"
                        )
                        raise ValueError(
                            f"Error: Detected multiple imaginary frequencies in TS calculation for {file}. "
                            f"Only one imaginary frequency is allowed for a valid TS. "
                            f"Please re-optimize the geometry to locate a true TS."
                        )
                    else:
                        logger.warning(
                            f"!! Transition state detected: ignored 1 imaginary frequency for {file}.\n"
                        )
                else:
                    logger.warning(
                        f"!! Replaced {thermochemistry.num_replaced_frequencies} imaginary frequency(s) with the cutoff value ({thermochemistry.cutoff:.1f} cm-1) for {file}.\n"
                    )

            if quasi_rrho:
                log(
                    "{:2} {:39} {:13.6f} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
                        index,
                        structure,
                        energy,
                        zero_point_energy,
                        enthalpy,
                        qrrho_enthalpy,
                        entropy_times_temperature,
                        qrrho_entropy_times_temperature,
                        gibbs_free_energy,
                        qrrho_gibbs_free_energy,
                    )
                )
            elif quasi_rrho_entropy:
                log(
                    "{:2} {:39} {:13.6f} {:10.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
                        index,
                        structure,
                        energy,
                        zero_point_energy,
                        enthalpy,
                        entropy_times_temperature,
                        qrrho_entropy_times_temperature,
                        gibbs_free_energy,
                        qrrho_gibbs_free_energy,
                    )
                )
            elif quasi_rrho_enthalpy:
                log(
                    "{:2} {:39} {:13.6f} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
                        index,
                        structure,
                        energy,
                        zero_point_energy,
                        enthalpy,
                        qrrho_enthalpy,
                        entropy_times_temperature,
                        gibbs_free_energy,
                        qrrho_gibbs_free_energy,
                    )
                )
            else:
                log(
                    "{:2} {:39} {:13.6f} {:10.6f} {:13.6f} {:10.6f} {:13.6f}\n".format(
                        index,
                        structure,
                        energy,
                        zero_point_energy,
                        enthalpy,
                        entropy_times_temperature,
                        gibbs_free_energy,
                    )
                )
            index += 1
        except (ValueError, TypeError, IndexError, AttributeError):
            logger.warning(
                f"!! Frequency information not found, skipped {file}.\n"
            )
    log("\n")
    logger.info(" * Done. Results saved to 'thermochemistry.dat'.")


if __name__ == "__main__":
    entry_point()
