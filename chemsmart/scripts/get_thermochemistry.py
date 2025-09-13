#!/usr/bin/env python
"""
Thermochemistry data extraction script.

This script extracts thermochemical properties (enthalpies, entropies,
Gibbs free energies) from computational chemistry output files and
performs thermochemical analysis and calculations.
"""

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
    "-d",
    "--directory",
    default=None,
    help="Directory in which to analyze files.",
)
@click.option(
    "-ft",
    "--filetype",
    default=None,
    type=click.Choice(["log", "out"], case_sensitive=False),
    help="Type of file to be analyzed, if directory is specified.",
)
@click.option(
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
    default=None,
    type=float,
    help="Concentration in mol/L (solution).",
)
@click.option(
    "-p",
    "--pressure",
    default=1.0,
    type=float,
    show_default=True,
    help="Pressure in Standard atmosphere (gas phase). "
         "Ignored if -c/--concentration is provided.",
)
@click.option(
    "-t",
    "--temperature",
    default=298.15,
    type=float,
    show_default=True,
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
    help="Use natural abundance weighted masses (True) or use most "
         "abundant masses (False).\nDefault to False, i.e., use single "
         "isotopic mass.",
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
def entry_point(
    filenames,
    directory,
    filetype,
    cutoff,
    entropy_cutoff,
    enthalpy_cutoff,
    concentration,
    pressure,
    temperature,
    alpha,
    weighted,
    quasi_rrho,
    quasi_rrho_entropy,
    quasi_rrho_enthalpy,
    units,
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
            pattern = "*.log"
            filenames = glob.glob(os.path.join(directory, pattern))
        elif filetype == "out":
            pattern = "*.out"
            filenames = glob.glob(os.path.join(directory, pattern))
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
                f"Unsupported file extension for '{file}'. "
                f"Only .log or .out files are accepted."
            )
            logger.error("Try 'get_thermochemistry.py --help' for help.")
            return

    # Display ChemSmart ASCII banner
    logger.info("\n")
    logger.info(
        "   "
        + " " * 25
        + "  ____ _   _ _____ __  __ ____  __  __    _    ____ _____ "
    )
    logger.info(
        "   "
        + " " * 25
        + " / ___| | | | ____|  \/  / ___||  \/  |  / \  |  _ \_   _|"
    )
    logger.info(
        "   "
        + " " * 25
        + "| |   | |_| |  _| | |\/| \___ \| |\/| | / _ \ | |_) || | "
    )
    logger.info(
        "   "
        + " " * 25
        + "| |___|  _  | |___| |  | |___) | |  | |/ ___ \|  _ < | | "
    )
    logger.info(
        "   "
        + " " * 25
        + " \____|_| |_|_____|_|  |_|____/|_|  |_/_/   \_\_| \_\|_| \n"
    )

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
    if concentration is not None:
        log(
            "   "
            + f"Concentration              : {concentration:.1f} mol/L"
            + "\n"
        )
    else:
        log("   " + f"Pressure                   : {pressure:.1f} atm" + "\n")
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
                pressure=pressure,
                use_weighted_mass=weighted,
                alpha=alpha,
                s_freq_cutoff=entropy_cutoff,
                h_freq_cutoff=enthalpy_cutoff,
            )
            structure = os.path.splitext(os.path.basename(file))[0]
            energy = thermochemistry.electronic_energy * unit_conversion
            zero_point_energy = (
                thermochemistry.zero_point_energy * unit_conversion
            )
            enthalpy = thermochemistry.enthalpy * unit_conversion
            qrrho_enthalpy = thermochemistry.qrrho_enthalpy * unit_conversion
            entropy_times_temperature = (
                thermochemistry.entropy_times_temperature * unit_conversion
            )
            qrrho_entropy_times_temperature = (
                thermochemistry.qrrho_entropy_times_temperature
                * unit_conversion
            )
            gibbs_free_energy = (
                thermochemistry.gibbs_free_energy * unit_conversion
            )
            if quasi_rrho:
                qrrho_gibbs_free_energy = (
                    thermochemistry.qrrho_gibbs_free_energy * unit_conversion
                )
            elif quasi_rrho_enthalpy:
                qrrho_gibbs_free_energy = (
                    thermochemistry.qrrho_gibbs_free_energy_qh
                    * unit_conversion
                )
            elif quasi_rrho_entropy:
                qrrho_gibbs_free_energy = (
                    thermochemistry.qrrho_gibbs_free_energy_qs
                    * unit_conversion
                )
            # Warning for imaginary frequency.
            if thermochemistry.imaginary_frequencies:
                if (
                    thermochemistry.job_type == "ts"
                    and thermochemistry.vibrational_frequencies[0] < 0.0
                ):
                    if len(thermochemistry.imaginary_frequencies) > 1:
                        continue
                    else:
                        logger.warning(
                            f"!! Transition state detected: ignored 1 imaginary frequency for {file}.\n"
                        )
                else:
                    continue

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

        except TypeError:
            log(
                "{:2} {:39} {:13.6f} {:<50}\n".format(
                    " ×",
                    structure,
                    energy,
                    "  Warning! Frequency information not found ...",
                )
            )
            continue

        except Exception as e:
            logger.error(f"{e}\n")
            continue
    log("\n")
    logger.info(" * Done. Results saved to 'thermochemistry.dat'.")


if __name__ == "__main__":
    entry_point()
