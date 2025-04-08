#!/usr/bin/env python
import glob
import logging
import os

import click

from chemsmart.analysis.thermochemistry import qRRHOThermochemistry

logger = logging.getLogger(__name__)

os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-f",
    default=100.0,
    type=float,
    show_default=True,
    help="Cutoff frequency for both entropy and enthalpy in wavenumbers",
)
@click.option(
    "--fs",
    default=100.0,
    type=float,
    show_default=True,
    help="Cutoff frequency for entropy in wavenumbers",
)
@click.option(
    "--fh",
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
    "-i",
    "--isotope",
    is_flag=True,
    default=False,
    show_default=True,
    help="Use natural abundance weighted masses (False) or use single isotope masses (True).",
)
@click.option(
    "-q",
    is_flag=True,
    default=False,
    show_default=True,
    help="Quasi-RRHO approximation for both entropy and enthalpy.",
)
@click.option(
    "--qs",
    is_flag=True,
    default=False,
    show_default=True,
    help="Quasi-RRHO approximation for entropy.",
)
@click.option(
    "--qh",
    is_flag=True,
    default=False,
    show_default=True,
    help="Quasi-RRHO approximation for enthalpy.",
)
@click.option(
    "-u",
    "--unit",
    default="Eh",
    show_default=True,
    type=click.Choice(
        ["Eh", "hartree", "eV", "kcal/mol", "kJ/mol"], case_sensitive=False
    ),
    help="Unit of energetic values.",
)
@click.argument("filenames", nargs=-1, type=click.Path(exists=True))
def get_thermo(
    f,
    fs,
    fh,
    concentration,
    temperature,
    alpha,
    isotope,
    q,
    qs,
    qh,
    unit,
    filenames,
):
    """Thermochemistry calculation script using quasi-RRHO approximation."""

    def log(message, output="thermochemistry.dat"):
        logger.info(message)
        with open(output, "a") as f:
            f.write(message)

    if f != 100.0:
        fs = f
        fh = f

    if unit.lower() == "ev":
        energy_unit = "eV"
        unit_conversion = 27.211386024367243  # 1 Eh = 27.211386024367243 eV
    elif unit.lower() == "kcal/mol":
        energy_unit = "kcal/mol"
        unit_conversion = (
            627.5094738898777  # 1 Eh = 627.5094738898777 kcal/mol
        )
    elif unit.lower() == "kj/mol":
        energy_unit = "kJ/mol"
        unit_conversion = (
            2625.4996387552483  # 1 Eh = 2625.4996387552483 kJ/mol
        )
    else:
        energy_unit = "Hartree"
        unit_conversion = 1.0

    files = []
    for filename in filenames:
        files.extend(glob.glob(filename))
    for file in files:
        if not file.endswith((".log", ".out")):
            logger.error(
                f"Unsupported file extension for '{file}'. Only .log or .out files are accepted."
            )
            logger.error("Try 'get_thermochemistry.py --help' for help.")
            return
    if not files:
        logger.error(
            "Please provide calculation output files on the command line."
        )
        logger.error("Try 'get_thermochemistry.py --help' for help.")
        return

    log("\n")
    log(
        "   "
        + " " * 25
        + "  ____ _   _ _____ __  __ ____  __  __    _    ____ _____ \n"
    )
    log(
        "   "
        + " " * 25
        + " / ___| | | | ____|  \/  / ___||  \/  |  / \  |  _ \_   _|\n"
    )
    log(
        "   "
        + " " * 25
        + "| |   | |_| |  _| | |\/| \___ \| |\/| | / _ \ | |_) || |  \n"
    )
    log(
        "   "
        + " " * 25
        + "| |___|  _  | |___| |  | |___) | |  | |/ ___ \|  _ < | |  \n"
    )
    log(
        "   "
        + " " * 25
        + " \____|_| |_|_____|_|  |_|____/|_|  |_/_/   \_\_| \_\|_|  \n\n"
    )
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
    if q or qs:
        log("   " + f"Entropy Frequency Cut-off  : {fs:.1f} cm-1" + "\n")
    if q or qh:
        log("   " + f"Enthalpy Frequency Cut-off : {fh:.1f} cm-1" + "\n")
    if q or qs or qh:
        log("   " + f"Damping Function Exponent  : {alpha}" + "\n")
    log(
        "   "
        + f"Mass Weighted              : {'Single Isotope Masses' if isotope else 'Natural Abundance Weighted Masses'}"
        + "\n"
    )
    log("   " + f"Energy Unit                : {energy_unit}" + "\n\n")
    if q or qs or qh:
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
    if q or qs:
        log("   - Entropic quasi-harmonic treatment: Grimme\n")
        log("     REF: Grimme, S. Chem. Eur. J. 2012, 18, 9955-9964\n\n")
    if q or qh:
        log("   - Enthalpy quasi-harmonic treatment: Head-Gordon\n")
        log(
            "     REF: Li, Y.; Gomes, J.; Sharada, S. M.; Bell, A. T.; Head-Gordon, M. J. Phys. Chem. C 2015, 119, 1840-1850\n\n"
        )
    log("\n")

    log(f" * Thermochemistry Results (in {energy_unit}): \n\n")
    if q:
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
    elif qs:
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
    elif qh:
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
    for file in files:
        try:
            thermochemistry = qRRHOThermochemistry(
                file,
                temperature=temperature,
                concentration=concentration,
                weighted_atomic_mass=not isotope,
                alpha=alpha,
                s_freq_cutoff=fs,
                h_freq_cutoff=fh,
            )
            structure = os.path.splitext(os.path.basename(file))[0]
            energy = thermochemistry.energies * unit_conversion
            zero_point_energy = (
                thermochemistry.zero_point_energy_hartree * unit_conversion
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
            if q:
                qrrho_gibbs_free_energy = (
                    thermochemistry.qrrho_gibbs_free_energy * unit_conversion
                )
            elif qh:
                qrrho_gibbs_free_energy = (
                    thermochemistry.qrrho_gibbs_free_energy_qh
                    * unit_conversion
                )
            elif qs:
                qrrho_gibbs_free_energy = (
                    thermochemistry.qrrho_gibbs_free_energy_qs
                    * unit_conversion
                )

            if q:
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
            elif qs:
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
            elif qh:
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
            pass
    log("\n")
    logger.info(" * Done. Results saved to 'thermochemistry.dat'.")


if __name__ == "__main__":
    get_thermo()
