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
    "-w",
    "--weight",
    is_flag=True,
    default=True,
    show_default=True,
    help="Use natural abundance weighted masses (True) or use single isotope masses (False).",
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
@click.argument("filenames", nargs=-1, type=click.Path(exists=True))
def get_thermo(
    f, fs, fh, concentration, temperature, alpha, weight, q, qs, qh, filenames
):
    """Thermochemistry calculation script using quasi-RRHO approximation."""
    if f != 100.0:
        fs = f
        fh = f

    files = []
    for filename in filenames:
        files.extend(glob.glob(filename))
    if len(files) == 0:
        logger.info(
            "Please provide calculation output files on the command line."
        )
        logger.info("Try 'get_thermochemistry.py --help' for help.")
        return

    logger.info(
        f"   Temperature = {temperature:.2f} Kelvin   Concentration = {concentration:.1f} mol/L"
    )
    logger.info(
        "   All energetic values below shown in Hartree unless otherwise specified."
    )
    logger.info("")
    if q or qs or qh:
        logger.info(
            f"   Damping function: dimensionless interpolator exponent of {alpha} will be used in the quasi-RRHO scheme."
        )
        logger.info(
            "   Chai and Head-Gordon: Long-range corrected hybrid density functionals with damped atom–atom dispersion corrections."
        )
        logger.info(
            "   REF: Chai, J.-D.; Head-Gordon, M. Phys. Chem. Chem. Phys. 2008, 10, 6615–6620"
        )
        logger.info("")
    if q or qs:
        logger.info(
            f"   Entropic quasi-harmonic treatment: frequency cut-off value of {fs:.1f} wavenumbers will be applied."
        )
        logger.info(
            "   QS = Grimme: Using a mixture of RRHO and Free-rotor vibrational entropies."
        )
        logger.info("   REF: Grimme, S. Chem. Eur. J. 2012, 18, 9955-9964")
        logger.info("")
    if q or qh:
        logger.info(
            f"   Enthalpy quasi-harmonic treatment: frequency cut-off value of {fh:.1f} wavenumbers will be applied."
        )
        logger.info(
            "   QH = Head-Gordon: Using an RRHO treatement with an approximation term for vibrational energy."
        )
        logger.info(
            "   REF: Li, Y.; Gomes, J.; Sharada, S. M.; Bell, A. T.; Head-Gordon, M. J. Phys. Chem. C 2015, 119, 1840-1850"
        )
        logger.info("")

    if q:
        logger.info(
            "   {:<39} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} {:>13}".format(
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
    elif qs:
        logger.info(
            "   {:<39} {:>13} {:>10} {:>13} {:>10} {:>10} {:>13} {:>13}".format(
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
    elif qh:
        logger.info(
            "   {:<39} {:>13} {:>10} {:>13} {:>13} {:>10} {:>13} {:>13}".format(
                "Structure", "E", "ZPE", "H", "qh-H", "T.S", "G(T)", "qh-G(T)"
            )
        )
    else:
        logger.info(
            "   {:<39} {:>13} {:>10} {:>13} {:>10} {:>13}".format(
                "Structure", "E", "ZPE", "H", "T.S", "G(T)"
            )
        )

    for file in files:
        thermochemistry = qRRHOThermochemistry(
            file,
            temperature=temperature,
            concentration=concentration,
            weighted_atomic_mass=weight,
            alpha=alpha,
            s_freq_cutoff=fs,
            h_freq_cutoff=fh,
        )
        structure = os.path.basename(file).split("/")[-1].split(".")[0]
        energy = thermochemistry.energies
        zero_point_energy = thermochemistry.zero_point_energy_hartree
        enthalpy = thermochemistry.enthalpy
        qrrho_enthalpy = thermochemistry.qrrho_enthalpy
        entropy_times_temperature = thermochemistry.entropy_times_temperature
        qrrho_entropy_times_temperature = (
            thermochemistry.qrrho_entropy_times_temperature
        )
        gibbs_free_energy = thermochemistry.gibbs_free_energy
        if q:
            qrrho_gibbs_free_energy = thermochemistry.qrrho_gibbs_free_energy
        elif qh:
            qrrho_gibbs_free_energy = (
                thermochemistry.qrrho_gibbs_free_energy_qh
            )
        elif qs:
            qrrho_gibbs_free_energy = (
                thermochemistry.qrrho_gibbs_free_energy_qs
            )

        if q:
            logger.info(
                "{:2} {:39} {:13.6f} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}".format(
                    "o",
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
            logger.info(
                "{:2} {:39} {:13.6f} {:10.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}".format(
                    "o",
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
            logger.info(
                "{:2} {:39} {:13.6f} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:13.6f} {:13.6f}".format(
                    "o",
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
            logger.info(
                "{:2} {:39} {:13.6f} {:10.6f} {:13.6f} {:10.6f} {:13.6f}".format(
                    "o",
                    structure,
                    energy,
                    zero_point_energy,
                    enthalpy,
                    entropy_times_temperature,
                    gibbs_free_energy,
                )
            )


if __name__ == "__main__":
    get_thermo()
