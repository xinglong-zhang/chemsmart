import functools
import logging
import os

import click

from chemsmart.analysis.thermochemistry import Thermochemistry
from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


def click_thermochemistry_options(f):
    """Common click options for Thermochemistry."""

    @click.option(
        "-d",
        "--directory",
        default=None,
        help="Directory in which to compute thermochemistry for all files.",
    )
    @click.option(
        "-t",
        "--filetype",
        default=None,
        help="Type of file to calculate thermochemistry for, if directory is specified.",
    )
    @click.option(
        "-f",
        "--filenames",
        type=str,
        multiple=True,
        default=None,
        help="Gaussian or ORCA output files for parsing thermochemistry.",
    )
    @click.option(
        "-cs",
        "--cutoff-entropy",
        default=None,
        type=float,
        show_default=True,
        help="Cutoff frequency for entropy in wavenumbers",
    )
    @click.option(
        "-ch",
        "--cutoff-enthalpy",
        default=None,
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
        "-p",
        "--pressure",
        default=1.0,
        type=float,
        show_default=True,
        help="Pressure in atm.",
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
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_thermochemistry_options
@click.pass_context
def thermochemistry(
    ctx,
    directory,
    filetype,
    filenames,
    cutoff_entropy,
    cutoff_enthalpy,
    concentration,
    pressure,
    temperature,
    alpha,
    weighted,
    units,
):
    if directory is not None:
        directory = os.path.abspath(directory)
        logger.info(
            f"Obtaining thermochemistry of files in directory: {directory}"
        )
        assert (
            filetype is not None
        ), "Type of files to calculate thermochemistry for should be given!"
        assert filenames is None, (
            "Filenames cannot be specified with directory!\n"
            "All files will be converted."
        )

        if filetype == "log":
            logger.info(
                "Obtaining thermochemistry of files in Gaussian .log files."
            )
            from chemsmart.io.gaussian.folder import GaussianLogFolder

            gaussian_log_folder = GaussianLogFolder(directory)
            gaussian_log_files = gaussian_log_folder.all_logfiles
            for file in gaussian_log_files:
                thermochemistry = Thermochemistry(
                    filename=file,
                    temperature=temperature,
                    concentration=concentration,
                    pressure=pressure,
                    use_weighted_mass=weighted,
                    alpha=alpha,
                    s_freq_cutoff=cutoff_entropy,
                    h_freq_cutoff=cutoff_enthalpy,
                )
                energy = thermochemistry.energies
                print(energy)
        elif filetype == "out":
            logger.info(
                "Obtaining thermochemistry of files in ORCA .out files."
            )
            from chemsmart.io.orca.folder import ORCAOutFolder

            orca_out_folder = ORCAOutFolder(directory)
            orca_out_files = orca_out_folder.all_outfiles
            for file in orca_out_files:
                thermochemistry = Thermochemistry(
                    filename=file,
                    temperature=temperature,
                    concentration=concentration,
                    pressure=pressure,
                    use_weighted_mass=weighted,
                    alpha=alpha,
                    s_freq_cutoff=cutoff_entropy,
                    h_freq_cutoff=cutoff_enthalpy,
                )
                energy = thermochemistry.energies
                print(energy)
        else:
            raise ValueError(
                f"Unsupported file extension for '{filetype}'\n. "
                f"Only Gaussian .log or ORCA .out files are accepted."
            )


@thermochemistry.result_callback()
@click.pass_context
def thermochemistry_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
