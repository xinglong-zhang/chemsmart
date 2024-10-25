#!/usr/bin/env python
# TODO
import logging
import os

import click

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.utils.utils import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-d",
    "--directory",
    default=None,
    help="Parent folder containing all the Gaussian/ORCA output files.",
)
@click.option(
    "-f",
    "--filename",
    default=None,
    help="Single Gaussian/ORCA file to be converted to xyz.",
)
@click.option(
    "-a/", "--append/--no-append", default=False, help="append to already existing file"
)
def entry_point(directory, filename, append):
    """Script to write Gaussian .log files or ORCA .out files to .xyz files, while keeping the same folder structure."""
    create_logger()
    # writing all .log files in a given folder to .xyz
    if directory is None:
        directory = "."

    g16_folder = GaussianLogFolder(folder=directory)
    all_logpaths = g16_folder.get_all_logs()

    # if input file is given, write only that file to .xyz
    if filename is not None:
        input_filepath = os.path.join(directory, filename)
        all_logpaths = os.path.abspath(input_filepath)
        all_logpaths = [all_logpaths]

    # if input file is not given
    for logpath in all_logpaths:
        logger.info(f"Converting file: {logpath}")
        logpath = os.path.abspath(logpath)
        logdir, logfile = os.path.split(logpath)
        basename = os.path.splitext(logfile)[0]
        xyzfile = f"{basename}.xyz"
        xyzpath = os.path.join(logdir, xyzfile)

        logoutput = Gaussian16Output(logfile=logpath)
        log_final_atoms = logoutput.get_atoms(include_failed_logfile=True)

        if log_final_atoms is None:
            logger.info(f"No readable structure from file {logpath}")
            continue

        log_empirical_formula = log_final_atoms.get_chemical_formula(empirical=True)
        log_final_energy_in_Hartree = round(
            logoutput.final_energy, 6
        )  # this is in Hartree from pymatgen
        # print(log_final_atoms.results()['energy']) # this is in eV as expected
        xyz_info = f"{logfile}    Empirical formula: {log_empirical_formula}    Energy(Hartree): {log_final_energy_in_Hartree}"
        logger.info(f"Writing outputfile to {xyzpath}")
        log_final_atoms.write(
            filename=xyzpath,
            format="xyz",
            comment=xyz_info,
            append=append,
            fmt="%15.6f",
        )


if __name__ == "__main__":
    entry_point()
