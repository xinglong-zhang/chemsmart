#!/usr/bin/env python
import logging
import os

import click

from chemsmart.io.gaussian.output import Gaussian16Output, GaussianLogFolder
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


@click.command()
@click.option(
    "-d",
    "--directory",
    default=None,
    help="Directory in which to convert files.",
)
@click.option(
    "-t",
    "--type",
    default=None,
    help="Type of file to be converted, if direcotry is specified.",
)
@click.option(
    "-f", "--filename", default=None, help="Input filename to be converted."
)
@click.option(
    "-o",
    "--output-filetype",
    default="xyz",
    help="Type of files to convert to, defaults to .xzy",
)
def entry_point(directory, type, filename, output_filetype):
    """Script for converting structures in different formats."""
    create_logger()
    if directory is not None:
        logger.info(f"Converting files in directory: {directory}")
        assert (
            type is not None
        ), "Type of file to be converted must be specified."
        convert_files(directory, type, output_filetype)
    else:
        logger.info(f"Converting file: {filename}")
        convert_files(directory, filename, output_filetype)

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

        log_empirical_formula = log_final_atoms.get_chemical_formula(
            empirical=True
        )
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
