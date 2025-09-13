#!/usr/bin/env python
"""
DIAS (Density-based Index of Aromaticity) plotting script.

This script creates plots and visualizations for DIAS aromaticity analysis
from Gaussian and ORCA calculations, helping to analyze aromatic character
in molecular systems.
"""

import click

from chemsmart.analysis.dias import GaussianDIASLogFolder, ORCADIASOutFolder
from chemsmart.utils.logger import create_logger


@click.command()
@click.option("-f", "--folder", default=".")
@click.option(
    "-p",
    "--program",
    default="gaussian",
    help="Type of program: Gaussian or ORCA output for plotting.",
)
@click.option(
    "-z/",
    "--zero/--no-zero",
    default=False,
    help="Set reference point all to zero.",
)
@click.option(
    "-o",
    "--outputname",
    default="dias",
    help="Output file name for DIAS data to be written",
)
@click.option("-e", "--extrapolate", default=True)
@click.option("-r", "--reverse/--no-reverse", default=False)
@click.option("-n", "--new-length", default=1000)
@click.option(
    "-k",
    "--k-value",
    type=int,
    default=3,
    help="Degree of the smoothing spline. Must be 1 <= k <= 5. "
         "k = 3 is a cubic spline. Default is 3.",
)
@click.option(
    "-a",
    "--atom-number1",
    type=int,
    required=True,
    default=None,
    help="Atom number 1 for the bond distance to plot on x-axis.",
)
@click.option(
    "-b",
    "--atom-number2",
    type=int,
    required=True,
    default=None,
    help="Atom number 2 for the bond distance to plot on x-axis.",
)
@click.option(
    "-x",
    "--ref-file",
    type=str,
    default=None,
    help="Reference file as zero for DIAS plot.",
)
def entry_point(
    folder,
    program,
    zero,
    outputname,
    extrapolate,
    reverse,
    new_length,
    k_value,
    atom_number1,
    atom_number2,
    ref_file,
):
    """
    Analyze and plot DIAS data from quantum chemistry calculations.
    
    Example usage: plot_dias.py -p orca -z -a 5 -b 7 -r
    """
    create_logger(debug=True, stream=True)
    
    # Initialize DIAS analysis based on program type
    if program.lower() == "gaussian":
        dias_folder = GaussianDIASLogFolder(
            folder=folder,
            atom1=atom_number1,
            atom2=atom_number2,
            zero=zero,
            outputname=f"gaussian_{outputname}",
            ref_file=ref_file,
        )
    elif program.lower() == "orca":
        dias_folder = ORCADIASOutFolder(
            folder=folder,
            atom1=atom_number1,
            atom2=atom_number2,
            zero=zero,
            outputname=f"orca_{outputname}",
            ref_file=ref_file,
        )
    else:
        raise TypeError(
            "Unknown program output files to plot DIAS. "
            "Supported outputs are from Gaussian or ORCA."
        )
        
    # Generate DIAS data and plots
    dias_folder.write_data()
    dias_folder.plot_dias(
        extrapolate=extrapolate,
        reversed=reverse,
        new_length=new_length,
        k=k_value,
    )


if __name__ == "__main__":
    entry_point()
