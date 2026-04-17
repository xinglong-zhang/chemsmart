from chemsmart.cli.gaussian import gaussian
from chemsmart.cli.grouper import grouper
from chemsmart.cli.iterate import iterate
from chemsmart.cli.mol import mol
from chemsmart.cli.nciplot import nciplot
from chemsmart.cli.orca import orca
from chemsmart.cli.thermochemistry import thermochemistry

subcommands = [
    gaussian,
    grouper,
    orca,
    mol,
    nciplot,
    thermochemistry,
    iterate,
]
