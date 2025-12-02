from chemsmart.cli.gaussian import gaussian
from chemsmart.cli.mol import mol
from chemsmart.cli.nciplot import nciplot
from chemsmart.cli.orca import orca
from chemsmart.cli.thermochemistry import thermochemistry
from chemsmart.cli.iterate import iterate

subcommands = [
    gaussian,
    orca,
    mol,
    nciplot,
    thermochemistry,
    iterate,
]
