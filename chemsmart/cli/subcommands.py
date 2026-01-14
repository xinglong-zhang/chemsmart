from chemsmart.cli.assemble import assemble
from chemsmart.cli.gaussian import gaussian
from chemsmart.cli.mol import mol
from chemsmart.cli.nciplot import nciplot
from chemsmart.cli.orca import orca
from chemsmart.cli.thermochemistry import thermochemistry

subcommands = [
    gaussian,
    orca,
    mol,
    nciplot,
    thermochemistry,
    assemble,
]
