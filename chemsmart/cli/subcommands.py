from chemsmart.cli.gaussian import gaussian
from chemsmart.cli.mol import mol
from chemsmart.cli.nciplot import nciplot
from chemsmart.cli.orca import orca
from chemsmart.cli.thermochemistry import thermochemistry
from chemsmart.cli.iterate import iterate

# Subcommands available for both 'run' and 'sub'
subcommands = [
    gaussian,
    orca,
    mol,
    nciplot,
    thermochemistry,
]

# Subcommands only available for 'run' (not submittable to queue)
run_subcommands = [
    iterate,
]
