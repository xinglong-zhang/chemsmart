from chemsmart.analysis.fukui import ORCA_FUKUI_MODES
from chemsmart.cli.fukui import register_fukui_cli
from chemsmart.cli.orca.orca import orca
from chemsmart.jobs.orca.fukui import ORCAFukuiJob

fukui = register_fukui_cli(orca, ORCAFukuiJob, modes=ORCA_FUKUI_MODES)
