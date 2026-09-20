from chemsmart.analysis.fukui import FUKUI_MODES
from chemsmart.cli.fukui import register_fukui_cli
from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.jobs.gaussian.fukui import GaussianFukuiJob

fukui = register_fukui_cli(
    gaussian, GaussianFukuiJob, modes=FUKUI_MODES, nbo_uses_wbi=True
)
