from chemsmart.cli.orca.orca import orca
from chemsmart.cli.redox import register_redox_cli
from chemsmart.jobs.orca.redox import ORCARedoxJob, ORCARedoxJobSettings

redox = register_redox_cli(orca, ORCARedoxJob, ORCARedoxJobSettings)
