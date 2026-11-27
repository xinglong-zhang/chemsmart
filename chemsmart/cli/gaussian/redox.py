from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.redox import register_redox_cli
from chemsmart.jobs.gaussian.redox import (
    GaussianRedoxJob,
    GaussianRedoxJobSettings,
)

redox = register_redox_cli(
    gaussian, GaussianRedoxJob, GaussianRedoxJobSettings
)
