import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.io.orca.input import ORCAInput
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("inp", cls=MyCommand)
@click_job_options
@click.pass_context
def inp(ctx, skip_completed, **kwargs):
    """Run an ORCA input job as it is. Only requires the file that is to be run."""
    filename = ctx.obj["filename"]

    from chemsmart.jobs.orca.job import ORCAInpJob

    return ORCAInpJob.from_filename(filename=filename, skip_completed=skip_completed, **kwargs)