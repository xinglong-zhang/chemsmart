"""Submission of jobs to queuing system via cli."""

import glob
import logging
import os

import click

from pyatoms.cli.defaults import DEFAULTS
from chemsmart.cli.jobrunner import jobrunner_options
from chemsmart.cli.logger import logger_options
from chemsmart.cli.subcommands import subcommands
from chemsmart.utils.logger import create_logger
from chemsmart.utils.cli import MyGroup
# from chemsmart.utils.cli import CtxObjArguments, MyGroup, determine_folders

logger = logging.getLogger(__name__)


# def prepare_job_from_source_folder(folder):
#     from pyatoms.analysis.results.results import NebOptimizationResults, Results
#
#     try:
#         results = Results.from_folder(folder)
#     except Exception as e:
#         raise ValueError(f'Could not create job from {folder}') from e
#
#     if not results.is_complete:
#         raise ValueError(f'Job from {folder} is not complete.')
#
#     atoms = results.transition_state_atoms if isinstance(results, NebOptimizationResults) else results.optimized_atoms
#
#     return atoms.set_calculated_magmoms_as_initial(error_if_no_magmoms=False)


@click.group(name='sub', cls=MyGroup)
@click.pass_context
@jobrunner_options
@logger_options
@click.option('-S', '--node-for-runner/--no-node-for-runner', default=DEFAULTS['node-for-runner'])
@click.option('-t', '--num-hours', type=float, default=None)
@click.option('-q', '--queue', type=str, default=DEFAULTS['queue'], help='queue')
@click.option('-Q', '--qos', type=str, default=None, help='qos')
@click.option('-V', '--verbose/--no-verbose', default=DEFAULTS['verbose'], help='turns on logging')
@click.option('-f', '--folders', '--folder', multiple=True, type=str, default=DEFAULTS['folders'], help='folder')
@click.option('--folderlist-file', type=click.Path(exists=True), default=None, help='file containing a list of folders')
@click.option('-D', '--dryrun/--no-dryrun', default=True, help='turns on dryrun')
@click.option(
    '-N',
    '--queue-manager-num-processes',
    type=int,
    default=DEFAULTS['queue-manager-num-processes'],
    help='runs jobs in parallel with N processes',
)
@click.option('-l', '--label', type=str, default=None)
@click.option('--multi/--no-multi', default=False, help='Submit all folders as separate jobs.')
@click.option(
    '-F',
    '--source-folder',
    type=click.Path(exists=True, resolve_path=True),
    default=DEFAULTS['source-folder'],
    help='read atoms from source folder',
)
@click.option(
    '--print-command/--no-print-command', default=DEFAULTS['print-command'], help='print the command generated'
)
@click.option('--test/--no-test', default=False, help='If true, job will not be submitted')
def sub(
    ctx,
    folders,
    folderlist_file,
    num_nodes,
    num_gpus,
    servername,
    num_gpus_per_node_per_job,
    exclude_localhost,
    num_cpus,
    source_folder,
    verbose,
    multi,
    queue_manager_num_processes,
    memory_gigs,
    test,
    **kwargs,
):
    # Set up logging
    create_logger(debug=debug, stream=stream)
    logger.info("Entering main program")

    if verbose:
        create_logger(stream=True, debug=True)

    if multi and queue_manager_num_processes is not None:
        raise ValueError(
            'Cannot use --multi and --queue-manager-num-processes together. '
            'Please choose one of them.\n'
            '    --multi submits all specified folders as single jobs.\n'
            '    -N/--queue-manager-num-processes submits all specified folders as a single job that is '
            'run in parallel by N processes.'
        )

    if len(folders) == 0 and folderlist_file is None:
        folders = glob.glob('*/') if multi else [os.getcwd()]

    folders = determine_folders(folders=None, folder_regex=folders, folderlist_file=folderlist_file)

    if len(folders) == 0:
        raise ValueError(f'Invalid folders specified: folders: {folders}')

    jobrunner = create_jobrunner(
        servername=servername,
        num_nodes=num_nodes,
        num_cpus=num_cpus,
        num_gpus=num_gpus,
        num_gpus_per_node_per_job=num_gpus_per_node_per_job,
        exclude_localhost=exclude_localhost,
        use_host_queues=False,  # no need host queues as we are not running
        memory_gigs=memory_gigs,
    )

    ctx.obj['jobrunner'] = jobrunner
    ctx.obj['source_folder'] = source_folder
    ctx.obj['other_queue_managers'] = {}
    ctx.obj['folder'] = folders
    ctx.obj['folders'] = folders if isinstance(folders, list) else [folders]

    if source_folder:
        atoms = prepare_job_from_source_folder(folder=source_folder)
        ctx.obj['atoms'] = atoms


@sub.result_callback(replace=True)
@click.pass_context
def process_pipeline(ctx, *args, qos, multi, **kwargs):  # noqa: PLR0915
    from pyatoms.jobs.job import JobBatch, Procedure
    from pyatoms.settings.server.server import Server
    from pyatoms.utils.processors import GlobalProcessStore

    def _clean_command(ctx):
        """Remove keywords used in sub.py but not in run.py."""
        # Get "sub" command
        command = [subcommand for subcommand in ctx.obj['subcommand'] if subcommand['name'] == 'sub']
        assert len(command) == 1
        command = command[0]
        keywords_not_in_run = [
            'print_command',
            'queue',
            'num_hours',
            'verbose',
            'node_for_runner',
            'test',
            'qos',
            'dryrun',
            'multi',
            'label',
            'folderlist_file',
        ]

        for keyword in keywords_not_in_run:
            del command['kwargs'][keyword]
        return ctx

    def _dry_run(job, jobrunner):
        if isinstance(job, JobBatch):
            for j in job:
                _dry_run(j, jobrunner)
        elif isinstance(job, Procedure):
            for j in job.steps:
                _dry_run(j, jobrunner)
        else:
            jobrunner.dry_run = True
            jobrunner.run(job)

    def _reconstruct_cli_args(ctx, multi, job):
        """Get cli args that reconstruct the command line."""
        commands = ctx.obj['subcommand']

        if multi:
            # Jobs should be submitted within each folder, so we need to edit the
            # "folders" value for the "sub" subcommand
            command = [subcommand for subcommand in commands if subcommand['name'] == 'sub']
            assert len(command) == 1
            command = command[0]
            command['kwargs']['folders']['value'] = '.'

        args = CtxObjArguments(commands, entry_point='sub')
        cli_args = args.reconstruct_command_line()[1:]  # remove the first element 'sub'
        if kwargs.get('print_command'):
            print(cli_args)
        return cli_args

    def _process_single_job(job, jobrunner, multi):
        # Perform a dry run to catch any user input errors
        if kwargs['dryrun']:
            _dry_run(job=job, jobrunner=jobrunner)

        if kwargs.get('test'):
            logger.warning('Not submitting as "test" flag specified.')

        submitscript_kwargs = {}
        if qos is not None:
            submitscript_kwargs.update({'qos': qos})

        cli_args = _reconstruct_cli_args(ctx, multi, job)
        # If --multi flag is specified, submit each job in its own folder

        server = Server.from_servername(kwargs.get('servername'))
        server.submit(
            job=job,
            label=kwargs.get('label'),
            submit_in_job_folder=multi,
            num_nodes=kwargs.get('num_nodes'),
            num_hours=kwargs.get('num_hours'),
            num_gpus=kwargs.get('num_gpus'),
            num_cpus=kwargs.get('num_cpus'),
            memory_gigs=kwargs.get('memory_gigs'),
            node_for_runner=kwargs.get('node_for_runner'),
            exclude_localhost=kwargs.get('exclude_localhost'),
            queue=kwargs.get('queue'),
            queue_manager_num_processes=kwargs.get('queue_manager_num_processes'),
            variable_time_minimum=kwargs.get('variable_time_minimum'),
            test=kwargs.get('test'),
            cli_args=cli_args,
            **submitscript_kwargs,
        )

    ctx = _clean_command(ctx)
    jobrunner = ctx.obj['jobrunner']
    job = args[0]

    if not isinstance(job, list):
        job = [job]

    batch = len(job) > 1
    if batch and not multi:
        job = [JobBatch(jobs=job)]

    for j in job:
        _process_single_job(job=j, jobrunner=jobrunner, multi=multi)

    for qm, processes in ctx.obj['other_queue_managers'].items():
        qm.close_queue()
        qm.join_processes(processes)

    logger.debug('Closing global processes')
    store = GlobalProcessStore()
    store.join_processes()


for subcommand in subcommands:
    sub.add_command(subcommand)
