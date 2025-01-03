"""CLI options for all jobs that can be run in this package."""

def click_job_options(f):
    import functools

    import click

    @click.option(
        '-S/-R',
        '--skip-completed/--no-skip-completed',
        is_flag=True,
        default=True,
        type=bool,
        help='To run completed job again. Use -R to rerun completed job.',
    )
    @click.option('--timeout-seconds', type=float)
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options
