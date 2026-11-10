"""Shared MECP state and seam-analysis options for engine entry points."""

import click


def _apply_options(function, options):
    for option in reversed(options):
        function = option(function)
    return function


def click_mecp_state_options(
    function=None, *, required=False, include_charges=True
):
    """Declare numbered state options; ORCA requires both and one shared charge."""

    def decorate(function):
        multiplicity_type = click.IntRange(min=1) if required else int
        options = [
            click.option(
                "--multiplicity1",
                "-m1",
                type=multiplicity_type,
                default=None,
                required=required,
                help="Spin multiplicity for state 1 (A).",
            ),
            click.option(
                "--multiplicity2",
                "-m2",
                type=multiplicity_type,
                default=None,
                required=required,
                help="Spin multiplicity for state 2 (B)."
                + ("" if required else " Defaults to multiplicity1 + 2."),
            ),
        ]
        if include_charges:
            options.extend(
                [
                    click.option(
                        "--charge1",
                        "-c1",
                        type=int,
                        default=None,
                        help="Charge for state 1 (A).",
                    ),
                    click.option(
                        "--charge2",
                        "-c2",
                        type=int,
                        default=None,
                        help="Charge for state 2 (B). Defaults to charge1.",
                    ),
                ]
            )
        return _apply_options(function, options)

    return decorate if function is None else decorate(function)


def click_mecp_hessian_options(
    *, method_default, step_default, step_type=float
):
    """Share names and units while retaining each entry point's defaults."""

    def decorate(function):
        return _apply_options(
            function,
            (
                click.option(
                    "--hessian-method",
                    type=click.Choice(
                        ["numerical", "analytic"], case_sensitive=False
                    ),
                    default=method_default,
                    show_default=True,
                    help="Source of both state Hessians, including every followseam check.",
                ),
                click.option(
                    "--hess-step-size",
                    type=step_type,
                    default=step_default,
                    show_default=step_default is not None,
                    help="Cartesian finite-difference displacement (Bohr; default 1e-3). Numerical Hessians only.",
                ),
            ),
        )

    return decorate


def click_mecp_follow_controls(function):
    return _apply_options(
        function,
        (
            click.option(
                "--seam-mode-displacement",
                "--displacement",
                type=click.FloatRange(min=0, min_open=True),
                default=0.05,
                show_default=True,
                help="Seam-following progress step (Angstrom).",
            ),
            click.option(
                "--seam-mode-max-steps",
                "--max-follow-steps",
                type=click.IntRange(min=1),
                default=30,
                show_default=True,
                help="Maximum constrained seam-following macro steps per direction.",
            ),
        ),
    )


def click_mecp_only_option(function):
    return click.option(
        "--only",
        is_flag=True,
        help="With numfreq or anafreq, analyze the supplied stationary crossing without optimization.",
    )(function)
