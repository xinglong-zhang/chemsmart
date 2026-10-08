"""Shared Click options for MECP entry points."""

import click


def click_mecp_state_options(function):
    """Apply the common charge and multiplicity options for two MECP states."""
    options = (
        click.option(
            "--multiplicity1",
            "-m1",
            type=int,
            default=None,
            help="Spin multiplicity for state 1.",
        ),
        click.option(
            "--multiplicity2",
            "-m2",
            type=int,
            default=None,
            help="Spin multiplicity for state 2. Defaults to multiplicity1 + 2.",
        ),
        click.option(
            "--charge1",
            "-c1",
            type=int,
            default=None,
            help="Charge for state 1.",
        ),
        click.option(
            "--charge2",
            "-c2",
            type=int,
            default=None,
            help="Charge for state 2. Defaults to charge1.",
        ),
    )
    for option in reversed(options):
        function = option(function)
    return function


def add_mecp_method_suffix(label, method):
    """Return a MECP label ending in ``_<method>_mecp``."""
    suffix = f"_{method}_mecp"
    if label.lower().endswith(suffix.lower()):
        return label
    if label.lower().endswith("_mecp"):
        return f"{label[:-5]}{suffix}"
    return f"{label}{suffix}"


def click_mecp_step_size_method_option(function):
    """Apply the common MECP optimizer-selection option."""
    return click.option(
        "--step-size-method",
        type=click.Choice(
            ["bb", "grow_shrink", "harvey"], case_sensitive=False
        ),
        default=None,
        help=(
            "MECP optimizer. 'harvey' uses the Harvey inverse-BFGS optimizer "
            "(default); 'bb' uses the Barzilai-Borwein secant rule; "
            "'grow_shrink' uses a merit-based grow/shrink rule."
        ),
    )(function)


def click_mecp_restart_option(function):
    """Apply the common interrupted-optimization restart option."""
    return click.option(
        "--restart/--no-restart",
        default=True,
        show_default=True,
        help="Resume an interrupted MECP optimization from its saved state.",
    )(function)


def click_mecp_frequency_options(function):
    """Apply the common post-MECP Hessian-analysis options."""
    options = (
        click.option(
            "--mecp-numfreq/--no-mecp-numfreq",
            default=False,
            show_default=True,
            help=(
                "After convergence, verify the seam minimum and write "
                "mass-weighted projected MECP frequencies and normal modes "
                "to <label>_mecp_freq.log."
            ),
        ),
        click.option(
            "--hess-step-size",
            type=float,
            default=None,
            help=(
                "Finite-difference step size (Bohr) for MECP Hessian "
                "analysis (default: 1e-3)."
            ),
        ),
        click.option(
            "--follow-seam-imaginary-mode/--no-follow-seam-imaginary-mode",
            default=False,
            show_default=True,
            help=(
                "If the projected MECP Hessian has a significant imaginary "
                "mode, displace in both directions, reoptimize both MECPs, "
                "and retain the lower seam minimum. Implies --mecp-numfreq."
            ),
        ),
        click.option(
            "--seam-mode-displacement",
            type=click.FloatRange(min=0.0, min_open=True),
            default=0.05,
            show_default=True,
            help=(
                "Cartesian norm (Angstrom) of each +/- projected-mode "
                "displacement used by --follow-seam-imaginary-mode."
            ),
        ),
        click.option(
            "--seam-mode-max-steps",
            type=click.IntRange(min=1),
            default=30,
            show_default=True,
            help=(
                "Maximum constrained mode-following steps in each +/- "
                "direction."
            ),
        ),
    )
    for option in reversed(options):
        function = option(function)
    return function
