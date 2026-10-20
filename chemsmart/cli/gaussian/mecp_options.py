"""Gaussian MECP options and compatibility exports for shared state options."""

import click

from chemsmart.cli.mecp import (
    click_mecp_follow_controls,
    click_mecp_hessian_options,
    click_mecp_state_options,
)

__all__ = [
    "add_mecp_method_suffix",
    "click_mecp_frequency_options",
    "click_mecp_restart_option",
    "click_mecp_state_options",
    "click_mecp_step_size_method_option",
]


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
    function = click_mecp_hessian_options(
        method_default="numerical", step_default=None
    )(function)
    function = click_mecp_follow_controls(function)
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
            "--follow-seam-imaginary-mode/--no-follow-seam-imaginary-mode",
            "--followseam",
            default=False,
            show_default=True,
            help=(
                "If the projected MECP Hessian has a significant imaginary "
                "mode, displace in both directions, reoptimize both MECPs, "
                "and retain the lower seam minimum. Implies --mecp-numfreq."
            ),
        ),
    )
    for option in reversed(options):
        function = option(function)
    return function
