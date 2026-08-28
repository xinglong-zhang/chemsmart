"""Shared Click options for MECP entry points."""

import click


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
