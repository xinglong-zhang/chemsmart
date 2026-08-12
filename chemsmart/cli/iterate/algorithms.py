"""Shared Iterate algorithm CLI subcommands."""

import click

from chemsmart.utils.cli import MyCommand

ALGORITHM_OVERRIDE_HELP = """\
The algorithm subcommand selects the final algorithm. If the input config
selects the same algorithm, CLI options that are not passed keep the input
config values. If the CLI switches to a different algorithm, options from the
other algorithm are discarded and this algorithm starts from its defaults.
Explicit CLI options always have the highest priority.
"""

JLGO_HELP = f"""\
Run Joint Lagrange Geometry Optimization (JLGO).

JLGO attaches one or more substituents in a single joint optimization. It is
useful when substituent placements should be optimized together instead of
independently sampled by ETKDG.

{ALGORITHM_OVERRIDE_HELP}
Example:

\b
chemsmart run iterate yaml -f config.yaml jlgo --max-starts 16
chemsmart run iterate direct -skf core.gjf -skg '[1]' -subf Me.gjf -subi 1 -subg '[1]' jlgo --max-starts 16
"""

ETKDG_HELP = f"""\
Optimize substituent positions with the RDKit ETKDGv3 algorithm.

ETKDG embeds candidate substituent geometries with RDKit distance geometry. The
default local mode keeps skeleton atoms fixed and re-embeds the substituent;
``--global`` re-embeds the whole molecule.

{ALGORITHM_OVERRIDE_HELP}
Example:

\b
chemsmart run iterate yaml -f config.yaml etkdg --num-conformers 50
chemsmart run iterate direct -skf core.gjf -skg '[1]' -subf Me.gjf -subi 1 -subg '[1]' etkdg --num-conformers 50
"""


def _collect_cli_options(options: dict) -> dict:
    """Drop unset CLI options before merging them over input settings."""
    return {
        name: value for name, value in options.items() if value is not None
    }


def register_iterate_algorithm_commands(parent, build_job):
    """Register shared Iterate algorithm subcommands on an input group."""

    @parent.command(name="jlgo", cls=MyCommand, help=JLGO_HELP)
    @click.option(
        "--adaptive-sampling/--no-adaptive-sampling",
        default=None,
        help="Run a fixed coarse sampling stage first; the six full-stage "
        "sampling/pruning options (--link-sphere-samples, "
        "--orientation-sphere-samples, --axial-samples, --candidate-pool-size, "
        "--preselect, --beam-width) only take effect when the coarse stage does "
        "not produce an acceptable optimized structure. Enabled by default. "
        "--max-starts and "
        "--slsqp-maxiter always apply. Use --no-adaptive-sampling to always "
        "apply the full sampling parameters.",
    )
    @click.option(
        "--link-sphere-samples",
        default=None,
        type=int,
        help="Full-stage number of linking-atom bond-sphere position samples "
        "(default: 48).",
    )
    @click.option(
        "--orientation-sphere-samples",
        default=None,
        type=int,
        help="Full-stage number of substituent principal-axis direction samples "
        "(default: 24).",
    )
    @click.option(
        "--axial-samples",
        default=None,
        type=int,
        help="Number of axial rotations per orientation direction (default: 4).",
    )
    @click.option(
        "--candidate-pool-size",
        default=None,
        type=int,
        help="Per-substituent candidate pool size kept after region exclusion "
        "(default: 20).",
    )
    @click.option(
        "--preselect",
        default=None,
        type=int,
        help="Top joint combinations fed into greedy start selection "
        "(default: 48).",
    )
    @click.option(
        "--beam-width",
        default=None,
        type=int,
        help="Beam width retained per layer during feasible-domain pruning "
        "(default: 4096).",
    )
    @click.option(
        "--max-starts",
        default=None,
        type=int,
        help="Maximum number of 6K-dimensional joint starts handed to SLSQP "
        "(default: 8).",
    )
    @click.option(
        "--slsqp-maxiter",
        default=None,
        type=int,
        help="Maximum SLSQP iterations per start (default: 200).",
    )
    @click.pass_context
    def jlgo(
        ctx,
        adaptive_sampling,
        link_sphere_samples,
        orientation_sphere_samples,
        axial_samples,
        candidate_pool_size,
        preselect,
        beam_width,
        max_starts,
        slsqp_maxiter,
    ):
        """Run Joint Lagrange Geometry Optimization (JLGO)."""
        cli_options = _collect_cli_options(
            {
                "use_adaptive_sampling": adaptive_sampling,
                "n_link_sphere": link_sphere_samples,
                "n_orientation_sphere": orientation_sphere_samples,
                "n_axial": axial_samples,
                "candidate_pool_size": candidate_pool_size,
                "preselect": preselect,
                "beam_width": beam_width,
                "max_starts": max_starts,
                "slsqp_maxiter": slsqp_maxiter,
            },
        )
        return build_job(
            ctx,
            cli_algorithm_name="jlgo",
            cli_options=cli_options,
        )

    @parent.command(name="etkdg", cls=MyCommand, help=ETKDG_HELP)
    @click.option(
        "--global/--local",
        default=None,
        help="Embedding mode. 'local' (default) keeps the skeleton fixed and "
        "only re-embeds the substituent; 'global' re-embeds every atom.",
    )
    @click.option(
        "--num-conformers",
        default=None,
        type=click.IntRange(min=1),
        help="Number of ETKDG conformers to try per attachment; the "
        "lowest-energy one is kept (default: 10).",
    )
    @click.option(
        "--random-seed",
        default=None,
        type=int,
        help="Base RDKit random seed (-1 for a non-reproducible random seed; "
        "default: 42).",
    )
    @click.option(
        "--max-iterations",
        default=None,
        type=click.IntRange(min=0),
        help="Maximum ETKDG embedding iterations (0 uses the RDKit default; "
        "default: 2000).",
    )
    @click.option(
        "--random-coords/--no-random-coords",
        default=None,
        help="Start embedding from random coordinates (usually more robust; "
        "enabled by default).",
    )
    @click.option(
        "--enforce-chirality/--no-enforce-chirality",
        default=None,
        help="Enforce the input chirality during embedding (disabled by default).",
    )
    @click.option(
        "--force-field",
        default=None,
        type=click.Choice(
            ["none", "uff", "mmff94", "mmff94s"],
            case_sensitive=False,
        ),
        help="Optional force-field post-optimization after embedding "
        "(default: none).",
    )
    @click.pass_context
    def etkdg(
        ctx,
        num_conformers,
        random_seed,
        max_iterations,
        random_coords,
        enforce_chirality,
        force_field,
        **cli_values,
    ):
        """Optimize substituent positions with the RDKit ETKDGv3 algorithm."""
        cli_options = _collect_cli_options(
            {
                "use_global_optimization": cli_values["global"],
                "num_conformers": num_conformers,
                "random_seed": random_seed,
                "max_iterations": max_iterations,
                "use_random_coordinates": random_coords,
                "enforce_chirality": enforce_chirality,
                "force_field": force_field,
            },
        )
        return build_job(
            ctx,
            cli_algorithm_name="etkdg",
            cli_options=cli_options,
        )
