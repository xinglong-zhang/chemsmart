from __future__ import annotations

import copy
import logging
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Callable, List, Optional, Protocol, Tuple

if TYPE_CHECKING:
    from chemsmart.io.molecules.structure import Molecule

    # One substituent attachment as passed to an analyzer builder:
    # (substituent, skeleton_link_index, substituent_link_index), with
    # 1-based link indices.
    SubstituentAssignment = Tuple[Molecule, int, int]

    class AnalyzerProtocol(Protocol):
        """Structural type implemented by every iterate analyzer.

        An analyzer joins a skeleton and its substituents into one molecule
        and returns the combined structure (or ``None`` on failure).
        """

        def run(self) -> Optional[Molecule]: ...

    # An analyzer builder turns a resolved options dict plus the skeleton and
    # its substituent assignments into a ready-to-run analyzer.
    AnalyzerBuilder = Callable[
        [dict, Molecule, List[SubstituentAssignment]], AnalyzerProtocol
    ]

logger = logging.getLogger(__name__)

# Default algorithm used when neither YAML nor CLI specify one.
DEFAULT_ALGORITHM_NAME = "etkdg"

# Valid slot-combination strategies for IterateJobSettings.
_VALID_COMBINATION_MODES = ("independent", "global")


@dataclass
class IterateAlgorithmConfig:
    """Lightweight algorithm configuration for iterate jobs.

    Attributes
    ----------
    name : str
        Canonical algorithm name (e.g. ``etkdg``).
    options : dict
        Algorithm-specific options (e.g. ``n_link_sphere``).
    """

    name: str = DEFAULT_ALGORITHM_NAME
    options: dict = field(default_factory=dict)

    def copy(self) -> "IterateAlgorithmConfig":
        """Return a deep copy of this configuration."""
        return copy.deepcopy(self)


@dataclass(frozen=True)
class AlgorithmSpec:
    """Specification describing a single iterate algorithm.

    Attributes
    ----------
    canonical_name : str
        Canonical (normalized) name of the algorithm.
    aliases : tuple[str, ...]
        Accepted aliases (including the canonical name) used to look up
        the spec. Matching is case-insensitive.
    default_options : dict
        Mapping of option name to default value. Also defines the set of
        options that are valid for the algorithm.
    option_validators : dict
        Mapping of option name to a callable ``validator(value)`` that
        raises :class:`ValueError` when the value is invalid.
    analyzer_builder : AnalyzerBuilder
        Factory ``builder(options, skeleton, substituents)`` returning the
        analyzer for a combination, where ``substituents`` is a list of
        ``(substituent, skeleton_link, substituent_link)`` tuples (1-based
        link indices). The analyzer joins the skeleton and every substituent
        (one or more) into one molecule and produces the combined structure
        in a single ``run()``. Required: every registered algorithm must
        provide a builder (there are no placeholder algorithms).
    """

    canonical_name: str
    aliases: tuple
    default_options: dict
    analyzer_builder: AnalyzerBuilder
    option_validators: dict = field(default_factory=dict)

    def __post_init__(self):
        # Freeze the mutable mappings so a registered spec cannot be
        # mutated in place (which would leak into every later resolution).
        object.__setattr__(
            self,
            "default_options",
            MappingProxyType(dict(self.default_options)),
        )
        object.__setattr__(
            self,
            "option_validators",
            MappingProxyType(dict(self.option_validators)),
        )
        # Every executable algorithm must ship a builder.
        if self.analyzer_builder is None:
            raise ValueError(
                f"AlgorithmSpec '{self.canonical_name}': analyzer_builder is "
                f"required."
            )
        # The canonical name must be one of its own aliases so lookups by
        # canonical name always resolve.
        if self.canonical_name.lower() not in {
            alias.lower() for alias in self.aliases
        }:
            raise ValueError(
                f"AlgorithmSpec '{self.canonical_name}': canonical name must "
                f"be included in its aliases {self.aliases}."
            )
        # Every default option must have a validator so no option can be set
        # without a defined validation policy.
        options_without_validator = set(self.default_options) - set(
            self.option_validators
        )
        if options_without_validator:
            raise ValueError(
                f"AlgorithmSpec '{self.canonical_name}': options "
                f"{sorted(options_without_validator)} have no validator."
            )


# --- Option validators (shared across algorithms) ------------------------


def _validate_positive_int(value) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError("expected an integer >= 1")


def _validate_non_negative_int(value) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("expected an integer >= 0")


# RDKit stores the ETKDG random seed in a C++ 32-bit signed int; values
# outside this range raise OverflowError when handed to RDKit.
_RANDOM_SEED_MIN = -(2**31)
_RANDOM_SEED_MAX = 2**31 - 1


def _validate_random_seed(value) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError("expected an integer")
    if not (_RANDOM_SEED_MIN <= value <= _RANDOM_SEED_MAX):
        raise ValueError(
            f"expected an integer within "
            f"[{_RANDOM_SEED_MIN}, {_RANDOM_SEED_MAX}]"
        )


def _validate_bool(value) -> None:
    if not isinstance(value, bool):
        raise ValueError("expected a boolean")


def _validate_force_field(value) -> None:
    choices = ("none", "uff", "mmff94", "mmff94s")
    if not isinstance(value, str) or value.lower() not in choices:
        raise ValueError(f"expected one of {list(choices)}")


# --- Analyzer builders (lazy import to keep this module lightweight) ------


def _build_lagrange_analyzer(
    options,
    skeleton,
    substituents,
):
    from chemsmart.jobs.iterate.jlgo import (
        IterateJointLagrangeAnalyzer,
    )

    # The Joint-Lagrange optimizer consumes the full substituent list in one
    # 6K-dimensional optimization; K=1 and K>1 share this single path.
    return IterateJointLagrangeAnalyzer(
        skeleton=skeleton,
        substituents=substituents,
        options=options,
    )


def _build_etkdg_analyzer(
    options,
    skeleton,
    substituents,
):
    from chemsmart.jobs.iterate.etkdg import IterateETKDGAnalyzer

    return IterateETKDGAnalyzer(
        skeleton=skeleton,
        substituents=substituents,
        options=options,
    )


def _build_alias_map(specs):
    """Build the lowercased alias -> spec map, rejecting inconsistencies.

    Raises
    ------
    ValueError
        If two specs share a canonical name, or if the same alias is
        declared by more than one spec.
    """
    alias_map = {}
    seen_canonical = set()
    for spec in specs:
        canonical = spec.canonical_name.lower()
        if canonical in seen_canonical:
            raise ValueError(
                f"Duplicate canonical algorithm name: '{spec.canonical_name}'."
            )
        seen_canonical.add(canonical)
        for alias in spec.aliases:
            key = alias.lower()
            if key in alias_map:
                raise ValueError(
                    f"Duplicate algorithm alias '{alias}' (already registered "
                    f"by '{alias_map[key].canonical_name}')."
                )
            alias_map[key] = spec
    return alias_map


# Registry of supported algorithms. Adding a new algorithm requires
# registering a new AlgorithmSpec here (with its option validators and
# analyzer builder) and, for CLI exposure, a subcommand in
# chemsmart/cli/iterate/yaml_cmd.py.
_ALGORITHM_SPECS = (
    AlgorithmSpec(
        canonical_name="jlgo",
        aliases=("jlgo", "lagrange", "lagrange_multipliers"),
        default_options={
            # Adaptive sampling first runs a fixed coarse stage; the six
            # full-stage sampling/pruning parameters (n_link_sphere,
            # n_orientation_sphere, n_axial, candidate_pool_size, preselect,
            # beam_width) only take effect when the coarse stage does not
            # produce an acceptable optimized structure. max_starts and
            # slsqp_maxiter always apply, in either stage.
            "use_adaptive_sampling": True,
            # Full-stage sampling of linking-atom bond-sphere positions.
            "n_link_sphere": 48,
            # Full-stage sampling of substituent principal-axis directions.
            "n_orientation_sphere": 24,
            # Axial rotations per orientation direction.
            "n_axial": 4,
            # Per-substituent candidate pool size after region exclusion.
            "candidate_pool_size": 20,
            # Top joint combinations fed into greedy start selection.
            "preselect": 48,
            # Beam width retained per layer during feasible pruning.
            "beam_width": 4096,
            # Maximum 6K-dimensional joint starts handed to SLSQP.
            "max_starts": 8,
            # Maximum SLSQP iterations per start.
            "slsqp_maxiter": 200,
        },
        option_validators={
            "use_adaptive_sampling": _validate_bool,
            "n_link_sphere": _validate_positive_int,
            "n_orientation_sphere": _validate_positive_int,
            "n_axial": _validate_positive_int,
            "candidate_pool_size": _validate_positive_int,
            "preselect": _validate_positive_int,
            "beam_width": _validate_positive_int,
            "max_starts": _validate_positive_int,
            "slsqp_maxiter": _validate_positive_int,
        },
        analyzer_builder=_build_lagrange_analyzer,
    ),
    AlgorithmSpec(
        canonical_name="etkdg",
        aliases=("etkdg",),
        default_options={
            # Embedding mode: local (default) fixes the skeleton atoms and
            # only regenerates the substituent; global re-embeds every atom.
            "use_global_optimization": False,
            # Number of ETKDG conformers embedded per attachment; the
            # lowest-energy one is kept.
            "num_conformers": 10,
            # RDKit random seed (a fixed value keeps results reproducible).
            "random_seed": 42,
            # Maximum ETKDG embedding iterations per conformer.
            "max_iterations": 2000,
            # Start embedding from random coordinates (more robust).
            "use_random_coordinates": True,
            # Enforce input chirality during embedding.
            "enforce_chirality": False,
            # Optional force-field refinement: none / uff / mmff94 / mmff94s.
            "force_field": "none",
        },
        option_validators={
            "use_global_optimization": _validate_bool,
            "num_conformers": _validate_positive_int,
            "random_seed": _validate_random_seed,
            "max_iterations": _validate_non_negative_int,
            "use_random_coordinates": _validate_bool,
            "enforce_chirality": _validate_bool,
            "force_field": _validate_force_field,
        },
        analyzer_builder=_build_etkdg_analyzer,
    ),
)

# Lowercased alias -> AlgorithmSpec.
_ALGORITHM_BY_ALIAS = _build_alias_map(_ALGORITHM_SPECS)


def available_algorithm_names() -> list:
    """Return the list of canonical algorithm names."""
    return [spec.canonical_name for spec in _ALGORITHM_SPECS]


def get_algorithm_spec(name: str) -> AlgorithmSpec:
    """Look up an :class:`AlgorithmSpec` by canonical name or alias.

    Parameters
    ----------
    name : str
        Algorithm name or alias (case-insensitive).

    Returns
    -------
    AlgorithmSpec
        The matching spec.

    Raises
    ------
    ValueError
        If the name does not match any registered algorithm.
    """
    if name is None:
        raise ValueError("Algorithm name cannot be None.")
    spec = _ALGORITHM_BY_ALIAS.get(str(name).strip().lower())
    if spec is None:
        raise ValueError(
            f"Unknown algorithm '{name}'. "
            f"Available algorithms: {available_algorithm_names()}."
        )
    return spec


def normalize_algorithm_name(name: str) -> str:
    """Normalize an algorithm name/alias to its canonical form.

    E.g. ``lagrange_multipliers`` -> ``jlgo``, ``etkdg`` -> ``etkdg``.
    """
    return get_algorithm_spec(name).canonical_name


def validate_algorithm_options(name: str, options: dict) -> dict:
    """Validate option keys and values for an algorithm.

    Parameters
    ----------
    name : str
        Algorithm name or alias.
    options : dict
        Options to validate.

    Returns
    -------
    dict
        A shallow copy of ``options`` (all keys recognized and values valid).

    Raises
    ------
    ValueError
        If ``options`` contains an unknown key, or a value that fails the
        algorithm's type/range/enum validator.
    """
    spec = get_algorithm_spec(name)
    if not options:
        return {}
    unknown = set(options) - set(spec.default_options)
    if unknown:
        raise ValueError(
            f"Unknown option(s) for algorithm '{spec.canonical_name}': "
            f"{sorted(unknown)}. "
            f"Allowed options: {sorted(spec.default_options)}."
        )
    for key, value in options.items():
        validator = spec.option_validators.get(key)
        if validator is None:
            continue
        try:
            validator(value)
        except ValueError as exc:
            raise ValueError(
                f"Invalid value for option '{key}' of algorithm "
                f"'{spec.canonical_name}': {exc} (got {value!r})."
            )
    return dict(options)


def resolve_algorithm_config(
    yaml_algorithm: dict = None,
    cli_algorithm_name: str = None,
    cli_options: dict = None,
) -> IterateAlgorithmConfig:
    """Merge built-in defaults, YAML config, and CLI overrides.

    Priority (lowest to highest):

    1. Built-in default (``etkdg`` + spec default options).
    2. YAML ``algorithm`` block (name and options).
    3. CLI algorithm subcommand (name only).
    4. CLI explicitly-provided options.

    If the CLI selects a different algorithm than the one declared in YAML,
    the YAML options (which belong to the other algorithm) are discarded
    rather than mixed into the newly selected algorithm.

    Parameters
    ----------
    yaml_algorithm : dict, optional
        Normalized YAML algorithm block with keys ``name`` and ``options``.
    cli_algorithm_name : str, optional
        Algorithm name selected via the CLI subcommand.
    cli_options : dict, optional
        Options explicitly provided on the command line (already filtered
        to only those the user actually passed).

    Returns
    -------
    IterateAlgorithmConfig
        The fully resolved algorithm configuration.
    """
    yaml_name = None
    yaml_options: dict = {}
    if yaml_algorithm:
        yaml_name = yaml_algorithm.get("name")
        yaml_options = yaml_algorithm.get("options") or {}

    # 1. Resolve the final algorithm name (CLI > YAML > default).
    if cli_algorithm_name is not None:
        final_name = normalize_algorithm_name(cli_algorithm_name)
    elif yaml_name is not None:
        final_name = normalize_algorithm_name(yaml_name)
    else:
        final_name = DEFAULT_ALGORITHM_NAME

    spec = get_algorithm_spec(final_name)

    # 2. Start from spec defaults.
    options = dict(spec.default_options)

    # 3. Merge YAML options only if YAML targets the same algorithm.
    if (
        yaml_name is not None
        and normalize_algorithm_name(yaml_name) == final_name
    ):
        options.update(validate_algorithm_options(final_name, yaml_options))

    # 4. Merge explicit CLI options (highest priority).
    if cli_options:
        options.update(validate_algorithm_options(final_name, cli_options))

    return IterateAlgorithmConfig(name=final_name, options=options)


class IterateJobSettings:

    def __init__(
        self,
        config_file=None,
        algorithm_config=None,
        combination_mode="independent",
    ):
        """
        Initialize iterate job settings.

        Parameters
        ----------
        config_file : str, optional
            Path to the YAML configuration file.
        algorithm_config : IterateAlgorithmConfig, optional
            Resolved algorithm configuration (name + options). When omitted,
            it falls back to the built-in default (``etkdg``
            with its default options). This is the single source of truth for
            the algorithm and its parameters.
        combination_mode : str, optional
            Combination strategy for skeleton slots.
            'independent' (default): each group generates combinations
            independently, results are merged.
            'global': all groups merge into one pool of position
            options, then a single Cartesian product is taken.
        """
        logger.debug("Initialized iterate job settings.")

        combination_mode = str(combination_mode).lower()
        if combination_mode not in _VALID_COMBINATION_MODES:
            raise ValueError(
                f"Invalid combination_mode '{combination_mode}'. "
                f"Expected one of {list(_VALID_COMBINATION_MODES)}."
            )

        # Path to the source configuration file; kept for diagnostics /
        # provenance only. It does not affect the algorithm or the output.
        self.config_file = config_file
        self.skeleton_list: list[dict] = []
        self.substituent_list: list[dict] = []
        self.combination_mode = combination_mode

        # Funnel every construction path (default, YAML, CLI, direct Python)
        # through resolve_algorithm_config so the stored config always has a
        # canonical algorithm name, validated options and complete defaults,
        # and cannot carry an invalid or partial state into the run.
        if algorithm_config is None:
            self.algorithm_config = resolve_algorithm_config()
        else:
            self.algorithm_config = resolve_algorithm_config(
                cli_algorithm_name=algorithm_config.name,
                cli_options=(algorithm_config.options or None),
            )

    def copy(self):
        """
        Create a deep copy of the current settings.
        """
        return copy.deepcopy(self)
