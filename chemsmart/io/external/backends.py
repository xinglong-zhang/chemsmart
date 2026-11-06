"""
Registry of external backends: the levels of theory chemsmart can hand
to a host program's optimiser.

A backend is a named ``build(options, calculator) -> ExternalEngine``
function. Built-in backends are registered below; others are added with
:func:`register_backend`, either in chemsmart or in a separate package
that declares an entry point in the ``chemsmart.external_backends``
group::

    # my_package/chemsmart_backends.py
    from chemsmart.io.external import ASEEngine, register_backend

    @register_backend("my-mlp", description="My ML potential")
    def build_my_mlp(options, calculator):
        def make(charge, multiplicity):
            from my_package import MyCalculator

            return MyCalculator(charge=charge, **options)

        return ASEEngine(make)

    # pyproject.toml of my_package
    [project.entry-points."chemsmart.external_backends"]
    my-mlp = "my_package.chemsmart_backends"

The backend is then available everywhere at once: ``-b my-mlp`` on
``chemsmart run gaussian/orca ... external`` and in the launchers.
"""

import importlib
import logging
from dataclasses import dataclass
from typing import Callable, Optional

from chemsmart.io.external.engines import ASEEngine, XTBEngine

logger = logging.getLogger(__name__)

ENTRY_POINT_GROUP = "chemsmart.external_backends"


@dataclass(frozen=True)
class ExternalBackend:
    """A registered backend.

    Attributes:
        name (str): Lower-case name used on the command line.
        build (callable): ``build(options, calculator)`` returning an
            :class:`~chemsmart.io.external.engines.ExternalEngine`. It runs
            inside the launcher, so heavy imports belong in the calculator
            factory it creates, not at module level.
        requires_calculator (bool): Whether ``--calculator`` is mandatory.
        validate (callable, optional): ``validate(options, calculator)``
            raising ``ValueError`` for bad input. Runs when the job is set
            up, so mistakes surface before the host program starts.
        description (str): One-line summary.
    """

    name: str
    build: Callable
    requires_calculator: bool = False
    validate: Optional[Callable] = None
    description: str = ""


_REGISTRY = {}
_plugins_loaded = False


def register_backend(
    name,
    *,
    requires_calculator=False,
    validate=None,
    description="",
    replace=False,
):
    """Decorator registering ``build(options, calculator)`` as a backend.

    Raises ``ValueError`` if ``name`` is taken, unless ``replace=True``.
    """
    key = str(name).lower()

    def decorator(build):
        if key in _REGISTRY and not replace:
            raise ValueError(f"External backend {key!r} already registered.")
        _REGISTRY[key] = ExternalBackend(
            name=key,
            build=build,
            requires_calculator=requires_calculator,
            validate=validate,
            description=description,
        )
        return build

    return decorator


def load_backend_plugins():
    """Import modules declared in the ``chemsmart.external_backends``
    entry-point group (once); importing them registers their backends."""
    global _plugins_loaded
    if _plugins_loaded:
        return
    _plugins_loaded = True
    from importlib.metadata import entry_points

    for entry_point in entry_points(group=ENTRY_POINT_GROUP):
        try:
            entry_point.load()
        except Exception as exc:
            logger.warning(
                f"Could not load external backend plugin "
                f"{entry_point.name!r} ({entry_point.value}): {exc}"
            )


def available_backends():
    """Names of all registered backends, plugins included."""
    load_backend_plugins()
    return tuple(_REGISTRY)


def get_backend(name):
    """Return the :class:`ExternalBackend` called ``name``."""
    load_backend_plugins()
    key = str(name).lower()
    if key not in _REGISTRY:
        raise ValueError(
            f"Unknown External backend {name!r}; choose from "
            f"{', '.join(_REGISTRY)}."
        )
    return _REGISTRY[key]


def validate_backend(name, options=None, calculator=None):
    """Check a backend configuration without building the engine.

    Returns the :class:`ExternalBackend`; raises ``ValueError``.
    """
    backend = get_backend(name)
    if backend.requires_calculator and not calculator:
        raise ValueError(
            f"backend {backend.name!r} needs a calculator dotted path, e.g. "
            "'ase.calculators.emt.EMT'."
        )
    if backend.validate is not None:
        backend.validate(dict(options or {}), calculator)
    return backend


def build_engine(config):
    """Create the engine described by a launcher configuration dict."""
    backend = validate_backend(
        config["backend"], config.get("options"), config.get("calculator")
    )
    return backend.build(
        dict(config.get("options") or {}), config.get("calculator")
    )


def import_object(dotted_path):
    """Import ``package.module.attribute`` and return the attribute."""
    module_name, _, attr = str(dotted_path).rpartition(".")
    if not module_name:
        raise ValueError(
            f"{dotted_path!r} is not a dotted path like "
            "'ase.calculators.emt.EMT'."
        )
    return getattr(importlib.import_module(module_name), attr)


# ── Built-in backends ─────────────────────────────────────────────────────


@register_backend(
    "ase",
    requires_calculator=True,
    description="Any ASE calculator, given by --calculator as a dotted path",
)
def _build_ase(options, calculator):
    target = import_object(calculator)

    def make(charge, multiplicity):
        return target(**options)

    return ASEEngine(make)


def _register_mace(name, function_name, description):
    @register_backend(name, description=description)
    def build(options, calculator):
        kwargs = {
            "model": "medium",
            "device": "cpu",
            "default_dtype": "float64",
            **options,
        }

        def make(charge, multiplicity):
            from mace import calculators

            return getattr(calculators, function_name)(**kwargs)

        return ASEEngine(make)


_register_mace("mace-off", "mace_off", "MACE-OFF organic force field")
_register_mace("mace-mp", "mace_mp", "MACE-MP foundation model")


@register_backend("aimnet2", description="AIMNet2 neural network potential")
def _build_aimnet2(options, calculator):
    options = dict(options)
    model = options.pop("model", "aimnet2")

    def make(charge, multiplicity):
        from aimnet.calculators import AIMNet2ASE

        return AIMNet2ASE(model, charge=charge, mult=multiplicity, **options)

    return ASEEngine(make)


@register_backend("uma", description="FAIRChem UMA universal model")
def _build_uma(options, calculator):
    options = dict(options)
    model = options.pop("model", "uma-s-1")
    device = options.pop("device", "cpu")
    task_name = options.pop("task_name", "omol")

    def make(charge, multiplicity):
        from fairchem.core import FAIRChemCalculator, pretrained_mlip

        predictor = pretrained_mlip.get_predict_unit(model, device=device)
        return FAIRChemCalculator(predictor, task_name=task_name, **options)

    return ASEEngine(make)


def _validate_xtb(options, calculator):
    XTBEngine(**options)  # raises ValueError for bad methods/solvents


@register_backend(
    "xtb",
    validate=_validate_xtb,
    description="GFN-xTB via the xtb executable",
)
def _build_xtb(options, calculator):
    return XTBEngine(**options)
