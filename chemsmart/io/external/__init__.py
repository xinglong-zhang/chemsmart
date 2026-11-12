"""
Program-neutral runtime for external levels of theory.

Host programs with an external-gradient interface (Gaussian
``External=``, ORCA ``ExtOpt``) keep their own optimisers and ask a
chemsmart launcher for energies, gradients and, for Gaussian, Hessians.
This package holds everything that does not depend on the host:

- :mod:`~chemsmart.io.external.engines`: :class:`ExternalEngine` and the
  ASE- and xtb-based engines (Angstrom in, atomic units out);
- :mod:`~chemsmart.io.external.backends`: the backend registry
  (:func:`register_backend`, :func:`build_engine`) and built-in backends;
- :mod:`~chemsmart.io.external.hessian`: analytic/finite-difference
  Hessians;
- :mod:`~chemsmart.io.external.protocol`: :class:`ExternalProtocol`, the
  base each host's file exchange implements, and the shared driver;
- :mod:`~chemsmart.io.external.runlog`: the per-job log
  ``<label>_external.log`` (one line per request: step, derivatives,
  energy, max/RMS gradient, Hessian source, wall time; plus backend
  warnings and, optionally, backend program output).

Adding a backend
    Register a ``build(options, calculator)`` function with
    :func:`register_backend`, in :mod:`~chemsmart.io.external.backends`
    or in another package through the ``chemsmart.external_backends``
    entry point (see that module). Nothing else changes: the CLI,
    settings and launchers of every host accept it.

Adding a host program
    1. ``chemsmart/io/<host>/external.py``: an :class:`ExternalProtocol`
       subclass for the host's files, and ``main(config, argv=None)``
       calling its :meth:`~ExternalProtocol.main`.
    2. ``chemsmart/jobs/<host>/external.py``: a
       :class:`~chemsmart.jobs.external.ExternalCalculator` subclass
       setting ``host``, ``runtime_module`` (step 1) and
       ``script_prefix``, plus how the host's input refers to the
       launcher.
    3. Host hooks: an ``external`` attribute on the host's job settings,
       the route/input keyword in its writer, and
       ``external.write_script(running_directory, server=...)`` in its
       job runner's ``_write_input``. The per-job log needs nothing:
       :class:`~chemsmart.jobs.runner.JobRunner` sets
       :data:`EXTERNAL_LOG_ENV` for every external job, and the shared
       driver writes the records.
    4. ``chemsmart/cli/<host>/external.py``: a group using
       :func:`chemsmart.cli.external.click_external_backend_options`,
       :func:`~chemsmart.cli.external.build_external_calculator` and
       :func:`~chemsmart.cli.external.use_external_level_of_theory`,
       wrapping the host's job subcommands.
    Tests mirror ``tests/test_gaussian_external.py`` and
    ``tests/test_orca_external.py``.
"""

from chemsmart.io.external.backends import (
    ENTRY_POINT_GROUP,
    ExternalBackend,
    available_backends,
    build_engine,
    get_backend,
    import_object,
    load_backend_plugins,
    register_backend,
    validate_backend,
)
from chemsmart.io.external.engines import (
    ASEEngine,
    ExternalEngine,
    XTBEngine,
    read_xtb_gradient,
    read_xtb_hessian,
)
from chemsmart.io.external.hessian import (
    DEFAULT_FD_STEP,
    HESSIAN_MODES,
    compute_hessian,
    finite_difference_hessian,
)
from chemsmart.io.external.protocol import (
    ExternalProtocol,
    ExternalRequest,
    ExternalResult,
    evaluate,
)
from chemsmart.io.external.runlog import EXTERNAL_LOG_ENV, ExternalRunLog

__all__ = [
    "ASEEngine",
    "DEFAULT_FD_STEP",
    "ENTRY_POINT_GROUP",
    "EXTERNAL_LOG_ENV",
    "ExternalBackend",
    "ExternalEngine",
    "ExternalProtocol",
    "ExternalRequest",
    "ExternalResult",
    "ExternalRunLog",
    "HESSIAN_MODES",
    "XTBEngine",
    "available_backends",
    "build_engine",
    "compute_hessian",
    "evaluate",
    "finite_difference_hessian",
    "get_backend",
    "import_object",
    "load_backend_plugins",
    "read_xtb_gradient",
    "read_xtb_hessian",
    "register_backend",
    "validate_backend",
]
