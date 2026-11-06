"""
Run ORCA optimisers (Opt, OptTS/ScanTS, IRC, constrained optimisations
and scans) on an external level of theory: machine-learning potentials,
xTB, or any ASE calculator.

:class:`ORCAExternalCalculator` is set as the ``external`` attribute of
:class:`~chemsmart.jobs.orca.settings.ORCAJobSettings`. When present, the
route uses ORCA's ``ExtOpt`` keyword in place of the QM method and basis,
the input gets a ``%method ProgExt "<launcher>" end`` block, and the ORCA
job runner writes the launcher next to the ``.inp`` file. The launcher
hands each request to :func:`chemsmart.io.orca.external.main`.

ORCA builds every Hessian it needs (TS searches, ``NumFreq``, IRC
initial Hessians) numerically from external gradients, so there is no
Hessian option here.
"""

import os

from chemsmart.jobs.external import EXTERNAL_BACKENDS, ExternalCalculator

__all__ = ["EXTERNAL_BACKENDS", "ORCAExternalCalculator"]


class ORCAExternalCalculator(ExternalCalculator):
    """Backend configuration for ORCA's ``ExtOpt`` interface.

    Takes the arguments of
    :class:`~chemsmart.jobs.external.ExternalCalculator`.
    """

    host = "ORCA"
    runtime_module = "chemsmart.io.orca.external"
    script_prefix = "chemsmart_orca_external"

    def progext(self, directory):
        """``ProgExt`` value: absolute launcher path in ``directory``.

        ORCA runs ``ProgExt`` directly (it must be executable and takes
        no interpreter), so the launcher's shebang selects the Python.
        """
        return os.path.abspath(os.path.join(directory, self.script_filename))
