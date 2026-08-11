"""PySCF CLI package.

The leaf imports attach the sp/opt/hess/td subcommands to the ``pyscf`` group,
so they are load-bearing rather than re-exports.
"""

from .hess import hess
from .opt import opt
from .pyscf import pyscf
from .singlepoint import sp
from .td import td

__all__ = ["hess", "opt", "pyscf", "sp", "td"]
