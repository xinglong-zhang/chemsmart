"""
ORCA command-line interface subcommands.

This module provides CLI subcommands for various ORCA quantum chemistry calculations,
including geometry optimizations, transition state searches, IRC calculations,
single point calculations, and constrained optimizations.
"""

from .inp import inp
from .irc import irc
from .modred import modred
from .opt import opt
from .orca import orca

# from .qmmm import qmmm
from .qrc import qrc
from .scan import scan
from .singlepoint import sp
from .ts import ts

__all__ = [
    "inp",
    "irc",
    "modred",
    "opt",
    "orca",
    "qrc",
    "scan",
    "sp",
    "ts",
    "qmmm",
]
__all__ = ["inp", "irc", "modred", "opt", "orca", "scan", "sp", "ts", "qmmm"]
