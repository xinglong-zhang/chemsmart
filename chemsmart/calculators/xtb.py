"""
XTB calculator interface for external program coupling.

This module provides the XTBCalculator class for interfacing the
semi-empirical tight-binding program xtb with external optimization
engines like Gaussian and ORCA. The calculator handles command-line
argument generation and configuration for xtb calculations.

The implementation is based on the xtb-gaussian interface from:
https://github.com/aspuru-guzik-group/xtb-gaussian

References:
    - xtb: https://github.com/grimme-lab/xtb
    - gau_xtb: http://sobereva.com/soft/gau_xtb/
"""

import logging
import shutil
from typing import Optional

logger = logging.getLogger(__name__)


class XTBCalculator:
    """
    XTB calculator configuration for external program interfaces.

    Manages the configuration and command-line generation for xtb
    calculations when used as an external calculator with optimization
    engines like Gaussian's TS optimizer. Handles parallelization,
    solvation, and other xtb-specific options.

    Attributes:
        method (str): XTB method level (e.g., 'gfn2', 'gfn1', 'gfn0').
        num_threads (int): Number of parallel threads for xtb.
        electronic_temperature (float): Electronic temperature in Kelvin.
        solvation (str): Implicit solvation model (e.g., 'gbsa', 'alpb').
        solvent (str): Solvent name for solvation calculations.
        accuracy (float): Numerical accuracy setting.
        max_iterations (int): Maximum SCF iterations.
        uhf (int): Number of unpaired electrons for spin state.
        xtb_executable (str): Path or name of xtb executable.
        wrapper_script (str): Path to the xtb-gaussian wrapper script.
        additional_args (str): Additional command-line arguments.

    Example:
        >>> calc = XTBCalculator(num_threads=4, solvent="water")
        >>> calc.command_string
        '-P 4 --alpb water'
    """

    # Supported solvation models
    SOLVATION_MODELS = ["gbsa", "alpb"]

    # Supported method levels
    XTB_METHODS = ["gfn2", "gfn1", "gfn0", "gfnff"]

    def __init__(
        self,
        method: str = "gfn2",
        num_threads: int = 1,
        electronic_temperature: Optional[float] = None,
        solvation: Optional[str] = None,
        solvent: Optional[str] = None,
        accuracy: Optional[float] = None,
        max_iterations: Optional[int] = None,
        uhf: Optional[int] = None,
        xtb_executable: str = "xtb",
        wrapper_script: str = "xtb-gaussian",
        additional_args: Optional[str] = None,
        log_all: bool = False,
    ):
        """
        Initialize XTB calculator with specified configuration.

        Args:
            method: XTB method level. Defaults to 'gfn2'.
            num_threads: Number of parallel threads. Defaults to 1.
            electronic_temperature: Electronic temperature in Kelvin.
            solvation: Solvation model ('gbsa' or 'alpb').
            solvent: Solvent name for implicit solvation.
            accuracy: Numerical accuracy (smaller = more accurate).
            max_iterations: Maximum SCF iterations.
            uhf: Number of unpaired electrons for spin state.
            xtb_executable: Path or name of xtb executable.
            wrapper_script: Path to xtb-gaussian wrapper script.
            additional_args: Additional xtb command-line arguments.
            log_all: If True, include full xtb output in Gaussian log.

        Raises:
            ValueError: If method or solvation model is not supported.
        """
        if method.lower() not in self.XTB_METHODS:
            raise ValueError(
                f"XTB method '{method}' not supported. "
                f"Supported methods: {self.XTB_METHODS}"
            )

        if solvation is not None and solvation.lower() not in self.SOLVATION_MODELS:
            raise ValueError(
                f"Solvation model '{solvation}' not supported. "
                f"Supported models: {self.SOLVATION_MODELS}"
            )

        self.method = method.lower()
        self.num_threads = num_threads
        self.electronic_temperature = electronic_temperature
        self.solvation = solvation.lower() if solvation else None
        self.solvent = solvent
        self.accuracy = accuracy
        self.max_iterations = max_iterations
        self.uhf = uhf
        self.xtb_executable = xtb_executable
        self.wrapper_script = wrapper_script
        self.additional_args = additional_args
        self.log_all = log_all

        logger.debug(
            f"XTBCalculator initialized with method={method}, "
            f"threads={num_threads}, solvent={solvent}"
        )

    @property
    def command_args(self) -> list:
        """
        Generate list of command-line arguments for xtb.

        Constructs the argument list based on calculator configuration,
        excluding the executable name itself.

        Returns:
            list: Command-line arguments for xtb calculation.
        """
        args = []

        # Log all flag must come first for xtb-gaussian wrapper
        if self.log_all:
            args.append("--log-all")

        # Parallelization
        if self.num_threads > 1:
            args.extend(["-P", str(self.num_threads)])

        # Method level (only add if not default gfn2)
        if self.method != "gfn2":
            args.append(f"--{self.method}")

        # Electronic temperature
        if self.electronic_temperature is not None:
            args.extend(["--etemp", str(self.electronic_temperature)])

        # Solvation
        if self.solvent is not None:
            solvation_model = self.solvation or "alpb"
            args.append(f"--{solvation_model}")
            args.append(self.solvent)

        # Accuracy
        if self.accuracy is not None:
            args.extend(["--acc", str(self.accuracy)])

        # Max iterations
        if self.max_iterations is not None:
            args.extend(["--iterations", str(self.max_iterations)])

        # Spin state (unpaired electrons)
        if self.uhf is not None:
            args.extend(["--uhf", str(self.uhf)])

        # Additional user arguments
        if self.additional_args is not None:
            args.extend(self.additional_args.split())

        return args

    @property
    def command_string(self) -> str:
        """
        Generate complete command string for xtb arguments.

        Returns:
            str: Space-separated string of xtb arguments.
        """
        return " ".join(self.command_args)

    @property
    def external_command(self) -> str:
        """
        Generate the complete external command for Gaussian.

        Creates the command string suitable for use in Gaussian's
        external keyword, including the wrapper script path.

        Returns:
            str: Complete command string for Gaussian external keyword.
        """
        if self.command_string:
            return f'{self.wrapper_script} {self.command_string}'
        return self.wrapper_script

    def is_available(self) -> bool:
        """
        Check if xtb and wrapper script are available on the system.

        Returns:
            bool: True if both xtb and wrapper are found in PATH.
        """
        xtb_available = shutil.which(self.xtb_executable) is not None
        wrapper_available = shutil.which(self.wrapper_script) is not None

        if not xtb_available:
            logger.warning(f"xtb executable '{self.xtb_executable}' not found in PATH")
        if not wrapper_available:
            logger.warning(
                f"Wrapper script '{self.wrapper_script}' not found in PATH"
            )

        return xtb_available and wrapper_available

    def copy(self) -> "XTBCalculator":
        """
        Create a deep copy of the calculator.

        Returns:
            XTBCalculator: Independent copy of this calculator.
        """
        return XTBCalculator(
            method=self.method,
            num_threads=self.num_threads,
            electronic_temperature=self.electronic_temperature,
            solvation=self.solvation,
            solvent=self.solvent,
            accuracy=self.accuracy,
            max_iterations=self.max_iterations,
            uhf=self.uhf,
            xtb_executable=self.xtb_executable,
            wrapper_script=self.wrapper_script,
            additional_args=self.additional_args,
            log_all=self.log_all,
        )

    def __repr__(self) -> str:
        return (
            f"XTBCalculator(method='{self.method}', "
            f"num_threads={self.num_threads}, "
            f"solvent={self.solvent!r})"
        )

    @classmethod
    def default(cls) -> "XTBCalculator":
        """
        Create a default XTB calculator with standard settings.

        Returns:
            XTBCalculator: Calculator with default configuration.
        """
        return cls()

    @classmethod
    def from_dict(cls, config: dict) -> "XTBCalculator":
        """
        Create XTB calculator from a configuration dictionary.

        Args:
            config: Dictionary with calculator parameters.

        Returns:
            XTBCalculator: Configured calculator instance.
        """
        return cls(**config)

    def to_dict(self) -> dict:
        """
        Convert calculator configuration to dictionary.

        Returns:
            dict: Calculator configuration parameters.
        """
        return {
            "method": self.method,
            "num_threads": self.num_threads,
            "electronic_temperature": self.electronic_temperature,
            "solvation": self.solvation,
            "solvent": self.solvent,
            "accuracy": self.accuracy,
            "max_iterations": self.max_iterations,
            "uhf": self.uhf,
            "xtb_executable": self.xtb_executable,
            "wrapper_script": self.wrapper_script,
            "additional_args": self.additional_args,
            "log_all": self.log_all,
        }
