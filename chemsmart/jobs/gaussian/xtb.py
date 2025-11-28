"""
Gaussian job with xtb external calculator for TS optimization.

This module provides job classes for using Gaussian's optimization
algorithms with xtb as the energy/gradient calculator. This is
particularly useful for transition state optimization where xtb
provides fast semi-empirical energies while Gaussian provides
sophisticated TS optimization methods like QST2 and Berny.

Based on the xtb-gaussian interface:
https://github.com/aspuru-guzik-group/xtb-gaussian
"""

import logging

from chemsmart.calculators.xtb import XTBCalculator
from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings

logger = logging.getLogger(__name__)


class GaussianXTBJobSettings(GaussianJobSettings):
    """
    Job settings for Gaussian calculations using xtb as external calculator.

    Extends GaussianJobSettings with configuration for the xtb calculator
    which is called via Gaussian's external keyword. This enables using
    Gaussian's optimization algorithms (especially TS optimization) with
    xtb's fast semi-empirical energy and gradient calculations.

    Attributes:
        calculator (XTBCalculator): XTB calculator configuration.
        All other attributes inherited from GaussianJobSettings.

    Note:
        When using external calculators, the functional and basis set
        are not used for energy/gradient calculations - they are
        replaced by the external calculator. However, they may still
        appear in the input file for compatibility.
    """

    def __init__(
        self,
        calculator=None,
        **kwargs,
    ):
        """
        Initialize Gaussian XTB job settings.

        Args:
            calculator (XTBCalculator, optional): XTB calculator config.
                If None, a default XTBCalculator is created.
            **kwargs: Additional arguments for GaussianJobSettings.
        """
        # Set defaults for external calculator jobs if not provided
        if "functional" not in kwargs or kwargs.get("functional") is None:
            kwargs["functional"] = "HF"
        if "basis" not in kwargs or kwargs.get("basis") is None:
            kwargs["basis"] = "STO-3G"

        super().__init__(**kwargs)

        if calculator is None:
            calculator = XTBCalculator.default()
        elif isinstance(calculator, dict):
            calculator = XTBCalculator.from_dict(calculator)

        self.calculator = calculator
        logger.debug(f"GaussianXTBJobSettings initialized with {calculator}")

    @property
    def external_command(self) -> str:
        """
        Get the external command string for Gaussian input.

        Returns:
            str: Command string for Gaussian's external keyword.
        """
        return self.calculator.external_command

    @property
    def route_string(self):
        """
        Generate the Gaussian route string with external keyword.

        Constructs the route line for Gaussian input, replacing the
        standard method/basis with the external calculator command.

        Returns:
            str: Complete route string with external keyword.
        """
        route_string = self._get_route_string_for_external()
        logger.debug(f"Route for XTB external settings: {route_string}")
        return route_string

    def _get_route_string_for_external(self):
        """
        Generate route string for external calculator jobs.

        Creates the route string using Gaussian's external keyword
        syntax. The external command invokes xtb for energy and
        gradient calculations.

        Returns:
            str: Route string with external calculator specification.
        """
        route_string = ""

        # Add #-tag prefix
        if self.dieze_tag is not None:
            route_string += f"#{self.dieze_tag}"
        else:
            route_string += "#"

        # Add external keyword with calculator command
        route_string += f' external="{self.external_command}"'

        # Add optimization keywords based on job type
        if self.additional_opt_options_in_route is not None:
            if self.job_type == "opt":
                route_string += (
                    f" opt=({self.additional_opt_options_in_route})"
                )
            elif self.job_type == "ts":
                if "calcall" in self.additional_opt_options_in_route:
                    route_string += (
                        f" opt=(ts,noeigentest,"
                        f"{self.additional_opt_options_in_route})"
                    )
                else:
                    route_string += (
                        f" opt=(ts,calcfc,noeigentest,"
                        f"{self.additional_opt_options_in_route})"
                    )
            elif self.job_type == "qst2":
                route_string += (
                    f" opt=(qst2,{self.additional_opt_options_in_route})"
                )
            elif self.job_type == "qst3":
                route_string += (
                    f" opt=(qst3,{self.additional_opt_options_in_route})"
                )
        else:
            if self.job_type == "opt":
                route_string += " opt"
            elif self.job_type == "ts":
                route_string += " opt=(ts,calcfc,noeigentest)"
            elif self.job_type == "qst2":
                route_string += " opt=qst2"
            elif self.job_type == "qst3":
                route_string += " opt=qst3"

        # Add frequency calculations if enabled
        if self.freq and not self.numfreq:
            route_string += " freq"
        elif not self.freq and self.numfreq:
            route_string += " freq=numer"

        # Add additional route parameters
        if self.additional_route_parameters is not None:
            route_string += f" {self.additional_route_parameters}"

        return route_string

    def copy(self):
        """
        Create a deep copy of the settings object.

        Returns:
            GaussianXTBJobSettings: Deep copy of settings.
        """
        import copy
        return copy.deepcopy(self)

    @classmethod
    def default(cls):
        """
        Create default Gaussian XTB job settings.

        Returns:
            GaussianXTBJobSettings: Default settings for XTB jobs.
        """
        return cls(
            calculator=XTBCalculator.default(),
            charge=None,
            multiplicity=None,
            chk=True,
            job_type="opt",
            title="Gaussian job with XTB external calculator",
            freq=False,
            numfreq=False,
        )

    @classmethod
    def for_ts_optimization(
        cls,
        num_threads=4,
        solvent=None,
        electronic_temperature=None,
        additional_opt_options=None,
        **kwargs,
    ):
        """
        Create settings for transition state optimization with xtb.

        Factory method for creating settings optimized for TS searches
        using Gaussian's Berny optimizer with xtb gradients.

        Args:
            num_threads: Number of xtb parallel threads.
            solvent: Solvent name for implicit solvation.
            electronic_temperature: Electronic temperature in Kelvin.
            additional_opt_options: Extra optimization keywords.
            **kwargs: Additional settings parameters.

        Returns:
            GaussianXTBJobSettings: Settings for TS optimization.
        """
        calculator = XTBCalculator(
            num_threads=num_threads,
            solvent=solvent,
            electronic_temperature=electronic_temperature,
        )

        # Default TS optimization options
        opt_options = additional_opt_options or "nomicro"

        return cls(
            calculator=calculator,
            job_type="ts",
            freq=True,
            additional_opt_options_in_route=opt_options,
            title="Transition state optimization with XTB",
            **kwargs,
        )

    @classmethod
    def for_qst2(
        cls,
        num_threads=4,
        solvent=None,
        additional_opt_options=None,
        **kwargs,
    ):
        """
        Create settings for QST2 reaction path search with xtb.

        Factory method for QST2 (Quasi-Newton Synchronous Transit)
        calculations using xtb for energy and gradient evaluation.

        Args:
            num_threads: Number of xtb parallel threads.
            solvent: Solvent name for implicit solvation.
            additional_opt_options: Extra optimization keywords.
            **kwargs: Additional settings parameters.

        Returns:
            GaussianXTBJobSettings: Settings for QST2 calculation.
        """
        calculator = XTBCalculator(
            num_threads=num_threads,
            solvent=solvent,
        )

        opt_options = additional_opt_options or "nomicro,bimolecular,recalcfc=5"

        return cls(
            calculator=calculator,
            job_type="qst2",
            freq=False,
            additional_opt_options_in_route=opt_options,
            title="QST2 reaction path search with XTB",
            **kwargs,
        )


class GaussianXTBJob(GaussianJob):
    """
    Gaussian job using xtb as external energy/gradient calculator.

    This job class enables using Gaussian's optimization algorithms
    (particularly TS optimization, QST2, and QST3) with xtb's fast
    semi-empirical calculations. The xtb-gaussian wrapper script
    handles the communication between Gaussian and xtb.

    Attributes:
        TYPE (str): Job type identifier ('g16xtb').
        molecule (Molecule): Molecular structure for calculation.
        settings (GaussianXTBJobSettings): Job configuration.
        label (str): Job identifier for file naming.

    Example:
        >>> from chemsmart.calculators.xtb import XTBCalculator
        >>> settings = GaussianXTBJobSettings.for_ts_optimization(
        ...     num_threads=4, solvent="methanol"
        ... )
        >>> job = GaussianXTBJob(molecule=mol, settings=settings, label="ts")
        >>> job.run()
    """

    TYPE = "g16xtb"

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Initialize a Gaussian XTB job.

        Args:
            molecule (Molecule): Molecular structure for calculation.
            settings (GaussianXTBJobSettings): Job configuration.
            label (str, optional): Job identifier.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional arguments for parent class.
        """
        if settings is None:
            settings = GaussianXTBJobSettings.default()

        # Ensure settings is the correct type
        if not isinstance(settings, GaussianXTBJobSettings):
            if isinstance(settings, GaussianJobSettings):
                # Convert regular settings to XTB settings
                logger.warning(
                    "Converting GaussianJobSettings to GaussianXTBJobSettings"
                )
                settings = GaussianXTBJobSettings(
                    calculator=XTBCalculator.default(),
                    **settings.__dict__,
                )
            else:
                raise ValueError(
                    f"Settings must be GaussianXTBJobSettings, "
                    f"got {type(settings)}"
                )

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @property
    def calculator(self):
        """
        Get the XTB calculator from settings.

        Returns:
            XTBCalculator: The configured XTB calculator.
        """
        return self.settings.calculator


class GaussianXTBTSJob(GaussianXTBJob):
    """
    Specialized job for TS optimization using xtb with Gaussian.

    Convenience class for transition state optimization that
    pre-configures appropriate settings for TS searches.

    Attributes:
        TYPE (str): Job type identifier ('g16xtbts').
    """

    TYPE = "g16xtbts"

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Initialize a Gaussian XTB TS optimization job.

        Args:
            molecule (Molecule): TS guess structure.
            settings (GaussianXTBJobSettings, optional): Job config.
            label (str, optional): Job identifier.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional arguments.
        """
        if settings is None:
            settings = GaussianXTBJobSettings.for_ts_optimization(**kwargs)

        # Ensure job_type is ts
        settings.job_type = "ts"

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )


class GaussianXTBQST2Job(GaussianXTBJob):
    """
    Job for QST2 reaction path search using xtb with Gaussian.

    Performs QST2 (Quasi-Newton Synchronous Transit) calculations
    to find transition states from reactant and product structures.
    Uses xtb for fast energy and gradient evaluation.

    Attributes:
        TYPE (str): Job type identifier ('g16xtbqst2').
        reactant (Molecule): Reactant structure.
        product (Molecule): Product structure.
    """

    TYPE = "g16xtbqst2"

    def __init__(
        self,
        reactant,
        product,
        settings=None,
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        """
        Initialize a Gaussian XTB QST2 job.

        Args:
            reactant (Molecule): Reactant molecular structure.
            product (Molecule): Product molecular structure.
            settings (GaussianXTBJobSettings, optional): Job config.
            label (str, optional): Job identifier.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional arguments.
        """
        if settings is None:
            settings = GaussianXTBJobSettings.for_qst2(**kwargs)

        # Store both structures
        self.reactant = reactant.copy()
        self.product = product.copy()

        # Use reactant as primary molecule
        super().__init__(
            molecule=reactant,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
