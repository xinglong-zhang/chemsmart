"""
Settings configuration for thermochemistry jobs.

This module provides the ThermochemistryJobSettings class for configuring
thermochemical calculations, including temperature, pressure, frequency
cutoffs, and various correction methods for entropy and enthalpy.
"""

import logging

logger = logging.getLogger(__name__)


class ThermochemistryJobSettings:
    """
    Configuration settings for thermochemistry calculations.

    This class manages all parameters needed for thermochemical property
    calculations, including temperature conditions, frequency corrections,
    and output formatting options.

    Attributes:
        temperature (float | None): Temperature in Kelvin.
        concentration (float | None): Concentration for solution-phase corrections.
        pressure (float): Pressure in atm (default 1.0).
        use_weighted_mass (bool): Use mass-weighted scheme when True.
        alpha (int): Scaling factor for QRRHO corrections.
        s_freq_cutoff (float | None): Entropy frequency cutoff (cm^-1).
        entropy_method (str | None): Entropy method (e.g., 'qrrho').
        h_freq_cutoff (float | None): Enthalpy frequency cutoff (cm^-1).
        energy_units (str): Energy unit label (e.g., 'hartree', 'kcal/mol').
        outputfile (str | None): Path to write thermochemistry results.
        overwrite (bool): Overwrite existing output when True.
        check_imaginary_frequencies (bool): Validate absence of imaginary modes.
    """

    def __init__(
        self,
        temperature=None,
        concentration=None,
        pressure=1.0,
        use_weighted_mass=False,
        alpha=4,
        s_freq_cutoff=None,
        entropy_method=None,
        h_freq_cutoff=None,
        energy_units="hartree",
        outputfile=None,
        overwrite=False,
        check_imaginary_frequencies=True,
    ):
        """
        Initialize thermochemistry job settings.

        Args:
            temperature (float, optional): Temperature for calculations (K)
            concentration (float, optional): Concentration for calculations
            pressure (float): Pressure for calculations (atm, default: 1.0)
            use_weighted_mass (bool): Whether to use weighted atomic masses
            alpha (int): Scaling factor for QRRHO corrections (default: 4)
            s_freq_cutoff (float, optional): Frequency cutoff for entropy
            entropy_method (str, optional): Method for entropy calculation
            h_freq_cutoff (float, optional): Frequency cutoff for enthalpy
            energy_units (str): Units for energy values (default: 'hartree')
            outputfile (str, optional): Path to output file
            overwrite (bool): Whether to overwrite existing output files
            check_imaginary_frequencies (bool): Check for imaginary frequencies
        """
        logger.debug("Initializing ThermochemistryJobSettings")
        self.temperature = temperature
        self.concentration = concentration
        self.pressure = pressure
        self.use_weighted_mass = use_weighted_mass
        self.alpha = alpha
        self.s_freq_cutoff = s_freq_cutoff
        self.entropy_method = entropy_method
        self.h_freq_cutoff = h_freq_cutoff
        self.energy_units = energy_units.lower()
        self.outputfile = outputfile
        self.overwrite = overwrite
        self.check_imaginary_frequencies = check_imaginary_frequencies

    def copy(self):
        """
        Create a deep copy of the settings instance.

        Returns:
            ThermochemistryJobSettings: New instance with identical settings
        """
        return ThermochemistryJobSettings(
            temperature=self.temperature,
            concentration=self.concentration,
            pressure=self.pressure,
            use_weighted_mass=self.use_weighted_mass,
            alpha=self.alpha,
            s_freq_cutoff=self.s_freq_cutoff,
            entropy_method=self.entropy_method,
            h_freq_cutoff=self.h_freq_cutoff,
            energy_units=self.energy_units,
            outputfile=self.outputfile,
            overwrite=self.overwrite,
            check_imaginary_frequencies=self.check_imaginary_frequencies,
        )

    @classmethod
    def from_dict(cls, settings_dict):
        """
        Create settings instance from a dictionary.

        Factory method that creates a ThermochemistryJobSettings instance
        from a dictionary containing configuration parameters.

        Args:
            settings_dict (dict): Dictionary containing settings parameters

        Returns:
            ThermochemistryJobSettings: Configured settings instance
        """
        return cls(**settings_dict)
