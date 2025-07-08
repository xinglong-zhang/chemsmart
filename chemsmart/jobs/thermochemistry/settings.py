import logging

logger = logging.getLogger(__name__)


class ThermochemistryJobSettings:
    """Settings for ThermochemistryJob."""

    def __init__(
        self,
        temperature=None,
        concentration=None,
        pressure=1.0,
        use_weighted_mass=False,
        alpha=4,
        s_freq_cutoff=None,
        h_freq_cutoff=None,
        energy_units="hartree",
    ):
        self.temperature = temperature
        self.concentration = concentration
        self.pressure = pressure
        self.use_weighted_mass = use_weighted_mass
        self.alpha = alpha
        self.s_freq_cutoff = s_freq_cutoff
        self.h_freq_cutoff = h_freq_cutoff
        self.energy_units = energy_units.lower()

    def copy(self):
        return ThermochemistryJobSettings(
            temperature=self.temperature,
            concentration=self.concentration,
            pressure=self.pressure,
            use_weighted_mass=self.use_weighted_mass,
            alpha=self.alpha,
            s_freq_cutoff=self.s_freq_cutoff,
            h_freq_cutoff=self.h_freq_cutoff,
            energy_units=self.energy_units,
        )

    @classmethod
    def from_dict(cls, settings_dict):
        return cls(**settings_dict)
