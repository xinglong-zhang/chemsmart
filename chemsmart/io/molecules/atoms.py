from ase import Atoms

class AtomsChargeMultiplicity(Atoms):
    """Simple ASE Atoms subclass with charge and spin multiplicity."""

    def __init__(self, charge, multiplicity, **kwargs):
        super().__init__(**kwargs)
        self.charge = charge
        self.spin_multiplicity = multiplicity
