from ase import Atoms


class AtomsChargeMultiplicity(Atoms):
    """Simple ASE Atoms subclass with charge and spin multiplicity."""

    # see pymatgen.io.ase.py ASEAtomsAdaptor.get_molecule()

    def __init__(self, charge, multiplicity, **kwargs):
        super().__init__(**kwargs)
        self.charge = charge
        self.spin_multiplicity = multiplicity

    def to_molecule(self):
        """Convert to Molecule object."""
        from chemsmart.io.molecules.structure import Molecule

        return Molecule(
            symbols=self.symbols,
            positions=self.positions,
            charge=self.charge,
            multiplicity=self.spin_multiplicity,
            frozen_atoms=None,
            pbc_conditions=self.pbc,
            translation_vectors=self.cell,
            info=self.info,
        )

    @classmethod
    def from_atoms(cls, atoms, charge=None, multiplicity=None):
        """Create AtomsChargeMultiplicity from ASE Atoms."""
        return cls(
            charge=charge,
            multiplicity=multiplicity,
            symbols=atoms.get_chemical_symbols(),
            positions=atoms.get_positions(),
            cell=atoms.get_cell(),
            pbc=atoms.get_pbc(),
            info=atoms.info,
        )
