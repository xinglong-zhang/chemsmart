from ase import Atoms


class AtomsChargeMultiplicity(Atoms):
    """Modified ASE Atoms subclass with charge and spin multiplicity.
    Main purpose for this is to be able to convert Atoms stored in
    ASE format to Molecule objects with charge and multiplicity."""

    def __init__(
        self, charge, multiplicity, frozen_atoms, energy, forces, **kwargs
    ):
        super().__init__(**kwargs)
        self._charge = charge
        self._multiplicity = multiplicity
        self.frozen_atoms = frozen_atoms
        self.energy = energy
        self.forces = forces

    @property
    def charge(self) -> int:
        return self._charge

    @charge.setter
    def charge(self, value):
        """Set the molecular charge.

        Args:
            value: The charge of the molecule (must be an integer).

        Raises:
            TypeError: If value is not an integer or cannot be converted to one.
            ValueError: If value is not a valid molecular charge (optional, depending on constraints).
        """
        if not isinstance(value, (int, float)) or (
            isinstance(value, float) and not value.is_integer()
        ):
            raise TypeError("Charge must be an integer")
        self._charge = int(value)

    @property
    def multiplicity(self) -> int:
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, value):
        """Set the molecular spin multiplicity.

        Args:
            value: The spin multiplicity of the molecule (must be an integer).

        Raises:
            TypeError: If value is not an integer or cannot be converted to one.
            ValueError: If value is not a valid spin multiplicity (optional, depending on constraints).
        """
        if not isinstance(value, (int, float)) or (
            isinstance(value, float) and not value.is_integer()
        ):
            raise TypeError("Spin multiplicity must be an integer")
        self._multiplicity = int(value)

    def to_molecule(self, charge=None, multiplicity=None):
        """Convert to Molecule object."""
        from chemsmart.io.molecules.structure import Molecule

        if charge is None:
            charge = self.charge
        if multiplicity is None:
            multiplicity = self.multiplicity

        return Molecule(
            symbols=self.symbols,
            positions=self.positions,
            charge=charge,
            multiplicity=multiplicity,
            frozen_atoms=self.frozen_atoms,
            pbc_conditions=self.pbc,
            translation_vectors=self.cell,
            energy=self.energy,
            forces=self.forces,
            velocities=self.get_velocities(),
            info=self.info,
        )

    @classmethod
    def from_atoms(cls, atoms, charge=None, multiplicity=None):
        """Create AtomsChargeMultiplicity from ASE Atoms."""
        from ase.calculators.calculator import (
            CalculatorError,
            PropertyNotImplementedError,
        )
        from ase.constraints import FixAtoms

        # check if the Atoms object has any constraints, if yes, convert to 1-indexed list
        frozen_atoms = []
        if atoms.constraints:
            # get indices of fix atoms --> frozen elements in Molecule
            if isinstance(atoms.constraints, list):
                for i, constraint in enumerate(atoms.constraints):
                    if isinstance(constraint, FixAtoms):
                        indices = FixAtoms.todict(constraint)["kwargs"][
                            "indices"
                        ]
                        # convert to 1-indexed list
                        for index in indices:
                            frozen_atoms.append(index + 1)
            elif isinstance(atoms.constraints, FixAtoms):
                indices = FixAtoms.todict(atoms.constraints)["kwargs"][
                    "indices"
                ]
                # convert to 1-indexed list
                for index in indices:
                    frozen_atoms.append(index + 1)
        if len(frozen_atoms) == 0:
            frozen_atoms = None
        else:
            # convert the 1-indexed list to masks where -1 means frozen and 0 means not frozen
            frozen_atoms = [
                -1 if i + 1 in frozen_atoms else 0 for i in range(len(atoms))
            ]

        # if no pbc, then should set translation vectors to None
        # cannot be set in ASE Atoms, as it defaults to zero Cell (still not None)

        # TODO: may need to convert velocities from ASE Atoms to Molecule
        # although this is not used in the current implementation

        # Check for valid calculator before calling get_potential_energy
        if atoms.calc is not None:
            try:
                energy = atoms.get_potential_energy()
                forces = atoms.get_forces()
                velocities = atoms.get_velocities()
            except (CalculatorError, PropertyNotImplementedError) as e:
                raise RuntimeError(
                    f"Failed to obtain energy or forces from {atoms.calc}. "
                    "Ensure that the Atoms object has a valid calculator."
                ) from e
        else:
            energy = None
            forces = None
            velocities = None

        return cls(
            charge=charge,
            multiplicity=multiplicity,
            frozen_atoms=frozen_atoms,
            energy=energy,
            forces=forces,
            symbols=atoms.get_chemical_symbols(),
            positions=atoms.get_positions(),
            pbc=atoms.get_pbc(),
            cell=atoms.get_cell(),  # translation vectors
            velocities=velocities,  # velocities
            info=atoms.info,
        )
