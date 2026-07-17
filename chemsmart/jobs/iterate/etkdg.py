"""RDKit ETKDG analyzer for Iterate structure generation."""

import logging
from typing import Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)


class IterateETKDGAnalyzer:
    """Attach one or more substituents to a skeleton using RDKit ETKDGv3.

    This analyzer attaches every substituent of a combination in a single
    pass. The skeleton and all substituents are joined into one RDKit
    molecule (skeleton atoms first, then each substituent in order); every
    skeleton-substituent bond is added explicitly, and ETKDGv3 then embeds
    the whole assembly exactly once. Substituents are never inserted
    sequentially through repeated ``Molecule`` round-trips, so the explicit
    connection bonds are preserved and the assembly never fragments.

    In the default ``local`` mode only the original processed skeleton atoms
    are held fixed through ``coordMap`` and pinned back to their exact
    reference coordinates; every substituent atom is free to move. In
    ``global`` mode every atom is re-embedded. ``num_conformers`` conformers
    are embedded and the lowest-energy one (by UFF energy) is kept.

    Supported ``options`` keys (see the ``etkdg`` entry in the algorithm
    registry): ``use_global_optimization``, ``num_conformers``,
    ``random_seed``, ``max_iterations``, ``use_random_coordinates``,
    ``enforce_chirality`` and ``force_field``.
    """

    #: Supported force fields for the optional post-embedding minimization.
    _SUPPORTED_FORCE_FIELDS = ("none", "uff", "mmff94", "mmff94s")
    #: Fixed settings for the optional force-field minimization.
    _FORCE_FIELD_MAX_ITERATIONS = 200
    _FORCE_FIELD_FORCE_TOL = 1e-4
    _FORCE_FIELD_ENERGY_TOL = 1e-6
    _NONBONDED_THRESHOLD = 100.0

    def __init__(
        self,
        skeleton: Molecule,
        substituents: list,
        options: Optional[dict] = None,
    ):
        """
        Parameters
        ----------
        skeleton : Molecule
            Preprocessed skeleton molecule. All old groups are removed in a
            single batch beforehand, so every substituent attaches to this
            same skeleton.
        substituents : list of tuple
            One ``(substituent, skeleton_link_index, substituent_link_index)``
            tuple per attachment, where ``substituent`` is a preprocessed
            :class:`Molecule` and both link indices are 1-based. A single
            attachment is simply a list holding one tuple.
        options : dict, optional
            ETKDG algorithm options (see the ``etkdg`` entry in the
            algorithm registry).
        """
        self.skeleton = skeleton
        # Store each attachment with 0-based link indices for RDKit.
        self.substituents = [
            (substituent, skeleton_link_index - 1, substituent_link_index - 1)
            for (
                substituent,
                skeleton_link_index,
                substituent_link_index,
            ) in substituents
        ]
        self.options = dict(options or {})

    @staticmethod
    def _prepare_molecule(molecule: Chem.Mol) -> Chem.Mol:
        """Return an embeddable copy of the combined molecule.

        Clears aromatic flags on atoms/bonds that are not in a ring,
        Kekulizes the remaining aromatic rings (falling back to single bonds
        on failure), sanitizes while skipping strict valence checks (so
        metals and unusual valences do not fail early), and removes any
        conformers. The input molecule is not modified.
        """
        mol = Chem.Mol(molecule)
        editable = Chem.RWMol(mol)
        for bond in editable.GetBonds():
            if (
                bond.GetBondType() == Chem.BondType.AROMATIC
                and not bond.IsInRing()
            ):
                bond.SetBondType(Chem.BondType.SINGLE)
                bond.SetIsAromatic(False)
        for atom in editable.GetAtoms():
            if atom.GetIsAromatic() and not atom.IsInRing():
                atom.SetIsAromatic(False)
        mol = editable.GetMol()
        try:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception:
            editable = Chem.RWMol(mol)
            for bond in editable.GetBonds():
                if bond.GetBondType() == Chem.BondType.AROMATIC:
                    bond.SetBondType(Chem.BondType.SINGLE)
                    bond.SetIsAromatic(False)
            for atom in editable.GetAtoms():
                atom.SetIsAromatic(False)
            mol = editable.GetMol()
        mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(
                mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES
            )
        except Exception:
            pass
        mol.RemoveAllConformers()
        return mol

    @staticmethod
    def _fixed_reference_positions(
        molecule: Chem.Mol, n_skeleton: int
    ) -> dict:
        """Read the skeleton atom coordinates (indices 0..n_skeleton-1)."""
        conformer = molecule.GetConformer()
        return {
            index: np.array(
                [
                    conformer.GetAtomPosition(index).x,
                    conformer.GetAtomPosition(index).y,
                    conformer.GetAtomPosition(index).z,
                ],
                dtype=float,
            )
            for index in range(n_skeleton)
        }

    def _embedding_parameters(self, reference_positions: dict):
        """Build ETKDGv3 parameters from the options and skeleton coordMap."""
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = int(self.options.get("random_seed", 42))
        parameters.maxIterations = int(
            self.options.get("max_iterations", 2000)
        )
        parameters.useRandomCoords = bool(
            self.options.get("use_random_coordinates", True)
        )
        parameters.enforceChirality = bool(
            self.options.get("enforce_chirality", False)
        )
        parameters.numThreads = 1
        if reference_positions:
            parameters.SetCoordMap(
                {
                    index: Point3D(
                        float(position[0]),
                        float(position[1]),
                        float(position[2]),
                    )
                    for index, position in reference_positions.items()
                }
            )
        return parameters

    @staticmethod
    def _set_fixed_coordinates(
        molecule: Chem.Mol, conformer_id: int, reference_positions: dict
    ) -> None:
        """Pin the fixed atoms back to their exact reference coordinates."""
        if not reference_positions:
            return
        conformer = molecule.GetConformer(conformer_id)
        for index, position in reference_positions.items():
            conformer.SetAtomPosition(
                index,
                Point3D(
                    float(position[0]),
                    float(position[1]),
                    float(position[2]),
                ),
            )

    def _force_field_optimize(
        self,
        molecule: Chem.Mol,
        conformer_id: int,
        fixed_indices: tuple,
        force_field: str,
    ) -> None:
        """Minimize one conformer with the requested force field.

        The fixed atoms are added as force-field fixed points so the
        skeleton does not move during minimization.
        """
        if force_field == "uff":
            field = AllChem.UFFGetMoleculeForceField(
                molecule,
                confId=conformer_id,
                vdwThresh=self._NONBONDED_THRESHOLD,
            )
        else:
            variant = "MMFF94s" if force_field == "mmff94s" else "MMFF94"
            properties = AllChem.MMFFGetMoleculeProperties(
                molecule, mmffVariant=variant
            )
            if properties is None:
                raise ValueError(f"{variant} parameters are unavailable")
            field = AllChem.MMFFGetMoleculeForceField(
                molecule,
                properties,
                confId=conformer_id,
                nonBondedThresh=self._NONBONDED_THRESHOLD,
            )
        if field is None:
            raise ValueError(f"could not construct {force_field} force field")
        for index in fixed_indices:
            field.AddFixedPoint(index)
        field.Initialize()
        field.Minimize(
            maxIts=self._FORCE_FIELD_MAX_ITERATIONS,
            forceTol=self._FORCE_FIELD_FORCE_TOL,
            energyTol=self._FORCE_FIELD_ENERGY_TOL,
        )

    def _rough_place_substituents(self) -> list:
        """Rigidly translate each substituent near its skeleton link atom.

        The exact placement is unimportant because ETKDG regenerates the
        substituent atoms; this only provides a complete, non-degenerate
        starting conformer and a sensible connecting-bond distance for every
        attachment. Each substituent is placed independently.
        """
        skel_positions = self.skeleton.positions
        centroid = skel_positions.mean(axis=0)

        placed = []
        for substituent, skeleton_link, substituent_link in self.substituents:
            sub_positions = substituent.positions
            skel_link_pos = skel_positions[skeleton_link]
            direction = skel_link_pos - centroid
            norm = float(np.linalg.norm(direction))
            if norm < 1e-6:
                direction = np.array([1.0, 0.0, 0.0])
            else:
                direction = direction / norm

            target = skel_link_pos + 1.5 * direction
            translation = target - sub_positions[substituent_link]
            new_positions = sub_positions + translation

            placed.append(
                Molecule(
                    symbols=substituent.chemical_symbols,
                    positions=new_positions,
                    charge=substituent.charge,
                    multiplicity=substituent.multiplicity,
                    frozen_atoms=substituent.frozen_atoms,
                )
            )
        return placed

    def _build_combined_rdmol(self) -> tuple[Chem.Mol, int]:
        """Join the skeleton and every substituent into one RDKit molecule.

        Each fragment is converted to RDKit separately, so its bonds come
        from its own geometry and no spurious cross-fragment bonds are
        inferred from overlapping coordinates. The fragments are merged in
        order (skeleton first, then each substituent), and one explicit
        single bond connects each skeleton link atom to its substituent link
        atom. Every substituent is therefore built into the same molecule
        before a single ETKDG embedding runs, keeping the topology correct
        regardless of the rough substituent placement.
        """
        n_skeleton = len(self.skeleton)
        placed_subs = self._rough_place_substituents()

        combined = self.skeleton.to_rdkit()
        offset = n_skeleton
        bond_pairs = []
        for (substituent, skeleton_link, substituent_link), placed in zip(
            self.substituents, placed_subs
        ):
            substituent_rdmol = placed.to_rdkit()
            bond_pairs.append((skeleton_link, offset + substituent_link))
            combined = Chem.CombineMols(combined, substituent_rdmol)
            offset += len(substituent)

        editable = Chem.RWMol(combined)
        for skeleton_atom, substituent_atom in bond_pairs:
            editable.AddBond(
                skeleton_atom, substituent_atom, Chem.BondType.SINGLE
            )
        combined_mol = editable.GetMol()
        combined_mol.UpdatePropertyCache(strict=False)
        return combined_mol, n_skeleton

    @staticmethod
    def _score(molecule: Chem.Mol, conformer_id: int) -> float:
        """Score a conformer by UFF energy (lower is better).

        Returns ``inf`` if UFF is unavailable, so that a successful
        embedding is still preferred over failure.
        """
        try:
            force_field = AllChem.UFFGetMoleculeForceField(
                molecule, confId=conformer_id
            )
            if force_field is None:
                return float("inf")
            return float(force_field.CalcEnergy())
        except Exception:
            return float("inf")

    def _merge_frozen(self) -> Optional[list]:
        """Merge skeleton and substituent frozen-atom flags (skeleton first)."""
        skeleton_frozen = self.skeleton.frozen_atoms
        substituent_frozens = [
            substituent.frozen_atoms for substituent, _, _ in self.substituents
        ]
        if skeleton_frozen is None and all(
            frozen is None for frozen in substituent_frozens
        ):
            return None
        merged = list(
            skeleton_frozen
            if skeleton_frozen is not None
            else [0] * len(self.skeleton)
        )
        for (substituent, _, _), frozen in zip(
            self.substituents, substituent_frozens
        ):
            merged.extend(
                frozen if frozen is not None else [0] * len(substituent)
            )
        return merged

    def _rdmol_to_molecule(
        self, molecule: Chem.Mol, conformer_id: int
    ) -> Molecule:
        """Convert an ETKDG result RDKit molecule back to a ``Molecule``."""
        conformer = molecule.GetConformer(conformer_id)
        n_atoms = molecule.GetNumAtoms()
        positions = np.array(
            [
                [
                    conformer.GetAtomPosition(i).x,
                    conformer.GetAtomPosition(i).y,
                    conformer.GetAtomPosition(i).z,
                ]
                for i in range(n_atoms)
            ],
            dtype=float,
        )
        symbols = list(self.skeleton.chemical_symbols)
        for substituent, _, _ in self.substituents:
            symbols.extend(substituent.chemical_symbols)
        return Molecule(
            symbols=symbols,
            positions=positions,
            charge=self.skeleton.charge,
            multiplicity=self.skeleton.multiplicity,
            frozen_atoms=self._merge_frozen(),
        )

    def run(self) -> Optional[Molecule]:
        """Embed the combined molecule with ETKDGv3 and keep the best conformer.

        Returns
        -------
        Molecule or None
            Combined skeleton + substituent molecule for the lowest-energy
            conformer, or ``None`` if embedding fails.
        """
        use_global = bool(self.options.get("use_global_optimization", False))
        force_field = str(self.options.get("force_field", "none")).lower()
        if force_field not in self._SUPPORTED_FORCE_FIELDS:
            raise ValueError(
                f"Unsupported force_field '{force_field}'; expected one of "
                f"{self._SUPPORTED_FORCE_FIELDS}."
            )
        num_conformers = max(1, int(self.options.get("num_conformers", 1)))

        combined, n_skeleton = self._build_combined_rdmol()
        molecule = self._prepare_molecule(combined)

        # In local mode the skeleton atoms are constrained to the coordinates
        # taken from the combined starting geometry; global mode frees them.
        reference_positions = (
            {}
            if use_global
            else self._fixed_reference_positions(combined, n_skeleton)
        )
        parameters = self._embedding_parameters(reference_positions)

        try:
            conformer_ids = list(
                AllChem.EmbedMultipleConfs(
                    molecule, numConfs=num_conformers, params=parameters
                )
            )
        except Exception as error:
            logger.warning(f"ETKDG embedding failed: {error}")
            conformer_ids = []

        if not conformer_ids:
            logger.error("ETKDG produced no conformers.")
            return None

        fixed_indices = () if use_global else tuple(range(n_skeleton))

        best_conformer_id = None
        best_score = None
        for conformer_id in conformer_ids:
            # Local mode: restore the skeleton to its exact coordinates.
            self._set_fixed_coordinates(
                molecule, conformer_id, reference_positions
            )
            if force_field != "none":
                try:
                    self._force_field_optimize(
                        molecule, conformer_id, fixed_indices, force_field
                    )
                    self._set_fixed_coordinates(
                        molecule, conformer_id, reference_positions
                    )
                except Exception as error:
                    logger.warning(
                        f"Force field '{force_field}' skipped: {error}"
                    )
            score = self._score(molecule, conformer_id)
            if best_score is None or score < best_score:
                best_score = score
                best_conformer_id = conformer_id

        return self._rdmol_to_molecule(molecule, best_conformer_id)
