import logging
import multiprocessing
import typing
from functools import partial
import numpy as np
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.structure_matcher import StructureMatcher
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
from chemsmart.io.molecules.structure import Molecule

logger = logging.getLogger(__name__)


def equality_func(x, y):
    return x == y


class Grouper:
    def __init__(self, objects, equality_comparison_func, num_procs=1):
        if equality_comparison_func is None:
            equality_comparison_func = equality_func

        self.objects = objects
        self.equality_comparison_func = equality_comparison_func
        self.num_procs = num_procs

    def group(self):
        if len(self.objects) == 0:
            return [], []

        if self.num_procs == 1:
            groups, idx_groups = self._group_no_processes()
        else:
            groups, idx_groups = self._group_with_processes()

        groups, idx_groups = self._post_processing(groups, idx_groups)

        return groups, idx_groups

    def unique(self):
        groups, _ = self.group()
        return [group[0] for group in groups]

    def unique_indices(self):
        _, idx_groups = self.group()
        return [idx_group[0] for idx_group in idx_groups]

    def _group_no_processes(self):
        raise NotImplementedError

    def _group_with_processes(self):
        raise NotImplementedError

    def _post_processing(self, groups, idx):
        # Subclasses can implement
        return groups, idx


class SelfConsistentGrouper(Grouper):
    def __init__(self, objects, equality_comparison_func, num_procs=1):
        super().__init__(
            objects=objects,
            equality_comparison_func=equality_comparison_func,
            num_procs=num_procs,
        )

    def _group_no_processes(self):
        return self._group_with_processes()

    def _group_with_processes(self):
        if len(self.objects) == 0:
            return [], []

        is_similar = self.similarity_matrix()

        sparse_is_similar = csr_matrix(is_similar)
        perms = reverse_cuthill_mckee(sparse_is_similar, symmetric_mode=True)

        reordered = sparse_is_similar[perms, :][:, perms]
        reordered = reordered.toarray()

        block_indices_reordered = self.determine_group_indices(reordered)
        block_indices_original = [
            sorted(perms[indices]) for indices in block_indices_reordered
        ]
        groups = [
            [self.objects[i] for i in indices]
            for indices in block_indices_original
        ]
        return groups, block_indices_original

    def similarity_matrix(self):
        ngroups = len(self.objects)
        is_similar = np.zeros([ngroups, ngroups])
        num_procs = min(self.num_procs, len(self.objects))

        if num_procs == 1:
            for i, o1 in enumerate(self.objects):
                for j, o2 in enumerate(self.objects):
                    if j >= i:
                        is_similar[i, j] = self.equality_comparison_func(
                            o1, o2
                        )
        else:
            with multiprocessing.Pool(num_procs) as p:
                logger.debug(
                    f"Setting up pool of {num_procs} processes to group {len(self.objects)} items"
                )
                for i, o1 in enumerate(self.objects):
                    p_func = partial(self.equality_comparison_func, o1)
                    is_similar[i, i:] = p.map(p_func, self.objects[i:])

        return is_similar + is_similar.T - np.eye(ngroups)

    @staticmethod
    def determine_group_indices(matrix):
        """Groups the rows and columns of the symmetric similarly matrix.

        The similarity matrix, M, consists of 1s and 0s where M_ij == 1 if
        i similar to j and M_ij == 0 is i is not similar to j.

        For a column or row (interchangeable) with index i, i belongs to a group if it
        is similar with all other members of the group. For example:

        >>> matrix = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1]])
        >>> SelfConsistentGrouper.determine_group_indices(matrix)
        [[0, 1], [2, 3]]

        >>> matrix = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        >>> SelfConsistentGrouper.determine_group_indices(matrix)
        [[0, 1], [2], [3]]
        """
        matrix = np.array(matrix)
        assert len(matrix.shape) == 2
        nrows = matrix.shape[0]
        groups = [[]]
        for i in range(nrows):
            found = False
            for group in groups:
                if all(matrix[i, group]):
                    group += [i]  # noqa: PLW2901
                    found = True
                    break

            if not found:
                groups += [[i]]

        return groups


class SequentialGrouper(Grouper):
    def __init__(self, objects, equality_comparison_func=None, num_procs=1):
        super().__init__(
            objects=objects,
            equality_comparison_func=equality_comparison_func,
            num_procs=num_procs,
        )

    def _group_with_processes(self):
        idx_ungrouped = range(len(self.objects))
        groups, idx_grouped, _unique_idx = [], [], []
        num_procs = min(self.num_procs, len(self.objects))
        with multiprocessing.Pool(num_procs) as p:
            logger.debug(
                f"Setting up pool of {num_procs} processes to group {len(self.objects)} items"
            )
            while idx_ungrouped:
                ungrouped_objects = [self.objects[i] for i in idx_ungrouped]
                logger.debug(f"{len(ungrouped_objects)} ungrouped objects")

                p_func = partial(
                    self.equality_comparison_func, ungrouped_objects[0]
                )
                matches = p.map(p_func, ungrouped_objects)

                if not matches[0]:
                    raise ValueError("Object is not equal to self.")

                groups.append(
                    [
                        o
                        for o, tf in zip(
                            ungrouped_objects, matches, strict=False
                        )
                        if tf
                    ]
                )
                idx_grouped.append(
                    [
                        i
                        for i, tf in zip(idx_ungrouped, matches, strict=False)
                        if tf
                    ]
                )

                idx_ungrouped = [
                    i
                    for i, tf in zip(idx_ungrouped, matches, strict=False)
                    if not tf
                ]
        logger.debug(f"{len(groups)} groups found")
        return groups, idx_grouped

    def _group_no_processes(self):
        idx_ungrouped = range(len(self.objects))
        groups, idx_grouped, _unique_idx = [], [], []
        while idx_ungrouped:
            ungrouped_objects = [self.objects[i] for i in idx_ungrouped]

            p_func = partial(
                self.equality_comparison_func, ungrouped_objects[0]
            )
            matches = [p_func(o) for o in ungrouped_objects]

            if not matches[0]:
                raise ValueError("Object is not equal to self.")

            groups.append(
                [
                    o
                    for o, tf in zip(ungrouped_objects, matches, strict=False)
                    if tf
                ]
            )
            idx_grouped.append(
                [
                    i
                    for i, tf in zip(idx_ungrouped, matches, strict=False)
                    if tf
                ]
            )

            idx_ungrouped = [
                i
                for i, tf in zip(idx_ungrouped, matches, strict=False)
                if not tf
            ]

        return groups, idx_grouped


def matching_function(x, y, scale, ltol=0.1, angle_tol=1, stol=0.18):
    return StructureMatcherWrapper.match_structures(
        ltol=ltol,
        angle_tol=angle_tol,
        stol=stol,
        structure=x,
        other_structure=y,
        scale=scale,
    )


class StructuralSelfConsistentGrouper(SelfConsistentGrouper):
    def __init__(
        self,
        atoms: [typing.Iterable[Atoms]],
        num_procs=1,
        scale=True,
        stol=0.18,
    ):
        equality_comparison_func = partial(
            matching_function, scale=scale, stol=stol
        )
        structures = [AseAtomsAdaptor.get_structure(atoms=a) for a in atoms]
        super().__init__(
            objects=structures,
            equality_comparison_func=equality_comparison_func,
            num_procs=num_procs,
        )
        self.atoms = atoms

    def _post_processing(self, groups, idx):
        """Return Atoms instead of Structures."""
        groups = [[self.atoms[i] for i in group] for group in idx]
        return groups, idx

    @classmethod
    def from_atoms(cls, atoms, **kwargs):
        return cls(atoms=atoms, **kwargs)

    @classmethod
    def from_structures(cls, structures, **kwargs):
        atoms = [AseAtomsAdaptor.get_atoms(structure=s) for s in structures]
        return cls.from_atoms(atoms=atoms, **kwargs)


class StructuralSequentialGrouper(SequentialGrouper):
    def __init__(self, atoms, num_procs=1, scale=True, stol=0.18):
        equality_comparison_func = partial(
            matching_function, scale=scale, stol=stol
        )
        structures = [AseAtomsAdaptor.get_structure(atoms=a) for a in atoms]
        super().__init__(
            objects=structures,
            equality_comparison_func=equality_comparison_func,
            num_procs=num_procs,
        )
        self.atoms = atoms

    def _post_processing(self, groups, idx):
        """Return Atoms instead of Structures."""
        groups = [[self.atoms[i] for i in group] for group in idx]
        return groups, idx

    @classmethod
    def from_atoms(cls, atoms, **kwargs):
        return cls(atoms=atoms, **kwargs)

    @classmethod
    def from_structures(cls, structures, **kwargs):
        atoms = [AseAtomsAdaptor.get_atoms(structure=s) for s in structures]
        return cls.from_atoms(atoms=atoms, **kwargs)


class StructureMatcherWrapper:
    @staticmethod
    def match_structures(
        structure,
        other_structure,
        ltol=0.1,
        stol=0.3,
        angle_tol=1,
        primitive_cell=False,
        scale=False,
        **kwargs,
    ):
        """Match two structures using pymatgen StructureMatcher.

        Default values from Tibi:
        I played with the script using LiBS2O72_mp-1020060_2181. Main finding:
        ltol and angle_tol are not useful parameters. I can go down to 0.005 and 0.1
        and that gave the same result as default 0.2 and 5 in every case I tried.
        On the other hand, stol is extremely sensitive. Small adjustment of default can immediately
        unmatch even the perfectly same TSs. Therefore, using default stol is basically
        mandatory while the other two parameters do not matter.
        """
        matcher = StructureMatcher(
            ltol=ltol,
            stol=stol,
            angle_tol=angle_tol,
            primitive_cell=primitive_cell,
            scale=scale,
            **kwargs,
        )

        return bool(matcher.fit(structure, other_structure))

    @staticmethod
    def match_molecules(molecule, other_molecule, **kwargs):
        molecule = Molecule.from_molecule(molecule)
        other_molecule = Molecule.from_molecule(other_molecule)

        if sorted(molecule.chemical_symbols) != sorted(
            other_molecule.chemical_symbols
        ):
            return False

        return StructureMatcherWrapper.match_structures(
            molecule, other_molecule, **kwargs
        )
