import copy
import datetime
import logging
import os
import re
import shutil
import sqlite3
import time
from collections.abc import Sequence
from contextlib import suppress
from functools import cached_property
from itertools import chain, pairwise
from random import shuffle
from tempfile import TemporaryDirectory

import ase.db
import ase.io
import h5py
import numpy as np
import pymysql
from ase import Atoms
from ase.calculators.calculator import PropertyNotImplementedError
from ase.calculators.singlepoint import SinglePointCalculator
from ase.db import connect
from ase.db.sqlite import SQLite3Database
from pymatgen.io.ase import AseAtomsAdaptor

from pyatoms.database.mysql import MySQLDatabase
from pyatoms.io.ase.atoms import AtomsWrapper, IAtomsWrapper
from pyatoms.io.ase.patches import _get_external_table_names
from pyatoms.utils.filelock import FileLock
from pyatoms.utils.mixins import IterableMixin
from pyatoms.utils.multiprocess import parallelizable
from pyatoms.utils.utils import file_cache, isoutlier, string2index, timed

logger = logging.getLogger(__name__)
SQLite3Database._get_external_table_names = _get_external_table_names


class Hashabledict(dict):
    def __hash__(self):
        return hash(frozenset(self))


class DatabaseFile:
    def __init__(self, filename):
        if not os.path.exists(filename):
            raise FileNotFoundError(f'{filename} could not be found')
        self.filename = filename

    def delete_entries(self, ids):
        if any(id < 1 for id in ids):
            raise ValueError('ids must be greater than 0. (1-indexed)')

        if any(id > self.num_images for id in ids):
            raise ValueError(
                f'ids must be less than or equal to the number of images in the database ({self.num_images})'
            )

        with connect(self.filename) as db:
            db.delete(ids)
            logger.info(f'id={ids} deleted from {self.filename}')

    @property
    def _lockfile(self):
        dirname = os.path.dirname(self.filename)
        basename = os.path.basename(self.filename)
        return os.path.join(dirname, f'.{basename}.lock')

    def to_dataset(self):
        return Dataset._from_database_file(self.filename)

    @property
    def num_images(self):
        with connect(self.filename) as db:
            return len(db)

    def backup(self, filename, timestamp=False):
        if timestamp:
            ts = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
            filename = f'{filename}.{ts}'

        with FileLock(self._lockfile):
            shutil.copy(src=self.filename, dst=filename)

    @property
    def is_corrupted(self):
        try:
            self.num_images  # noqa: B018 - needed to check if reading database triggers exception
            return False
        except (sqlite3.DatabaseError, sqlite3.OperationalError):
            logger.info('Encountered sqlite3 error while reading dataset. Dataset is likely corrupted')
            return True
        except (TypeError, OSError):
            logger.exception(
                'Encountered TypeError/OSError while reading dataset. Dataset is likely corrupted '
                '- check carefully the error message!'
            )
            return True
        except Exception as e:
            raise e


class Dataset(Sequence, IterableMixin):
    def __init__(self, images=None, key_value_pairs=None, datas=None, metadata=None, ase_db_attributes=None):
        if images is None:
            images = []

        if key_value_pairs is None:
            key_value_pairs = [{}] * len(images)

        if datas is None:
            datas = [{}] * len(images)

        if metadata is None:
            metadata = {}

        if ase_db_attributes is None:
            ase_db_attributes = [{}] * len(images)

        if len(images) != len(key_value_pairs):
            raise ValueError(
                f'Number of images and key_value_pairs not equal [{len(images)} != {len(key_value_pairs)}]'
            )

        if len(images) != len(datas):
            raise ValueError('Length of images and datas have to be equal')

        if len(images) != len(ase_db_attributes):
            raise ValueError('Length of images and ase_db_attributes have to be equal')

        images = [self._copy_atoms(image) for image in images]

        self._images = list(images)
        self._key_value_pairs = list(key_value_pairs)
        self._datas = list(datas)
        self._ase_db_attributes = list(ase_db_attributes)

        self.metadata = metadata

    def _iterables(self):
        return {
            'images': self._images,
            'key_value_pairs': self._key_value_pairs,
            'datas': self._datas,
            'ase_db_attributes': self._ase_db_attributes,
        }

    def __repr__(self):
        return f'<{self.__class__.__qualname__}(num_images={len(self)})>'

    def __delitem__(self, key):
        del self._images[key]
        del self._key_value_pairs[key]
        del self._datas[key]
        del self._ase_db_attributes[key]

    def __add__(self, other):
        if self.__class__ != other.__class__:
            raise NotImplementedError(f'Cannot add {self.__class__} to {other.__class__}')

        return self.from_images(
            images=self._images + other._images,
            key_value_pairs=self._key_value_pairs + other._key_value_pairs,
            datas=self._datas + other._datas,
            ase_db_attributes=self._ase_db_attributes + other._ase_db_attributes,
        )

    def __len__(self):
        assert len(self._images) == len(self._key_value_pairs)
        return len(self._images)

    def __getitem__(self, key):
        if self._is_int(key):
            return self._images[int(key)]

        new_iterables = {k: self._index(v, key) for k, v in self._iterables().items()}
        return type(self)(**new_iterables, metadata=self.metadata)

    def rattle(self, **kwargs):
        new = self.copy()
        for image in new._images:
            image.rattle(**kwargs)
        return new

    @property
    def ase_ids(self):
        return [ase_db_attribute['id'] for ase_db_attribute in self.ase_db_attributes]

    def average_num_atoms(self):
        return np.mean([len(image) for image in self._images])

    def evaluate_errors(self, calculator, **kwargs):
        from pyatoms.analysis.dataset import DatasetErrors

        return DatasetErrors.from_calculator_and_dataset(calculator=calculator, dataset=self, **kwargs)

    def remove_calculators(self):
        new_images = [image.copy() for image in self._images]  # copying removes calculators
        return Dataset.from_images(
            images=new_images,
            key_value_pairs=self.key_value_pairs,
            datas=self.datas,
            ase_db_attributes=self.ase_db_attributes,
        )

    def remove_images_without_calc(self):
        return self[[image.calc is not None for image in self.images]]

    def _immutable_images(self):
        return [IAtomsWrapper.from_atoms(image) for image in self._images]

    def remove_duplicates(self):
        hashable_kvps = [Hashabledict(pair) for pair in self._key_value_pairs]
        hashable_datas = [Hashabledict(pair) for pair in self._datas]
        tups = list(zip(self._immutable_images(), hashable_kvps, hashable_datas, strict=False))

        unique_indices = [tups.index(i) for i in set(tups)]

        return self[unique_indices]

    def has_forces(self):
        """Returns True if all images have forces."""
        try:
            [image.get_forces() for image in self._images]
        except PropertyNotImplementedError:
            return False
        except RuntimeError as e:
            if 'Atoms object has no calculator' in str(e):
                return False
        return True

    def has_energy(self, force_consistent=False):
        """Returns True if all images have energies."""
        try:
            [image.get_potential_energy(force_consistent=force_consistent) for image in self._images]
        except PropertyNotImplementedError:
            return False
        except RuntimeError as e:
            if 'Atoms object has no calculator' in str(e):
                return False
        return True

    def has_stress(self):
        """Returns True if all images have stresses."""
        try:
            [image.get_stress() for image in self._images]
        except PropertyNotImplementedError:
            return False
        except RuntimeError as e:
            if 'Atoms object has no calculator' in str(e):
                return False
        return True

    def write_to_poscars(self, folder):
        if os.path.exists(folder):
            raise ValueError(f'{folder} already exists')

        os.mkdir(folder)
        for idx, image in enumerate(self._images):
            subfolder = os.path.join(folder, f'{idx}')
            poscar = os.path.join(subfolder, 'POSCAR')

            os.mkdir(subfolder)
            image.write(poscar, vasp5=True, direct=True)

    def write(self, filename, append=False, format=None, overwrite=False, **kwargs) -> None:
        """Writes the dataset to a file.

        Args:
            filename (str): Filename to write to.
            append (bool): Append to the database file if it already exists.
            format (str): Format of the file to write to. If None, then the format is inferred from the filename.
            overwrite (bool): Overwrite the database file if it already exists.
            **kwargs (dict): Keyword arguments for ase.io.write.
        """
        if (not append and not overwrite) and os.path.exists(filename):
            raise FileExistsError(f'{filename} already exists!')

        # Note: If we reach here, then the database file will be overwritten. Therefore,
        # we do not need to explicitly remove the database even if overwrite=True.
        if filename.endswith('.h5'):
            self._write_to_hdf(filename=filename, force_consistent=False, **kwargs)
        elif filename.endswith('.db'):
            self._write_to_sqlite3_database(filename=filename, append=append, **kwargs)
        elif filename.startswith('mysql://'):
            self._write_to_mysql_database(filename=filename, append=append, **kwargs)
        else:
            ase.io.write(filename=filename, images=self._images, append=append, format=format)

    def _write_to_mysql_database(self, filename, append=True, write_metadata=True):
        database = MySQLDatabase.from_path(filename)

        if database.exists and not append:
            raise ValueError('append cannot be False as database already exists')

        if not database.exists:
            database = database.create()

        metadata = self.metadata if write_metadata else None

        writer = DatabaseWriter(
            images=self._images, key_value_pairs=self._key_value_pairs, datas=self._datas, metadata=metadata
        )
        writer.write(filename=filename, append=append)

    def _write_to_sqlite3_database(self, filename, append=True, write_metadata=True):
        metadata = self.metadata if write_metadata else None

        writer = DatabaseWriter(
            images=self._images, key_value_pairs=self._key_value_pairs, datas=self._datas, metadata=metadata
        )
        writer.write(filename=filename, append=append)

    def _write_to_hdf(self, filename, force_consistent=False):
        writer = HDFWriter(images=self._images, force_consistent=force_consistent)
        writer.write(filename=filename)

    def unique_chemical_symbols(self):
        elements = set()
        for image in self._images:
            elements.update(tuple(image.get_chemical_symbols()))
        return sorted(elements)

    def unique_elements(self):
        return self.unique_chemical_symbols()

    @property
    def images(self):
        return self._copied_images

    @property
    def _copied_images(self):
        return [image.copy(copy_calculator=True) for image in self._images]

    @property
    def pymatgen_structures(self):
        return [AseAtomsAdaptor.get_structure(atoms=image) for image in self._images]

    @property
    def key_value_pairs(self):
        return copy.deepcopy(self._key_value_pairs)

    @key_value_pairs.setter
    def key_value_pairs(self, key_value_pairs):
        if (num_kvps := len(key_value_pairs)) != len(self):
            raise ValueError(f'Length of key value pairs ({num_kvps}) not equal to length of {self}')

        self._key_value_pairs = key_value_pairs

    @property
    def datas(self):
        return self._datas.copy()

    @property
    def ase_db_attributes(self):
        return self._ase_db_attributes.copy()

    def energies(self, force_consistent=False):
        return np.array([a.get_potential_energy(force_consistent=force_consistent) for a in self._images])

    @timed(log_threshold=1)
    def find(self, atoms, return_none=True):
        """Returns the image that matches the input atoms."""
        try:
            idx = self._images.index(atoms)
            return self._images[idx]
        except ValueError as e:
            if 'is not in list' in str(e) and return_none:
                return None
            raise ValueError('Could not find atoms in dataset') from e

    def with_new_calculator(self, calculator):
        """Apply a new calculator to all images in the dataset."""
        images = []
        for _image in self._images:
            image = _image.copy()
            image.calc = calculator
            images.append(image)
        return Dataset.from_images(
            images, key_value_pairs=self.key_value_pairs, datas=self.datas, ase_db_attributes=self.ase_db_attributes
        )

    def energetic_outliers(self, calculator=None, value_threshold=None, zscore_threshold=3.5, force_consistent=False):
        if calculator is not None:
            new_images = []
            for _image in self._images:
                image = _image.copy()
                image.set_calculator(calculator)
                new_images += [image]

            calc_energies = [image.get_potential_energy(force_consistent=force_consistent) for image in new_images]
            calc_energies = np.array(calc_energies)

            saved_energies = np.array(self.energies(force_consistent=force_consistent))
            energetic_diff = np.abs(calc_energies - saved_energies)
            logger.info(f'Diff: {energetic_diff}')
        else:
            energetic_diff = np.array(self.energies(force_consistent=force_consistent))

        is_outlier = isoutlier(energetic_diff, value_threshold=value_threshold, zscore_threshold=zscore_threshold)
        idx_outliers = np.where(is_outlier)[0]
        logger.info(f'Images {idx_outliers} are outliers (energetic_diff: {energetic_diff[is_outlier]})')
        return is_outlier

    def finalize_results(self):
        """If image calculators are not SinglePointCalculators, then the results have to be finalized.

        This is useful for writing the results to a database, as the energies and forces will only
        be written if the calculator is a SinglePointCalculator.
        """
        images = self._copied_images
        for image in images:
            results = image.results()
            calc = SinglePointCalculator(atoms=image, **results)
            image.calc = calc
        return Dataset.from_images(
            images, key_value_pairs=self.key_value_pairs, datas=self.datas, ase_db_attributes=self.ase_db_attributes
        )

    def remove_energetic_outliers(self, calculator=None, value_threshold=None, zscore_threshold=3.5):
        is_outlier = self.energetic_outliers(
            calculator=calculator, value_threshold=value_threshold, zscore_threshold=zscore_threshold
        )
        return self[~is_outlier]

    def plot_energies_histogram(self, bins=None, binwidth=None, density=True, ax=None, label=None):
        from pyatoms.utils.plot import histogram

        histogram(
            values=self.energies(),
            bins=bins,
            binwidth=binwidth,
            density=density,
            ax=ax,
            label=label,
            xlabel='Energies / eV',
        )

    def forces(self, apply_constraint=False, ravel=False):
        forces = [a.get_forces(apply_constraint=apply_constraint) for a in self._images]
        return np.hstack([f.ravel() for f in forces]) if ravel else np.array(forces, dtype=object)

    def stress(self, apply_constraint=False, voigt=True, suppress_not_implemented_error=False):
        """Returns the stress for each image in the dataset.

        Args:
            apply_constraint (bool): Whether to apply the constraints to the stress.
            suppress_not_implemented_error (bool): Whether to ignore PropertyNotImplementedError.
                If property is not implemented for an image, then the stress will be set to NaN.
            voigt (bool): Whether to return the stress in Voigt notation.
        """
        try:
            stress = np.array([a.get_stress(apply_constraint=apply_constraint, voigt=voigt) for a in self._images])
        except ase.calculators.calculator.PropertyNotImplementedError as e:
            if suppress_not_implemented_error:
                logger.info('Calculator does not support stress')
                stress = np.array([[np.nan] * 6] * len(self))
            else:
                raise e
        return stress

    def max_forces(self, apply_constraint=False):
        return np.array([a.max_force(apply_constraint=apply_constraint) for a in self._images])

    def copy(self):
        return self.from_images(
            images=self._images,
            key_value_pairs=self._key_value_pairs.copy(),
            datas=self._datas.copy(),
            ase_db_attributes=self._ase_db_attributes.copy(),
        )

    def sorted_by_energies(self, force_consistent=False):
        indices = np.argsort(self.energies(force_consistent=force_consistent))
        return self[indices]

    def shuffle(self):
        """Perform in-place shuffling of the dataset."""
        indices = list(range(len(self)))
        shuffle(indices)

        # Atoms object cannot be passed into np.array,
        # they will get flattened into individual Atom objects!
        self._images = [self._images[i] for i in indices]
        self._key_value_pairs = list(np.array(self._key_value_pairs)[indices])
        self._datas = list(np.array(self._datas)[indices])
        self._ase_db_attributes = list(np.array(self._ase_db_attributes)[indices])

    def extend(self, images, key_value_pairs=None, datas=None, ase_db_attributes=None):
        """Perform in-place extension of the dataset."""
        if key_value_pairs is None:
            key_value_pairs = [{}] * len(images)

        if datas is None:
            datas = [{}] * len(images)

        if ase_db_attributes is None:
            ase_db_attributes = [{}] * len(images)

        if len(images) != len(key_value_pairs):
            raise ValueError('Length of images and key_value_pairs have to be equal')

        if len(images) != len(datas):
            raise ValueError('Length of images and datas have to be equal')

        if len(images) != len(ase_db_attributes):
            raise ValueError('Length of images and ase_db_attributes have to be equal')

        images = [self._copy_atoms(image) for image in images]
        logger.info(f'Extending {self} by {len(images)} images')
        self._images.extend(images)
        self._key_value_pairs.extend(key_value_pairs)
        self._datas.extend(datas)
        self._ase_db_attributes.extend(ase_db_attributes)

    def _copy_atoms(self, atoms):
        new_atoms = AtomsWrapper.from_atoms(atoms)
        return new_atoms.copy(copy_calculator=True)

    def remove_images(self, indices):
        """Removes images from the dataset."""
        new = self.copy()
        for idx in sorted(indices, reverse=True):
            del new[idx]
        return new

    def chunks(self, chunk_size=None, num_chunks=None):
        if chunk_size is not None and num_chunks is not None:
            raise ValueError('Can only choose chunk_size or num_chunks')

        if num_chunks is None:
            num_chunks = int(np.ceil(len(self) / chunk_size))

        if chunk_size is None:
            chunk_size = int(np.ceil(len(self) / num_chunks))

        chunks = []
        for i in range(num_chunks):
            indices = range(i * chunk_size, min((i + 1) * chunk_size, len(self)))
            chunk = self[indices]
            chunks.append(chunk)
        return chunks

    def split_into_two(self, fraction=None, num_datapoints=None):
        """Split into two or datasets."""
        if fraction is not None and num_datapoints is not None:
            raise ValueError('Specify either fraction or num_datapoints, not both.')

        if num_datapoints is not None:
            num_datapoints = [num_datapoints, len(self) - num_datapoints]

        if fraction is not None:
            fraction = [fraction, 1 - fraction]

        datasets = self.split(fraction=fraction, num_datapoints=num_datapoints)
        return datasets[0], datasets[1]

    def split(self, fraction=None, num_datapoints=None):
        """Split into two or more datasets."""
        if fraction is not None and num_datapoints is not None:
            raise ValueError('Specify either fraction or num_datapoints, not both.')

        if fraction is not None and (sum_fractions := sum(fraction)) != 1:
            raise ValueError(f'Sum of fractions {sum_fractions} is not equal to 1.')

        if fraction is not None:
            num_datapoints = [int(f * len(self)) for f in fraction]

            # Take care of any rounding errors in the last set
            num_datapoints[-1] += len(self) - sum(num_datapoints)

        if (sum_datapoints := sum(num_datapoints)) != len(self):
            raise ValueError(
                f'Sum of num_datapoints {sum_datapoints} not equal to the number of images in the dataset {len(self)}.'
            )

        def _split_once(start, stop):
            return self[range(start, stop)]

        cumsum = [0, *np.cumsum(num_datapoints).tolist()]
        return [_split_once(start, stop) for start, stop in pairwise(cumsum)]

    def split_into_groups_by_chemical_formula(self):
        chemical_formulas = [image.get_chemical_formula() for image in self._images]
        unique_chemical_formulas = list(set(chemical_formulas))

        subdatasets = []
        for uf in unique_chemical_formulas:
            indices = [i for i, cf in enumerate(chemical_formulas) if cf == uf]
            subdatasets.append(self[indices])
        return subdatasets

    def split_by_cells(self):
        cells = [str(image.cell) for image in self._images]
        unique_cells = list(set(cells))
        subdatasets = []
        for ucell in unique_cells:
            indices = [i for i, cell in enumerate(cells) if cell == ucell]
            subdatasets.append(self[indices])
        return subdatasets

    @classmethod
    def combine(cls, datasets):
        images = list(chain.from_iterable([dataset.images for dataset in datasets]))
        kvps = list(chain.from_iterable([dataset.key_value_pairs for dataset in datasets]))
        datas = list(chain.from_iterable([dataset.datas for dataset in datasets]))
        ase_db_attributes = list(chain.from_iterable([dataset.ase_db_attributes for dataset in datasets]))
        return cls(images=images, key_value_pairs=kvps, datas=datas, ase_db_attributes=ase_db_attributes)

    @classmethod
    def empty_dataset(cls, metadata=None):
        return cls(images=[], metadata=metadata)

    @classmethod
    def from_images(cls, images, **kwargs):
        if isinstance(images, Dataset):
            return images

        if not isinstance(images, list):
            images = [images]

        return cls(images=images, **kwargs)

    @classmethod
    def _from_ase_string(cls, string, **kwargs):
        """Reads images from a string in the format [filename]@[index].

        Args:
            string (str): String with format [filename]@[index].
            kwargs (dict): Keyword arguments for Dataset.from_files method.
        """
        if '@' not in string:
            raise ValueError('String has to be in the format [filename]@[index]')

        filename, index = string.split('@')
        return cls.from_file(filename=filename, index=index, **kwargs)

    @classmethod
    def from_regex_filenames(cls, regex_filenames, parallel=True, **kwargs):
        from pyatoms.utils.cli import determine_paths

        dataset_files = determine_paths(
            regex=regex_filenames,
            return_list=True,
            is_dir=False,
        )

        if parallel:
            datasets = cls.from_file(dataset_files, parallel=parallel, **kwargs)
        else:
            datasets = [cls.from_file(filename, **kwargs) for filename in dataset_files]

        dataset = cls.combine(datasets)

        if len(dataset) == 0:
            raise AssertionError(f'Dataset from {regex_filenames} has no images!')
        return dataset

    @parallelizable(progress_bar=True)
    @classmethod
    def from_file(cls, filename, **kwargs):
        """Reads images from a file."""
        # TODO: change to use AtomsWrapper.from_filepath for all
        if '@' in filename and not filename.startswith('mysql://'):
            return cls._from_ase_string(filename, **kwargs)
        if '.db' in filename:
            return cls._from_database_file(filename, **kwargs)
        if '.cfg' in filename:
            return cls._from_mtp_config_file(filename, **kwargs)
        if '.h5' in filename:
            return cls._from_h5_file(filename, **kwargs)
        if 'mysql://' in filename:
            return cls.from_mysql_path(filename, **kwargs)

        return cls._from_ase_readable_file(filename, **kwargs)

    @classmethod
    def _from_ase_readable_file(cls, filename, index=':', **kwargs):
        images = ase.io.read(filename, index=index)
        if isinstance(images, Atoms):
            images = [images]
        return cls.from_images(images, **kwargs)

    @classmethod
    def from_mysql_path(cls, filename, **kwargs):
        return Dataset._from_database_file_cached(filename, **kwargs)

    @classmethod
    def _from_mtp_config_file(cls, filename, **kwargs):
        from atomistic_ml.io.mtp import MTPConfiguration

        configs = MTPConfiguration.list_from_file(filename)
        images = [config.to_ase_atoms() for config in configs]
        return cls.from_images(images, **kwargs)

    @classmethod
    def _from_h5_file(cls, filename, **kwargs):
        reader = HDFReader(filename)
        return reader.read_database()

    @classmethod
    def _from_database_file(cls, filename, **kwargs):
        if os.path.getsize(filename) == 0:
            return cls()
        return Dataset._from_database_file_cached(filename, **kwargs)

    @file_cache(maxsize=64)
    @staticmethod
    def _from_database_file_cached(filename, **kwargs):
        reader = DatabaseReader(filename)
        return reader.read(**kwargs)

    def max_occurrence_key(self, inlist):
        # show results of `ase db final.db --show-values functional/basis`
        # and find out which functional/basis has max occurrences in the db
        values, counts = np.unique(inlist, return_counts=True)
        max_counts = max(counts)
        for i, count in enumerate(counts):
            if count == max_counts:
                return values[i]
        return None

    def clean_db(self, filter_keywords=None):
        """Use to clean the dataset based on the filter_keywords supplied.

        For example, remove data points where different functional/basis set is used
        :param filter_keyword: str or list containing the keywords to filter the db by
        :return: a cleaned dataset.
        """
        # get all images and kv_pairs
        db_images = self.images
        db_kv_pairs = self.key_value_pairs
        db_datas = self.datas
        db_ase_db_attributes = self.ase_db_attributes

        # get all keys for the db; same key for all entries (but may have different values)
        db_keywords = db_kv_pairs[0].keys()

        # filter_keyword is supplied
        if filter_keywords is None:
            logger.info(
                f'No keyword for filtering is supplied. Dataset is not filtered.\n'
                f'Available keywords in the database are {db_keywords} '
            )
            return None

        # check if filter_keywords is a string (e.g. 'functional, basis') or a list (e.g., ['functional', 'basis']
        if isinstance(filter_keywords, str):
            keywords_for_filtering = filter_keywords.split(',')
            # 'functional, basis' --> ['functional', ' basis'] --> ['functional', 'basis']
            keywords_for_filtering = [string.strip() for string in keywords_for_filtering]
        elif isinstance(filter_keywords, list):
            keywords_for_filtering = filter_keywords

        # obtain all corresponding filter values from filter_keywords
        filter_values = []
        for keyword in keywords_for_filtering:
            if keyword not in db_keywords:
                raise ValueError(
                    f'Keyword for filtering "{keyword}" is not one of the available keys in the db: {db_keywords}'
                )

            db_filter = [kv_pair[keyword] for kv_pair in db_kv_pairs]
            # find the max occurrence of the values according to the filter key
            filter_value = self.max_occurrence_key(db_filter)
            filter_values.append(filter_value)

        # keep only Atoms where their kv_pairs have values matching the filter_value
        filtered_images = []
        filtered_kv_pairs = []
        filtered_datas = []
        filtered_ase_db_attributes = []
        for num in range(len(db_kv_pairs)):
            # looping through each member of db_kv_pairs, where db_kv_pairs[num] corresponds to each kv_pair dict
            # returns a list of booleans to see if the filter_value match with db_kv_pairs value
            keep_data = [
                db_kv_pairs[num][keywords_for_filtering[i]] == filter_values[i]
                for i in range(len(keywords_for_filtering))
            ]
            # keep_data is True if all keep_data members are True
            if all(keep_data):
                filtered_images.append(db_images[num])
                filtered_kv_pairs.append(db_kv_pairs[num])
                filtered_datas.append(db_datas[num])
                filtered_ase_db_attributes.append(db_ase_db_attributes[num])
        return Dataset(
            images=filtered_images,
            key_value_pairs=filtered_kv_pairs,
            datas=filtered_datas,
            ase_db_attributes=filtered_ase_db_attributes,
        )

    def to_schnet_dataset(self, force_consistent=False):
        try:
            from schnetpack.data import AtomsData as SchnetAtomsData
        except ImportError as e:
            raise ImportError('Schnetpack needed.') from e
        # makes use of AtomsData()._deprecation_update() from schnetpack.data.atoms.py to read ase db directly
        property_list = []
        available_properties = []
        for at in self.images:
            # All properties need to be stored as numpy arrays.
            # Note: The shape for scalars should be (1,), not ()
            # Note: GPUs work best with float32 data
            energy = np.array([float(at.get_potential_energy(force_consistent=force_consistent))], dtype=np.float32)
            try:
                forces = at.get_forces(apply_constraint=False)
                property_list.append({'energy': energy, 'forces': forces})
                available_properties = ['energy', 'forces']
            except Exception:
                property_list.append({'energy': energy})
                available_properties = ['energy']

        import tempfile

        new_db_temp_dir = tempfile.mkdtemp()
        new_db_path = os.path.join(new_db_temp_dir, 'temp.db')

        new_dataset = SchnetAtomsData(new_db_path, available_properties=available_properties)
        new_dataset.add_systems(self.images, property_list)
        return new_dataset

    @classmethod
    def from_lasp_files(cls, structure_file, forces_file):
        from atomistic_ml.io.lasp import LASPAtomsBuilder

        builder = LASPAtomsBuilder(structure_file=structure_file, forces_file=forces_file)
        images = builder.build()
        return cls.from_images(images)

    def all_cells_equal(self, atol=1e-6):
        cells = [image.cell for image in self._images]
        return np.allclose(cells[0], cells, atol=atol)

    def to_torchani_data(self):
        import tempfile

        import torchani

        h5_dir = tempfile.mkdtemp()
        h5_path = os.path.join(h5_dir, 'temp.h5')
        self.write(h5_path)

        additional_properties = []
        if self.has_forces():
            additional_properties += ['forces']

        return torchani.data.load(h5_path, additional_properties=additional_properties)


class LazyDataset:
    def __init__(self, filename, **kwargs):
        if not os.path.exists(filename):
            raise FileNotFoundError(f'{filename} could not be found')

        filename = os.path.abspath(filename)

        self.filename = filename
        self.kwargs = kwargs

    @cached_property
    def _dataset(self):
        return Dataset.from_file(self.filename, **self.kwargs)

    def copy(self, return_lazy=True):
        """Copy the LazyDataset.

        Args:
            return_lazy (bool): If False, the underlying Dataset will be copied and returned as a Dataset.
                If True, the LazyDataset properties will be copied instead.
        """
        if return_lazy:
            return copy.copy(self)

        return self._dataset.copy()

    def __len__(self):
        """Refined here as getattr doesn't work on built-in methods."""
        return len(self._dataset)

    def __getattr__(self, attr):
        if attr == '__setstate__':
            raise AttributeError(attr)
        return self._dataset.__getattribute__(attr)

    @classmethod
    def from_file(cls, filename, **kwargs):
        return cls(filename, **kwargs)


class HDFWriter:
    def __init__(self, images, force_consistent=False):
        images = [AtomsWrapper.from_atoms(image) for image in images]
        self.images = images
        self.force_consistent = force_consistent

    def write(self, filename):
        if os.path.exists(filename):
            raise FileExistsError(f'Filename {filename} already exists')

        import h5py

        # in h5py version >3 need to explicitly specify "w"
        with h5py.File(filename, 'w') as f:
            chemical_formulas = [image.get_chemical_formula() for image in self.images]
            unique_chemical_formulas = list(set(chemical_formulas))
            for uf in unique_chemical_formulas:
                self._write_group(f=f, formula=uf)

    def _write_group(self, f, formula):
        images = [image for image in self.images if image.get_chemical_formula() == formula]
        images = [image.sorted_by_elements() for image in images]
        coordinates = np.array([image.get_positions() for image in images])
        energies = np.array([image.get_potential_energy(force_consistent=self.force_consistent) for image in images])
        forces = np.array([image.get_forces(apply_constraint=False) for image in images])
        cells = np.array([image.cell for image in images])
        # may need to get volume too?
        # volumes = np.array([image.get_volume() for image in images])
        from h5py import string_dtype

        dt = string_dtype()
        # species = np.array(images[0].get_chemical_symbols(), dtype=dt)
        species = [np.array(image.get_chemical_symbols(), dtype=dt) for image in images]

        group = f.create_group(name=formula)
        group.create_dataset(name='coordinates', data=coordinates)
        group.create_dataset(name='energies', data=energies)
        group.create_dataset(name='forces', data=forces)
        group.create_dataset(name='species', data=species)
        group.create_dataset(name='cells', data=cells)
        # added volumes from the data images
        # group.create_dataset(name='volumes', data=volumes)
        # may need to write cell and pbc too?


class HDFReader:
    """HDFReader to read .h5 files that are created specifically using the HDFWriter above ONLY.

    *may* need functionalities to convert dataset.write_to_hdf() data back to dataset.
    """

    def __init__(self, filename):
        self.filename = filename

    @property
    def group_keys(self):
        """:return: group keys of .h5 file, keys here are formula of atoms from HDFWriter above"""
        with h5py.File(self.filename, 'r') as f:
            return list(f.keys())

    @property
    def base_items(self):
        with h5py.File(self.filename, 'r') as f:
            return list(f.items())

    def get_all_atoms(self):
        images = []
        with h5py.File(self.filename, 'r') as f:
            for key in self.group_keys:
                key_item = list(f.get(key).items())
                images += self._item_to_atoms(key_item)
        return images

    def _item_to_dict(self, item):
        """Convert each item in .h5 database into a dictionary.

        Dictionary containis key-value pair where key is the str of property name e.g., 'species', 'coordinates'
        and value is the value of the property.
        """
        all_data = {}
        for each_tuple in item:
            all_data[each_tuple[0]] = each_tuple[1][()]  # each_tuple[1].value is deprecated
        return all_data

    def _item_to_atoms(self, item):
        cells = None
        all_data = self._item_to_dict(item)
        all_data_keys = list(all_data.keys())
        if 'species' in all_data_keys:
            species = all_data['species']
        if 'coordinates' in all_data_keys:
            coordinates = all_data['coordinates']
        if 'energies' in all_data_keys:
            energies = all_data['energies']
            len_images = len(energies)
        if 'forces' in all_data_keys:
            forces = all_data['forces']
        if 'cells' in all_data_keys:
            cells = all_data['cells']

        # assemble ase Atoms object
        all_atoms_in_item = []
        for i in range(len_images):
            specie = [s.decode() for s in species[i]]
            if cells is not None:
                atoms = Atoms(symbols=specie, positions=coordinates[i, :, :], cell=cells[i, :, :])
            else:
                atoms = Atoms(symbols=specie, positions=coordinates[i, :, :])
            atoms.calc = SinglePointCalculator(
                atoms=atoms,
                energy=energies[i],
                forces=forces[i, :, :],
            )
            all_atoms_in_item.append(atoms)
        return all_atoms_in_item

    def read_database(self):
        images = self.get_all_atoms()
        return Dataset.from_images(images)


class DatabaseReader:
    def __init__(self, filename):
        if '.db' in filename:
            filename = os.path.abspath(filename)
        elif 'mysql://' in filename:
            pass
        else:
            raise ValueError(f'DatabaseReader can only handle .db or mysql:// paths. Got {filename}')
        self.filename = filename

    def read(self, index=':', copy_to_tmp_if_fail=True):
        try:
            return self._read(index=index)
        except sqlite3.OperationalError as e:
            # File system may not support writing/reading of sqlite files, so try different filesystem
            if not self.copy_to_tmp_if_fail:
                raise ValueError('Could not read database file due to sqlite3 error') from e

            logger.info(f'sqlite3 disk I/O error in reading {self.filename}. Copying to tmp and read it there.')
            with TemporaryDirectory() as tempdir:
                tempdb = os.path.join(tempdir, 'tmp.db')
                shutil.copy(src=self.filename, dst=tempdb)
                reader = DatabaseReader(tempdb)
                return reader.read(index=index, copy_to_tmp_if_fail=False)

    def _read(self, index=':'):
        metadata = self._read_metadata()
        ids = self.get_ids(index)
        rows = self._read_rows(ids)
        atoms, key_value_pairs, datas, ase_db_attributes = self._rows_to_atoms(rows)

        return Dataset(
            images=atoms,
            key_value_pairs=key_value_pairs,
            datas=datas,
            metadata=metadata,
            ase_db_attributes=ase_db_attributes,
        )

    def get_ids(self, index=':'):
        """Get the ids of the images to read from the database.

        Cannot assume that the ids are sequential (holes can be left behind via deletion etc.),
        so need to get all ids and then select the ones to read.
        """
        with ase.db.connect(self.filename) as db, db.managed_connection() as con:
            cur = con.cursor()
            cur.execute('SELECT id FROM systems')
            ids = [result[0] for result in cur.fetchall()]

        return ids[string2index(index, return_iterable=True)]

    def _read_metadata(self):
        with ase.db.connect(self.filename) as db:
            return db.metadata

    @property
    def is_mysql_database(self):
        return 'mysql://' in self.filename

    @property
    def num_images(self):
        with ase.db.connect(self.filename) as db:
            return len(db)

    def _read_rows(self, ids):
        """Reads rows from the database given the ids."""
        if len(ids) == 0:
            return []

        db = ase.db.connect(self.filename)
        with db.managed_connection() as con:
            cur = con.cursor()

            # placeholder is '?' for sqlite3 and '%s' for mysql
            placeholder = '?' if not self.is_mysql_database else '%s'
            placeholders = ','.join([placeholder] * len(ids))
            cur.execute(f'SELECT * FROM systems WHERE id IN ({placeholders})', ids)
            values = cur.fetchall()

        return [db._convert_tuple_to_row(v) for v in values]

    def _rows_to_atoms(self, rows):
        atoms = []
        kvps = []
        datas = []
        ase_db_attributes = []

        for atom_row in rows:
            atoms.append(atom_row.toatoms())
            kvps.append(atom_row.key_value_pairs)
            datas.append(atom_row.data)
            ase_db_attributes.append(
                {
                    'id': atom_row.id,
                    'unique_id': atom_row.unique_id,
                    'ctime': atom_row.ctime,
                    'mtime': atom_row.mtime,
                    'user': atom_row.user,
                }
            )
        return atoms, kvps, datas, ase_db_attributes


class DatabaseWriter:
    """Generalized writer for ASE databases.

    Args:
        images (list): List of Atoms
        metadata (dict): Metadata of the database. If metadata is None, it will not be written to the database.
        key_value_pairs (list[dict]): List of dictionaries corresponding to the key_value_pairs of each image
        datas (list[dict]): List of dictionaries corresponding to the datas of each image
    """

    def __init__(self, images=None, key_value_pairs=None, datas=None, metadata=None):
        if images is None:
            images = []

        if key_value_pairs is None:
            key_value_pairs = [{}] * len(images)

        if datas is None:
            datas = [{}] * len(images)

        self.images = images
        self.key_value_pairs = key_value_pairs
        self.datas = datas
        self.metadata = metadata

    def write(self, filename, append=True, tmp_staging=True, overwrite_metadata=False) -> list:
        def _write_metadata(db, overwrite_metadata):
            if self.metadata is None:
                # None indicates no writing of metadata
                return

            if db.metadata not in ({}, self.metadata) and not overwrite_metadata:
                raise ValueError(
                    f'self.metadata ({self.metadata} does not match that in database ({db.metadata}). '
                    f'Set overwrite_metadata=True to overwrite'
                )

            db.metadata = self.metadata

        def _write_images(db):
            ids = []
            for image, pairs, data in zip(self.images, self.key_value_pairs, self.datas, strict=False):
                self._remove_pyatoms_constraints(image)
                id = db.write(image, key_value_pairs=pairs, data=data)
                ids.append(id)
            return ids

        if append is False:
            with suppress(FileNotFoundError):
                os.remove(filename)

        self._remove_old_lock_files(filename)
        filelock = self._get_filelock_path(filename)
        try:
            # Use own filelock as it is much better in terms of performance for concurrent writes.
            # ASE DB filelock sometimes gets stuck

            with ase.db.connect(filename) as db, FileLock(filelock):
                ids = _write_images(db)
                _write_metadata(db, overwrite_metadata)

        except sqlite3.OperationalError:
            # Error with some filesystems. Use tmp staging and copy back.
            if tmp_staging:
                ids = self._stage_using_tmp_and_copy(filename, append)
            else:
                logger.exception(f'Failed to write to {filename}')
                raise
        except pymysql.err.OperationalError as e:
            is_unknown_db_error = re.search(r'1049, "Unknown database.*?', str(e))
            if not is_unknown_db_error:
                raise

            # Database does not exist, so create one and try again
            MySQLDatabase.create_from_path(filename)
            ids = self.write(filename=filename, append=append, tmp_staging=tmp_staging)
        return ids

    def _get_filelock_path(self, filename):
        if 'mysql://' in filename:
            return None

        dirname = os.path.dirname(filename)
        basename = os.path.basename(filename)
        return os.path.join(dirname, f'.{basename}.filelock')

    def _remove_pyatoms_constraints(self, image):
        """Remove pyatoms constraints from an image. Because pyatoms constraints cannot be written to db."""
        from pyatoms.io.ase.constraints import HookeanPyatoms

        constraints = image.constraints
        constraints = [c for c in constraints if not isinstance(c, HookeanPyatoms)]
        image.set_constraint(constraints)

    def _stage_using_tmp_and_copy(self, filename, append):
        """Try to write to tempdir and then copy the file back to the original path.

        This is under the assumption that the tempdir is likely to have a different filesystem.
        Useful if filesystem of original path may not support writing/reading of sqlite files.
        """
        logger.info('Encountered sqlite3 disk I/O error. Trying to write to tmp and copy back')
        with TemporaryDirectory() as tempdir:
            tmp_filename = os.path.join(tempdir, 'tmp.db')
            if filename == tmp_filename:
                raise ValueError

            if append and os.path.exists(filename):
                shutil.copy(src=filename, dst=tmp_filename)

            ids = self.write(filename=tmp_filename, tmp_staging=False, append=append)
            shutil.move(tmp_filename, filename)
        return ids

    def _remove_old_lock_files(self, filename, age_threshold=300):
        """Remove old lock files that might be left behind after database is deleted."""
        lockfile = filename + '.lock'
        journal_file = filename + '-journal'
        files = [lockfile, journal_file, filename]
        if not os.path.exists(lockfile) and not os.path.exists(journal_file):
            return

        # Remove lockfiles if db is not present.
        if not os.path.exists(filename):
            for f in files:
                with suppress(FileNotFoundError):
                    os.remove(f)
            return

        # If db is present, only remove lockfiles if they are old
        while True:
            try:
                st = os.stat(lockfile)
                logger.info(f'Lock file {lockfile} exists.')
            except FileNotFoundError:
                # Could be that we caught the DB when it was writing something. If so, lockfile will
                # disappear, and we don't need to do anything
                logger.info('Lock file is gone')
                break

            age = time.time() - st.st_mtime
            logger.info(f'Lock file age {age}, threshold: {age_threshold}')
            if age > age_threshold:
                logger.info(f'Lock file age {age} s > {age_threshold} s. Removing lock and other files.')

                for f in files:
                    with suppress(FileNotFoundError):
                        os.remove(f)
                break
            time.sleep(1)
