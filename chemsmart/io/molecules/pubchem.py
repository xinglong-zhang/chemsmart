"""Functions for querying pubchem for structures, adapted from ASE's pubchem.py module."""

import logging
import os
import urllib.request
import warnings
from contextlib import suppress
from functools import lru_cache
from io import StringIO
from tempfile import TemporaryDirectory
from urllib.error import HTTPError, URLError

from ase.data.pubchem import analyze_input, base_url
from chemsmart.io.molecules.structure import Molecule

try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem.AllChem import EmbedMolecule
except ImportError:
    rdkit = None

logger = logging.getLogger(__name__)


def search_pubchem_raw(search, field, suffix: str = "3d"):
    """Search pubchem for structure.

    Changelog:
    1. Added suffix attribute.
    2. Added silent flag to print statements in except blocks.
    3. Changed default len of conformer_id to 1.
    """
    suffix = "sdf?record_type=" + suffix

    if field == "conformers":
        # We don't use the "compound" flag when looking for conformers
        url = f"{base_url}/{field}/{search}/{suffix}"
    else:
        url = f"{base_url}/compound/{field}/{search}/{suffix}"

    r = urllib.request.urlopen(url)

    # Check if there are confomers and warn them if there are
    if field != "conformers":
        # Set default conformer_ids to len==1, as there are pubchem structures with no conformer information
        # which would return "HTTPError:PUGREST.NotFound"
        conformer_ids = ["conformer information not found."]

        with suppress(HTTPError, URLError):
            # Test if there is any conformer information of the structure
            conformer_ids = search_pubchem_raw(
                search, field, suffix="conformers/JSON"
            )

        if len(conformer_ids) > 1:
            warnings.warn(
                f'The structure "{search}" has more than one '
                "conformer in PubChem. By default, the "
                "first conformer is returned, please ensure"
                " you are using the structure you intend to"
                " or use the "
                "`ase.data.pubchem.pubchem_conformer_search`"
                " function",
                stacklevel=2,
            )

    return r.read().decode("utf-8")


@lru_cache(maxsize=128)
def pubchem_search(*args, fail_silently=True, **kwargs):
    """Search pubchem for structure.

    Changelog:
    1. Try suffix =='3d' and '2d'.
    2. Reads 2d sdf files by rdkit to create more accurate structures.
    3. Simplify parse_pubchem_raw() to read atom structure.
    """

    def _pubchem_search():
        search, field = analyze_input(*args, **kwargs)

        try:
            raw_pubchem = search_pubchem_raw(search, field, suffix="3d")
        except (HTTPError, URLError) as e:
            error_substrings = ["HTTP Error 400", "HTTP Error 404"]
            if not any(substring in str(e) for substring in error_substrings):
                raise e
            logger.info(
                f"Error getting 3D {field} structure from pubchem for {search}, trying 2D"
            )

            try:
                raw_pubchem = search_pubchem_raw(search, field, suffix="2d")
            except (HTTPError, URLError) as e:
                raise ValueError(
                    "Error getting structure from pubchem!"
                ) from e

            raw_pubchem = _pubchem_2d_to_3d(raw_pubchem)
            logger.info(f"Successfully got 2D {field} structure for {search}")

        f_like = StringIO(raw_pubchem)

        from chemsmart.utils.utils import sdf2molecule

        sdflines = f_like.readlines()
        molecule = sdf2molecule(sdflines)
        assert isinstance(molecule, Molecule)
        return molecule

    try:
        return _pubchem_search()
    except ValueError as e:
        if fail_silently:
            logger.info(
                "Error getting structure from pubchem. Returning None as fail_silently=True"
            )
            return None
        raise e


def _pubchem_2d_to_3d(data):
    """Transform 2D SDF into 3D. Assumes only one molecule in data."""
    if rdkit is None:
        raise ImportError(
            "rdkit package needed to convert PubChem 2d structures into 3d structures. "
            "Please install rdkit via `pip install rdkit to continue using StructureBuilder."
        )

    def _embed_molecule(mol, use_random_coords=False):
        return EmbedMolecule(
            mol2,
            randomSeed=0xF00D,
            maxAttempts=100000,
            useRandomCoords=use_random_coords,
        )

    with TemporaryDirectory() as tempdir:
        raw_pubchem_sdf = os.path.join(tempdir, "raw.sdf")
        with open(raw_pubchem_sdf, "w") as f:
            f.write(data)

        suppl = list(Chem.SDMolSupplier(raw_pubchem_sdf))
        assert len(suppl) == 1
        mol = suppl[0]

        # Check if the molecule has multiple components. If so, embedding will fail to give a reasonable result.
        props = mol.GetPropsAsDict()
        if props["PUBCHEM_COMPONENT_COUNT"] > 1:
            raise RuntimeError(
                f'Cannot convert multi-component molecule: "{mol.GetProp("_Name")}"'
            )

        mol2 = Chem.AddHs(mol)

        # RDKit can fail to embed molecules, which would return a returncode of -1.
        # Try with random coords if it fails, which is supposed to help. See e.g.:
        # https://github.com/rdkit/rdkit/issues/2996#issuecomment-606464769
        returncode = _embed_molecule(mol2, use_random_coords=False)
        if returncode != 0:
            returncode = _embed_molecule(mol2, use_random_coords=True)

        if returncode != 0:
            raise RuntimeError(
                f'Could not embed molecule {mol2.GetProp("_Name")}'
            )

        return str(Chem.MolToMolBlock(mol2))
