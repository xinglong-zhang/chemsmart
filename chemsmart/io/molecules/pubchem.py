"""Functions for querying PubChem for structures, adapted from ASE's pubchem.py module."""

import logging
import os
import warnings
from contextlib import suppress
from functools import lru_cache
from io import StringIO
from tempfile import TemporaryDirectory

import requests
from ase.data.pubchem import analyze_input, base_url
from requests.exceptions import HTTPError as RequestsHTTPError
from requests.exceptions import RequestException, Timeout
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from chemsmart.io.molecules.structure import Molecule

try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem.AllChem import EmbedMolecule
except ImportError:
    rdkit = None

logger = logging.getLogger(__name__)


def search_pubchem_raw(search, field, suffix: str = "3d", timeout: int = 10):
    """Search PubChem for structure.

    Changelog:
    1. Added suffix attribute.
    2. Added silent flag to print statements in except blocks.
    3. Changed default len of conformer_id to 1.

    Args:
        search (str): The search term (e.g., CID, SMILES, name).
        field (str): The field to search (e.g., 'cid', 'smiles', 'name', 'conformers').
        suffix (str): The suffix for the request (e.g., '3d', '2d'). Default is '3d'.
        timeout (int): Request timeout in seconds. Default is 10.

    Returns:
        str: The raw SDF or JSON response decoded as UTF-8.

    Raises:
        requests.exceptions.Timeout: If the request times out.
        requests.exceptions.HTTPError: For HTTP-related errors (e.g., 404, 400).
        requests.exceptions.RequestException: For other network issues.
    """
    suffix = "sdf?record_type=" + suffix

    if field == "conformers":
        # We don't use the "compound" flag when looking for conformers
        url = f"{base_url}/{field}/{search}/{suffix}"
    else:
        url = f"{base_url}/compound/{field}/{search}/{suffix}"

    # Use requests with timeout
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()  # Raises HTTPError for 400, 404, etc.

    # Check if there are conformers and warn them if there are
    if field != "conformers":
        # Set default conformer_ids to len==1, as there are PubChem structures with no conformer information
        # which would return "HTTPError:PUGREST.NotFound"
        conformer_ids = ["conformer information not found."]

        with suppress(RequestsHTTPError, RequestException):
            # Test if there is any conformer information of the structure
            conformer_json = search_pubchem_raw(
                search, field, suffix="conformers/JSON", timeout=timeout
            )
            # Parse JSON to check for conformers (assuming a list or dict response)
            import json

            conformer_data = json.loads(conformer_json)
            conformer_ids = conformer_data.get("Conformers", conformer_ids)

        if len(conformer_ids) > 1:
            warnings.warn(
                f'The structure "{search}" has more than one '
                "conformer in PubChem. By default, the "
                "first conformer is returned, please ensure "
                "you are using the structure you intend to "
                "or use the "
                "`ase.data.pubchem.pubchem_conformer_search` "
                "function",
                stacklevel=2,
            )

    return response.text  # Already UTF-8 decoded by requests


@lru_cache(maxsize=128)
@retry(
    stop=stop_after_attempt(3),  # Retry up to 3 times
    wait=wait_exponential(
        multiplier=1, min=2, max=10
    ),  # Wait 2, 4, then 8 seconds
    retry=retry_if_exception_type(Timeout),  # Retry on timeout
    before_sleep=lambda retry_state: logger.debug(
        f"Retrying PubChem search (attempt {retry_state.attempt_number}) after {retry_state.idle_for}s..."
    ),
)
def pubchem_search(*args, fail_silently=True, **kwargs):
    """Search PubChem for structure.

    Changelog:
    1. Try suffix =='3d' and '2d'.
    2. Reads 2d sdf files by rdkit to create more accurate structures.
    3. Simplify parse_pubchem_raw() to read atom structure.

    Args:
        *args: Variable positional arguments for search.
        fail_silently (bool): If True, return None on failure instead of raising an error.
        **kwargs: Keyword arguments for search (e.g., cid, smiles, name).

    Returns:
        Molecule or None: The molecule object or None if retrieval fails and fail_silently=True.

    Raises:
        ValueError: If structure cannot be retrieved and fail_silently=False.
    """

    def _pubchem_search():
        search, field = analyze_input(*args, **kwargs)

        try:
            raw_pubchem = search_pubchem_raw(
                search, field, suffix="3d", timeout=10
            )
        except RequestsHTTPError as e:
            error_substrings = ["400", "404"]  # Match status codes as strings
            if not any(substring in str(e) for substring in error_substrings):
                logger.warning(
                    f"Unexpected HTTP error fetching 3D {field} for {search}: {str(e)}"
                )
                raise  # Retry on unexpected HTTP errors via outer handler
            logger.info(
                f"Error getting 3D {field} structure from PubChem for {search}, trying 2D"
            )

            try:
                raw_pubchem = search_pubchem_raw(
                    search, field, suffix="2d", timeout=10
                )
                raw_pubchem = _pubchem_2d_to_3d(raw_pubchem)
                logger.info(
                    f"Successfully got 2D {field} structure for {search}"
                )
            except RequestsHTTPError as e:
                logger.error(
                    f"Failed to get 2D {field} structure for {search}: {str(e)}"
                )
                raise ValueError(
                    "Error getting structure from PubChem!"
                ) from e

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
                "Error getting structure from PubChem. Returning None as fail_silently=True"
            )
            return None
        raise e
    except Timeout as e:
        logger.warning(
            f"Timeout fetching structure for {kwargs.get('search', args)}: {str(e)}"
        )
        raise  # Trigger retry
    except RequestException as e:
        logger.warning(
            f"Network error fetching structure for {kwargs.get('search', args)}: {str(e)}"
        )
        raise  # Trigger retry


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
