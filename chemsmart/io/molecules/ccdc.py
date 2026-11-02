"""
Functions for querying CCDC (Cambridge Crystallographic Data Centre) for structures.

This module provides functionality to download CIF files from the CCDC database
using deposition numbers.
"""

import logging
import os
from functools import lru_cache
from tempfile import gettempdir

import requests
from requests.exceptions import HTTPError as RequestsHTTPError
from requests.exceptions import RequestException, Timeout
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

logger = logging.getLogger(__name__)

# CCDC base URL for structure retrieval
# NOTE: This is a placeholder URL. The actual CCDC API endpoint should be configured
# based on the specific CCDC web service being used. Users may need to update this
# or provide their own endpoint via environment variable or configuration file.
CCDC_BASE_URL = os.environ.get(
    "CCDC_BASE_URL", "https://www.ccdc.cam.ac.uk/structures"
)


def get_ccdc_cache_dir():
    """
    Get the cache directory for CCDC CIF files.

    Returns:
        str: Path to cache directory
    """
    cache_dir = os.path.join(gettempdir(), "chemsmart_ccdc_cache")
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


@lru_cache(maxsize=128)
@retry(
    stop=stop_after_attempt(3),  # Retry up to 3 times
    wait=wait_exponential(
        multiplier=1, min=2, max=10
    ),  # Wait 2, 4, then 8 seconds
    retry=retry_if_exception_type(Timeout),  # Retry on timeout
    before_sleep=lambda retry_state: logger.debug(
        f"Retrying CCDC download (attempt {retry_state.attempt_number}) "
        f"after {retry_state.idle_for}s..."
    ),
)
def fetch_cif_from_ccdc(deposition_number, cache=True, timeout=30):
    """
    Download a CIF file from CCDC database by deposition number.

    Args:
        deposition_number (str or int): CCDC deposition number (e.g., 1428476)
        cache (bool): Whether to cache the downloaded file. Default is True.
        timeout (int): Request timeout in seconds. Default is 30.

    Returns:
        str: Path to the downloaded CIF file

    Raises:
        ValueError: If the deposition number is invalid or structure not found
        requests.exceptions.Timeout: If the request times out
        requests.exceptions.RequestException: For other network issues

    Example:
        >>> cif_path = fetch_cif_from_ccdc(1428476)
        >>> mol = Molecule.from_cif_file(cif_path)
    """
    # Convert to string and validate
    dep_num = str(deposition_number).strip()
    if not dep_num.isdigit():
        raise ValueError(
            f"Invalid CCDC deposition number: {deposition_number}. "
            "Must be a positive integer."
        )

    # Check cache first
    cache_dir = get_ccdc_cache_dir()
    cached_file = os.path.join(cache_dir, f"{dep_num}.cif")

    if cache and os.path.exists(cached_file):
        logger.info(f"Using cached CIF file for CCDC {dep_num}")
        return cached_file

    # Construct download URL
    # Note: The actual CCDC API endpoint may vary. This is a placeholder.
    # The real implementation would need to use the actual CCDC API or scraping logic.
    url = f"{CCDC_BASE_URL}/Search?Ccdcid={dep_num}&DatabaseToSearch=Published"

    logger.info(f"Downloading CIF file for CCDC deposition {dep_num}")

    try:
        # Make the request
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()

        # Check if we got a valid CIF file
        content = response.text
        if not content or len(content.strip()) == 0:
            raise ValueError(
                f"Empty response received for deposition {dep_num}. "
                "The structure may not exist or may not be publicly available."
            )

        # Validate CIF format by attempting to parse it with ASE
        try:
            import tempfile

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".cif", delete=False
            ) as tmp:
                tmp.write(content)
                tmp_path = tmp.name

            # Try to read with ASE to validate format
            from ase.io import read as ase_read

            try:
                ase_read(tmp_path, format="cif")
            finally:
                os.unlink(tmp_path)

        except Exception as e:
            raise ValueError(
                f"Invalid CIF data received for deposition {dep_num}. "
                f"CIF validation failed: {str(e)}"
            ) from e

        # Save to cache
        if cache:
            with open(cached_file, "w") as f:
                f.write(content)
            logger.info(f"Cached CIF file for CCDC {dep_num} at {cached_file}")
            return cached_file
        else:
            # Save to temporary file
            import tempfile

            temp_file = tempfile.NamedTemporaryFile(
                mode="w", suffix=".cif", delete=False
            )
            temp_file.write(content)
            temp_file.close()
            return temp_file.name

    except RequestsHTTPError as e:
        if e.response.status_code == 404:
            raise ValueError(
                f"CCDC deposition {dep_num} not found. "
                "Please verify the deposition number is correct."
            ) from e
        elif e.response.status_code in [400, 403]:
            raise ValueError(
                f"Access denied or invalid request for CCDC deposition {dep_num}. "
                "The structure may require authentication or may not be publicly available."
            ) from e
        else:
            logger.error(
                f"HTTP error {e.response.status_code} fetching CCDC {dep_num}: {e}"
            )
            raise
    except Timeout as e:
        logger.warning(f"Timeout fetching CCDC deposition {dep_num}: {e}")
        raise
    except RequestException as e:
        logger.error(f"Network error fetching CCDC deposition {dep_num}: {e}")
        raise


def ccdc_search(deposition_number, fail_silently=False, **kwargs):
    """
    Search CCDC for structure and return a Molecule object.

    Args:
        deposition_number (str or int): CCDC deposition number
        fail_silently (bool): If True, return None on failure instead of raising.
            Default is False.
        **kwargs: Additional arguments passed to fetch_cif_from_ccdc

    Returns:
        Molecule or None: The molecule object or None if retrieval fails
            and fail_silently=True.

    Raises:
        ValueError: If structure cannot be retrieved and fail_silently=False.

    Example:
        >>> mol = ccdc_search(1428476)
        >>> if mol:
        ...     print(mol.chemical_formula)
    """
    from chemsmart.io.molecules.structure import Molecule

    try:
        cif_path = fetch_cif_from_ccdc(deposition_number, **kwargs)
        molecule = Molecule.from_cif_file(cif_path)

        if molecule is not None:
            logger.info(
                f"Successfully created molecule from CCDC deposition {deposition_number}"
            )
        return molecule

    except (ValueError, RequestException) as e:
        if fail_silently:
            logger.info(
                f"Error getting structure from CCDC {deposition_number}. "
                "Returning None as fail_silently=True"
            )
            return None
        raise ValueError(
            f"Failed to retrieve structure from CCDC deposition {deposition_number}: {e}"
        ) from e
