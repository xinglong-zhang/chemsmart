"""
Functions for retrieving structures from the Cambridge Crystallographic
Data Centre (CCDC) by deposition number.

Two retrieval routes are supported, tried in order:

1. Local CSD Python API (preferred): if the licensed ``ccdc`` package
   (CSD Python API) is installed and activated, the structure is
   retrieved directly from the local CSD installation. In this case
   ``ccdc_search()`` returns the deposited *discrete molecule*
   (e.g., one benzene molecule for CCDC 725244), which is what is
   normally needed as input for molecular quantum chemistry jobs.
2. HTTP download (fallback): otherwise the CIF file is downloaded from
   the CCDC web service (base URL configurable via the ``CCDC_BASE_URL``
   environment variable) and parsed with ASE. Note that this route
   returns the crystallographic structure read from the CIF (i.e., the
   expanded unit cell), not a discrete molecule.
"""

import logging
import os
import tempfile
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

# CCDC base URL for structure retrieval via HTTP (fallback route).
# NOTE: This is a placeholder URL. The actual CCDC web endpoint should be
# configured based on the specific CCDC web service being used, via the
# CCDC_BASE_URL environment variable. Users with a CSD licence are
# encouraged to use the local CSD Python API instead (preferred route).
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


def validate_deposition_number(deposition_number):
    """
    Validate and normalize a CCDC deposition number.

    Args:
        deposition_number (str or int): CCDC deposition number
            (e.g., 1428476).

    Returns:
        str: Normalized deposition number string.

    Raises:
        ValueError: If the deposition number is not a positive integer.
    """
    dep_num = str(deposition_number).strip()
    if not dep_num.isdigit() or int(dep_num) <= 0:
        raise ValueError(
            f"Invalid CCDC deposition number: {deposition_number}. "
            "Must be a positive integer."
        )
    return dep_num


def csd_api_available():
    """
    Check whether the local CSD Python API is installed and licensed.

    Returns:
        bool: True if the ``ccdc`` package can be imported and its
            licence check passes, False otherwise.
    """
    try:
        import ccdc  # noqa: F401
    except (ImportError, RuntimeError) as e:
        # ImportError: csd-python-api not installed
        # RuntimeError: CSD licence not activated
        logger.debug(f"Local CSD Python API unavailable: {e}")
        return False
    return True


def _search_csd_entries(deposition_number):
    """
    Search the local CSD for entries matching a deposition number.

    Args:
        deposition_number (str or int): CCDC deposition number.

    Returns:
        list: List of ``ccdc.entry.Entry`` objects matching the
            deposition number.

    Raises:
        ValueError: If no structure in the CSD matches the deposition
            number.
    """
    from ccdc.search import TextNumericSearch

    dep_num = validate_deposition_number(deposition_number)
    searcher = TextNumericSearch()
    searcher.add_ccdc_number(int(dep_num))
    hits = searcher.search()

    if not hits:
        raise ValueError(
            f"No structure found in the CSD for CCDC deposition "
            f"number {dep_num}. Please verify the deposition number "
            "is correct."
        )

    if len(hits) > 1:
        logger.info(
            f"CCDC deposition {dep_num} matched {len(hits)} CSD entries; "
            f"using the first entry: {hits[0].entry.identifier}"
        )
    return [hit.entry for hit in hits]


def _molecule_from_csd_entry(entry):
    """
    Build a Molecule object from the discrete molecule of a CSD entry.

    Unlike parsing the exported CIF (which yields the crystallographic
    unit cell), this uses ``entry.molecule`` and therefore returns the
    deposited discrete molecule with its full set of atoms.

    Args:
        entry (ccdc.entry.Entry): A CSD entry.

    Returns:
        Molecule: Molecule object with symbols and 3D positions.

    Raises:
        ValueError: If the entry has no 3D coordinates.
    """
    from ase import Atoms

    from chemsmart.io.molecules.atoms import AtomsChargeMultiplicity

    csd_molecule = entry.molecule
    symbols = [atom.atomic_symbol for atom in csd_molecule.atoms]

    positions = []
    for atom in csd_molecule.atoms:
        if atom.coordinates is None:
            raise ValueError(
                f"CSD entry {entry.identifier} has no 3D coordinates "
                f"for atom {atom.label}; cannot build a 3D structure."
            )
        positions.append(tuple(atom.coordinates))

    logger.debug(
        f"Building molecule from CSD entry {entry.identifier} "
        f"({len(symbols)} atoms)"
    )
    ase_atoms = Atoms(symbols=symbols, positions=positions)
    return AtomsChargeMultiplicity.from_atoms(ase_atoms).to_molecule()


def _export_cif_via_csd_api(deposition_number, cached_file):
    """
    Export the crystal structure of a CSD entry to a CIF file using the
    local CSD Python API.

    Args:
        deposition_number (str): Validated CCDC deposition number.
        cached_file (str): Destination path for the exported CIF file.

    Returns:
        str: Path to the exported CIF file.
    """
    from ccdc.io import CrystalWriter

    entries = _search_csd_entries(deposition_number)
    entry = entries[0]

    with CrystalWriter(cached_file) as writer:
        writer.write(entry.crystal)

    logger.info(
        f"Exported CIF file for CCDC {deposition_number} "
        f"({entry.identifier}) to {cached_file}"
    )
    return cached_file


def _validate_cif_content(content, deposition_number):
    """
    Validate that downloaded content is parseable as a CIF file by ASE.

    Args:
        content (str): Text content to validate.
        deposition_number (str): CCDC deposition number (for error
            messages).

    Raises:
        ValueError: If the content cannot be parsed as CIF.
    """
    from ase.io import read as ase_read

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".cif", delete=False
    ) as tmp:
        tmp.write(content)
        tmp_path = tmp.name

    try:
        ase_read(tmp_path, format="cif")
    except Exception as e:
        raise ValueError(
            f"Invalid CIF data received for deposition "
            f"{deposition_number}. CIF validation failed: {str(e)}"
        ) from e
    finally:
        os.unlink(tmp_path)


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
def _download_cif_via_http(deposition_number, cached_file, cache, timeout):
    """
    Download a CIF file from the CCDC web service via HTTP.

    Args:
        deposition_number (str): Validated CCDC deposition number.
        cached_file (str): Destination path for the cached CIF file.
        cache (bool): Whether to save the file to the cache directory;
            if False, a temporary file is used instead.
        timeout (int): Request timeout in seconds.

    Returns:
        str: Path to the downloaded CIF file.

    Raises:
        ValueError: If the structure is not found, access is denied, or
            the downloaded content is not a valid CIF file.
        requests.exceptions.Timeout: If the request times out
            (after retries).
        requests.exceptions.RequestException: For other network issues.
    """
    url = (
        f"{CCDC_BASE_URL}/Search?Ccdcid={deposition_number}"
        "&DatabaseToSearch=Published"
    )

    logger.info(
        f"Downloading CIF file for CCDC deposition {deposition_number}"
    )

    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()

        content = response.text
        if not content or len(content.strip()) == 0:
            raise ValueError(
                f"Empty response received for deposition "
                f"{deposition_number}. The structure may not exist or "
                "may not be publicly available."
            )

        # Validate CIF format by attempting to parse it with ASE
        _validate_cif_content(content, deposition_number)

        if cache:
            with open(cached_file, "w") as f:
                f.write(content)
            logger.info(
                f"Cached CIF file for CCDC {deposition_number} "
                f"at {cached_file}"
            )
            return cached_file

        # Save to temporary file
        temp_file = tempfile.NamedTemporaryFile(
            mode="w", suffix=".cif", delete=False
        )
        temp_file.write(content)
        temp_file.close()
        return temp_file.name

    except RequestsHTTPError as e:
        if e.response.status_code == 404:
            raise ValueError(
                f"CCDC deposition {deposition_number} not found. "
                "Please verify the deposition number is correct."
            ) from e
        elif e.response.status_code in [400, 403]:
            raise ValueError(
                f"Access denied or invalid request for CCDC deposition "
                f"{deposition_number}. The structure may require "
                "authentication or may not be publicly available."
            ) from e
        else:
            logger.error(
                f"HTTP error {e.response.status_code} fetching CCDC "
                f"{deposition_number}: {e}"
            )
            raise
    except Timeout as e:
        logger.warning(
            f"Timeout fetching CCDC deposition {deposition_number}: {e}"
        )
        raise
    except RequestException as e:
        logger.error(
            f"Network error fetching CCDC deposition "
            f"{deposition_number}: {e}"
        )
        raise


def fetch_cif_from_ccdc(deposition_number, cache=True, timeout=30):
    """
    Download a CIF file from CCDC database by deposition number.

    Uses the local CSD Python API if it is installed and licensed;
    otherwise falls back to downloading from the CCDC web service.
    Downloaded files are cached on disk so that repeated lookups do not
    repeat the query.

    Args:
        deposition_number (str or int): CCDC deposition number
            (e.g., 1428476).
        cache (bool): Whether to cache the downloaded file.
            Default is True.
        timeout (int): Request timeout in seconds (HTTP route only).
            Default is 30.

    Returns:
        str: Path to the downloaded CIF file.

    Raises:
        ValueError: If the deposition number is invalid, the structure
            is not found, or the retrieved data is not a valid CIF file.
        requests.exceptions.Timeout: If the request times out.
        requests.exceptions.RequestException: For other network issues.

    Example:
        >>> cif_path = fetch_cif_from_ccdc(1428476)
        >>> mol = Molecule.from_cif_file(cif_path)
    """
    dep_num = validate_deposition_number(deposition_number)

    # Check cache first
    cached_file = os.path.join(get_ccdc_cache_dir(), f"{dep_num}.cif")
    if cache and os.path.exists(cached_file):
        logger.info(f"Using cached CIF file for CCDC {dep_num}")
        return cached_file

    # Preferred route: local licensed CSD Python API
    if csd_api_available():
        return _export_cif_via_csd_api(dep_num, cached_file)

    # Fallback route: HTTP download from the CCDC web service
    logger.debug(
        "Local CSD Python API not available; falling back to HTTP "
        "download for CCDC structure retrieval."
    )
    return _download_cif_via_http(
        dep_num, cached_file, cache=cache, timeout=timeout
    )


def ccdc_search(deposition_number, fail_silently=False, **kwargs):
    """
    Search CCDC for structure and return a Molecule object.

    If the local CSD Python API is available, the deposited discrete
    molecule is returned (preferred). Otherwise the CIF file is
    downloaded via HTTP and parsed with ASE, which returns the
    crystallographic structure (unit cell) contained in the CIF.

    Args:
        deposition_number (str or int): CCDC deposition number.
        fail_silently (bool): If True, return None on failure instead of
            raising. Default is False.
        **kwargs: Additional arguments passed to fetch_cif_from_ccdc
            (HTTP route only).

    Returns:
        Molecule or None: The molecule object or None if retrieval fails
            and fail_silently=True.

    Raises:
        ValueError: If structure cannot be retrieved and
            fail_silently=False.

    Example:
        >>> mol = ccdc_search(725244)  # benzene
        >>> if mol:
        ...     print(mol.chemical_formula)  # 'C6H6'
    """
    from chemsmart.io.molecules.structure import Molecule

    try:
        if csd_api_available():
            # Preferred route: discrete molecule from local CSD
            entries = _search_csd_entries(deposition_number)
            molecule = _molecule_from_csd_entry(entries[0])
        else:
            # Fallback route: CIF download + ASE parsing
            logger.debug(
                "Local CSD Python API not available; retrieving CIF "
                "from CCDC web service. Note that the structure read "
                "from the CIF is the crystallographic unit cell."
            )
            cif_path = fetch_cif_from_ccdc(deposition_number, **kwargs)
            molecule = Molecule.from_cif_file(cif_path)

        if molecule is not None:
            logger.info(
                f"Successfully created molecule from CCDC deposition "
                f"{deposition_number}"
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
            f"Failed to retrieve structure from CCDC deposition "
            f"{deposition_number}: {e}"
        ) from e
