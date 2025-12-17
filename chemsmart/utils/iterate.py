"""Utility functions for iterate module."""

import logging
import os

logger = logging.getLogger(__name__)

ITERATE_TEMPLATE = """# Chemsmart Iterate Configuration Template
# =========================================
# This template defines skeletons and substituents for iterative structure generation.
#
# Structure Sources (mutually exclusive - use only one per entry):
#   - file_path: Path to molecular structure file (.xyz, .com, .log, .sdf, etc.)
#   - smiles: SMILES string for the molecule
#   - pubchem: PubChem identifier (CID, name, or SMILES)
#
# Required Fields:
#   - label: Unique identifier for this structure
#   - link_index: Atom index (1-based) where substitution occurs
#
# Optional Fields:
#   - skeleton_indices: Atom indices to keep, format: "1, 3-10, 15" (1-based)
#                       If not specified, all atoms except link_index are kept
#
# Note: Leave empty for null/None values. Strings don\'t need quotes.

# ==============================================================================
# SKELETONS
# ==============================================================================
# Define base molecular scaffolds that will be modified with substituents.
# Each skeleton is a list item (starts with "-")

skeletons:
  # Example 1: Skeleton from file
  - file_path: /path/to/skeleton1.xyz
    smiles:
    pubchem:
    label: benzene_core
    link_index: 6,7,8
    skeleton_indices: 1-5

  # Example 2: Skeleton from SMILES
  - file_path:
    smiles: c1ccccc1
    pubchem:
    label: phenyl_ring
    link_index: 1
    skeleton_indices:

  # Add more skeletons as needed...


# ==============================================================================
# SUBSTITUENTS
# ==============================================================================
# Define functional groups or fragments to attach to skeletons.
# Each substituent is a list item (starts with "-")

substituents:
  # Example 1: Substituent from file
  - file_path: /path/to/methyl.xyz
    smiles:
    pubchem:
    label: methyl
    link_index: 1

  # Example 2: Substituent from SMILES
  - file_path:
    smiles: C
    pubchem:
    label: CH3
    link_index: 1

  # Example 3: Substituent from PubChem
  - file_path:
    smiles:
    pubchem: methane
    label: methyl_pubchem
    link_index: 1

  # Add more substituents as needed...


# ==============================================================================
# USAGE NOTES
# ==============================================================================
# 1. Atom indices are 1-based (first atom is index 1)
# 2. link_index specifies the atom that will form the new bond
# 3. For skeletons, the atom at link_index is typically replaced/removed
# 4. For substituents, the atom at link_index bonds to the skeleton
# 5. skeleton_indices format: "1, 3-10, 15, 20-25" (ranges and individual indices)
# 6. Only one source (file_path, smiles, pubchem) should be non-null per entry
# 7. Use ~ for null/None values
"""


def generate_template(
    output_path: str = "iterate_template.yaml", overwrite: bool = False
) -> str:
    """
    Generate a template YAML configuration file for iterate jobs.

    Parameters
    ----------
    output_path : str
        Path to write the template file. Default is 'iterate_template.yaml'.
    overwrite : bool
        If True, overwrite existing file. Default is False.

    Returns
    -------
    str
        Path to the generated template file.

    Raises
    ------
    FileExistsError
        If file exists and overwrite is False.
    """
    # Add .yaml extension if not present
    if not output_path.endswith((".yaml", ".yml")):
        output_path = f"{output_path}.yaml"

    # Check if file exists
    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError(
            f"File '{output_path}' already exists. "
            f"Use --overwrite or specify a different filename."
        )

    # Write template
    with open(output_path, "w") as f:
        f.write(ITERATE_TEMPLATE)

    logger.info(f"Generated template: {output_path}")
    return output_path
