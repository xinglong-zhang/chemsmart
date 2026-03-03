import logging
import re
import xml.etree.ElementTree as ET
from collections import Counter

import numpy as np

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class SDFFile(FileMixin):
    """
    SDF file object.
    """

    def __init__(self, filename):
        self.filename = filename

    @property
    def molecule(self):
        return self.get_molecule()

    def get_molecule(self):
        list_of_symbols = []
        cart_coords = []
        # sdf line pattern containing coordinates and element type
        from chemsmart.utils.repattern import sdf_pattern

        for line in self.contents:
            match = re.match(sdf_pattern, line)
            if match:
                x = float(match.group(1))
                y = float(match.group(2))
                z = float(match.group(3))
                atom_type = str(match.group(4))
                list_of_symbols.append(atom_type)
                cart_coords.append((x, y, z))

        cart_coords = np.array(cart_coords)

        if len(list_of_symbols) == 0 or len(cart_coords) == 0:
            raise ValueError("No coordinates found in the SDF file!")

        return Molecule.from_symbols_and_positions_and_pbc_conditions(
            list_of_symbols=list_of_symbols, positions=cart_coords
        )


class CDXFile(FileMixin):
    """
    ChemDraw file object for reading .cdx and .cdxml files.

    Supports both binary (.cdx) and XML-based (.cdxml) ChemDraw formats.
    Uses RDKit for parsing and generates 3D coordinates using EmbedMolecule.

    Args:
        filename (str or pathlib.Path): Path to the ChemDraw file.
    """

    def __init__(self, filename):
        from pathlib import Path

        # Accept both str and pathlib.Path
        self.filename = (
            str(filename) if isinstance(filename, Path) else filename
        )

    @property
    def molecules(self):
        """
        Return all molecules from the ChemDraw file.
        """
        return self._parse_chemdraw_file()

    def get_molecules(self, index="-1", return_list=False):
        """Return molecule(s) for the requested 1-based ChemDraw index."""
        molecules = self.molecules

        if index == ":":
            selection = molecules
        else:
            parsed_index = (
                index
                if isinstance(index, slice)
                else string2index_1based(str(index))
            )
            if isinstance(parsed_index, slice):
                selection = molecules[parsed_index]
            else:
                selection = molecules[parsed_index]

        if return_list:
            return selection if isinstance(selection, list) else [selection]
        return selection

    def _parse_chemdraw_file(self):
        """
        Parse the ChemDraw file and return a list of Molecule objects.

        Uses RDKit to parse the file and generate 3D coordinates.
        Falls back to Open Babel (via obtain_mols_from_cdx_via_obabel) for .cdx
        if RDKit cannot read it.

        Returns:
            list[Molecule]: List of Molecule objects with 3D coordinates.

        Raises:
            ValueError: If no molecules can be read or no valid 3D structures
            can be generated from the ChemDraw file.
        """
        from pathlib import Path

        from rdkit import Chem
        from rdkit.Chem import AllChem

        suffix = Path(self.filename).suffix.lower()

        rdkit_mols = []
        try:
            # NOTE: RDKit's MolsFromCDXMLFile always supports CDXML.
            # CDX files are only supported if RDKit was built with
            # ChemDraw CDX support.
            rdkit_mols = list(
                Chem.MolsFromCDXMLFile(self.filename, removeHs=False)
            )
        except Exception as e:
            logger.debug(
                f"RDKit MolsFromCDXMLFile failed for {self.filename}: {e}"
            )

        # Fallback for .cdx: use Open Babel helper if RDKit gave nothing.
        if not rdkit_mols and suffix == ".cdx":
            logger.debug(
                f"RDKit did not return any molecules for {self.filename}; "
                "falling back to Open Babel CDX reader."
            )
            try:
                from chemsmart.utils.io import obtain_mols_from_cdx_via_obabel

                rdkit_mols = obtain_mols_from_cdx_via_obabel(self.filename)
            except Exception as e:
                logger.debug(
                    f"Open Babel CDX fallback failed for {self.filename}: {e}"
                )

        if not rdkit_mols:
            raise ValueError(
                f"No molecules could be read from ChemDraw file: {self.filename}"
            )

        molecules = []
        for rdkit_mol in rdkit_mols:
            if rdkit_mol is None:
                continue

            # Add explicit hydrogens for proper structure
            rdkit_mol = Chem.AddHs(rdkit_mol)

            # Generate 3D coordinates
            try:
                # Try to embed the molecule to get 3D coordinates
                result = AllChem.EmbedMolecule(rdkit_mol, randomSeed=42)
                if result == -1:
                    # Embedding failed, try with random coordinates
                    result = AllChem.EmbedMolecule(
                        rdkit_mol,
                        useRandomCoords=True,
                        randomSeed=42,
                    )
                    if result == -1:
                        logger.warning(
                            f"Could not generate 3D coordinates for a molecule "
                            f"in {self.filename}. Skipping this molecule."
                        )
                        continue

                # Optimize the geometry (may fail for exotic atom types)
                try:
                    AllChem.MMFFOptimizeMolecule(rdkit_mol)
                except Exception as e:
                    logger.debug(
                        f"MMFF optimization failed for a molecule in {self.filename}: {e}"
                    )

            except Exception as e:
                logger.warning(
                    f"Error generating 3D coordinates for molecule in {self.filename}: {e}"
                )
                continue

            # Convert RDKit Mol to Molecule
            mol = Molecule.from_rdkit_mol(rdkit_mol)
            molecules.append(mol)

        if not molecules:
            raise ValueError(
                "No valid molecules with 3D coordinates could be generated "
                f"from ChemDraw file: {self.filename}"
            )

        return molecules

    # ------------------------------------------------------------------
    # CDXML atom-colour helpers
    # ------------------------------------------------------------------


class PKaCDXFile(CDXFile):
    """Specialized CDXFile subclass for pKa calculations that need to identify"""

    def parse_cdxml_element_colors(self):
        """Parse the CDXML file and return per-atom colour information.

        Walks every ``<fragment>`` element in document order and collects
        ``<n>`` (node / atom) elements.  Atoms without an ``Element``
        attribute are implicit carbons and are still tracked for correct
        index mapping.

        For each atom the method also inspects the ``<t>``/``<s>`` text
        spans that make up the displayed label.  When a label contains
        an "H" character rendered in a **different** colour from the
        heavy-atom symbol (e.g. ``<s color="0">O</s><s color="4">H</s>``
        for a phenol –OH), the colour of that "H" span is recorded in
        ``implicit_h_color``.  This enables detection of functional-group
        protons (–OH, –NH, …) that are *not* represented as separate
        ``<n>`` nodes but instead as implicit hydrogens on the heavy atom.

        Returns:
            list[dict]: One entry per CDXML atom, in document order, with
                keys:

                * ``cdxml_id`` – the ``id`` attribute of the ``<n>`` element
                * ``element`` – integer atomic number (6 for C when omitted)
                * ``color`` – integer colour-table index (``0`` when absent)
                * ``symbol`` – element symbol string (``"C"``, ``"H"``, …)
                * ``num_hydrogens`` – number of implicit hydrogens attached
                  (from the ``NumHydrogens`` CDXML attribute, or ``None``)
                * ``implicit_h_color`` – colour-table index of an "H" span
                  in the label that differs from the atom's own colour, or
                  ``None`` when no such span exists

        Raises:
            ValueError: If the file cannot be parsed as valid CDXML.

        Note:
            The atom order matches the order RDKit uses when reading the
            same CDXML through ``Chem.MolsFromCDXMLFile``, **before**
            ``Chem.AddHs`` adds implicit hydrogens.
        """
        from rdkit.Chem import GetPeriodicTable

        pt = GetPeriodicTable()

        try:
            tree = ET.parse(self.filename)
        except ET.ParseError as exc:
            raise ValueError(
                f"Failed to parse CDXML file {self.filename}: {exc}"
            ) from exc

        root = tree.getroot()
        list_of_elements = []

        for fragment in root.iter("fragment"):
            # Skip nested fragments (Nicknames / abbreviations) that are
            # children of another <n> element.
            parent = self._find_parent(root, fragment)
            if parent is not None and parent.tag == "n":
                continue

            for node in fragment.findall("n"):
                # Skip ExternalConnectionPoint nodes
                if node.get("NodeType") == "ExternalConnectionPoint":
                    continue

                cdxml_id = node.get("id")
                element_num = int(node.get("Element", "6"))  # default C
                color = int(node.get("color", "0"))
                num_h_attr = node.get("NumHydrogens")
                num_hydrogens = (
                    int(num_h_attr) if num_h_attr is not None else None
                )

                # ---- inspect <t>/<s> text spans ----
                # Collect (text, color) pairs from all <s> children.
                spans = []  # list of (text_content, int_color)
                for t_elem in node.iter("t"):
                    for s_elem in t_elem.iter("s"):
                        s_text = (s_elem.text or "").strip()
                        s_color = int(s_elem.get("color", "0"))
                        if s_text:
                            spans.append((s_text, s_color))

                # Determine the atom's own display colour.
                # Priority: node-level color attribute > first <s> span colour.
                if color == 0 and spans:
                    # Use the colour of the first non-"H" span (the heavy-atom
                    # symbol) as the atom colour, falling back to the first span.
                    for s_text, s_color in spans:
                        if s_text.upper() not in ("H",):
                            color = s_color
                            break
                    else:
                        color = spans[0][1]

                # Detect a functional-group "H" rendered in a *different*
                # colour from the atom's own colour.  For example the phenol
                # label ``<s color="0">O</s><s color="4">H</s>`` yields
                # implicit_h_color = 4.
                implicit_h_color = None
                for s_text, s_color in spans:
                    if "H" in s_text and s_color != color:
                        implicit_h_color = s_color
                        break

                try:
                    symbol = pt.GetElementSymbol(element_num)
                except Exception:
                    symbol = "?"

                list_of_elements.append(
                    {
                        "cdxml_id": cdxml_id,
                        "element": element_num,
                        "color": color,
                        "symbol": symbol,
                        "num_hydrogens": num_hydrogens,
                        "implicit_h_color": implicit_h_color,
                    }
                )

        return list_of_elements

    @staticmethod
    def _find_parent(root, target):
        """Return the parent element of *target* inside *root*, or None."""
        for parent in root.iter():
            for child in parent:
                if child is target:
                    return parent
        return None

    def get_colored_proton_index(self, color_code=None):
        """Identify the 1-based atom index of a proton marked by colour.

        Supports two kinds of coloured protons:

        1. **Explicit H nodes** – a standalone ``<n Element="1">`` atom
           drawn in a different colour (e.g. a separate –H fragment).
        2. **Functional-group H** – a heavy-atom node whose label contains
           an "H" character in a *different* colour from the atom symbol
           (e.g. ``<s color="0">O</s><s color="4">H</s>`` for phenol –OH).
           The hydrogen is implicit (``NumHydrogens >= 1``) and only
           becomes an explicit atom after RDKit's ``Chem.AddHs``.

        Two detection modes are available:

        **Default mode** (``color_code=None``):
            Auto-detect: find the atom (explicit H *or* functional-group H)
            whose colour differs from the majority.

        **User-specified mode** (``color_code=<int>``):
            Find atoms / functional-group H spans matching that colour.

        Args:
            color_code (int | None): Colour-table index to look for.
                ``None`` triggers automatic detection.

        Returns:
            int: 1-based index of the hydrogen atom in the Molecule
                produced by ``CDXFile.get_molecules`` (i.e. after
                ``Chem.AddHs`` and 3-D embedding).

        Raises:
            ValueError: When the proton cannot be identified unambiguously.
        """
        atoms = self.parse_cdxml_element_colors()

        if not atoms:
            raise ValueError(f"No atoms found in CDXML file: {self.filename}")

        if color_code is not None:
            return self._proton_index_by_color(atoms, color_code)
        return self._proton_index_auto(atoms)

    # ----- private helpers for get_colored_proton_index -----

    def _proton_index_by_color(self, atoms, color_code):
        """Return the 1-based proton index for a user-specified colour."""
        # Case 1: explicit H node with that colour
        explicit_h = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["color"] == color_code and a["symbol"] == "H"
        ]
        if len(explicit_h) == 1:
            idx, atom = explicit_h[0]
            proton_index = idx + 1
            logger.info(
                f"User-specified colour: explicit H at 1-based index "
                f"{proton_index} (CDXML id={atom['cdxml_id']}, "
                f"color={color_code})."
            )
            return proton_index
        if len(explicit_h) > 1:
            raise ValueError(
                f"Multiple explicit H atoms ({len(explicit_h)}) have "
                f"color code {color_code}."
            )

        # Case 2: functional-group implicit H with that colour
        fg_matches = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["implicit_h_color"] == color_code
        ]
        if len(fg_matches) == 0:
            # Also check if any non-H atoms have that node colour (give
            # a clear message distinguishing "no match" from "not H")
            any_match = [a for a in atoms if a["color"] == color_code]
            if any_match:
                raise ValueError(
                    f"Atom(s) with color code {color_code} found, but "
                    f"none are hydrogen. Symbols: "
                    f"{[a['symbol'] for a in any_match]}."
                )
            raise ValueError(
                f"No atoms with color code {color_code} found in "
                f"{self.filename}."
            )
        if len(fg_matches) > 1:
            raise ValueError(
                f"Multiple functional-group H spans ({len(fg_matches)}) "
                f"have color code {color_code}."
            )

        cdxml_idx, atom = fg_matches[0]
        return self._resolve_implicit_h_index(atoms, cdxml_idx, atom)

    def _proton_index_auto(self, atoms):
        """Auto-detect the uniquely coloured proton."""
        # Gather *all* colours that appear – both node colours and
        # implicit-H span colours.
        all_colors = []
        for a in atoms:
            all_colors.append(a["color"])
            if a["implicit_h_color"] is not None:
                all_colors.append(a["implicit_h_color"])

        color_counts = Counter(all_colors)
        if len(color_counts) < 2:
            raise ValueError(
                "All atoms in the CDXML file share the same colour. "
                "Cannot auto-detect the proton to remove. "
                "Use -cl/--color-code to specify the colour explicitly, "
                "or use -pi/--proton-index."
            )

        majority_color = color_counts.most_common(1)[0][0]

        # Candidates: explicit H with a non-majority colour
        explicit_h = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["symbol"] == "H" and a["color"] != majority_color
        ]

        # Candidates: functional-group implicit H whose span colour
        # differs from the majority
        fg_h = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["implicit_h_color"] is not None
            and a["implicit_h_color"] != majority_color
        ]

        total_candidates = len(explicit_h) + len(fg_h)

        if total_candidates == 0:
            # Collect all uniquely coloured atoms for the error message
            unique_atoms = [a for a in atoms if a["color"] != majority_color]
            if unique_atoms:
                symbols = [a["symbol"] for a in unique_atoms]
                raise ValueError(
                    f"Uniquely coloured atom(s) found but none are "
                    f"hydrogen (found: {symbols}). Only hydrogen atoms "
                    f"can be removed for pKa calculations."
                )
            raise ValueError(
                "No uniquely coloured atom found in the CDXML file."
            )

        if total_candidates > 1:
            raise ValueError(
                f"Multiple uniquely coloured hydrogen atoms found "
                f"({total_candidates}). Cannot determine which proton "
                f"to remove. Use -cl/--color-code or -pi/--proton-index "
                f"to specify explicitly."
            )

        # Exactly one candidate
        if explicit_h:
            idx, atom = explicit_h[0]
            proton_index = idx + 1
            logger.info(
                f"Auto-detected explicit H at 1-based index "
                f"{proton_index} (CDXML id={atom['cdxml_id']}, "
                f"color={atom['color']})."
            )
            return proton_index

        cdxml_idx, atom = fg_h[0]
        return self._resolve_implicit_h_index(atoms, cdxml_idx, atom)

    def _resolve_implicit_h_index(self, atoms, cdxml_idx, atom):
        """Map a functional-group implicit H to a 1-based Molecule index.

        After ``Chem.AddHs`` the implicit hydrogens are appended to the
        RDKit molecule.  This helper parses the CDXML with RDKit (without
        removing Hs), calls ``AddHs``, and finds the H atom bonded to the
        heavy atom at *cdxml_idx*.

        Args:
            atoms: Full atom list from ``parse_cdxml_atom_colors``.
            cdxml_idx: 0-based position of the heavy atom in *atoms*.
            atom: The dict entry for that heavy atom.

        Returns:
            int: 1-based index in the final Molecule.
        """
        from rdkit import Chem

        num_h = atom.get("num_hydrogens")
        if num_h is None or num_h < 1:
            raise ValueError(
                f"Atom at CDXML position {cdxml_idx + 1} "
                f"(id={atom['cdxml_id']}, {atom['symbol']}) has a "
                f"coloured H in its label but NumHydrogens is "
                f"{num_h!r}. Cannot identify an implicit hydrogen."
            )

        # Read the molecule with RDKit (same pipeline as _parse_chemdraw_file)
        rdkit_mols = list(
            Chem.MolsFromCDXMLFile(self.filename, removeHs=False)
        )
        if not rdkit_mols or rdkit_mols[0] is None:
            raise ValueError(
                f"RDKit could not read {self.filename} for implicit-H "
                f"resolution."
            )
        rdkit_mol = rdkit_mols[0]
        rdkit_mol = Chem.AddHs(rdkit_mol)

        # The heavy atom keeps its original 0-based index after AddHs.
        heavy_atom = rdkit_mol.GetAtomWithIdx(cdxml_idx)
        h_indices = []
        for bond in heavy_atom.GetBonds():
            other = bond.GetOtherAtom(heavy_atom)
            if other.GetSymbol() == "H":
                h_indices.append(other.GetIdx())

        if not h_indices:
            raise ValueError(
                f"No hydrogen bonded to {atom['symbol']} "
                f"(CDXML id={atom['cdxml_id']}) after AddHs."
            )

        # Use the *first* H bonded to this heavy atom (for OH / NH there
        # is typically exactly one).
        proton_0based = h_indices[0]
        proton_index = proton_0based + 1  # 1-based

        logger.info(
            f"Functional-group H detected on {atom['symbol']} "
            f"(CDXML id={atom['cdxml_id']}): 1-based Molecule index "
            f"{proton_index} (implicit_h_color={atom['implicit_h_color']})."
        )
        return proton_index

    # ------------------------------------------------------------------
    # PKaMolecule factory
    # ------------------------------------------------------------------

    def get_pka_molecule(
        self,
        proton_index=None,
        color_code=None,
        index="-1",
    ):
        """Return a :class:`PKaMolecule` with the acidic proton resolved.

        The proton is identified in priority order:

        1. **User-supplied atom index** – ``proton_index`` (1-based).
        2. **User-supplied colour code** – ``color_code``; the matching
           coloured hydrogen is looked up via
           :meth:`get_colored_proton_index`.
        3. **Auto-detection** – the uniquely coloured hydrogen is found
           automatically from the CDXML colour table.

        Args:
            proton_index (int | None): Explicit 1-based proton index.
                Skips all colour-based detection when provided.
            color_code (int | None): CDXML colour-table index.  Used
                only when *proton_index* is ``None``.
            index (str): Which molecule to return when the file
                contains several fragments (``"-1"`` = last, default).

        Returns:
            PKaMolecule: Molecule object with ``proton_index`` attribute.

        Raises:
            ValueError: If the proton cannot be identified or validated.
        """
        from chemsmart.io.molecules.structure import PKaMolecule

        mol = self.get_molecules(index=index)

        if proton_index is None:
            proton_index = self.get_colored_proton_index(
                color_code=color_code,
            )

        return PKaMolecule(molecule=mol, proton_index=proton_index)

    def get_pka_molecules(
        self,
        proton_index=None,
        color_code=None,
    ):
        """Return all molecules as :class:`PKaMolecule` instances.

        Each molecule gets the same ``proton_index`` /
        ``color_code`` resolution applied.  Useful for multi-fragment
        CDXML files where every fragment shares the same colouring
        convention.

        Args:
            proton_index (int | None): Explicit 1-based proton index.
            color_code (int | None): CDXML colour-table index.

        Returns:
            list[PKaMolecule]: One ``PKaMolecule`` per fragment.
        """
        from chemsmart.io.molecules.structure import PKaMolecule

        if proton_index is None:
            proton_index = self.get_colored_proton_index(
                color_code=color_code,
            )

        molecules = self.molecules
        return [
            PKaMolecule(molecule=mol, proton_index=proton_index)
            for mol in molecules
        ]
