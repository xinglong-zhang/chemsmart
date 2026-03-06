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
    """Specialized ``CDXFile`` subclass used in pKa workflows.

    This class adds utilities on top of :class:`CDXFile` for:

    * parsing ChemDraw CDXML atom colours and element labels,
    * distinguishing explicitly drawn hydrogens from implicit/functional-group
      protons based on label structure and colour differences, and
    * providing structured atom and proton information that can be consumed
      by higher-level pKa logic (for example, to build :class:`PKaMolecule`
      instances or similar pKa-aware molecule objects).

    The class does **not** perform the pKa calculation itself; instead it
    focuses on converting CDXML drawings into a representation where
    ionizable sites and their associated protons can be identified reliably.
    """

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
        root = self._parse_cdxml_root()

        list_of_elements = []
        for fragment in self._iter_top_level_fragments(root):
            list_of_elements.extend(self._parse_fragment_nodes(fragment, pt))

        return list_of_elements

    @staticmethod
    def _find_parent(root, target):
        """Return the parent element of *target* inside *root*, or None."""
        for parent in root.iter():
            for child in parent:
                if child is target:
                    return parent
        return None

    def _parse_cdxml_root(self):
        """Parse CDXML and return the XML root element."""
        try:
            tree = ET.parse(self.filename)
        except ET.ParseError as exc:
            raise ValueError(
                f"Failed to parse CDXML file {self.filename}: {exc}"
            ) from exc
        return tree.getroot()

    def _iter_top_level_fragments(self, root):
        """Yield top-level fragment elements, excluding nested node fragments."""
        for fragment in root.iter("fragment"):
            parent = self._find_parent(root, fragment)
            if parent is not None and parent.tag == "n":
                continue
            yield fragment

    def _parse_fragment_nodes(self, fragment, periodic_table):
        """Convert all valid CDXML nodes in one fragment into atom dicts."""
        fragment_atoms = []
        for node in fragment.findall("n"):
            if node.get("NodeType") == "ExternalConnectionPoint":
                continue

            cdxml_id = node.get("id")
            element_num = int(node.get("Element", "6"))
            color = int(node.get("color", "0"))
            num_h_attr = node.get("NumHydrogens")
            num_hydrogens = int(num_h_attr) if num_h_attr is not None else None

            spans = []
            for t_elem in node.iter("t"):
                for s_elem in t_elem.iter("s"):
                    s_text = (s_elem.text or "").strip()
                    s_color = int(s_elem.get("color", "0"))
                    if s_text:
                        spans.append((s_text, s_color))

            if color == 0 and spans:
                for s_text, s_color in spans:
                    if s_text.upper() not in ("H",):
                        color = s_color
                        break
                else:
                    color = spans[0][1]

            implicit_h_color = None
            for s_text, s_color in spans:
                if "H" in s_text and s_color != color:
                    implicit_h_color = s_color
                    break

            try:
                symbol = periodic_table.GetElementSymbol(element_num)
            except Exception:
                symbol = "?"

            fragment_atoms.append(
                {
                    "cdxml_id": cdxml_id,
                    "element": element_num,
                    "color": color,
                    "symbol": symbol,
                    "num_hydrogens": num_hydrogens,
                    "implicit_h_color": implicit_h_color,
                }
            )

        return fragment_atoms

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

        When ``proton_index`` is ``None`` and ``color_code`` is ``None``,
        per-fragment automatic proton detection is used: each fragment's
        dominant colour is computed independently, and the uniquely
        coloured hydrogen in that fragment is selected.

        When ``proton_index`` or ``color_code`` is supplied, the same
        value is applied to every fragment (legacy behaviour).

        Args:
            proton_index (int | None): Explicit 1-based proton index.
            color_code (int | None): CDXML colour-table index.

        Returns:
            list[PKaMolecule]: One ``PKaMolecule`` per fragment.
        """
        if proton_index is None and color_code is None:
            return self.get_pka_molecules_auto()

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

    # ------------------------------------------------------------------
    # Per-fragment colour parsing and auto-detection
    # ------------------------------------------------------------------

    def parse_cdxml_fragment_colors(self):
        """Parse per-fragment atom colour information from a CDXML file.

        Unlike :meth:`parse_cdxml_element_colors`, which returns a flat
        list across all fragments, this method returns a **list of lists**
        where each inner list contains the atom-colour dicts for one
        top-level ``<fragment>`` element.

        The per-fragment grouping preserves the 1-to-1 correspondence
        between fragments and the ``Molecule`` objects returned by
        :meth:`CDXFile.molecules` / ``Chem.MolsFromCDXMLFile``.

        Returns:
            list[list[dict]]: One sub-list per top-level fragment.
                Each dict has the same keys as
                :meth:`parse_cdxml_element_colors`.

        Raises:
            ValueError: If the file cannot be parsed as valid CDXML.
        """
        from rdkit.Chem import GetPeriodicTable

        pt = GetPeriodicTable()
        root = self._parse_cdxml_root()

        fragments_atoms = []
        for fragment in self._iter_top_level_fragments(root):
            fragment_atoms = self._parse_fragment_nodes(fragment, pt)
            if fragment_atoms:
                fragments_atoms.append(fragment_atoms)

        return fragments_atoms

    def _detect_proton_in_fragment(self, atoms, fragment_index=None):
        """Auto-detect the uniquely coloured proton within a single fragment.

        The logic mirrors :meth:`_proton_index_auto` but operates on a
        fragment-local atom list and returns a **fragment-local** 0-based
        index (the offset within the fragment's atom list before
        ``AddHs``).

        For implicit / functional-group hydrogens the method returns a
        tuple ``(heavy_atom_local_idx, atom_dict)`` so the caller can
        resolve the final Molecule index via RDKit ``AddHs``.

        Args:
            atoms: List of atom dicts for one fragment (from
                :meth:`parse_cdxml_fragment_colors`).
            fragment_index: Optional fragment number for error messages.

        Returns:
            dict: With keys:
                * ``type`` – ``"explicit"`` or ``"implicit"``
                * ``local_idx`` – 0-based position in the fragment's
                  atom list (before ``AddHs``)
                * ``atom`` – the atom dict

        Raises:
            ValueError: When the proton cannot be identified.
        """
        frag_info = (
            f" (fragment {fragment_index})"
            if fragment_index is not None
            else ""
        )

        all_colors = []
        for a in atoms:
            all_colors.append(a["color"])
            if a["implicit_h_color"] is not None:
                all_colors.append(a["implicit_h_color"])

        color_counts = Counter(all_colors)
        if len(color_counts) < 2:
            raise ValueError(
                f"All atoms in fragment{frag_info} share the same colour. "
                "Cannot auto-detect the proton to remove."
            )

        majority_color = color_counts.most_common(1)[0][0]

        explicit_h = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["symbol"] == "H" and a["color"] != majority_color
        ]
        # Functional-group implicit H on *heavy* atoms only — skip
        # explicit H nodes to avoid double-counting.
        fg_h = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["symbol"] != "H"
            and a["implicit_h_color"] is not None
            and a["implicit_h_color"] != majority_color
        ]

        total = len(explicit_h) + len(fg_h)

        if total == 0:
            unique = [a for a in atoms if a["color"] != majority_color]
            if unique:
                raise ValueError(
                    f"Uniquely coloured atom(s) found in fragment{frag_info} "
                    f"but none are hydrogen (found: "
                    f"{[a['symbol'] for a in unique]})."
                )
            raise ValueError(
                f"No uniquely coloured atom found in fragment{frag_info}."
            )

        if total > 1:
            raise ValueError(
                f"Multiple uniquely coloured H atoms ({total}) in "
                f"fragment{frag_info}. Cannot determine which proton "
                f"to remove."
            )

        if explicit_h:
            idx, atom = explicit_h[0]
            return {"type": "explicit", "local_idx": idx, "atom": atom}

        idx, atom = fg_h[0]
        return {"type": "implicit", "local_idx": idx, "atom": atom}

    def get_pka_molecules_auto(self):
        """Per-fragment proton auto-detection → list of PKaMolecule.

        For each top-level ``<fragment>`` in the CDXML file:

        1. Parse atom colours within that fragment.
        2. Identify the dominant colour of the fragment.
        3. Find the uniquely coloured hydrogen (explicit or
           functional-group implicit).
        4. Map the proton back to the 1-based index in the
           ``Molecule`` produced after ``Chem.AddHs`` + embedding.
        5. Wrap the molecule as a :class:`PKaMolecule`.

        Falls back to :meth:`get_colored_proton_index` (file-global
        detection) if per-fragment parsing produces only one fragment.

        Returns:
            list[PKaMolecule]: One ``PKaMolecule`` per fragment with
                the per-fragment ``proton_index`` attached.

        Raises:
            ValueError: If detection fails for any fragment.
        """
        from rdkit import Chem

        from chemsmart.io.molecules.structure import PKaMolecule

        fragments_atoms = self.parse_cdxml_fragment_colors()
        molecules = self.molecules  # list[Molecule], one per fragment

        if len(fragments_atoms) != len(molecules):
            logger.warning(
                f"Fragment count ({len(fragments_atoms)}) differs from "
                f"molecule count ({len(molecules)}). Falling back to "
                f"global proton detection."
            )
            proton_index = self.get_colored_proton_index()
            return [
                PKaMolecule(molecule=mol, proton_index=proton_index)
                for mol in molecules
            ]

        # Read RDKit mols once for implicit-H resolution
        rdkit_mols = list(
            Chem.MolsFromCDXMLFile(self.filename, removeHs=False)
        )
        rdkit_mols_h = []
        for rm in rdkit_mols:
            if rm is not None:
                rdkit_mols_h.append(Chem.AddHs(rm))
            else:
                rdkit_mols_h.append(None)

        pka_molecules = []
        for frag_idx, (frag_atoms, mol, rdkit_mol_h) in enumerate(
            zip(fragments_atoms, molecules, rdkit_mols_h)
        ):
            detection = self._detect_proton_in_fragment(
                frag_atoms, fragment_index=frag_idx + 1
            )

            if detection["type"] == "explicit":
                # Explicit H: local_idx is 0-based in the pre-AddHs atom
                # list.  After AddHs the original atoms keep their
                # indices, so 1-based proton_index = local_idx + 1.
                proton_index = detection["local_idx"] + 1
            else:
                # Implicit / functional-group H: need to find the H
                # bonded to the heavy atom via RDKit after AddHs.
                if rdkit_mol_h is None:
                    raise ValueError(
                        f"RDKit molecule for fragment {frag_idx + 1} is "
                        f"None; cannot resolve implicit H index."
                    )
                heavy_idx = detection["local_idx"]
                heavy_atom = rdkit_mol_h.GetAtomWithIdx(heavy_idx)
                h_indices = [
                    bond.GetOtherAtom(heavy_atom).GetIdx()
                    for bond in heavy_atom.GetBonds()
                    if bond.GetOtherAtom(heavy_atom).GetSymbol() == "H"
                ]
                if not h_indices:
                    raise ValueError(
                        f"No hydrogen bonded to "
                        f"{detection['atom']['symbol']} "
                        f"(CDXML id={detection['atom']['cdxml_id']}) "
                        f"in fragment {frag_idx + 1} after AddHs."
                    )
                proton_index = h_indices[0] + 1  # 1-based

            logger.info(
                f"Fragment {frag_idx + 1}: detected proton at 1-based "
                f"index {proton_index} "
                f"(type={detection['type']}, "
                f"cdxml_id={detection['atom']['cdxml_id']})."
            )
            pka_molecules.append(
                PKaMolecule(molecule=mol, proton_index=proton_index)
            )

        return pka_molecules
