import re
from functools import cached_property

import numpy as np

from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.periodictable import PeriodicTable as pt
from chemsmart.utils.utils import string2index_1based

p = pt()


class PDBFile(FileMixin):
    """
    Parser for PDB coordinate files.

    This class handles PDB files containing molecular coordinates, supporting
    both single-model files and multi-model files (e.g., NMR ensembles or
    MD trajectories). It can extract molecular geometries including residue
    metadata such as atom names, residue names, chain IDs, and record types.

    PDB files use a fixed-width column format (PDB v3.3), so this parser
    reads raw lines preserving whitespace for accurate column slicing.
    """

    def __init__(self, filename):
        """
        Initialize PDB file parser.

        Args:
            filename (str): Path to the PDB file to parse
        """
        self.filename = filename

    def __repr__(self):
        return f"PDBFile({self.filename})"

    def __str__(self):
        return f"PDBFile object with filename: {self.filename}"

    # ------------------------------------------------------------------
    # Raw line access (PDB is fixed-width; must not strip leading spaces)
    # ------------------------------------------------------------------

    @cached_property
    def raw_lines(self):
        """
        Read and cache file lines preserving column formatting.

        Unlike the inherited ``contents`` property which strips whitespace,
        this preserves leading spaces required for PDB column parsing.

        Returns:
            list[str]: Lines with only the trailing newline removed.
        """
        with open(self.filepath, "r") as f:
            return [line.rstrip("\n") for line in f.readlines()]

    # ------------------------------------------------------------------
    # High-level convenience properties
    # ------------------------------------------------------------------

    @cached_property
    def num_atoms(self):
        """
        Get number of atoms from the last model in the PDB file.

        Returns:
            int: Number of atoms in the molecular structure
        """
        return self.molecule.num_atoms

    @cached_property
    def molecule(self):
        """
        Get the last molecule (model) from the PDB file.

        Returns:
            Molecule: Last molecular structure in the file
        """
        return self.get_molecules(index="-1")

    # ------------------------------------------------------------------
    # Molecule extraction
    # ------------------------------------------------------------------

    def get_molecules(self, index=":", return_list=False):
        """
        Extract molecular structures from PDB file.

        Parses all MODEL/ENDMDL blocks (or the single implicit model)
        and returns the requested subset of molecules.

        Args:
            index (str): Index specification for molecule selection.
                Uses 1-based indexing (e.g., ``"1"`` for first,
                ``"-1"`` for last, ``":"`` for all).
            return_list (bool): Whether to always return a list.

        Returns:
            Molecule or list[Molecule]: Single molecule or list of
            molecules depending on *index* and *return_list*.

        Raises:
            ValueError: If no ATOM/HETATM records are found in the file.
        """
        models = self._parse_models()
        if not models:
            raise ValueError(
                f"No ATOM/HETATM records found in PDB file: {self.filepath}"
            )

        parsed_index = string2index_1based(index)
        molecules = models[parsed_index]

        if return_list and not isinstance(molecules, list):
            return [molecules]
        return molecules

    # ------------------------------------------------------------------
    # Internal parsing
    # ------------------------------------------------------------------

    def _parse_models(self):
        """
        Parse all PDB models from the cached raw lines into
        ``Molecule`` objects.

        Returns:
            list[Molecule]: One ``Molecule`` per MODEL block (or one
            for the entire file when no MODEL records are present).
        """
        lines = self.raw_lines
        models = []
        current_atom_lines = []
        seen_model_records = False

        def flush_current_model():
            if current_atom_lines:
                models.append(
                    self._molecule_from_atom_lines(current_atom_lines)
                )

        for line in lines:
            record_name = line[:6].strip()
            if record_name == "MODEL":
                # Flush any accumulated atoms before starting a new model
                if current_atom_lines:
                    flush_current_model()
                seen_model_records = True
                current_atom_lines = []
                continue

            if record_name == "ENDMDL":
                flush_current_model()
                current_atom_lines = []
                continue

            if record_name in {"ATOM", "HETATM"}:
                current_atom_lines.append(line.rstrip("\n"))

        if seen_model_records:
            flush_current_model()
        elif current_atom_lines:
            models.append(self._molecule_from_atom_lines(current_atom_lines))

        return models

    @staticmethod
    def _molecule_from_atom_lines(atom_lines):
        """
        Build a ``Molecule`` from parsed PDB ATOM/HETATM lines.

        Extracts element symbols, 3-D coordinates, and residue-level
        metadata (atom names, residue names/numbers, chain IDs, record
        types) from fixed-width PDB columns.

        Args:
            atom_lines (list[str]): Raw PDB lines (ATOM or HETATM records).

        Returns:
            Molecule: A new ``Molecule`` populated with coordinates and
            PDB metadata stored in both direct attributes and the
            ``info`` dictionary.
        """
        from chemsmart.io.molecules.structure import Molecule

        symbols = []
        positions = []
        atom_names = []
        residue_names = []
        residue_numbers = []
        chain_ids = []
        record_types = []

        for line in atom_lines:
            record_types.append(line[:6].strip() or "HETATM")
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip() or "MOL"
            chain_id = line[21:22].strip() or ""
            residue_number_text = line[22:26].strip()
            residue_number = (
                int(residue_number_text) if residue_number_text else 1
            )

            element_text = line[76:78].strip()
            if element_text:
                symbol = p.to_element(element_text)
            else:
                symbol = PDBFile._infer_element_from_atom_name(atom_name)

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            symbols.append(symbol)
            positions.append([x, y, z])
            atom_names.append(atom_name or symbol)
            residue_names.append(residue_name)
            residue_numbers.append(residue_number)
            chain_ids.append(chain_id)

        molecule = Molecule(
            symbols=symbols,
            positions=np.array(positions, dtype=float),
            info={
                "atom_name": atom_names,
                "residue_name": residue_names,
                "residue_number": residue_numbers,
                "chain_id": chain_ids,
                "record_type": record_types,
            },
        )
        molecule.atom_names = atom_names
        molecule.residue_names = residue_names
        molecule.residue_numbers = residue_numbers
        molecule.chain_ids = chain_ids
        molecule.record_type = record_types
        return molecule

    @staticmethod
    def _infer_element_from_atom_name(atom_name):
        """
        Infer element symbol from PDB atom-name field when columns 77-78
        are blank.

        This uses PDB-specific heuristics instead of interpreting the
        full atom name as an element symbol. The rules are:

        - Strip leading digits (e.g. ``1H`` → ``H``).
        - Keep only letters from the remaining name.
        - Try a normalized two-letter element candidate first (e.g.
          ``FE`` → ``Fe``), unless the atom name matches a common
          biomolecular atom-label pattern such as ``CA`` for C-alpha.
        - Otherwise, use only the first letter as the element.

        Args:
            atom_name (str): The 4-character PDB atom name field.

        Returns:
            str: Canonical element symbol.

        Raises:
            ValueError: If the element cannot be inferred.
        """
        if not atom_name:
            raise ValueError(
                f"Unable to infer element from atom name '{atom_name}'"
            )

        # Remove leading digits (e.g. "1H", "2CA")
        name = re.sub(r"^[0-9]+", "", atom_name.strip())
        # Keep only letters thereafter
        cleaned = re.sub(r"[^A-Za-z]", "", name)

        if not cleaned:
            raise ValueError(
                f"Unable to infer element from atom name '{atom_name}'"
            )

        if len(cleaned) == 1:
            try:
                return p.to_element(cleaned[0].upper())
            except Exception:
                raise ValueError(
                    f"Unable to infer element from atom name '{atom_name}'"
                )

        # Common biomolecular atom-label prefixes that should remain
        # single-letter element assignments, e.g. CA/CB/CG for
        # carbon alpha/beta/gamma.
        ambiguous_biomolecular_names = {
            "CA",
            "CB",
            "CG",
            "CD",
            "CE",
            "CZ",
            "CH",
            "HA",
            "HB",
            "HG",
            "HD",
            "HE",
            "HZ",
            "HH",
            "HN",
            "ND",
            "NE",
            "NH",
            "NZ",
            "OD",
            "OE",
            "OG",
            "OH",
            "SD",
            "SG",
        }

        normalized_two_letter = f"{cleaned[0].upper()}{cleaned[1].lower()}"
        if cleaned[:2].upper() not in ambiguous_biomolecular_names:
            try:
                return p.to_element(normalized_two_letter)
            except Exception:
                pass

        candidate = cleaned[0].upper()

        try:
            return p.to_element(candidate)
        except Exception:
            raise ValueError(
                f"Unable to infer element from atom name '{atom_name}'"
            )

    # ------------------------------------------------------------------
    # PDB v3.3 output formatting
    # ------------------------------------------------------------------

    @staticmethod
    def format_pdb_block(molecule, pdb_block):
        """
        Normalize a raw RDKit PDB block to strict PDB v3.3 atom-line
        formatting, honouring any residue/chain metadata stored on
        *molecule*.

        This is the authoritative formatter used by
        ``Molecule.to_pdb()``.  It replaces every ATOM/HETATM line
        produced by RDKit with a column-exact record built from the
        molecule's own coordinate and metadata arrays, so that chain
        IDs, residue names, occupancy, etc. are always taken from the
        ``Molecule`` object rather than from whatever RDKit inferred.

        Args:
            molecule (Molecule): Source of coordinates and metadata.
            pdb_block (str): Raw PDB text returned by
                ``Chem.MolToPDBBlock``.

        Returns:
            str: PDB v3.3-formatted text ending with an ``END`` record.
        """
        lines = pdb_block.splitlines()
        formatted_lines = []
        atom_index = 0

        for line in lines:
            record_name = line[:6].strip()
            if record_name in {"ATOM", "HETATM"}:
                formatted_lines.append(
                    PDBFile._format_pdb_atom_line(
                        molecule,
                        atom_index=atom_index,
                        default_record_name=record_name,
                    )
                )
                atom_index += 1
                continue

            if record_name == "END":
                continue

            if line.strip():
                formatted_lines.append(line.rstrip())

        # Append any atoms that RDKit may not have emitted
        while atom_index < molecule.num_atoms:
            formatted_lines.append(
                PDBFile._format_pdb_atom_line(
                    molecule,
                    atom_index=atom_index,
                    default_record_name="HETATM",
                )
            )
            atom_index += 1

        formatted_lines.append("END")
        return "\n".join(formatted_lines) + "\n"

    @staticmethod
    def _format_pdb_atom_line(
        molecule, atom_index, default_record_name="HETATM"
    ):
        """Create one strict-width PDB v3.3 ATOM/HETATM record for
        *atom_index* of *molecule*."""
        atom = PDBFile._pdb_atom_object(molecule, atom_index)
        raw_atom_type = (
            PDBFile._pdb_atom_field(
                molecule,
                atom,
                atom_index,
                ("atom_type", "type", "symbol", "element"),
            )
            or molecule.chemical_symbols[atom_index]
        )
        element = PDBFile._pdb_element_symbol(raw_atom_type)

        record_name = PDBFile._pdb_atom_field(
            molecule,
            atom,
            atom_index,
            ("record_type", "pdb_record_type", "record_name"),
        )
        if record_name is None:
            hetero = PDBFile._pdb_atom_field(
                molecule,
                atom,
                atom_index,
                ("hetero", "hetatm", "is_hetero"),
            )
            if hetero is not None:
                record_name = "HETATM" if bool(hetero) else "ATOM"
        record_name = (record_name or default_record_name or "HETATM").upper()
        if record_name not in {"ATOM", "HETATM"}:
            record_name = default_record_name

        atom_name = PDBFile._pdb_atom_field(
            molecule,
            atom,
            atom_index,
            ("atom_name", "atom_names", "name", "atom_label", "label"),
        ) or str(raw_atom_type)
        atom_name_field = PDBFile._pdb_atom_name_field(atom_name, element)

        residue_name = (
            PDBFile._pdb_atom_field(
                molecule,
                atom,
                atom_index,
                (
                    "residue_name",
                    "residue_names",
                    "resname",
                    "residue",
                    "pdb_residue_name",
                ),
            )
            or "MOL"
        )
        residue_name = str(residue_name).strip()[:3].rjust(3)

        residue_number = PDBFile._pdb_atom_field(
            molecule,
            atom,
            atom_index,
            (
                "residue_number",
                "residue_numbers",
                "resseq",
                "resseqs",
                "pdb_residue_number",
            ),
        )
        residue_number = 1 if residue_number is None else int(residue_number)

        chain_id = PDBFile._pdb_atom_field(
            molecule,
            atom,
            atom_index,
            ("chain_id", "chain", "chain_ids", "pdb_chain_id"),
        )
        chain_id = "" if chain_id is None else str(chain_id).strip()[:1]

        occupancy = PDBFile._pdb_atom_field(
            molecule,
            atom,
            atom_index,
            ("occupancy", "pdb_occupancy"),
        )
        occupancy = 1.00 if occupancy is None else float(occupancy)

        temp_factor = PDBFile._pdb_atom_field(
            molecule,
            atom,
            atom_index,
            ("temp_factor", "bfactor", "b_factor", "pdb_temp_factor"),
        )
        temp_factor = 0.00 if temp_factor is None else float(temp_factor)

        x, y, z = molecule.positions[atom_index]
        serial = atom_index + 1

        return (
            f"{record_name:<6}{serial:>5} "
            f"{atom_name_field}"
            f" {residue_name} {chain_id:1}"
            f"{residue_number:>4} "
            f"   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{occupancy:>6.2f}{temp_factor:>6.2f}"
            f"          {element:>2}"
        )

    @staticmethod
    def _pdb_atom_object(molecule, atom_index):
        """Return an atom-like object from ``molecule.atoms`` when available."""
        atoms = getattr(molecule, "atoms", None)
        if atoms is None:
            return None
        try:
            if len(atoms) == molecule.num_atoms:
                return atoms[atom_index]
        except TypeError:
            return None
        return None

    @staticmethod
    def _pdb_atom_field(molecule, atom, atom_index, keys):
        """Resolve atom metadata from *atom* object, *molecule* attrs, or
        ``molecule.info`` dict."""
        for key in keys:
            value = PDBFile._pdb_lookup(atom, key)
            if value is not None:
                return value

            value = PDBFile._pdb_lookup(molecule, key)
            resolved = PDBFile._pdb_indexed_value(value, atom_index)
            if resolved is not None:
                return resolved

            value = PDBFile._pdb_lookup(molecule.info, key)
            resolved = PDBFile._pdb_indexed_value(value, atom_index)
            if resolved is not None:
                return resolved
        return None

    @staticmethod
    def _pdb_lookup(container, key):
        """Return ``container[key]`` or ``getattr(container, key)``."""
        if container is None:
            return None
        if isinstance(container, dict):
            return container.get(key)
        return getattr(container, key, None)

    @staticmethod
    def _pdb_indexed_value(value, atom_index):
        """Extract the per-atom value at *atom_index* from a scalar,
        list, tuple, or numpy array."""
        if value is None:
            return None
        if isinstance(value, np.ndarray):
            if value.ndim == 0:
                return value.item()
            if atom_index < len(value):
                return value[atom_index]
            return None
        if isinstance(value, (list, tuple)):
            if atom_index < len(value):
                return value[atom_index]
            return None
        return value

    @staticmethod
    def _pdb_element_symbol(atom_type):
        """Map atom type/symbol text to a canonical PDB element symbol."""
        return p.to_element(str(atom_type))

    @staticmethod
    def _pdb_atom_name_field(atom_name, element):
        """Format *atom_name* into the 4-character PDB columns 13-16."""
        cleaned = "".join(str(atom_name).strip().split())[:4] or element
        if len(element.strip()) == 1 and not cleaned[:1].isdigit():
            return cleaned.rjust(4)
        return cleaned.ljust(4)

    # ------------------------------------------------------------------
    # Writing helpers
    # ------------------------------------------------------------------

    @staticmethod
    def write(molecule, filename, mode="w", **kwargs):
        """
        Write a ``Molecule`` to PDB format file.

        This is a convenience wrapper around ``Molecule.write_pdb``.

        Args:
            molecule (Molecule): The molecule to write.
            filename (str): Output PDB file path.
            mode (str): File write mode. Default ``'w'``.
            **kwargs: Forwarded to ``Molecule.write_pdb`` (e.g.
                ``flavor``, ``add_bonds``, ``bond_cutoff_buffer``,
                ``adjust_H``).
        """
        molecule.write_pdb(filename, mode=mode, **kwargs)
