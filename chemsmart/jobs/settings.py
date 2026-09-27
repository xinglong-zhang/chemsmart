import importlib
import logging
import os
from collections.abc import Mapping
from dataclasses import dataclass

import yaml

from chemsmart.utils.utils import update_dict_with_existing_keys

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class NativeWord:
    """A program's own word a project would write into its input.

    ``field`` is the settings field that carries it and ``word`` the word
    (empty when the whole field is the finding: a field that replaces the
    input or the route the host writes).  ``states`` says what the word
    asks the program for, and ``route`` names the typed setting that
    asks for the same thing -- or, where none exists, what the host
    writes instead.
    """

    field: str
    word: str
    states: str
    route: str


def native_input_fields(program):
    """The settings fields of *program* whose values reach its input verbatim.

    Declared by each program's settings module (``NATIVE_INPUT_FIELDS``);
    a program that declares none has none a project can set.
    """

    module = _program_settings_module(program)
    return tuple(getattr(module, "NATIVE_INPUT_FIELDS", ()) or ())


def agent_refused_fields(program):
    """The native fields *program* refuses whenever an Agent sets them.

    Declared by each program's settings module (``AGENT_REFUSED_FIELDS``),
    beside the ``native_words`` table that refuses them, so the capability
    an Agent reads can leave them out instead of offering a field its
    project tool refuses.
    """

    module = _program_settings_module(program)
    return tuple(getattr(module, "AGENT_REFUSED_FIELDS", ()) or ())


def project_native_words(program, sections):
    """Every native word project *sections* carry that has a typed route.

    Returns ``((section, NativeWord), ...)`` from the program's own table
    (``native_words`` in its settings module).  A word the table does not
    route is not returned: it stays the reviewer's to read, verbatim.
    Used where an Agent authors a project; a person's project keeps every
    field the loader accepts.
    """

    describe = getattr(_program_settings_module(program), "native_words", None)
    if describe is None or not isinstance(sections, Mapping):
        return ()
    found = []
    for section, settings in sections.items():
        if isinstance(settings, Mapping):
            found.extend((str(section), item) for item in describe(settings))
    return tuple(found)


def _program_settings_module(program):
    try:
        return importlib.import_module(f"chemsmart.jobs.{program}.settings")
    except ImportError:
        return None


def native_route_words(value):
    """The words of a route-channel value: a string or a list of strings.

    Words are split on whitespace outside parentheses, so
    ``IRC=(MaxPoints=80, StepSize=10)`` stays one word.
    """

    if value is None:
        return ()
    if isinstance(value, (list, tuple)):
        text = " ".join(str(item) for item in value)
    else:
        text = str(value)
    words, depth, current = [], 0, []
    for char in text:
        if char == "(":
            depth += 1
        elif char == ")":
            depth = max(0, depth - 1)
        if char.isspace() and depth == 0:
            if current:
                words.append("".join(current))
                current = []
            continue
        current.append(char)
    if current:
        words.append("".join(current))
    return tuple(words)


def native_field_is_set(value):
    """Whether a native field holds anything a writer would write."""

    if value is None or value is False:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (list, tuple, dict)):
        return bool(value)
    return True


#: What a ChemSmart functional literal means, whichever program runs it.
#:
#: A project literal named here is one functional in every program that
#: accepts it.  Each program's settings class spells it in that program's
#: own vocabulary (or refuses it where the program has no spelling), and
#: each result reader reads the applied form back from the program's own
#: record; this table is the meaning both sides answer to, so no program
#: table restates it.  A literal absent from the table passes through to
#: each program unchanged, which is safe only where the programs agree on
#: the name.
#:
#: Measured, not recalled (CUHK Slurm 2149277/2149278, water and seven
#: other species at tight matched numerics): Gaussian ``B3LYP``, ORCA
#: ``B3LYP/G`` (which prints ``LDAOpt .... VWN-3``) and PySCF ``b3lypg``
#: (libxc 402) agree to 2.2e-6 Eh; ORCA's bare ``B3LYP`` (``VWN-5``) and
#: PySCF ``b3lyp5`` agree to 8.6e-7 Eh and lie 0.032-0.154 Eh above the
#: first form, 2.34 kcal/mol apart in the vertical IP of water and 0.38
#: in the C-Cl homolysis of CH3Cl.  Before this table ``b3lyp`` meant the
#: first form in Gaussian and PySCF and the second in ORCA.
FUNCTIONAL_IDENTITIES = {
    "b3lyp": {
        "functional_family": "b3lyp",
        "correlation_convention": "vwn3_gaussian",
        "definition": (
            "B3LYP of Stephens et al. (1994) with the VWN RPA local "
            "correlation: Gaussian's B3LYP (VWN functional III), ORCA's "
            "B3LYP/G, libxc HYB_GGA_XC_B3LYP"
        ),
    },
    "b3lyp5": {
        "functional_family": "b3lyp",
        "correlation_convention": "vwn5",
        "definition": (
            "B3LYP with the VWN functional V local correlation: ORCA's and "
            "TURBOMOLE's bare B3LYP, libxc HYB_GGA_XC_B3LYP5"
        ),
    },
    # Gaussian's spelling differs, and Gaussian completes a route word
    # that prefixes one of its keywords: ``pbe0`` ran PBE0-DH, a double
    # hybrid, 0.0365 Eh above PBE0 on water (CUHK Slurm 2149277). ORCA
    # PBE0 and PySCF pbe0 agree to 1.0e-5 Eh there, PBE to 1.4e-5 Eh --
    # the precision of the PW92 constants, well inside the numerics of a
    # relative energy.
    "pbe0": {
        "functional_family": "pbe0",
        "correlation_convention": "pw92",
        "definition": (
            "PBE0 (25% exact exchange on PBE): Gaussian's PBE1PBE, ORCA's "
            "PBE0, libxc HYB_GGA_XC_PBEH"
        ),
    },
    "pbe": {
        "functional_family": "pbe",
        "correlation_convention": "pw92",
        "definition": "the PBE GGA: Gaussian's PBEPBE, ORCA's PBE",
    },
    # ORCA's BP86 puts P86 on the Perdew-Wang 92 local correlation (it
    # prints ``LDAOpt .... PW91-LDA``); Gaussian's BP86 and libxc's B88+P86
    # put it on Perdew-Zunger 81. At tight numerics ORCA sits 0.96-2.67 mEh
    # below Gaussian on seven species, 0.74 kcal/mol apart in the vertical
    # IP of water and 1.06 in the C-Cl homolysis of CH3Cl, while Gaussian
    # and PySCF agree in those to 0.003 kcal/mol (their totals differ by
    # 6e-5-2.3e-4 Eh inside P86) -- CUHK Slurm 2149487.
    "bp86": {
        "functional_family": "bp86",
        "correlation_convention": "pz81",
        "definition": (
            "Becke 88 exchange with Perdew 86 correlation on the "
            "Perdew-Zunger 81 local correlation: Gaussian's BP86, libxc "
            "B88+P86"
        ),
    },
    "bp86-pw92": {
        "functional_family": "bp86",
        "correlation_convention": "pw92",
        "definition": (
            "Becke 88 exchange with Perdew 86 correlation on the Perdew-Wang "
            "92 local correlation: ORCA's BP86"
        ),
    },
}

#: Other spellings of a literal in ``FUNCTIONAL_IDENTITIES``.
FUNCTIONAL_LITERAL_SYNONYMS = {
    "b3lypg": "b3lyp",
    "b3lyp/g": "b3lyp",
    "b3lyp-g": "b3lyp",
    "b3lyp-vwn5": "b3lyp5",
    "pbe1pbe": "pbe0",
    "pbepbe": "pbe",
}


def canonical_functional_literal(value):
    """Return the ChemSmart literal a functional spelling means.

    A spelling in ``FUNCTIONAL_IDENTITIES`` or its synonym table answers
    the literal it names; anything else answers ``None``, which means the
    literal carries no cross-program definition here and passes through.
    """

    if value is None:
        return None
    key = str(value).strip().lower()
    key = FUNCTIONAL_LITERAL_SYNONYMS.get(key, key)
    return key if key in FUNCTIONAL_IDENTITIES else None


def functional_identity(value):
    """Return the program-neutral identity of a functional spelling.

    ``{"literal", "functional_family", "correlation_convention"}`` for a
    literal the table defines, ``None`` otherwise.
    """

    literal = canonical_functional_literal(value)
    if literal is None:
        return None
    record = FUNCTIONAL_IDENTITIES[literal]
    return {
        "literal": literal,
        "functional_family": record["functional_family"],
        "correlation_convention": record["correlation_convention"],
    }


def functional_resolution_record(
    *, program, functional, ab_initio, native, source
):
    """What one program is told for a project functional, as a host record.

    Every program's settings module answers its receipt through this one
    shape: the literal the project named, the literal it means, the native
    spelling the program's own writer produces (``native``, passed in so
    the record states the writer's output rather than a copy of its table)
    and the program-neutral identity.  ``literal_preserved`` says the name
    has no cross-program definition here, which is not a claim that the
    programs agree on it.
    """

    base = {
        "schema_version": "chemsmart.functional-resolution.v2",
        "program": str(program),
        "source": str(source),
    }
    if ab_initio is not None and str(ab_initio).strip():
        return {
            **base,
            "status": "not_applicable",
            "requested_method_kind": "ab_initio",
            "requested_literal": None,
            "canonical_literal": "",
            "applied_native": None,
            "functional_family": "wavefunction",
            "correlation_convention": "not_applicable",
            "rule_id": f"{program}.functional.not_applicable_ab_initio",
        }
    if functional is None or not str(functional).strip():
        return {
            **base,
            "status": "missing",
            "requested_method_kind": "dft",
            "requested_literal": None,
            "canonical_literal": "",
            "applied_native": None,
            "functional_family": "",
            "correlation_convention": "unresolved",
            "rule_id": f"{program}.functional.missing",
        }
    requested = str(functional).strip()
    identity = functional_identity(requested)
    if identity is None:
        return {
            **base,
            "status": "literal_preserved",
            "requested_method_kind": "dft",
            "requested_literal": requested,
            "canonical_literal": requested.lower(),
            "applied_native": None if native is None else str(native),
            "functional_family": "no_cross_program_definition",
            "correlation_convention": "not_declared",
            "rule_id": f"{program}.functional.literal_preserved",
        }
    return {
        **base,
        "status": "canonical_literal",
        "requested_method_kind": "dft",
        "requested_literal": requested,
        "canonical_literal": identity["literal"],
        "applied_native": None if native is None else str(native),
        "functional_family": identity["functional_family"],
        "correlation_convention": identity["correlation_convention"],
        "rule_id": f"{program}.functional.{identity['literal']}",
    }


#: Which excited-state manifolds a reference has, whichever program runs
#: the response.  A closed-shell (singlet) reference has spin-adapted
#: singlet and triplet excitations, asked for alone or together; an
#: open-shell reference has one spin-conserving manifold, ``unrestricted``,
#: whose roots are not spin eigenfunctions.  This is a fact about the
#: reference and not about a program, so every program's settings ask this
#: one rule rather than restating it.
TD_CLOSED_SHELL_MANIFOLDS = ("singlet", "singlet_triplet", "triplet")
TD_OPEN_SHELL_MANIFOLD = "unrestricted"


def td_manifold_reference_refusal(manifold, multiplicity):
    """Why *manifold* cannot belong to a reference of *multiplicity*.

    Returns None when the pair is admissible or either is unknown, and
    otherwise the sentence a settings class refuses with, naming the
    manifold the reference does have.
    """

    if manifold is None or multiplicity is None:
        return None
    word = str(manifold).strip().lower()
    closed_shell = int(multiplicity) == 1
    if word in TD_CLOSED_SHELL_MANIFOLDS and not closed_shell:
        return (
            f"state_manifold={word!r} is spin-adapted and needs a "
            f"closed-shell (singlet) reference; multiplicity {multiplicity} "
            "has one manifold, state_manifold: "
            f"{TD_OPEN_SHELL_MANIFOLD}."
        )
    if word == TD_OPEN_SHELL_MANIFOLD and closed_shell:
        return (
            "A closed-shell (singlet) reference asks for singlet, triplet "
            "or singlet_triplet excitations; state_manifold: "
            f"{TD_OPEN_SHELL_MANIFOLD} is the one manifold of an open-shell "
            "reference."
        )
    return None


#: What ``broken_symmetry: true`` asks for, whichever program runs it: the
#: lowest spin-polarised (broken-symmetry, open-shell singlet) solution of
#: the unrestricted equations at Ms = 0 -- the bound multiplicity 1 -- where
#: one lies below the spin-symmetric solution, and the spin-symmetric one
#: where none does.  Each program's settings write it with the mechanism
#: measured to reach that solution there on one-electron sites, and say
#: which (the measurements below).  A mixing guess mixes one orbital pair:
#: on two S = 1 Ni(II) centres ORCA's reached a state 33.40 mEh above the
#: one a flip of one centre from the quintet reaches, with the spin on the
#: bridges (R11 truth-2, CUHK 2157086), which ORCA's ``site_spin_flip``
#: writes; Gaussian's and PySCF's mechanisms were not measured there:
#:
#: - Gaussian: the method unrestricted (``U``) from ``guess=mix``, the alpha
#:   HOMO and LUMO of Gaussian's guess mixed 50:50;
#: - ORCA: ``HFTyp UHF`` with ``GuessMix`` at this angle, the same mixing of
#:   ORCA's guess;
#: - PySCF: the restricted solution's own RHF/RKS -> UHF/UKS instability
#:   followed by PySCF's stability analysis into the unrestricted space, and
#:   internal instabilities followed until stable.
#:
#: Whether the symmetry actually broke is a fact about the result, read
#: from its <S**2>, never assumed from the request.  Measured, not recalled
#: (R10 Q18 oracles O0 and O0b, CUHK Slurm 2153330 and 2153375, fixed
#: geometries, def2-SVP, matched tight numerics): one solution in all three
#: programs for H2 at 2.00 A (UB3LYP -1.0184866 Eh, <S**2> 0.7053, within
#: 1.2e-8 Eh), ethylene twisted 85 degrees (-78.4223308, within 6e-8) and
#: p-benzyne (-230.704724, <S**2> 0.9703, within 1.1e-6), and the restricted
#: energy where none lies lower (H2 at 0.74 A, planar ethylene under B3LYP).
#: A mixing guess is not one mechanism across programs: which orbitals are a
#: guess's HOMO and LUMO is a program fact -- PySCF's mix of its own guess
#: collapsed on p-benzyne, and PySCF's default ``init_guess_breaksym`` left
#: H2 at 2.00 A restricted -- and at an exactly degenerate geometry (the
#: 90-degree twist) Gaussian's mixing reached a solution 26 kcal/mol above
#: the one ORCA and PySCF reach, with <S**2> 1.001 that does not show it.
#: The particular solution is not portable; the request and its evidence
#: are.
BROKEN_SYMMETRY_MULTIPLICITY = 1
BROKEN_SYMMETRY_GUESS_MIX_DEGREES = 45


def broken_symmetry_request(value):
    """The boolean a ``broken_symmetry`` setting holds, or a refusal.

    ``None`` and ``False`` ask for nothing; ``True`` asks for the
    broken-symmetry open-shell singlet.  Any other value is refused rather
    than read as truthy: a string such as ``"mix"`` would otherwise ask for
    a program-specific guess the typed field exists to replace.
    """

    if value is None:
        return False
    if isinstance(value, bool):
        return value
    raise ValueError(
        f"broken_symmetry takes true or false, got {value!r}: it asks for an "
        "unrestricted open-shell singlet started from a symmetry-breaking "
        "guess, and each program's own mechanism is written by the host."
    )


#: The sentence tail every program's translation ends with: where the
#: answer to the request is read.  One text, so the three programs'
#: reviews point at one set of evidence.
BROKEN_SYMMETRY_EVIDENCE_SENTENCE = (
    "Whether the spin symmetry broke is read from the result, not assumed: "
    "its level states the reference that ran (rks, uks, ...) and "
    "broken_symmetry, spin_square gives <S**2> against 0 for the singlet, "
    "and a request whose solution stayed spin-symmetric raises the anomaly "
    "spin.broken_symmetry_request_unbroken with its numbers."
)


def broken_symmetry_refusal(requested, multiplicity):
    """Why a broken-symmetry request cannot be translated for *multiplicity*.

    Returns None when nothing was requested, the state is unknown (a
    project validated before a molecule is bound) or the state is the
    singlet the request defines; otherwise the sentence a writer refuses
    with, naming the routes that remain.
    """

    if not requested or multiplicity is None:
        return None
    if int(multiplicity) == BROKEN_SYMMETRY_MULTIPLICITY:
        return None
    return (
        f"broken_symmetry asks for the open-shell singlet (Ms = 0, "
        f"multiplicity {BROKEN_SYMMETRY_MULTIPLICITY}) reached from a "
        "symmetry-breaking guess; this node binds multiplicity "
        f"{multiplicity}, which is already an open-shell reference and has "
        "no broken-symmetry translation here. Bind multiplicity 1 for the "
        "broken-symmetry singlet, or remove broken_symmetry for the "
        "high-spin state."
    )


#: Settings whose values are each program's own words for its numerics:
#: ORCA's TightSCF is not Gaussian's scf=tight, DefGrid2 is no Gaussian
#: grid, and neither program's optimiser counts cycles as the other's does.
#: A settings object converted from another program's file does not carry
#: them across; it keeps its own program's defaults, as it always did.
PROGRAM_OWN_NUMERICS = ("scf_convergence", "defgrid", "geom_maxiter")


def without_program_own_numerics(settings):
    """A settings mapping without another program's numerics words."""

    values = dict(settings if isinstance(settings, dict) else vars(settings))
    for name in PROGRAM_OWN_NUMERICS:
        values.pop(name, None)
    return values


# Public top-level vocabulary owned by the loader below.  These are section
# names, not the Click job inventory: the two sets intentionally differ.
MOLECULAR_GAS_PHASE_JOB_SECTIONS = (
    "crest",
    "dias",
    "irc",
    "modred",
    "nci",
    "neb",
    "opt",
    "resp",
    "scan",
    "set",
    "traj",
    "ts",
    "uvvis",
    "wbi",
)
MOLECULAR_COMMON_DIRECT_SECTIONS = (
    *MOLECULAR_GAS_PHASE_JOB_SECTIONS,
    "qmmm",
    "sp",
    "td",
)


def molecular_project_section_names(program):
    """Return the exact top-level section vocabulary this loader accepts."""

    if program not in {"gaussian", "orca"}:
        raise ValueError(f"Unsupported molecular project program: {program}")
    names = {"gas", "solv", *MOLECULAR_COMMON_DIRECT_SECTIONS}
    if program == "gaussian":
        names.add("link")
    return tuple(sorted(names))


def molecular_project_section_sources(project_config, *, program, jobtype):
    """Return loader-order YAML sections for a Gaussian/ORCA job.

    This is the section-selection rule used by the molecular project loader:
    the historical phase section is applied first and an explicit job section
    overrides it.  Keeping this rule in the loader module lets capability and
    agent observations project the real semantics instead of maintaining a
    second phase/jobtype table.
    """

    available = {
        name
        for name, settings in project_config.items()
        if settings is not None
    }
    direct = jobtype if jobtype in available else None
    if jobtype in {"qmmm", "link"}:
        return (direct,) if direct is not None else ()
    # A ``td:`` section is read on its own: the loader seeds it from the
    # stage defaults, never from ``gas:`` or ``solv:``.  This answered
    # ``('solv', 'td')`` for a project with no ``gas:``, so the Agent's
    # project observation reported the solv level feeding a td stage whose
    # route the writer then built with no method at all.
    if jobtype == "td" and (direct is not None or "gas" in available):
        return (direct,) if direct is not None else ()
    if jobtype == "sp":
        phase = "solv" if "solv" in available else "gas"
    else:
        phase = "gas" if "gas" in available else "solv"
    return tuple(
        dict.fromkeys(
            name
            for name in (phase, direct)
            if name is not None and name in available
        )
    )


class MolecularJobSettings:
    """Common base job settings for molecular
    systems using Gaussian and ORCA jobs."""

    def __init__(
        self,
        ab_initio=None,
        functional=None,
        dispersion=None,
        basis=None,
        semiempirical=None,
        defgrid=None,
        charge=None,
        multiplicity=None,
        freq=True,
        numfreq=False,
        jobtype=None,
        title=None,
        solvent_model=None,
        solvent_id=None,
        additional_route_parameters=None,
        route_to_be_written=None,
        modred=None,
        gen_genecp_file=None,
        heavy_elements=None,
        heavy_elements_basis=None,
        light_elements_basis=None,
        custom_solvent=None,
        forces=False,
        input_string=None,
        **kwargs,
    ):
        self.ab_initio = ab_initio
        self.functional = functional
        self.dispersion = dispersion
        self.basis = basis
        self.semiempirical = semiempirical
        self.defgrid = defgrid
        self.charge = charge
        self.multiplicity = multiplicity
        self.freq = freq
        self.numfreq = numfreq
        self.jobtype = jobtype
        self.title = title
        self.solvent_model = solvent_model
        self.solvent_id = solvent_id
        self.additional_route_parameters = additional_route_parameters
        self.route_to_be_written = route_to_be_written
        self.modred = modred
        self.gen_genecp_file = gen_genecp_file
        self.heavy_elements = heavy_elements
        self.heavy_elements_basis = heavy_elements_basis
        self.light_elements_basis = light_elements_basis

        if custom_solvent is not None:
            if not isinstance(custom_solvent, str):
                raise ValueError(
                    "Custom solvent parameters must be a string! It can be either a string"
                    "giving the path of the custom solvent file or the custom solvent parameters"
                    "in free string format."
                )
            if os.path.exists(custom_solvent) and os.path.isfile(
                custom_solvent
            ):
                self.set_custom_solvent_via_file(custom_solvent)
            else:
                self.custom_solvent = custom_solvent
            # check that the last line of custom_solvent
            # is empty, if not, add an empty line
            if self.custom_solvent[-1] != "\n":
                self.custom_solvent += "\n"
            logger.debug(f"Custom solvent parameters: {self.custom_solvent}")
        else:
            self.custom_solvent = None

        self.forces = forces
        self.input_string = input_string

    def remove_solvent(self):
        self.solvent_model = None
        self.solvent_id = None
        self.custom_solvent = None

    def update_solvent(self, solvent_model=None, solvent_id=None):
        """Update solvent model and solvent identity for implicit solvation.

        Solvent models available: ['pcm', 'iefpcm',
        'cpcm', 'smd', 'dipole', 'ipcm', 'scipcm'].
        """
        # update only if not None; do not update to default value of None
        if solvent_model is not None:
            self._check_solvent(solvent_model)
            self.solvent_model = solvent_model

        if solvent_id is not None:
            self.solvent_id = solvent_id

    def _check_solvent(self, solvent_model):
        pass

    def modify_solvent(self, remove_solvent=False, **kwargs):
        if not remove_solvent:
            self.update_solvent(**kwargs)
        else:
            self.remove_solvent()

    def set_custom_solvent_via_file(self, filename):
        if not os.path.exists(os.path.expanduser(filename)):
            raise ValueError(f"File {filename} does not exist!")

        # path to the file for custom_solvent parameters
        with open(filename) as f:
            lines = f.readlines()

        lines = [line.strip() for line in lines]

        self.custom_solvent = "\n".join(lines)

    @classmethod
    def from_dict(cls, settings_dict):
        return cls(**settings_dict)

    @classmethod
    def from_database(
        cls,
        filepath,
        record_index=None,
        record_id=None,
        structure_index="-1",
        structure_id=None,
    ):
        """Create job settings from a chemsmart database file.

        With record selectors (record_index/record_id), this reconstructs
        source metadata from the selected record and charge/multiplicity from
        the selected structure.
        With a global structure selector (structure_id), this uses defaults
        and fills only charge/multiplicity from the selected structure.
        """
        from chemsmart.database.database import Database
        from chemsmart.database.utils import resolve_record
        from chemsmart.utils.utils import string2index_1based

        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"Database file not found: {filepath}")

        db = Database(filepath)
        record_selected = record_index is not None or record_id is not None
        if structure_id is not None and record_selected:
            raise ValueError(
                "Use either structure_id or record_index/record_id, not both."
            )

        settings = cls.default()

        if structure_id is not None:
            full_sid = db.get_structure_by_partial_id(structure_id)
            structure = db.get_structure(full_sid)
            if structure is None:
                raise ValueError(
                    f"No structure found with ID '{structure_id}'."
                )
            settings.charge = structure.get("charge")
            settings.multiplicity = structure.get("multiplicity")
            settings.title = (
                "Job prepared from chemsmart database "
                f"{os.path.basename(filepath)}"
            )
            logger.info(
                "Created JobSettings from database: "
                f"charge={settings.charge}, "
                f"multiplicity={settings.multiplicity}"
            )
            return settings

        record = resolve_record(
            db,
            record_index=record_index,
            record_id=record_id,
            return_list=False,
        )
        if record is None:
            return None

        meta = record.get("meta", {})
        molecules = record.get("molecules", [])
        selected_index = string2index_1based(str(structure_index))
        if isinstance(selected_index, slice):
            raise ValueError(
                "Database-aware jobs support one structure at a time."
            )
        try:
            structure = molecules[selected_index] if molecules else {}
        except IndexError as exc:
            raise ValueError(
                f"Structure index {structure_index} out of range for "
                "selected record."
            ) from exc
        settings.charge = structure.get("charge")
        settings.multiplicity = structure.get("multiplicity")
        settings.functional = meta.get("method")
        settings.basis = meta.get("basis")
        settings.jobtype = meta.get("jobtype")
        settings.solvent_model = meta.get("solvent_model")
        settings.solvent_id = meta.get("solvent_id")
        settings.custom_solvent = meta.get("custom_solvent")
        settings.title = (
            "Job prepared from chemsmart database "
            f"{os.path.basename(filepath)}"
        )
        logger.info(
            "Created JobSettings from database: "
            f"charge={settings.charge}, "
            f"multiplicity={settings.multiplicity}"
        )
        return settings


def read_molecular_job_yaml(filename, program="gaussian"):
    # read in defaults, if exists
    file_directory = os.path.dirname(filename)
    default_file = os.path.join(file_directory, "defaults.yaml")
    default_config = {}
    if os.path.exists(default_file):
        with open(default_file) as f:
            default_config = yaml.safe_load(f)
        logger.debug(
            f"Using the following pre-set defaults: \n{default_config}"
        )
    else:
        logger.warning("Default file settings does not exist.\n")
        if program == "gaussian":
            from chemsmart.settings.gaussian import GaussianJobSettings

            default_config = GaussianJobSettings.default().__dict__
        elif program == "orca":
            from chemsmart.settings.orca import ORCAJobSettings

            default_config = ORCAJobSettings.default().__dict__
        else:
            # other programs may be implemented in future
            pass
        logger.debug(
            f"Using the following pre-set defaults: \n{default_config}"
        )

    # job types
    gas_phase_jobs = list(MOLECULAR_GAS_PHASE_JOB_SECTIONS)
    sp_job = ["sp"]
    td_job = ["td"]
    link_job = ["link"] if program == "gaussian" else []
    qmmm_job = ["qmmm"]
    all_jobs = gas_phase_jobs + sp_job + td_job

    def stage_defaults(job):
        """Lift a project stage into the settings class that owns its keys."""

        if program == "gaussian" and job == "td":
            from chemsmart.jobs.gaussian.settings import (
                GaussianTDDFTJobSettings,
            )

            return GaussianTDDFTJobSettings(**default_config).__dict__.copy()
        if program == "gaussian" and job == "link":
            from chemsmart.jobs.gaussian.settings import (
                GaussianLinkJobSettings,
            )

            return GaussianLinkJobSettings(**default_config).__dict__.copy()
        return default_config.copy()

    # read in project config
    with open(filename) as f:
        project_config = yaml.safe_load(f)
        logger.debug(
            f"Project settings from yaml {filename}: \n{project_config}"
        )
    if not isinstance(project_config, dict):
        raise TypeError(f"{program} project YAML root must be a mapping")
    accepted_sections = set(molecular_project_section_names(program))
    unknown_sections = sorted(
        set(project_config).difference(accepted_sections)
    )
    unknown_stage_sections = sorted(
        name for name in unknown_sections if name not in default_config
    )
    if unknown_stage_sections:
        raise ValueError(
            f"Unknown {program} project section(s): "
            + ", ".join(unknown_stage_sections)
        )

    # populate job settings for different jobs
    all_project_configs = {}  # store all job settings in a dict

    # check if solv settings exist
    solv_config = project_config.get("solv", None)

    # check if separate gas phase settings exist
    gas_config = project_config.get("gas", None)
    qmmm_config = project_config.get("qmmm", None)

    direct_stage_present = any(
        project_config.get(job) is not None
        for job in all_jobs + qmmm_job + link_job
    )
    if gas_config is None and solv_config is None and not direct_stage_present:
        # Neither a phase section nor a direct calculation stage is present,
        # so there is nothing to merge.  A standalone ``td:`` project is a
        # complete fixed-geometry calculation and must not need a dummy
        # ``gas:`` section merely to satisfy the historical phase layout.
        found = (
            ", ".join(sorted(str(key) for key in project_config))
            if isinstance(project_config, dict) and project_config
            else "nothing"
        )
        raise ValueError(
            f"{program} project settings in {filename} define neither a "
            "'gas'/'solv' section nor a supported direct job section; "
            f"found: {found}."
        )
    if unknown_sections:
        raise ValueError(
            f"Unknown {program} project section(s): "
            + ", ".join(unknown_sections)
        )

    if gas_config is None:
        # no settings for gas phase; using
        # implicit solvation model for all jobs
        # (except td and qmmm, which will use their own configurations)
        for job in all_jobs:
            sources = molecular_project_section_sources(
                project_config, program=program, jobtype=job
            )
            phase_config = next(
                (
                    project_config[name]
                    for name in sources
                    if name in {"gas", "solv"}
                ),
                {},
            )
            all_project_configs[job] = (
                default_config.copy()
            )  # populate defaults
            if job == "sp":
                # A single-point project must not acquire the route-building
                # frequency job merely because only the solv phase was
                # declared. An explicit ``freq`` in the section still
                # overrides this below.
                all_project_configs[job]["freq"] = False
            all_project_configs[job]["jobtype"] = job  # update jobtype
            all_project_configs[job] = update_dict_with_existing_keys(
                all_project_configs[job], phase_config
            )
            if job == "td":
                # A td stage is a vertical spectrum at the supplied
                # geometry, as the td: branch below says; the solv phase
                # is borrowed here for its level, and neither the shared
                # default's ``freq: true`` nor that phase's own frequency
                # flag makes it a frequency job. Gaussian wrote
                # ``freq TD(...)`` -- an excited-state frequency
                # calculation -- for a solv-only project, where ORCA's td
                # settings refuse ``freq`` outright.
                all_project_configs[job]["freq"] = False
    else:
        # settings for gas phase exist - also solv settings exist
        for job in gas_phase_jobs:  # jobs using gas config
            sources = molecular_project_section_sources(
                project_config, program=program, jobtype=job
            )
            phase_config = next(
                (
                    project_config[name]
                    for name in sources
                    if name in {"gas", "solv"}
                ),
                {},
            )
            all_project_configs[job] = (
                default_config.copy()
            )  # populate defaults
            all_project_configs[job]["jobtype"] = job  # update jobtype
            all_project_configs[job] = update_dict_with_existing_keys(
                all_project_configs[job], phase_config
            )
            try:
                # Try updating with gas_config first
                all_project_configs[job] = update_dict_with_existing_keys(
                    all_project_configs[job], phase_config
                )
            except Exception as e:
                logger.warning(
                    f"Updating job '{job}' with gas_config failed ({e}). "
                    f"Falling back to qmmm_config."
                )
                # Fallback: try updating with qmmm_config
                all_project_configs[job] = update_dict_with_existing_keys(
                    all_project_configs[job], qmmm_config
                )
        # ``sp`` reads ``solv:`` because the canonical workflow optimises in
        # gas phase and takes the single point in solvent.  When a project
        # describes *only* gas phase, the phase that single point belongs to is
        # not ambiguous -- it is gas phase, and it inherits ``gas:``.  Without
        # this the fallback is one-sided: a ``solv:``-only project feeds every
        # job type, while a ``gas:``-only project feeds ``sp`` nothing, so the
        # settings reach the writer carrying no level of theory and the run
        # dies with "neither ab initio nor DFT is specified" and a zero-byte
        # input.  No working project can depend on that, because a
        # ``gas:``-only ``sp`` is an unconditional failure today.
        sp_sources = molecular_project_section_sources(
            project_config, program=program, jobtype="sp"
        )
        sp_phase = next(
            (name for name in sp_sources if name in {"gas", "solv"}),
            None,
        )
        sp_inherits_gas = sp_phase == "gas"
        sp_config = project_config.get(sp_phase, {})
        for job in sp_job:  # jobs using solv config
            all_project_configs[job] = (
                default_config.copy()
            )  # populate defaults
            # turn off freq calculation for single point calculations
            all_project_configs[job]["freq"] = False
            all_project_configs[job]["jobtype"] = job  # update jobtype
            all_project_configs[job] = update_dict_with_existing_keys(
                all_project_configs[job], sp_config or {}
            )
            if sp_inherits_gas:
                # An explicit ``freq:`` under ``solv:`` overrides sp's default,
                # because ``solv:`` is sp's own section and an author writing
                # it there means it.  ``gas:`` is not sp's section -- it is
                # being borrowed for its level of theory -- and it carries
                # ``freq: true`` in most projects because it describes the
                # optimisation.  Inheriting that would turn a requested single
                # point into an opt+freq, so freq stays off on this path.
                all_project_configs[job]["freq"] = False

    if program == "orca" and "neb" in all_project_configs:
        # ``gas:`` supplies the electronic-structure settings shared by the
        # path calculation, but a phase-level ``freq: true`` belongs to the
        # optimization described by that phase. A NEB job should not silently
        # acquire it from the borrowed phase section. A chemist can still
        # request one explicitly under the job's own ``neb:`` section.
        all_project_configs["neb"]["freq"] = False

    # check if td settings exist (optional)
    if "td" in project_config:
        td_config = project_config["td"]
        for job in td_job:  # jobs using td config s
            all_project_configs[job] = stage_defaults(job)
            if program in {"gaussian", "orca"}:
                # A ``td`` stage is a vertical spectrum at the supplied
                # geometry. The shared molecular-settings default predates
                # that stage and requests a frequency calculation. Do not
                # silently add one; an explicit ``td.freq`` below still wins.
                all_project_configs[job]["freq"] = False
            all_project_configs[job]["jobtype"] = job  # update jobtype
            all_project_configs[job] = update_dict_with_existing_keys(
                all_project_configs[job], td_config
            )

    # check if qmmm settings exist (optional)
    if "qmmm" in project_config:
        qmmm_config = project_config["qmmm"]
        for job in qmmm_job:  # jobs using qmmm config
            all_project_configs[job] = (
                default_config.copy()
            )  # populate defaults
            all_project_configs[job]["jobtype"] = job  # update jobtype
            logger.debug(
                f"Updating qmmm job settings: {all_project_configs[job]} with {qmmm_config}"
            )
            for k, v in qmmm_config.items():
                logger.debug(f"Updating qmmm job settings: {k} with {v}")
                all_project_configs[job][k] = v

    # Canonical ChemSmart projects may describe a stage directly (``opt:``,
    # ``sp:``, and so on) instead of routing every job through the historical
    # ``gas:``/``solv:`` pair.  Build the legacy defaults first, then let an
    # explicit stage override only that stage.  This keeps old projects
    # readable while allowing a model to express a paper workflow without
    # inventing settings for unrelated jobs.
    for job in all_jobs + qmmm_job + link_job:
        stage_config = project_config.get(job)
        if stage_config is None:
            continue
        if not isinstance(stage_config, dict):
            raise TypeError(f"Project section `{job}` must be a mapping")
        if job not in all_project_configs:
            all_project_configs[job] = stage_defaults(job)
            all_project_configs[job]["jobtype"] = job
        if program == "orca" and job in {"irc", "neb", "ts"}:
            # The shared ORCA defaults describe one-geometry jobs and do not
            # contain the settings owned by path calculations.  Lift an
            # explicit path stage into its real settings class before
            # applying it.  Otherwise ``irc: {direction: both}`` and
            # ``neb: {nimages: ...}`` are rejected as unknown keys even though
            # their CLI jobs and native writers support them.
            #
            # ``ts`` belongs here for the same reason and was missing. The
            # model-visible tool surface instructs that a transition-state
            # search seeded from a validated producer's Hessian "must set
            # inhess: true" in its project ts section -- and ``inhess`` lives
            # on ORCATSJobSettings, not on the shared defaults, so that
            # instruction raised `Keyword 'inhess' is not in list of
            # keywords`. The producer rule it belongs to is recorded in the
            # charter as admitted intent that has never executed; an
            # instruction the loader refuses is why.
            from chemsmart.jobs.orca.settings import (
                ORCAIRCJobSettings,
                ORCANEBJobSettings,
                ORCATSJobSettings,
            )

            settings_class = {
                "irc": ORCAIRCJobSettings,
                "neb": ORCANEBJobSettings,
                "ts": ORCATSJobSettings,
            }[job]
            all_project_configs[job] = settings_class(
                **all_project_configs[job]
            ).__dict__.copy()
        if program == "gaussian":
            # The same lift, for the same reason, on the program that
            # never got it. Gaussian's shared defaults describe a
            # one-geometry job, so every setting a path, a response or a
            # link job owns -- ``direction`` and ``maxpoints`` on an IRC,
            # ``nstates`` and ``root`` on a TD -- was refused as an
            # unknown key although the CLI takes it and the native writer
            # writes it. Observed live (CUHK r9g-g1, 2026-09-21): a
            # session asked for ``irc: {maxpoints: 50}`` to bound a
            # reaction path, was told the keyword was not in the list of
            # keywords, and fell back to ``additional_route_parameters:
            # maxpoints=50``, which appends a bare token beside
            # ``irc(...)`` instead of inside it -- a route keyword that
            # is not one.
            #
            # Which class owns which section is the Gaussian loader's own
            # declaration, read from there rather than written again.
            from chemsmart.settings.gaussian import (
                gaussian_jobtype_settings_classes,
            )

            gaussian_class = gaussian_jobtype_settings_classes().get(job)
            if gaussian_class is not None:
                all_project_configs[job] = gaussian_class(
                    **all_project_configs[job]
                ).__dict__.copy()
        all_project_configs[job] = update_dict_with_existing_keys(
            all_project_configs[job], stage_config
        )

    return all_project_configs
