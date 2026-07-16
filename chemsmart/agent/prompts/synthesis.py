"""Budgeted prompt construction for ChemSmart CLI command synthesis."""

from __future__ import annotations

import json
from typing import Any

from chemsmart.agent.compact_cli_schema import (
    compact_cli_signature,
    compact_signature_paths,
)

JsonDict = dict[str, Any]
SYNTHESIS_PROMPT_MAX_CHARS = 8192

_BASE_POLICY = """You synthesize one legal ChemSmart CLI command.

Return ONLY one JSON object:
{"status":"ready|needs_clarification|infeasible","command":"chemsmart ...|","explanation":"brief user rationale","confidence":"low|medium|high","missing_info":[],"alternatives":[]}

Rules:
- Use only paths/flags/values in the compact signature. Never emit shell operators, redirects, environment assignments, multiple commands, debug, or verbose flags.
- Signature args use `flags[:type][:{choices}][! required][[] repeatable]`; omitted type means a string value.
- `ready` requires one complete command starting with `chemsmart`; ask only for a truly missing required slot; use `infeasible` when the CLI cannot express the request.
- Remote server, submit, queue, scheduler, walltime, or HPC intent -> `chemsmart sub`; explicit local intent -> `chemsmart run`. Resource flags belong before the program.
- Explicit ORCA -> `orca`; explicit Gaussian/G16 or otherwise workspace Gaussian -> `gaussian`. Preserve the user's program, kind, path, project, server, charge, multiplicity, atom indices, and numeric settings.
- `-p/--project` selects project YAML; `-P/--pubchem` is only for an explicit PubChem fetch. Never combine `-f` and `-P`.
- `-f/--filename` is the input path. Program-level `-m` is multiplicity; run/sub-level `-m` is memory. Use long `--num-steps` and `--nstates` to avoid `-n` ambiguity.
- Runtime project YAML owns method/basis/solvent unless the user explicitly overrides them. Do not invent scientific values.
- Prefer validation/fake/test flags only when the user requested preview or testing; do not silently change requested execution semantics.
"""

_PACKS = {
    "td": """Gaussian TD-DFT uses `gaussian ... td` with `--states {singlets|triplets|50-50}`, `--nstates N`, `--root N`, and `--eqsolv {eqsolv|noneqsolv}`. Do not invent `--singlets` or `--solvent`; the selected project must have a `td` block.
Example: `chemsmart run gaussian -p pyridine -f pyridine.xyz -c 0 -m 1 td --states singlets --nstates 5 --root 1`. ORCA has no td subcommand.""",
    "scan_modred": """Preserve scan versus freeze intent. Gaussian relaxed scan: `gaussian ... scan --coordinates '[[1,2]]' --step-size '[0.1]' --num-steps '[10]'`; add `--constrained-coordinates '[[1,3]]'` for frozen coordinates. ORCA scan uses `--coordinates '[1,2]' --dist-start 1.0 --dist-end 2.0 --num-steps 10`. `gaussian|orca modred -c '[[1,2]]'` freezes at the input geometry with no stepping. Vary/range -> scan; freeze/fix/hold -> modred.
Example: `chemsmart run gaussian -p demo -f mol.xyz -c 0 -m 1 scan --coordinates '[[1,2]]' --step-size '[0.1]' --num-steps '[10]'`.""",
    "nci": """Raw NCIPLOT for `.wfn`/`.cube` -> `run|sub nciplot -f FILE`; NCI in a Gaussian workflow -> `run|sub gaussian ... nci`.
Example: `chemsmart run nciplot -f dimer.wfn`.""",
    "irc_qrc": """Do not confuse IRC and QRC. IRC follows the reaction path; QRC displaces along an imaginary mode. Preserve requested direction, predictor/recorrect, max points, source checkpoint/log, and mode settings using only exposed flags.
Example: `chemsmart run gaussian -p demo -f ts.log -c 0 -m 1 irc`.""",
    "ts": """Use the program's `ts` subcommand, not plain opt. Preserve requested Hessian strategy, maxstep, charge, multiplicity, and input path; do not duplicate TS route keywords in additional route parameters.
Example: `chemsmart sub -s chemnode1 -n 16 gaussian -p demo -f guess.xyz -c 0 -m 1 ts`.""",
    "qmmm": """QMMM atom regions are scientific intent. Preserve exact 1-based QM/high/low atom sets, charge, multiplicity, and low-level method. Never infer a missing partition or convert ranges to whitespace-only lists.
Example: `chemsmart run orca -p demo -f complex.xyz -c 0 -m 1 opt qmmm --high-level-atoms 1-8`.""",
    "neb": """ORCA NEB requires distinct reactant/product endpoints. Preserve endpoint ownership, `--nimages`, and requested NEB/NEB-CI/NEB-TS job option; ask rather than invent a missing endpoint.
Example: `chemsmart run orca -p demo neb -r reactant.xyz -p product.xyz --nimages 8`.""",
    "dias_wbi": """DIAS requires fragment definitions and fragment states; WBI is bond-index analysis. Do not substitute one for the other. Preserve flat fragment-1 indices and let runtime derive the complement when that is the CLI contract.""",
    "crest": """CREST/conformer requests must preserve the requested conformer count or selection strategy and job type. External CREST protocol context is not a Gaussian method override.""",
}


def build_synthesis_system_prompt(
    cli_schema: JsonDict,
    *,
    compact: bool = True,
    max_chars: int = SYNTHESIS_PROMPT_MAX_CHARS,
) -> str:
    """Build a request-scoped prompt and fail closed on budget overflow."""

    signature = compact_cli_signature(cli_schema)
    paths = compact_signature_paths(signature)
    packs = _context_packs(paths)
    schema_json = json.dumps(
        signature,
        sort_keys=True,
        separators=(",", ":") if compact else None,
        indent=None if compact else 2,
    )
    sections = [_BASE_POLICY.strip(), *packs, f"Compact CLI signature:\n{schema_json}"]
    prompt = "\n\n".join(section for section in sections if section)
    if len(prompt) > max_chars:
        raise ValueError(
            f"Synthesis prompt is {len(prompt)} chars; budget is {max_chars}. "
            "Use request-scoped schema pruning; full-schema fallback is review-only."
        )
    return prompt


def _context_packs(paths: set[str]) -> list[str]:
    tokens = {token for path in paths for token in path.split()}
    selected: list[str] = []
    if "td" in tokens:
        selected.append(_PACKS["td"])
    if {"scan", "modred"}.intersection(tokens):
        selected.append(_PACKS["scan_modred"])
    if "nci" in tokens or "nciplot" in tokens:
        selected.append(_PACKS["nci"])
    if {"irc", "qrc"}.intersection(tokens):
        selected.append(_PACKS["irc_qrc"])
    if "ts" in tokens:
        selected.append(_PACKS["ts"])
    if "qmmm" in tokens:
        selected.append(_PACKS["qmmm"])
    if "neb" in tokens:
        selected.append(_PACKS["neb"])
    if {"dias", "wbi"}.intersection(tokens):
        selected.append(_PACKS["dias_wbi"])
    if {"crest", "traj"}.intersection(tokens):
        selected.append(_PACKS["crest"])
    return selected
