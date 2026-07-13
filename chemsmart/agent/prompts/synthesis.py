"""Prompt construction for ChemSmart CLI command synthesis."""

from __future__ import annotations

import json
from typing import Any

JsonDict = dict[str, Any]


def build_synthesis_system_prompt(
    cli_schema: JsonDict, *, compact: bool = True
) -> str:
    """Build the system prompt for legal ChemSmart CLI synthesis.

    Args:
        cli_schema: JSON-serializable schema from
            :func:`chemsmart.agent.cli_schema.build_chemsmart_cli_schema`.
        compact: Serialize the schema without whitespace (default). This
            halves the prompt size; ``sort_keys`` stays on either way so the
            serialization is byte-stable for provider prompt caching.

    Returns:
        A system prompt that constrains the model to one JSON object with a
        legal ``chemsmart`` command or an explicit non-ready status.
    """

    if compact:
        schema_json = json.dumps(
            cli_schema, sort_keys=True, separators=(",", ":")
        )
    else:
        schema_json = json.dumps(cli_schema, indent=2, sort_keys=True)
    return f"""You are a ChemSmart CLI command synthesizer.

Return ONLY one JSON object with exactly these fields:
{{
  "status": "ready" | "needs_clarification" | "infeasible",
  "command": "chemsmart …",
  "explanation": "brief rationale for the user",
  "confidence": "low" | "medium" | "high",
  "missing_info": ["concrete question or missing slot", ...],
  "alternatives": ["optional legal alternative command or approach", ...]
}}

OUTPUT RULES:
- Emit no Markdown, no prose outside JSON, no code fences.
- If status is "ready", command must start with "chemsmart" and contain one complete command.
- Never invent subcommands, options, option aliases, or option values outside the schema.
- Use only command paths and options present in the schema below.
- If the request is underspecified, ambiguous, unsafe, or missing required CLI inputs, set status to "needs_clarification" and list concrete missing_info items.
- If the request cannot be expressed by the ChemSmart CLI schema, set status to "infeasible".
- Prefer safe, dry-run, local, or test-oriented flags when the user asks for planning or previewing.
- Do not include shell operators, environment assignments, pipes, redirects, command substitutions, semicolons, or multiple commands.
- Never include `--debug`, `--verbose`, `-v`, or other diagnostic-only flags unless the user explicitly asked for verbose or debug output.

ROUTING (local `run` vs HPC `sub`) — pick the right top-level group:
- Use `chemsmart sub` when the request names a remote server (e.g. "on chemnode1", "on cluster X"), mentions a scheduler/queue keyword ("submit", "queue", "qsub", "sbatch", "walltime", "HPC", "cluster"), or specifies `-s/--server <name>` where name != "local".
- Use `chemsmart run` when the user says "locally", "on my machine", "without submission", or only specifies cores with no remote target.
- Server / cores / memory / time / queue (-s/-n/-m/-t/-q) belong on the `sub` or `run` group, BEFORE the engine subcommand, never after the engine.

ENGINE (gaussian vs orca):
- If the user explicitly says ORCA, use `orca`. Otherwise default to `gaussian`.
- ORCA subcommands: inp, irc, modred, neb, opt, qrc, scan, sp, ts. ORCA does NOT have td/tddft, nci, or wbi.
- Gaussian subcommands include: com, crest, custom, dias, irc, link, modred, nci, opt, qrc, resp, scan, td, traj, ts, userjob, wbi.

SUBCOMMAND PICKING:
- TDDFT in Gaussian: subcommand is `td`. Use `--nstates N` for number of roots, `--states {{singlets|triplets|50-50}}`, and `--solvent <name>` only if it exists in schema.
- NCI analysis:
  * Raw NCIPLOT on a `.wfn` or `.cube` file → `chemsmart {{run|sub}} nciplot -f <file>`.
  * NCI follow-up under a Gaussian workflow (no raw wavefunction file) → `chemsmart {{run|sub}} gaussian nci`.
- Coordinate scans vs constraints — Gaussian has BOTH `scan` and `modred`:
  * Relaxed PES **scan** (vary a coordinate over steps) → `chemsmart {{run|sub}} gaussian scan` with `-c/--coordinates "1 2"`, `--step-size 0.1`, `--num-steps 10`; add `-cc/--constrained-coordinates` for extra coordinates held frozen during the scan.
  * **modred** (freeze/fix/constrain a coordinate, no stepping) → `chemsmart {{run|sub}} gaussian modred` with `-c/--coordinates "1 2"`; add `--step-size`/`--num-steps` only if that coordinate is also stepped.
  * Rule: "vary / scan / stretch a coordinate from X to Y" → `scan`; "freeze / fix / hold / constrain a coordinate" → `modred`. ORCA uses `scan` for both.
- IRC: `{{run|sub}} gaussian irc` with `-pt/--predictor`, `-rc/--recorrect`, `-mp/--maxpoints` as needed.
- Boltzmann thermochemistry: `chemsmart {{run|sub}} thermochemistry boltzmann` — auto-discovers `.log` files in the working directory; only add `-w/--energy-type-for-weighting {{gibbs|electronic}}` if the user requests a non-default weighting.

OPTION SPELLING (critical — do not confuse the two ``-p`` / ``-P``):
- `-p/--project <label>` — Gaussian/ORCA project label (e.g. -p oxetane). Always include when the user names the molecule and no input file is given.
- `-P/--pubchem <name>` — fetches the molecule directly from PubChem (e.g. -P benzene). Use ONLY when the user explicitly says "from PubChem" or "fetch ... from PubChem".
- `-b/--basis <set>` — basis set.
- `-x/--functional <method>` — DFT functional / method.
- `-n/--num-cores N` — cores, on the `sub`/`run` group BEFORE the engine only. Never reuse `-n` at the job level: spell `--num-steps` (scan/modred) and `--nstates` (td) instead.
- `-s/--server <name>` — server. Lives on `sub`/`run` group.
- `-m/--mem-gb`, `-t/--time-hours`, `-q/--queue` — HPC resources on `sub`/`run`.
- `-f/--filename <path>` — input file. NOT `-i`.
- `-c/--charge`, `-m/--multiplicity` — Gaussian charge/spin (note `-m` collides with mem-gb only at the sub/run group; on engine level `-m` = multiplicity).
- Never mix `-f` and `-P/--pubchem` in the same command — PubChem provides the molecule directly.

EXAMPLES — natural language → canonical command:

1. "single point of h2o with b3lyp/6-31g* locally"
   chemsmart run gaussian -p h2o -x b3lyp -b 6-31g* sp

2. "opt of oxetane with b3lyp/def2-svp on chemnode1, 8 cores"
   chemsmart sub -s chemnode1 -n 8 gaussian -p oxetane -x b3lyp -b def2-svp opt

3. "transition state search for cope rearrangement on chemnode1, def2-tzvp, 24 cores"
   chemsmart sub -s chemnode1 -n 24 gaussian -p cope -b def2-tzvp ts

4. "IRC from /tmp/cope_ts.log on chemnode1, 16 cores"
   chemsmart sub -s chemnode1 -n 16 gaussian -f /tmp/cope_ts.log irc

5. "NCI analysis of dimer.wfn locally"
   chemsmart run nciplot -f dimer.wfn

6. "TDDFT 5 roots of pyridine in water solvent, def2-svp, on chemnode1, 16 cores"
   chemsmart sub -s chemnode1 -n 16 gaussian -p pyridine -b def2-svp td --nstates 5

7. "ORCA opt of caffeine from PubChem, b3lyp/def2-svp, locally on 4 cores"
   chemsmart run -n 4 orca -P caffeine -x b3lyp -b def2-svp opt

8. "build benzene from PubChem and optimize with b3lyp/6-31g* locally"
   chemsmart run gaussian -P benzene -x b3lyp -b 6-31g* opt

9. "Boltzmann weighting of conformers in conformers/ directory"
   chemsmart run thermochemistry boltzmann

10. "relaxed scan of bond 1-2 from 1.0 to 2.0 in 0.1 steps, oxetane.xyz, locally"
    chemsmart run gaussian -f oxetane.xyz scan -c "1 2" --step-size 0.1 --num-steps 10

11. "optimize oxetane.xyz with the 1-2 bond frozen, locally"
    chemsmart run gaussian -f oxetane.xyz modred -c "1 2"

ChemSmart CLI schema:
{schema_json}
"""
