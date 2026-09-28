# Evidence graph

What this repository has learned, kept as links instead of prose. Nothing
under `chemsmart/` imports it.

| file | what it is |
|---|---|
| `graph.yaml` | the authored join: fundamentals, lessons pointing at the text that states them, negative results with their falsifiers, the commits that replaced what a sentence says, and per-rule annotations |
| `claims.yaml` | what the project may and may not yet claim, with the strongest counterexample and the missing experiment for each |
| `loop/graph.py` | the graph, derived on every call from the live registries (rules, gates, policies, guides, tests, charter sections and topics, and the capability registry's concepts) and joined to `graph.yaml`; nothing is stored, so it cannot go stale |
| `loop/census.py` | what each instruction surface costs, measured from the live tree |
| `loop/replay.py` | what archived Agent runs did, recomputed from their hash-chained event streams; `--transcripts ROOT` prints their public transcripts as readable turns |
| `loop/signed_words.py` | the classes of word the host signs, each with its one signer, the record it lands in, and how a census checks it |
| `loop/resign.py` | archived goals' settlement words re-signed on the imported tree |
| `loop/resign_producing.py` | each archived goal re-signed on the code that produced it, where that code still hashes to the digest its job printed |
| `loop/resign_stationarity.py` | archived stationarity words (a stationary-point order, a free energy) re-signed on the imported tree |
| `loop/word_reader.py` | an independent reader, importing nothing from `chemsmart`: is each signed word true of the records it cites? |
| `loop/refusal_census.py` | every refusal the Agent met, grouped by what refused it and how the goal ended |
| `loop/receipt_refusals.py` | whether a receipt the decision gate refused as "not one it minted" was in fact minted by the host |
| `loop/receipt_gate_replay.py` | every refusal the decision gate made, re-read through the imported tree's own gate functions |
| `loop/literal_claims.py` | every row a host-rendered analysis report shows, classed by whether it stands on a model-authored literal, and of what dimension |
| `loop/certificate_census.py` | every archived completion certificate that named claims on a failed criterion, re-read against the goal's other streams |
| `loop/composition_census.py` | every approved cycle's stages, programs and producer edges, deduplicated by bundle content, goal and plan digest; admitted vs realised edges and cross-approval lifts |
| `loop/task_code_audit.py` | whether a task-named operation, CLI or job type is a general conversion or task code, and whether the Agent can reach it |
| `loop/replay_composition.py` | whether an archived cycle's plan, producer rules and operations re-admit at a given tree |
| `loop/compare_words.py` | two re-signings of one population compared goal by goal: the word, the qualified rows it wrote, the executor's analysis word |
| `loop/literal_unique.py` | the literal census counted once per distinct report, since archived workspace copies carry the same rendered report more than once |
| `loop/matched_turns.py` | an archived goal prefix replayed through the tree `PYTHONPATH` names, then real provider turns, recording the model, the calls in view and the prefix's faithfulness |
| `loop/matched_outcomes.py` | matched-turn outcomes classified from the host's own records by rules fixed before the first sample |
| `loop/r11_behav/` | R11's two-model baseline and arm runners, analyses and blind classifiers, byte for byte as run |
| `loop/evidence_census.py` | where the Agent's work leaves the typed evidence layer (native words, unbound prose numbers, vocabulary refusals, reuse attempts), with denominators and a shuffled-session control |
| `loop/evidence_census_report.py` | the census's figures recomputed from its outputs |
| `loop/probe_mode_selection.py` | whether a plan written before a Hessian picks a mode by which atoms move, through the host's public tool surface, with each selection's runner-up |

```bash
PY=/opt/anaconda3/bin/python          # never pip install
$PY .agents/research/loop/graph.py find <word | rule id | commit>
$PY .agents/research/loop/graph.py why <node>
$PY .agents/research/loop/graph.py cost | orphans | stats | check
$PY .agents/research/loop/census.py
$PY .agents/research/loop/replay.py [ROOT ...]
```

A lesson or a result is added to `graph.yaml` with a pointer that `check` can
resolve; it is never restated in `AGENTS.md`. Edit `graph.yaml` as text: its
comment header is the file's own instructions, and loading and re-dumping the
file keeps every node and erases them. An anchor is a verbatim sentence and may
cross a line break in the prose it points at.

**Concepts arrive through their registries.** Two closed vocabularies are
nodes whole: `program_jobtype:<program>:<engine>:<jobtype>` with its ladder
rung (what can run, and how far it has been proven) and `signal:<id>` (what
the host may observe about a result and hand to a session). Open vocabularies
-- settings, selectors, operations, constants -- become nodes only when a
rule's text or one of its declared boundaries, a guide body, or an edge in
`graph.yaml` names them; the registry holds hundreds and the graph is not its
mirror (`CONCEPT_SOURCES` in `graph.py` is the whole policy). A sentence that
names a concept, or states its boundary, `explains` it: `why signal:<id>`
lists the sentences a session could have learned it from, and `orphans` lists
the signals no sentence explains and any live rule a commit supersedes.
