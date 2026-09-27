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
