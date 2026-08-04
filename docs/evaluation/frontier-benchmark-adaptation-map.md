# Frontier Benchmark Adaptation Map

## Purpose and decision vocabulary

This map records which benchmark constructs ChemSmart may reuse from three
computational-chemistry agent systems. It was source-checked on 2026-08-04. It
is an architecture and benchmark-design record, not evidence that ChemSmart has
run these benchmarks or matched the published systems.

- **Adopted**: the source construct can be retained with the same evaluation
  purpose.
- **Adapted**: the construct is retained only after translation into
  ChemSmart's typed project-YAML, command-DAG, evidence, and validator model.
- **Rejected**: the construct conflicts with ChemSmart's fixed authority
  boundary.
- **No analogue**: the source measures a capability outside the current
  ChemSmart surface; no substitute score is invented.

Except for a case whose four typed record digests are explicitly named, these
labels classify source constructs at the narrative design level only. They do
not instantiate `BenchmarkSourceLedgerV1`, `BenchmarkChallengeFactorsV1`,
`BenchmarkRubricProfileV1`, or `BenchmarkAdaptationRecordV1`, and they do not
authorize an episode.

Across every adaptation, a model proposes scientific intent. ChemSmart owns
project materialization, live-schema command compilation, generated program
scripts or inputs, preview, execution authorization, validation, and terminal
state. Model-authored native inputs, runtime shell authority, and runtime tool
forging are never imported as benchmark requirements.

## Evidence inventory

| System | Evidence status | Exact locators | Public implementation status |
|---|---|---|---|
| El Agente Gráfico | Preprint: [arXiv 2602.17902v1](https://arxiv.org/abs/2602.17902v1), PDF SHA-256 `789dab02298e5ec0e1a1a9be8a3684efd204bcb8ab1f5a1e57b4002194769767`; [AI4X-AC 2026 oral abstract](https://openreview.net/forum?id=YShauPAyTw), CC BY 4.0 | [Chat-transcript repository](https://github.com/jb2197/ElAgenteGrafico-ChatTranscript/tree/d2bb273875288254d6d3f5ac984b312063cbacc3), commit `d2bb273875288254d6d3f5ac984b312063cbacc3` | The repository identifies itself as chat transcripts and includes case-study files. It is not treated as a verified reusable release of the Gráfico framework. No separate replication dataset was verified. |
| El Agente Q | Peer-reviewed Matter article: [DOI 10.1016/j.matt.2025.102263](https://doi.org/10.1016/j.matt.2025.102263); [OSTI full text](https://www.osti.gov/pages/servlets/purl/3010680), retrieved PDF SHA-256 `3ce71f70a77dddefc0c054ae7b4f13fc45f53a1f3eb708a362d65c0020dc874f` | [Replication data DOI](https://doi.org/10.5683/SP3/JU2BQK), [released dataset version 1.1 metadata](https://borealisdata.ca/api/datasets/:persistentId/versions/1.1?persistentId=doi:10.5683/SP3/JU2BQK), CC BY-NC-SA 4.0 | The article and trace/artifact package are benchmarkable. No verified official reusable source-code repository was found; none is inferred from the data package. |
| El Agente Forjador | Preprint: [arXiv 2604.14609v1](https://arxiv.org/abs/2604.14609v1), CC BY 4.0, PDF SHA-256 `adb460a944ad9f6bff02be611efb2a4b82b4d1d8ad78a1131c843282813cd0b2` | [Replication data DOI](https://doi.org/10.5683/SP3/0YOMKL), [released dataset version 1.0 metadata](https://borealisdata.ca/api/datasets/:persistentId/versions/1.0?persistentId=doi:10.5683/SP3/0YOMKL), CC0 1.0; ChemSmart's sealed [Q03 preregistration](frontier-benchmark-preregistration-v1.md) | The public workspaces, questions, scorers, ground truths, and methodology rubrics are usable under their recorded roles. No separate official framework source repository was verified. |

## El Agente Gráfico — narrative mapping

### Task and scoring structure

The core benchmark contains six university-level quantum-chemistry families at
two difficulty levels: organic and inorganic molecular analysis,
hydrogen-abstraction/carbocation energetics, cycloalkane ring strain, pKa of
halogenated acetic acids, and TDDFT excited-state energies. Each task was run
ten times, for 120 runs per evaluated model. Conformer-ensemble spectroscopy
and metal-organic-framework design are extension studies, not members of the
core benchmark.

The source uses a deterministic evaluator over structured results and an LLM
judge for completeness, reasoning, and reporting. Numerical dimensions include
functional/basis/charge/multiplicity integrity, energy tolerance of 0.01 Ha,
geometry RMSD tolerance of 0.15 angstrom, HOMO--LUMO-gap tolerance of 0.1 Ha,
and requested dipole, symmetry, and frequency checks. It also records tokens,
cost, API requests, task duration, context saturation, error-recovery cost, and
carryover tokens. The source pass rule combines a perfect numerical score with
an LLM-judge score above 0.90; ChemSmart does not import that aggregation rule.

This is a narrative summary of the reported protocol, not a frozen scorer. No
Gráfico typed benchmark record exists until exact prompt/input, scorer,
ground-truth/methodology and rubric artifacts, criterion weights/tolerances,
licenses, and content hashes are bound.

| Source construct | Narrative disposition | ChemSmart analogue or boundary |
|---|---|---|
| Typed execution graph, typed state, and symbolic object references | Adopted | Typed scientific workflow state, immutable artifact bindings, and Runtime V2 events. |
| Persistent knowledge graph | Adapted | Canonical artifacts and replayable receipts externalize state; a knowledge store cannot override project, CLI, or validator authority. |
| Deterministic chemical-property dimensions and trace-efficiency metrics | Adopted | Report scientific correctness separately from provider tokens, latency, requests, repair cost, and context carryover. |
| LLM evaluation of completeness and narrative quality | Adapted | Supplementary evidence only; deterministic graders remain primary and source-comparable and ChemSmart-strict scores remain separate. |
| General code-execution tool and any native-input or shell-mediated path | Rejected | Project YAML plus the live-schema command compiler is the only executable preparation surface. |
| MOF design and conformer-ensemble extensions | No analogue | Outside this PySCF factorial campaign; no proxy score is assigned. |

## El Agente Q — narrative mapping

### Task and scoring structure

The six task families are organic analysis, inorganic analysis, carbocation
formation enthalpy/free energy, cycloalkane ring strain, aqueous pKa, and
TDDFT absorption. Difficulty is task-specific rather than a single monotonic
scale: Level 1 may supply fewer systems or more procedural hints, while Level 2
may supply larger XYZ batches, omit multiplicity or methodology hints, or
require calibration. Each identical prompt was run ten times. A run receives a
task-specific value from 0 to 100; the paper reports the mean and standard
deviation across the ten runs. Exact criterion weights are not imported unless
they are bound to a verified source-rubric artifact.

This is likewise a narrative source mapping. The replication package makes
later instantiation possible, but no Q typed record exists until the selected
task prompt, coordinates, scorer/rubric and ground-truth artifacts, exact
criterion weights/units, and content hashes are pinned.

| Source construct | Narrative disposition | ChemSmart analogue or boundary |
|---|---|---|
| Six chemistry task families and their per-task difficulty ladders | Adopted | Eligible for future source-comparable scoring only after a typed ledger binds the selected prompt, inputs, rubric/scorer artifacts, weights, units, and hashes. |
| Hierarchical expert roles and concise upward reports | Adapted | Complexity-gated, bounded specialists receive immutable packets; one coordinator owns molecular identity, final YAML, DAG integration, and readiness. |
| Task-specific rubric attainment over repeated runs | Adapted | Keep the original-comparable score distinct from ChemSmart's all-required identity, project, command, preview, and evidence gates. |
| Exported observable action traces and documented failure modes | Adopted | Persist visible actions, tool results, artifacts, repairs, and outcomes; provider-private reasoning is not evidence. |
| Agent-authored ORCA input, shell/SLURM control, and automatic resubmission or geometry mutation | Rejected | ChemSmart generates inputs and argv; execution and any changed scientific workflow require host validation and the applicable approval. |
| Success at freely authoring native syntax or using shell commands | No analogue | These are intentionally absent capabilities, so ChemSmart receives neither credit nor penalty for them. |
| SMILES-to-3D autonomy in tasks without approved exact coordinates | No analogue | Such cases remain ineligible or blocked until a separately approved deterministic coordinate contract exists. |

## El Agente Forjador

### Task and scoring structure

Forjador evaluates 24 tasks: 13 quantum-chemistry tasks and 11
quantum-dynamics tasks. Its quantum-chemistry portion uses six families shared
with the preceding benchmarks, with organic Levels 0--2 and Levels 1--2 for
inorganic analysis, carbocations, ring strain, pKa, and TDDFT. The study
compares zero-shot tool generation, curriculum-built tool reuse, and an
evaluator-only baseline. Each task is graded on accuracy and methodology, and
the source score is their mean; API cost and wall time are reported separately.
The versioned CC0 package supplies task prompts, deterministic scorers, ground
truths, and methodology rubrics.

| Source construct | Decision | ChemSmart analogue or boundary |
|---|---|---|
| Task/difficulty taxonomy and external deterministic scorers | Adopted | Preserve source tasks, scorer artifacts, units, and hidden-answer boundaries under a source ledger. |
| Accuracy and methodology channels | Adapted | Report both source channels separately, then add an independent ChemSmart-strict channel; do not collapse normalized or narrative success into strict success. |
| Tool analysis and reuse through progressive disclosure | Adapted | Retrieve existing live CLI capabilities, project schemas, validators, and scoped knowledge packs; composition produces a typed project-backed command DAG. |
| Runtime generation, editing, or execution of tools and scripts | Rejected | New capability is reviewed and integrated outside a benchmark episode; the model cannot forge its runtime surface. |
| Tool-library growth and reuse benefit | No analogue | ChemSmart does not dynamically grow a tool library. Provider cost, token use, and latency remain descriptive efficiency metrics. |
| Eleven quantum-dynamics tasks | No analogue | They are outside the current ChemSmart computational-program surface and receive no substitute score. |

`FORJADOR-Q03-ORGANIC-L2-V1` is already sealed in the linked preregistration:
six supplied XYZ systems, mixed neutral/anionic/cationic states, gas-phase
HF/def2-SVP, per-system `OPT -> HESS`, and property reporting. It remains held
out and is not authorized for provider or engine execution in the present
factorial campaign.

## Cross-system reporting rule

Any future result must name the source benchmark version, source-comparable
rubric, ChemSmart adaptation record, model/provider block, and execution state.
It must report source scores, ChemSmart-strict gates, and efficiency metrics as
separate channels. A pass on one adapted case cannot establish repeated
reliability, cross-domain generality, reproduction, provider superiority,
production readiness, or state-of-the-art performance.
