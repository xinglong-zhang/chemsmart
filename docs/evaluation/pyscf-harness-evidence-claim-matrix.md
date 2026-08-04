# PySCF Harness Evidence Claim Matrix

## Purpose and evidence boundary

This matrix is the claim baseline before the recovered provider campaign. It
separates code that exists from behavior checked by deterministic tests and
behavior observed from a live model. It does not infer scientific correctness,
factorial effects, or generality from implementation size, test presence, a
terminal label, or a provider transport response.

The committed baseline is `7c35e06487b0d62f67733dcbb167837f65e18214`.
The four baseline commits are:

- `c986180d2b3234aea8bb9b0d068f09362749615d` -- PySCF and
  GPU4PySCF program workflows;
- `728b1908239895311429271adc60879197dfd28e` -- scientific DAG and
  factorial harness;
- `40f1f8c15162636b61edb0dee4d89dc374d59fed` -- tests and fixtures;
- `7c35e06487b0d62f67733dcbb167837f65e18214` -- design, protocol, and
  maintenance documentation.

For local evidence pointers below, `RUN_ROOT` means:

```text
/Users/hongjiseung/.codex/runs/chemsmart-pyscf-agent-24h
```

Claim classes have the following meanings:

- **Implemented**: committed code exposes the stated behavior. This does not
  establish that the behavior passed a check or survived a live model session.
- **Deterministically verified**: a bounded deterministic receipt or exact
  byte comparison observed the stated behavior. This is not live-model or
  scientific validation.
- **Live-model verified**: a preserved provider episode observed the stated
  behavior under its exact case, arm, source freeze, and model configuration.
- **Unverified**: code or an experimental design exists, but no admissible
  observation establishes the claim.
- **Retired or invalid**: the named metric or interpretation must not be used
  as evidence in the recovered campaign.

## Claim matrix

| ID | Class | Bounded claim | Primary evidence | Explicit non-claim |
|---|---|---|---|---|
| IMP-01 | Implemented | The program layer contains project-YAML-driven PySCF CPU SP, OPT, and HESS workflows, GPU4PySCF planning/preview paths, result validation, environment/process observation, and bounded CPU TD/TDA preview. | Commit `c986180d`; principal paths `chemsmart/jobs/pyscf/`, `chemsmart/cli/pyscf/`, `chemsmart/settings/pyscf.py`, and `chemsmart/utils/process_observation.py`. | No GPU calculation was run. TD/TDA real execution is not enabled. Code presence is not proof of a converged or scientifically suitable calculation. |
| IMP-02 | Implemented | Runtime V2 contains typed scientific workflow/DAG, project, command, preview, capability, execution, and artifact-handoff contracts. | Commit `728b1908`; principal paths `chemsmart/agent/workflows.py`, `chemsmart/agent/execution.py`, `chemsmart/agent/tool_runtime.py`, and `chemsmart/agent/runtime/`. | No claim of paper-level planning, broad chemistry coverage, or production readiness. |
| IMP-03 | Implemented | The harness contains a coordinator-only path, bounded specialists, full and causal feedback projections, a fresh read-only critic, host-oracle grading, D/F/C configuration, and campaign analysis contracts. | Commit `728b1908`; `chemsmart/agent/specialists.py`, `live_specialists.py`, `feedback.py`, `reviews.py`, and `experiments/`. | Existence of all factor paths is not evidence that each factor has run or improves an outcome. |
| IMP-04 | Implemented | Alibaba Token Plan and DeepSeek adapters are connected to the provider-neutral session runner. | Commit `728b1908`; `chemsmart/agent/runtime/alibaba.py`, `runtime/deepseek.py`, `provider_config.py`, and `services/unified_session.py`. | Adapter implementation does not establish current quota, availability, or scientific performance. |
| DET-01 | Deterministically verified | One focused milestone rerun passed 394/394 checks after four stale fixture-contract failures were corrected without changing the production contract. The recorded scope includes PySCF settings/CLI/project YAML, TD preview, provenance, validation, resources, Runtime V2 replay/data edges, specialists, host oracle, critic disposition, campaign controller, and implementation freeze. | `RUN_ROOT/milestones/m1-program-dag-focused-suite.json` (SHA-256 `f5446025eccfdc92ff1fa012cfb34b96ce0bb0d1263bb16d1a51923d072f1295`). | This is a focused suite, not a full-suite, live-model, engine-execution, or scientific validation result. |
| DET-02 | Deterministically verified | The current committed baseline byte-matches all 420 worktree files captured by the v16 implementation snapshot. The v16 manifest records restoration verification and manifest digest `7a67561e1fb92a218a8646bf905f3c87a72c08b1cc0dc1dfccdc8c18f30e382a`. | `RUN_ROOT/campaign/frozen-oracle-integrity-gate-v16/implementation-freeze-manifest.json` (file SHA-256 `6a958bcbbf5a20a07d8767bb0c8f5fe3d943675e74c35dd67a85f9c60136625b`) plus the committed tree at `7c35e064`. | The manifest does not prove behavior for files outside its captured set and does not make v16 a factorial result. |
| DET-03 | Deterministically verified | The host oracle for the v16 producer-edge case graded SP and OPT as sibling nodes, HESS as dependent on optimized geometry, resolvable initial nodes as previewed, and the future HESS node as unresolved. | v16 episode envelope described in LIVE-02, `grade.checks` and `ledger`. | This verifies one fixture/oracle pair only; it is not evidence for all DAG shapes. |
| LIVE-01 | Live-model verified | A Qwen liveness call observed exact model ID `qwen3.8-max`, `xhigh` reasoning, HTTP 200, and the public response `PONG qwen3.8-max`. | `RUN_ROOT/transport/qwen-liveness-20260804T0243.json` (SHA-256 `19efdb0db2a2e300a1e66d0be0584621559ba656a925e556440f1cf48fe7d0b9`). | This proves provider/model transport liveness only, not ChemSmart session, tool-loop, factorial, or scientific behavior. |
| LIVE-02 | Live-model verified | The exact case/arm `QP-DEV-006.d1-ffull-c0.r0` reached factor-valid, safety-green, `previewed`, strict deterministic pass in v14, v15, and v16. The latest v16 observation used 4 provider sessions, 18 attempts, 435,319 input tokens, 48,848 output tokens, 29,723 reasoning tokens, 33 successful tool calls, and 3 rejected tool calls. | `RUN_ROOT/campaign/frozen-output-envelope-gate-v14/episode-public-evidence/QP-DEV-006.d1-ffull-c0.r0/episode-envelope.json` (SHA-256 `30d87b4ccb274c6246cf61617ddb0bfe59b65df0b0bc97248869775432cd7d8a`); corresponding v15 envelope (SHA-256 `4154aa575357c0a835c11a2bb442594772e12f78f80cb7a0d47a96a577fddd74`); `RUN_ROOT/campaign/frozen-oracle-integrity-gate-v16/episode-public-evidence/QP-DEV-006.d1-ffull-c0.r0/episode-envelope.json` (SHA-256 `0a994afebe9b5d56ea3e720621a334ac372081b2dfac45f6911fb00ce05a3c0e`). | These are three evolving integration gates for one case and one arm, not three independent factorial replicates. They do not establish D, F, or C effects or cross-case transfer. |
| UNV-01 | Unverified | Live coordinator-only behavior (`D=0`) remains unobserved by a gradeable model episode. | No admissible live envelope exists before the recovered campaign. | Do not infer lower cost or equal accuracy merely from the absence of specialist sessions. |
| UNV-02 | Unverified | Live causal-feedback behavior (`F=causal`) and total provider-token savings remain unobserved. | Projection code and deterministic tests exist in commits `728b1908` and `40f1f8c1`; no admissible live causal arm exists. | Serialized byte reduction is not token or episode-cost reduction. |
| UNV-03 | Unverified | A fresh critic (`C=1`) has not produced an admissible live factorial observation, and seeded critic recall/false-rejection/false-pass metrics have not been measured. | Critic and seeded-oracle contracts exist in `chemsmart/agent/reviews.py` and `experiments/qwen_pyscf_critic_efficacy.py`; no qualifying live receipt exists. | Do not claim critic efficacy, improved coordinator correctness, or reduced false passes. |
| UNV-04 | Unverified | QP-DEV-007 transfer, cross-case behavior, GPU execution, critic effectiveness, broad computational-chemistry generality, and production-default suitability remain open. | No qualifying artifact currently discharges these claims. | Implementation and one SP/OPT/HESS preview case cannot discharge them. |
| INV-01 | Retired or invalid | The twelve v17 public envelopes are transport-failure evidence, not factorial observations. Across 12 dispatches and 35 provider attempts, recorded input, output, and reasoning tokens and successful/failed tool calls were all zero. Twelve of 24 scheduled episodes were never dispatched. | `RUN_ROOT/campaign/frozen-primary-dfc-v17/implementation-freeze-manifest.json` (SHA-256 `6a958bcbbf5a20a07d8767bb0c8f5fe3d943675e74c35dd67a85f9c60136625b`), `freeze-bundle.json` (SHA-256 `943893b8a517c5ff5761f563408fc795a2471c2d1631d2be1d67e90ff6bc5a11`), `primary-ledgers/batch-001.json` through `batch-006.json`, and `episode-public-evidence/*/episode-envelope.json`. | The 11 `inconclusive`/factor-invalid labels and one `fail`/quota-exhausted label must not be analyzed as D/F/C outcomes. The v17 collapse does not contradict the v16 pass because both bind the same implementation manifest and v17 observed no model generation. |
| INV-02 | Retired or invalid | A session terminal label such as `complete` is not a scientific success metric. | The grade and scientific state are separate fields in the v14-v16 envelopes and Runtime V2 contracts. | Never substitute terminal state for deterministic grade, preview, validation, or reproduction. |
| INV-03 | Retired or invalid | Pilot and integration-gate runs v10-v16 are excluded from factorial-effect estimation. | `docs/evaluation/qwen-pyscf-dfc-24h-protocol.md` and the corresponding frozen campaign directories. | They may diagnose harness behavior but cannot estimate component effects. |
| INV-04 | Retired or invalid | Normalized or repaired structured output must not be reported as raw strict model success. | Output-envelope receipts and the raw-versus-normalized distinction in the committed protocol and campaign records. | A host normalization demonstrates recoverable intent, not raw wire compliance. |
| INV-05 | Retired or invalid | Presence of committed tests, test count, or source volume is not model-performance evidence. | Commit `40f1f8c1` versus deterministic receipt DET-01 and live receipt LIVE-02. | Keep implementation, deterministic verification, and live-model verification separate. |

## Interpretation of the C factor

In the primary D x F x C design, `C=1` is a **detector-and-cost factor only**.
The coordinator candidate and its raw deterministic grade are immutable before
the critic runs. The critic may emit findings against a detached candidate,
but it cannot repair that candidate, alter project YAML or the command DAG,
approve or execute work, change readiness, or replace the coordinator grade.

Consequently:

- C outcomes are critic realization, findings, detector agreement with a
  preregistered oracle when one exists, false findings, provider tokens,
  latency, and session/tool cost.
- A difference in coordinator strict grades between independently sampled C=0
  and C=1 episodes is not a critic-caused improvement. The critic runs after
  that grade's candidate is frozen.
- Critical recall, overall recall, false rejection, and false-pass reduction
  require the separate seeded-oracle study. No such live study has run.
- A future design in which critic findings trigger coordinator repair would be
  a different intervention and must receive a new factor name and freeze.

This interpretation is implemented in
`chemsmart/agent/experiments/qwen_pyscf_critic_efficacy.py`, which fixes
`primary_dfc_critic_interpretation` to `descriptive-only`, and in
`chemsmart/agent/reviews.py`, where critic evaluation is post-review analysis
rather than a corrected scientific grade.

## Admission rule for new evidence

A new provider episode may extend the live-model rows only when it has a unique
dispatch/attempt identity, nonzero model observation, a source/configuration/
oracle binding, valid factor realization for the factor being interpreted,
and a complete public evidence envelope. Pre-generation quota or rate-limit
failures remain transport observations. They cannot consume or grade a
scientific factorial cell.
