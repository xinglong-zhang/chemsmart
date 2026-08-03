# Runtime contract reference

ChemSmart already has `RuntimeV2Mode`, `TaskEnvelope`, `AgentDecision`,
`ToolReceipt`, `ArtifactRef`, `RuntimeEvent`, and a hash-chained event store.
Treat these as the extension point.

## Command-compiled extension

CommandWorkflowSpec v1 is the model-facing plan for calculation preparation.
It contains a workflow ID, task-spec ID, CLI-schema digest, ordered immutable
CommandNode records, and no executable shell or native-engine input text. Each
node declares a Click command path, canonical parameter-name/value map,
project reference and digest, input ArtifactBinding, charge, multiplicity,
execution intent, dependencies, and expected artifact classes. An
ArtifactBinding contains a stable artifact ID, hash, and producer-node ID; it
does not accept a model-selected path or placeholder.

The compiler emits CanonicalCommandInvocation and CommandPreflightReceipt
objects. Bind the rendered argv/display command to the schema, project, input,
environment, safe-preview artifact, independent parser observation, and
semantic round-trip result. A CommandCounterexample carries only a rule ID,
failed field, expected/observed values, and evidence reference for bounded
repair. It is not a prompt to regenerate a free-form command.

## Additive future contract

Introduce versioned payloads only after their fixture and replay requirements
are specified:

- `ProviderCapabilities`: protocol, structured-output support, continuation
  mode, context/tool limits, and supported parallelism;
- `ScientificTaskSpec`: molecule/artifact reference, electronic state,
  program, job kind, method settings, constraints, requested observable, units,
  assumptions, and expected evidence;
- `TaskNode` and `TaskGraph`: immutable inputs, dependencies, allowed tools,
  budget, approval scope, expected outputs, verifier, and deterministic join;
- `ResourceBudget`: token, tool-call, cost, wall-time, and compute ceilings;
- `ApprovalRequest` and `ApprovalResolution`: exact hashes, scope, actor,
  expiry, and one-shot decision;
- `EvidenceRef`, `ValidationReceipt`, `ClaimRecord`, `ReviewFinding`, and
  `ReportManifest`.

Version every new event payload and retain a registry from event kind to
payload model. Existing v1 events must replay unchanged. Opaque provider state
may support a continuation but is explicitly non-evidentiary.

## Provider-continuation boundary

Do not add a generic provider field such as `previous_response_id` to the core.
The current OpenAI adapter uses Chat Completions, where continuation is the
sanitized assistant tool-call message plus tool results in history. A
Responses-style continuation is a distinct adapter capability with a distinct
tool-output protocol.

If a provider needs resumable state, bind an opaque adapter-owned checkpoint to
session/turn, provider and wire protocol, resolved model, tool-schema/scope
digest, sanitized-history digest, remaining budgets, and approval/resource
envelope. Store raw state in a private adapter sidecar, expose only an opaque
reference to public runtime state, and invalidate it on any mismatch. Replay
must be idempotent and must never carry an approval into a changed invocation.

## Event lifecycle

Use the future lifecycle only after typed contracts exist:

`goal → scientific specification → task graph → preflight → approval →
execution → validation → review → evidence-bound report → terminal state`

Record terminal failure and blocked states explicitly. A report is not a
terminal success unless all mandatory receipts and validators are present.

Add command-workflow, command-preflight, command-counterexample, and
safe-preview payloads as versioned events rather than replacing older
RuntimeEvent payloads. Legacy logs must replay without a migration-induced
command or permission effect.
