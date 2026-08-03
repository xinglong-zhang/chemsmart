# Approval and evaluation reference

## Approval matrix

| Operation | Default | Required evidence |
| --- | --- | --- |
| Read-only inspection or schema generation | allowed | command and artifact receipt |
| Fixture or fake execution | allowed within task scope | deterministic validation result |
| Real local calculation | explicit approval | exact command, inputs, environment, cost/resource bound |
| Scheduler submission, cancellation, or retry | explicit approval | exact job artifact, scheduler target, resource bound |
| Paid API, remote execution, or publication | explicit approval | provider/target, budget, disclosure scope |

Invalidate approval whenever a bound input, project, executable, environment,
or command hash changes.

For a user-authorized API validation, borrow a standard Keychain secret only
for the request. Verify liveness and quota sufficiency without writing the
secret, raw header, or full response into a receipt. DeepSeek is the only
model-validation provider in this roadmap. Elsevier, SerpAPI, and Tavily are
literature-discovery or full-text-verification services, never chemistry
execution providers. Treat an Elsevier 403 as an entitlement result until
proven otherwise; do not label it an invalid key by default. Never top up,
change a plan, or exceed the user-owned quota. Once the user has authorized
the current development phase, calls within that recorded quota do not need
per-call approval; a new provider, target, quota expansion, or billing change
does.

## Bounded delegation

Dispatch only if subtasks have independent inputs and a typed merge operation.
The coordinator owns the task graph; workers own no shared mutable artifact.
Critics receive artifacts and declared assumptions, not persuasive self-reports.

## Evaluation rule

Keep a single-agent reference path. Compare it with any subagent or critic
configuration under fixed model, prompt, tool schema, task set, and budget.
Use deterministic outcome graders first. A component stays experimental until
it improves the preregistered metric without creating approval bypasses,
fabricated evidence, or false scientific passes.

Do not run pytest, Ruff, or broad checks after every edit. At a material
milestone run one focused suite, with at most one evidence-driven rerun. Run
the full test/lint/schema/link/citation/secret/diff gate only at the declared
freeze milestone.
