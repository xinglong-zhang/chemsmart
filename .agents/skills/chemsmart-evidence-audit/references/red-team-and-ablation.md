# Red-team and ablation reference

Use a preregistered factorial study before enabling optional agent components.
Hold model/provider revision, prompt, tool schema, tasks, budgets, and order
policy fixed. Use held-out chemistry families and repeated paired runs.

Before the factorial study, compare the direct-string command front end (A0)
with typed CommandWorkflowSpec plus deterministic compiler (A1). Require A1 to
have 100% schema-valid rendering, parser acceptance, and render determinism,
zero raw native-input authoring/hallucinated options/shell injection, no more
than a two-point semantic-preview-success regression, and either at most 1.25x
token/cost or a significant bounded-repair reduction. Compiler authority is a
safety boundary; efficacy language requires this paired evidence.

## Components

- task decomposition: off/on;
- structured evidence/report generation: off/on;
- independent read-only critique: off/on.

Measure task success, chemical validity, clean-environment reproducibility,
false-pass and unsupported-claim rate, approval violations, critic precision
and recall, tool/parser errors, time, token use, cost, and handoff loss.

Adopt a component only when its preregistered benefit is supported without a
safety regression. Required red lines are zero approval bypasses, fabricated
evidence, artifact mutation, and successful completion with a required failed
deterministic gate.

Use deterministic graders first, expert rubrics second, and LLM judges only as
supplementary analysis.

For literature evidence, use Elsevier, SerpAPI, and Tavily only through
existing user-owned quota and a Keychain lease. Store citation provenance and
access/error class, never credentials or raw authorization headers. A 403 from
Elsevier is an entitlement outcome unless independently shown to be a key
failure.
