# Wizard Vision Review

`/wizard` already does more than configuration scaffolding: it probes a real HPC environment, normalizes scheduler/software facts into a portable YAML contract, and refreshes a sidecar node cache. Read through that lens, v1 is the first control-plane primitive for a chemistry execution fabric, not just a setup helper.

## 1. Latent capabilities wizard unlocks beyond setup — automated cross-cluster benchmarking, queue capacity-aware scheduling, GPU vs CPU job routing.

The deepest unlock is that chemsmart now has a machine-readable description of **where** work can run, not just **how** to write an input file. The current pipeline already discovers scheduler family, queue defaults, scratch policy, software candidates, and per-site transport mode. That turns cluster setup into reusable state.

In the next 6-24 months, that state can support three high-value behaviors:

- **Automated cross-cluster benchmarking.** The agent can submit the same dry-run-validated workflow to multiple known sites, compare queue latency, runtime, failure rate, and cost proxies, then learn which cluster is best for Gaussian optimizations, ORCA DLPNO jobs, or NCI post-processing.
- **Queue capacity-aware scheduling.** The cache already stores partition summaries and freshness. That can evolve from static `QUEUE_NAME` autofill into “pick the fastest acceptable queue right now” for a given walltime/memory/GPU envelope.
- **GPU vs CPU routing.** Today `choose_queue()` explicitly prefers non-GPU queues unless GPU evidence is strong. Tomorrow that same normalized queue metadata can route ORCA or future GPU-aware engines to GPU partitions while keeping Gaussian CPU-bound by default.

The key point is that `/wizard` creates a **runtime topology layer**. Once the agent knows cluster capabilities in a normalized way, the planner can stop treating HPC as a single opaque endpoint and start reasoning about it as a portfolio of execution options.

## 2. Multi-agent extensions — specialized critique agents per scheduler family; method-recommendation agent reading cache + choosing functional/basis; cluster-doctor agent warning of misconfigurations.

The current agent loop and tool registry already imply a modular future: planning, critique, runtime validation, dry-run generation, and HPC submission are separate steps. `/wizard` gives each future specialist better ground truth.

Three concrete extensions stand out:

- **Scheduler-family critique agents.** A SLURM critic, PBS critic, LSF critic, or SGE critic could inspect wizard evidence, refresh output, and `validate_runtime` failures with scheduler-specific heuristics. That is materially better than one generic LLM trying to remember every scheduler dialect.
- **Method-recommendation agent with runtime awareness.** Today `recommend_method` is intentionally conservative and project-rule based. A stronger agent could combine molecular/task context with cache-derived execution facts: recommend a cheaper DFT screen on the busy CPU cluster, route a larger ORCA job to the cluster that actually exposes the right MPI/module stack, or prefer methods whose dependencies are verified on-path.
- **Cluster-doctor agent.** Because wizard records topology, module candidates, scratch paths, and later verification/refresh outcomes, a doctor agent could flag subtle issues early: stale cache, path-only installs that break in batch mode, scratch configured but not writable, mismatched ORCA+MPI modules, or a queue choice that silently strips GPU access.

What matters here is not “more agents” as a buzzword. It is that `/wizard` produces the structured environment evidence those specialists would need to be genuinely useful rather than decorative.

## 3. Data flywheel — telemetry from wizard probes → site-specific module-name corpus → fine-tune the LLM ranking step (research §3 future work).

The feasibility write-up was right to keep v1 deterministic-first and reserve LLM use for ambiguity around module naming and environment selection. That narrow ambiguity zone is exactly where a durable data flywheel can form.

A plausible sequence:

1. Wizard probes collect scheduler evidence, module-system flavor, candidate module names, PATH findings, scratch/project hints, and the final rendered YAML.
2. Later stages contribute labels: `wizard_verify`, `wizard_refresh`, `validate_runtime`, successful dry runs, successful submissions, and eventually completed chemistry jobs.
3. Those labels identify which module candidate or environment stanza actually worked for a given site fingerprint.
4. The resulting corpus becomes training data for the one part of the system where learned ranking is justified: choosing among multiple plausible site-local module stacks.

That is much more valuable than generic fine-tuning on chemistry prompts. It is **operational fine-tuning** on high-signal, low-volume, site-specific decisions that foundation models do not see on the public internet.

If handled carefully, chemsmart could accumulate an anonymized corpus like:

- scheduler family + queue names + hostname domain
- module system type/version
- candidate Gaussian/ORCA/NCIPLOT module lines
- chosen YAML stanza
- downstream success/failure signal

Over time, the LLM stops being a guesser and becomes a ranker trained on real HPC outcomes.

## 4. Cross-cluster federation — registry of known clusters; wizard learns from other users' (anonymized) successful YAML for same site.

The next step after per-user learning is federation. Many users will probe the same university or national-lab cluster repeatedly. Re-solving that configuration from scratch is wasted effort.

A cross-cluster registry could store anonymized cluster fingerprints and successful normalized YAML patterns for each site. Then a new wizard run would not start from zero; it would start from:

- known scheduler family for that hostname/domain,
- historically successful module names,
- common queue defaults,
- known scratch conventions,
- common project/account field patterns,
- common failure warnings for that site.

Done well, this does not replace local probing. It **primes** local probing. The deterministic survey still checks live facts, but the system begins with a prior instead of a blank slate.

That would make chemsmart progressively faster and more reliable at the exact part competitors will keep treating as manual HPC glue: the first-mile and last-mile site adaptation problem.

## 5. 10x extension — what makes this category-defining? (Q-Chem/Psi4/NWChem support? DFT method ranking from literature RAG? auto-MEP construction?)

The category-defining version is not “wizard, but for more YAML.” It is a system that jointly chooses:

- the **right chemistry method**,
- the **right engine**,
- the **right cluster**, and
- the **right workflow decomposition**.

A credible 10x roadmap would combine four expansions:

- **Broader engine coverage.** Q-Chem, Psi4, NWChem, xTB/CREST, and eventually mixed CPU/GPU backends. More engines increase the value of routing logic and make wizard’s environment model central rather than peripheral.
- **Literature-grounded method ranking.** A RAG layer over benchmark papers, method reviews, and project-local defaults could move `recommend_method` from static project rules toward evidence-backed shortlists for TS searches, weak interactions, heavy elements, or open-shell systems.
- **Automated benchmark ladders.** Instead of choosing one method, the agent could assemble a staged plan: cheap prescreen, intermediate refinement, high-level single-point, with cluster-aware placement for each stage.
- **Workflow synthesis beyond single jobs.** Auto-MEP construction, conformer funnels, ORCA-after-Gaussian refinement, or reaction-path campaigns become much more viable once the system knows both chemistry tools and execution venues.

The 10x outcome is a chemistry copilot that does not merely generate inputs. It builds and executes a **computational strategy**.

## 6. Strategic moat — why would this stay valuable as foundation models improve?

As foundation models improve, generic planning and explanation get cheaper. That weakens products whose only differentiator is “LLM in front of a CLI.” It strengthens products that own rare operational context.

`/wizard` points toward exactly that rarer asset:

- site-specific scheduler behavior,
- module-name distributions,
- scratch/account conventions,
- validation and submission outcomes,
- cluster performance history for chemistry workloads,
- mapping from chemistry intent to infrastructure reality.

A stronger model can read a manual better. It still does not automatically know that one lab’s ORCA stack only works with a particular MPI module, that a cluster’s “default” queue is terrible for DLPNO jobs, or that a scratch path passes a probe but fails in batch context. Those are empirical, local, and continuously updated facts.

So the moat is not the raw model. It is the **closed-loop execution memory** around chemistry HPC operations. Better models will likely increase the value of that memory because they can exploit it more effectively.

## 7. Risk — where "visionary" becomes anti-pattern (over-engineering, scope creep, replacing chemist judgment).

There are real failure modes.

- **Over-engineering the control plane.** If wizard becomes a giant autonomy layer before the chemistry loop is stable, the project risks building orchestration theater instead of reliable throughput.
- **Scope creep from setup intelligence to universal automation.** Not every cluster needs a federated memory graph; not every job needs adaptive routing. A lot of value will come from keeping the deterministic path boring and trustworthy.
- **Replacing chemist judgment where evidence is thin.** Method selection, benchmark interpretation, and reaction-path planning should remain advisory unless the confidence and provenance are explicit.
- **Optimizing infrastructure instead of science.** The “fastest queue” is not always the right choice if it pushes users toward a weaker method or hides numerical-risk tradeoffs.
- **Telemetry/privacy mistakes.** Operational logs are strategically valuable precisely because they are site-specific; that also makes them sensitive. Federation only works if anonymization, opt-in, and provenance are strong.
- **Conflating LLM polish with capability.** The research note’s AI-vs-deterministic boundary is a strength. If future work blurs that boundary and lets the model hallucinate scheduler or module facts, the system loses trust fast.

The healthy framing is: `/wizard` should become the trusted infrastructure memory beneath the chemist-facing agent, not a reason to sideline human review or to automate decisions that remain scientifically contestable.
