# /wizard chemist review

## 1. Pain points wizard claims to solve — cite README's manual baseline; quantify chemist's pain of hand-editing.

The README’s current baseline is still “configure, then hand-check and hand-edit.” It tells the user that `make configure` sets up `~/.chemsmart`, but then explicitly says the user should inspect the generated files, copy a scheduler template into a custom server YAML, edit modules / scratch / queue details, and even create the scratch symlink by hand.[^readme-manual] For a chemist, that is exactly the annoying part of fresh-cluster onboarding: not YAML syntax, but collecting the cluster facts that have to be correct all at once.

In practice, a manual server YAML means checking at least: scheduler family, queue/partition, walltime, memory, core count, project/account, scratch path, submit command, Gaussian path, ORCA path, module loads, and any login/env boilerplate. On a familiar cluster with decent docs, that is usually 30–60 minutes for an experienced user. On a new academic cluster with fuzzy documentation, it is commonly 2–4 hours spread across multiple “submit, fail, edit, resubmit” cycles. The real pain is not typing; it is the risk that a file looks valid but is operationally wrong.

## 2. Real-world workflow integration — fits how chem groups actually use HPC (lab-shared YAMLs, individual user dirs, multi-cluster researchers)?

Mostly yes, with an important caveat. The wizard’s target is the user-local `~/.chemsmart/server/<name>.yaml`, which matches how academic groups really work: the lab may share a starter config, but each researcher still has their own home directory, account string, scratch path, and shell environment. The explicit `HOST` handling plus `/wizard-verify` and `/wizard-refresh` are also a good fit for the very common case where one postdoc uses multiple clusters and sometimes works from a laptop, sometimes from the login node itself.[^wizard-workflow]

Where it is weaker is group standardization. Real groups usually converge on “the lab’s known-good YAML for Cluster X” and then copy that file forward. I do not see `/wizard` replacing that social workflow; I see it helping generate the first draft that a senior user then blesses. That matters because module choices, MPI compatibility, Gaussian login scripts, and scratch conventions are often lab knowledge or cluster-admin knowledge, not something a probe can fully infer.

Bottom line: good fit for individual onboarding and multi-cluster users; partial fit for shared lab maintenance; not a substitute for one vetted group reference file per cluster.

## 3. Failure modes that WASTE COMPUTE — list silent-wrong-YAML scenarios (wrong NUM_HOURS → walltime kill, wrong MEM_GB → OOM, wrong NUM_CORES → underutilization, wrong scratch dir, GPU vs CPU queue confusion).

The good news is that `/wizard` is deterministic-first. The bad news is that its current validation is mostly structural, not compute-safety validation.[^wizard-validate] From a chemist’s perspective, the expensive failures are:

- **Wrong `NUM_HOURS` → walltime kill.** The current normalization caps detected walltime at 24 hours even when the queue default or max is higher.[^wizard-queue] That is conservative in software terms, but for chemistry it is dangerous: TS searches, IRCs, large conformer batches, DLPNO single-points, and frequency jobs can die repeatedly at 24 h and waste both queue time and scratch I/O.
- **Wrong `MEM_GB` → OOM or severe thrashing.** Even if the parsed memory value is syntactically valid, the wizard does not know whether that number is appropriate for *your* workload. A Gaussian frequency on a modest organic system and an ORCA DLPNO job on a metal complex do not live in the same memory regime. Structural validation will not catch that.[^wizard-validate]
- **Wrong `NUM_CORES` / `NUM_THREADS` → underutilization or poor scaling.** The wizard maps threads directly to cores and uses whatever queue default it found.[^wizard-schema] That can be fine, but it can also leave ORCA badly tuned, or give you a “valid” job that burns allocation inefficiently for days.
- **Wrong scratch directory → failed or half-failed jobs after queue wait.** Scratch discovery does test candidates, but rendering still keeps program `SCRATCH: True` even when writability was not confirmed.[^wizard-schema] That is exactly the kind of silent wrongness that chemists hate: the submission looks polished, the job starts, then Gaussian or ORCA dies because scratch is unwritable, tiny, or purged unexpectedly.
- **GPU vs CPU queue confusion → wrong hardware class.** Queue choice prefers an explicit default, otherwise an enabled/started non-GPU queue with a short name.[^wizard-queue] That is sensible for safety, but it means a GPU-capable ORCA workflow may quietly land on CPU resources unless the cluster exposes GPU defaults clearly. The reverse problem also exists: a CPU chemistry job can get pushed toward a scarce accelerator queue and sit forever.
- **Program-stack ambiguity → “runs” on the wrong environment.** When software is not already on `PATH`, module candidates are inferred and the first candidate is loaded.[^wizard-schema] For ORCA especially, the wrong MPI/module pairing can turn into a job that starts, runs poorly, or fails late.

These are not hypothetical paper cuts. They are the standard ways a seemingly reasonable submission file burns a day or a week of cluster time.

## 4. Manual editing comparison — wall-clock time saved per setup; risk traded (what experienced chemists notice that wizard can't).

For a fresh cluster, I think `/wizard` can realistically cut setup from roughly 30–60 minutes down to 5–15 minutes for an experienced user, and from 2–4 hours down to perhaps 20–40 minutes for a junior user who would otherwise hunt through cluster docs. That is real value.

But the risk trade is obvious: the wizard saves clerical effort, while the chemist still has to supply judgment. An experienced user will notice things the wizard cannot reliably know from probes alone:

- whether the “default” queue is a short debug queue rather than the queue you actually want;
- whether walltime/memory defaults are acceptable for TS, IRC, DLPNO, or frequency workloads;
- whether Gaussian needs a site-specific login script;
- whether the ORCA module must match a specific MPI stack;
- whether scratch is fast local SSD, shared filesystem, quota-limited, or aggressively purged;
- whether the correct account is the user’s default, the group default, or a specific grant for this project.

So the right mental model is: `/wizard` saves drafting time, not review time. If a chemist skips review, the time savings can be paid back many times over in wasted compute.

## 5. Adoption blockers in a typical academic group (PI buy-in, "I trust my hand-edited YAML" inertia, fear of automation, training overhead).

The biggest blocker is trust. Most computational chemists have at least one story where “helpful automation” produced a plausible-looking file that was wrong in the one field that mattered. That makes people cling to hand-edited YAMLs they already trust.

Typical group-level blockers I would expect:

- **PI/postdoc skepticism:** no one wants to explain a wasted allocation because a wizard guessed wrong.
- **Known-good YAML inertia:** once a lab has one working server file, most people would rather copy it than regenerate anything.
- **Fear of silent mistakes:** automation that is 95% right is still scary when the remaining 5% can waste a week of queue time.
- **Training overhead:** users must learn not only `/wizard`, but also when to override it, when to run `/wizard-verify`, and when to refresh stale cluster info.
- **Ownership problem:** if scheduler queues, modules, or scratch policies change, who in the group is responsible for revalidating the generated YAMLs?

In other words, the barrier is not usability alone. It is organizational confidence.

## 6. Verdict — Would you use it for fresh cluster onboarding? Recommend to a postdoc? What ONE thing would you require before approving group-wide adoption?

**Would I use it for fresh cluster onboarding?** Yes — as a first-draft generator and checklist.

**Would I recommend it to a postdoc?** Yes, but only if they already know the cluster well enough to review the output against either cluster documentation or one known-good manual submission.

**Would I approve it for group-wide adoption today?** Not yet.

**The one thing I would require first:** a mandatory compute-safety confirmation gate that forces a human to explicitly sign off on queue, walltime, memory, cores, scratch, and program stack against a lab-approved reference YAML before the wizard output is treated as ready. The current `/wizard-verify` is useful, but it verifies transport wiring, not that the generated resource choices are chemically sane.[^wizard-verify]

That is the difference between “helpful setup automation” and “automation I would trust with expensive jobs.”

[^readme-manual]: `README.md:69-75`, `README.md:84-130`, `README.md:197-201`.
[^wizard-schema]: `chemsmart/agent/wizard/render.py:1-12`, `chemsmart/agent/wizard/render.py:75-123`.
[^wizard-validate]: `chemsmart/agent/wizard/validate.py:11-18`, `chemsmart/agent/wizard/validate.py:50-71`.
[^wizard-queue]: `chemsmart/agent/wizard/normalize.py:16-55`.
[^wizard-workflow]: `chemsmart/agent/cli.py:437-509`, `chemsmart/agent/wizard/orchestrator.py:32-58`, `chemsmart/agent/wizard/refresh.py:31-68`.
[^wizard-verify]: `chemsmart/agent/wizard/verify.py:30-126`.
