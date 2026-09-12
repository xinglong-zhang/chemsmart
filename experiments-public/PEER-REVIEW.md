# Referee report

**Manuscripts.** Two case reports reconstructing autonomous ChemSmart Agent runs:

- Case A — `experiments/ino3-r12/REVIEW.md`, nickel(II) bis(thiolate) oxidation
  (settlement `achieved_with_observations`).
- Case B — `experiments/po3-r19/REVIEW.md`, thermal azide–alkyne cycloaddition
  regiochemistry (settlement `unreachable_from_evidence`).

**What I did.** Read both reports and both `make_figures.py`; re-ran both figure
scripts; re-derived ~45 numbers from the raw ledgers, event streams, workspace
records, node outputs and artifacts using the repository's own readers
(`chemsmart.analysis.result_readers.reader_for("orca")`) and perceiver
(`chemsmart.io.molecules.perception.adjacency_matrix`); checked every blockquote
in both reports mechanically against the raw corpus; inspected all four figures
at full resolution. I did not modify either report or anything under `raw/`.

---

## Verdict

**Case A — `ino3-r12`: accept after revision (should-fix), with one blocking
correction.** The arithmetic is sound to machine precision, the quoting is
exact, and the narrative is genuinely followable. But the report presents an
uncertainty budget as fully evidence-backed while omitting the one number the
host itself displayed next to the constant that sets the zero of the potential
scale: the source's own stated 0.05–0.1 V accuracy on the 4.988 V absolute
ferrocene potential. The budget carries 0.008 V there — the *internal
consistency* of two constructions of the scale, which is not its accuracy. That
term propagates one-for-one into the delivered couple, and the report's
"the final error bar contains no magnitude the host cannot read off a receipt"
reads as completeness when it is not. `met` survives (0.195 V against a 0.200 V
tolerance if the source accuracy is used at its 0.1 V end), so the conclusion
holds — but the margin goes from comfortable to 3 %, and a referee cannot let
that stay unsaid. Four further omissions (a sixth cycle-1 finding, three
`declared_observable_misses` on the final completion receipt, the absent
functional-dependent-geometry term, the missing seal/harness pin) all point the
same way: the report is scrupulous about what the record says and less
scrupulous about what the record says is *missing*.

**Case B — `po3-r19`: accept after revision (should-fix), with two blocking
corrections.** This is the better of the two as a piece of scientific writing:
the negative result is not just admitted but made the centre of the case, and
figure B1 is a model of how to show 26 000 s spent on the wrong coordinate. Both
blocking items are the same species of defect — a load-bearing number the report
invites the reader to check, which does not check out as written. (1) The
arithmetic floor that *is* the settlement, `u-floor-without-method-kcal` =
0.845452, is a three-term quadrature; the report, its figure title and its figure
caption all name two terms, which combine to 0.830, not 0.845. (2) §6.4 is headed
"How the host graded the pre-registered prediction" and silently drops the host's
own `declared_after_evidence: true` flag, which sits on four of the five
diagnostics tabulated there — and one of those four says in its own
`expectation_basis` that the workspace record "already carries" the numbers it
predicts. The charter's stated reason for that field is that "a restatement is
never displayed as a pre-registration"; the table does precisely that. Neither
error changes a conclusion. Both must be fixed before a reader is asked to trust
the case's honesty as its headline virtue.

**Findings: 3 blocking, 18 should-fix, 10 minor** (31 in all, numbered A1–A15
and B1–B16). The seven chemistry concerns C1–C7 are listed separately at the
end and are not counted as findings against the writing, except where a
numbered finding cross-references them.

---

## Findings — Case A (`ino3-r12`)

### A1. BLOCKING — the reference constant's own stated accuracy is absent from the budget and from the report
`experiments/ino3-r12/REVIEW.md` §5 ("It pre-registered a refutation of its own
reference scale"), §6 (budget table, row "ferrocene reference-scale
consistency | 0.008000 | measured"), §8 ("The final error bar contains no
magnitude the host cannot read off a receipt"), and `figures/figA3-couple-and-budget.png`
panel (b), row "ferrocene reference scale 0.0080 V".

The host's own cycle-1 analysis report
(`raw/ino3-r12/.chemsmart-agent/goals/goal-ino3-r12/runs/cycle-1/analysis/partial-analysis-report.md`,
line 35) prints the constant with its provenance and this clause:

> … electron at rest in the gas phase; **the source states a likely accuracy of
> 0.05-0.1 V**

I confirmed the same clause in `chemsmart/analysis/literature_constants.py`
(`ferrocene_absolute_reduction_potential_acetonitrile_namazian2010`). The string
`0.05-0.1 V` appears four times in the cycle-1 event stream and in the cycle-1
public transcript. It appears **nowhere** in the report, and nothing of that
magnitude appears in the budget.

What the 0.008 V term actually is, in the record's own words on the claim:

> "Ferrocene reference scale: the gap between the registered namazian2010
> absolute Fc+/Fc potential in acetonitrile (4.988 V) and the same scale reached
> through the registered Fc-vs-SCE (0.38 V) plus SCE-absolute (4.60 V) pair."

That is agreement between two constructions of one scale. It is not the scale's
accuracy, and the two are not interchangeable: the delivered couple is
E_abs(complex, PBE0/CPCM) − E_abs(Fc, G3(MP2)-RAD/COSMO-RS), so an error in the
ferrocene leg enters volt for volt and does not cancel against the complex leg,
which was computed by a different method entirely. The settlement's
"leading systematic is common-mode to the subtraction" argument covers the
*solvent model*; it does not cover a method mismatch between the two legs.

**Correction requested.** In §6, add a row to the budget table naming the
source-stated 0.05–0.1 V accuracy of the reference constant as an absent term,
beside the CPCM row, and change §6's single "one term the agent named as
unmeasurable" to "two absent terms, one of which the case named". In §5, after
the 0.008 V measurement, add one sentence distinguishing scale *consistency*
from scale *accuracy*. In §8, state the sensitivity: with the reference term at
0.05 V the quadrature total is 0.174 V and at 0.1 V it is 0.195 V, so `met`
against ±0.2 V survives but with ~3 % headroom rather than 16 %. (I re-derived
those: 0.167348 → 0.174474 → 0.194785 V.) Add the same absent row to figure
A3(b).

### A2. SHOULD-FIX — the budget terms are labelled "1σ half-width"; the record calls them full spreads, and records no combination rule at all
§6 table header ("1σ half-width / V"), §6 prose ("Five separately measured axes
in quadrature"), and `figA3` panel (b) x-axis ("1σ half-width contributed / V").

The claim's own component meanings, read from
`raw/ino3-r12/.chemsmart-agent/runs/live-20260909T104643999449Z-b6464705b8-77f484d1/events.jsonl`
(`analysis_claims_recorded` → `claims[e-couple-vs-fc-pme3].uncertainty_components`), say:

- functional: "the **largest displacement** of the couple when … B3LYP
  (-0.129045 V) or TPSSh (-0.146959 V)…"
- basis: "the measured shift … **carried in full** as the scale of the uncomputed
  residual beyond triple-zeta"
- thermochemical model: "**the full spread** of the def2-SVP couple across
  electronic-only, harmonic RRHO and Grimme quasi-RRHO"

I verified each numerically: u-functional 0.146959 = max−min over
{−0.394149, −0.523195, −0.541109}; u-treatment 0.012906 = max−min over
{−0.477605, −0.472742, −0.485648}. These are full ranges, not half-widths, and
not σ of anything. Separately, the ino3 claim carries
`uncertainty_combination: null` — no `rule`, no `dependence`, no `coverage`.
Case B's claim carries all three, including the explicit phrase "each taken as a
1-sigma half-width", which is why the identical axis label is *correct* on figure
B4 and *not* correct on figure A3.

**Correction requested.** Retitle the column and the axis to "measured spread /
displacement (V)"; replace "Five separately measured axes in quadrature" with a
sentence saying the five magnitudes are full spreads combined in quadrature by
an expression node, and that the claim records no combination rule, dependence
or coverage statement — noting that Case B's does, since the two reports invite
comparison.

### A3. SHOULD-FIX — the final completion receipt is a delivery *with stated limitations*, and §0 does not say so
§0 ("Errors and limits first", four items) and §6 ("The spin answer and the
quartet").

The cycle-3 `analysis_completion_evaluated` payload carries:

```
limitation_output_ids = ["declared_observable:spin-ni-pme3",
                         "declared_observable:spin-s1-pme3",
                         "declared_observable:spin-s2-pme3"]
declared_observable_misses = [three lines, one per spin observable]
```

Half of the goal's six declared observables are on the final completion
receipt as limitations. The settlement resolves them at the goal grain
("delivered in an earlier cycle: spin-ni-pme3 (delivered in cycle 1) …"), so the
settlement is right — but the report's §0 lists four things "a reviewer should
know before the science" and this is not one of them. Case B flags the exactly
analogous cycle-6 miss line in its §0.5 *and* its §9 and says "I am not
smoothing it over"; Case A is silent. Related: no completion record in this goal
ever scored the three spin observables against their declared bands (0.60–0.95,
0.05–0.25) — both completion records show them `not_comparable` with empty
delivered values.

**Correction requested.** Add a fifth item to §0 naming the three
`limitation_output_ids` on the cycle-3 completion receipt and the settlement
line that resolves them, and add one sentence to §6 saying the declared spin
bands were never scored by a completion record because the delivery sat in an
earlier cycle.

### A4. SHOULD-FIX — one of the six cycle-1 findings is not reported, and it is the control on the smallest budget term
§3.5: "The analysis chain came back **partial** with six findings, and one of
them matters".

The cycle-1 partial-analysis report lists six. The report accounts for five
(`expr-modes-spinsq` skipped, plus its two cascaded validations and one cascaded
claim node in §3.4; `claims-redox` failed in §3.5). The sixth is not mentioned:

> `- val-svp-replication (scientific_validation): failed -- scientific
> validation input is not typed evidence from its planned producer`

From the cycle-1 review packet, that node carried rule
`r-sp-reproduces-registered-opt`: `maximum_absolute_less_equal`, threshold
0.002, unit eV, over output `repro-diff-ev` — i.e. the planned pass/fail control
that the single point on the 0.001 Å-promoted geometry reproduces the registered
optimisation's own energy. That is the control on precisely the
`geometry materialisation` term the §6 budget later carries at 0.000983 V. The
number was eventually computed as an expression at cycle 3, and it would have
passed the threshold — but the planned verdict never landed and was never
re-planned.

**Correction requested.** Name `val-svp-replication` in §3.5 alongside
`claims-redox`, state the rule and threshold it carried, and note in §6 that the
0.000983 V geometry term is an expression result and not a verdict, because its
planned validation failed at cycle 1.

### A5. SHOULD-FIX — functional-dependent geometry relaxation is absent from the budget and, unlike CPCM, is not named as absent
§6 (budget table), §8 ("It said what it could not do"), and the settlement's six
`decision_uncertainties`.

All three species (neutral, doublet cation, quartet cation) are PBE0-D3BJ/def2-SVP
optimisations; I confirmed the geometries and their zero imaginary-mode counts.
The functional term — 0.147 V, the largest in the budget — is measured with
B3LYP and TPSSh single points *at PBE0 geometries*. The agent's own recorded
rejection makes the point sharply:

> "PBE0 is the only level at which both partners are relaxed with their own
> frequencies, while B3LYP and TPSSh are **vertical single points on PBE0
> geometries**; averaging would weight an incomplete chain equally with the
> complete one."

The same objection bears on using that spread as the functional *error term*: it
bounds the vertical part and not the relaxation part. The case never considered
re-optimising at another functional — I read all 14 `approaches_recorded` entries
and the only re-optimisation considered was at def2-TZVP (basis) and for the
quartet gap. Unlike CPCM differential solvation, this term is not named as
unmeasured anywhere, so §6's "the one term the agent named as unmeasurable" and
§8's "It said what it could not do" overstate the coverage.

**Correction requested.** Add one row to the §6 budget as a second absent term
("functional-dependent geometry relaxation — not measured, not named by the
case"), and in §4.2, where the "vertical single points" rejection is quoted
approvingly, add the referee's observation that the same wording qualifies the
term's use as the dominant uncertainty.

### A6. SHOULD-FIX — figure A1 panel (c) is subtitled with the wrong coordination geometry
`figures/figA1-structures.png`; source at `make_figures.py:184`.

Panel (c), the doublet cation at S–Ni–S 151.0° / P–Ni–P 158.0° (I reproduced
both), is subtitled **"tetrahedrally distorted; 0 imaginary modes"**. Panel (d),
the quartet at 105.9° / 106.4° — the actual pseudo-tetrahedron — carries no
geometry descriptor at all, only "+0.219 eV above (c)". The report's own prose
(§4.2) has it the right way round: "the quartet is pseudo-tetrahedral … where the
doublet is only **partly flattened** out of square planar". A chemist reading the
figure alone gets the wrong picture of both cations.

**Correction requested.** Change panel (c)'s subtitle to "square plane buckled
to 151°/158°; 0 imaginary modes" and panel (d)'s to "pseudo-tetrahedral,
106°/106°; +0.219 eV above (c) at 298.15 K".

### A7. SHOULD-FIX — the figure script cannot run outside this host, and §0.1 does not say so
§0.1 and §10.

§0.1 correctly says the workspace is not self-contained. What it does not say is
that `make_figures.py` hard-codes the out-of-tree location
(`make_figures.py:36`, `UPSTREAM = /home/chemsmart/agent-campaigns/ax41-refine-100/novel-round-3/...`)
and reads three of the four figure-A1 panels from it (lines 154–156), so the
script exits on `FileNotFoundError` for anyone without that directory. The
couple, the gap, every thermochemical input and the one delivered spin claim
(0.80934, from `nicat-opt-pme3`) also come from there. Only figure A2 and panel
A1(a) are rebuildable from `raw/`.

**Correction requested.** Add to §0.1: which figure panels and which delivered
claims come from outside the tree, and that `make_figures.py` requires that
absolute path to run. §10's "upstream results, NOT in this tree" row should say
"required by make_figures.py".

### A8. SHOULD-FIX — no seal and no harness commit
§10 provenance table.

Case B names both its seal (`seals/ROUND17-PO3-R19-CLEAN-REPLAY.md`) and its
harness pin (`e08c34fe`). Case A names neither. From
`/home/chemsmart/agent-campaigns/ax41-refine-100/novel-round-7/seals/SUFFICIENCY-3-INO3-R12.md`:
this window's seal is that file and the harness was **`d692818c`**
(`feat/agent-accountable-precision`, 40 commits on `5bbade01`). This is not
bookkeeping: it is the reason finding A9 exists.

**Correction requested.** Add a seal row and a harness row to §10.

### A9. SHOULD-FIX — "pre-registered" needs qualifying for the cycle-3 diagnostic
§5 ("It pre-registered a refutation of its own reference scale. The first typed
act declared a *new* diagnostic, `u-reference-scale-pme3`, before reading any
constant").

Both halves are literally true within cycle 3 — I confirmed the declaration is
event sequence 10–12 and the first extraction is sequence 21. But the goal
already had two cycles of extraction, thermochemistry, expression and claim
receipts. The host did **not** stamp `declared_after_evidence` here (I grepped:
zero occurrences in all three ino3 sessions), whereas the later harness that ran
Case B stamps exactly this pattern (four occurrences there). So the absence of
the flag is a harness-version artefact of the earlier tree, not evidence of
pre-registration, and the report — which names no harness — leaves a reader to
conclude the wrong thing, especially reading the two cases side by side.

**Correction requested.** Qualify §5 to "pre-registered relative to its own
measurement, though declared at cycle 3 after two cycles of physics; this window
predates the goal-scoped `declared_after_evidence` flag that Case B's harness
applies", and add the harness pin per A8.

### A10. SHOULD-FIX — figure legibility
`figA1-structures.png`: the "P2" label is clipped to ".2" in panel (c) and cut
off in panel (d); panel (d) is drawn at a visibly smaller scale than (a)–(c), so
the Ni–S contraction / Ni–P elongation the caption asks the reader to see cannot
be compared by eye across panels; panel (a) uses unnumbered "S"/"P" where (b)–(d)
use S1/S2/P1/P2. `figA3` panel (a): the TPSSh line and its label sit in light blue
on the light grey interval band and are near-illegible in print.

**Correction requested.** Fix label placement so no label overlaps an atom;
share one axis scale across the four A1 panels; number the labels in panel (a);
darken the TPSSh series in A3(a) or move its label out of the band.

### A11. MINOR — the spin-sum tolerance is overstated
§7, figure A2: "each row sums to 1.000 e within 2 × 10⁻⁶". The B3LYP Mulliken row
sums to 0.9999969999999999, a deviation of 3.0 × 10⁻⁶. Correction: "within
3 × 10⁻⁶".

### A12. MINOR — an electronic energy called a free energy
§4: "moves **the ionisation free energy** from 4.5094 eV to 4.5880 eV". Those two
numbers are SCF electronic-energy differences (I reproduce 4.509411 and
4.588004 eV from the four total energies). The free-energy couple is the separate
−0.472742 V RRHO quantity. Correction: "the vertical-basis ionisation *energy*".

### A13. MINOR — two linear sums, no bridge
§4.2 quotes cycle 2's "linear sum of the **six** components (0.2649 V)"; §6 gives
0.247442 V over **five**. Both are right for their cycle (the settlement's fifth
`decision_uncertainty` explains the re-decomposition). Correction: one clause in
§6 saying the linear sum changed with the decomposition between cycles.

### A14. MINOR — an unmarked elision in a fenced block
§3.3 shows the TPSSh cation project YAML as a `solv:` block. The file
(`runs/live-20260909T0859*/projects/proj-orca-sp-svp-tpssh-cat.yaml`) also carries
a `gas:` section with the same seven fields. The `solv:` section is quoted
exactly; the elision is unmarked. Correction: add `# (the file also carries the
gas: section)` or an ellipsis.

### A15. MINOR — the delegated approval is disclosed only in the provenance table
§3.3 and §8 present "one displayed human decision" as a headline invariant; §10
records `granted_by: claude-owner-delegated-reviewer`. I verified the envelope
(`goal.json`). Correction: name the delegation where the strength is claimed,
not only where the digests are listed.

---

## Findings — Case B (`po3-r19`)

### B1. BLOCKING — §6.4 displays four host-flagged post-hoc restatements under the heading "pre-registered prediction"
`experiments/po3-r19/REVIEW.md` §6.4 ("How the host graded the pre-registered
prediction") and its six-row table.

The table's values are faithful — I matched every one against the final
`analysis_completion_evaluated` record in session `…-99b6acd6`. What the table
omits is a field the host wrote on the same rows:

| observable | agreement | `declared_after_evidence` |
|---|---|---|
| `dg-ts-esterc4-minus-esterc5-353k` (cycle 1) | indeterminate | *(none)* |
| `c4-restart-imaginary-mode-count` (cycle 3) | agreed | **true** |
| `ddg-basis-shift-tzvp-kcal` (cycle 5) | agreed | **true** |
| `ddg-functional-shift-tzvp-kcal` (cycle 5) | agreed | **true** |
| `ddg-regio-signed-spread-across-levels-353k` (cycle 6) | agreed | **true** |

Exactly one of the five is an unflagged pre-registration, and it is the only one
the host declined to score. The charter's reason for the field is explicit:
"a restatement is never displayed as a pre-registration." The table does that,
under a heading that says "pre-registered", and §6.4 closes with "That is the
cleanest instance of honesty-by-construction in either case."

The worst row is not merely post-hoc but self-declaredly so. From the cycle-6
declaration in the ledger, `ddg-regio-signed-spread-across-levels-353k`'s own
`expectation_basis`:

> "The host-written workspace record **already carries** the three levels' signed
> differences (-0.456, -0.221 and +0.675 kcal/mol …), so their range **cannot be
> smaller than about 1.13 kcal/mol**; this diagnostic is declared to make that
> level dependence a typed, scored quantity rather than prose."

It was then scored `agreed` for landing at 1.1308 inside a 0.5–1.5 band. The
stated motive (typing a fact rather than leaving it in prose) is legitimate; a
`agreed` grade presented as a successful prediction is not. `c4-restart-imaginary-mode-count`
is similar in kind though milder: its own `expectation_basis` notes the reached
structure's "parsed mode list already carries an imaginary mode".

**Correction requested.** Add a `declared_after_evidence` column to §6.4's table
with the host's values; retitle the section (e.g. "How the host graded the
predictions, and which one was actually pre-registered"); state that one
diagnostic was declared from numbers already on the record and quote its own
`expectation_basis` admitting it; and move the "cleanest instance of
honesty-by-construction" sentence so it attaches to the single cycle-1
diagnostic it is true of.

### B2. BLOCKING — the arithmetic floor that carries the settlement does not reproduce from the two terms the report names
§6.3 ("It converted the refusal from an argument into a receipt"),
`figures/figB4-level-ladder.png` panel (b) title, and the figure B4 caption in §7.

The quoted sentence, the report's endorsement of it, the figure's panel title and
the figure's caption all name two terms:

> "with the electronic-structure term set to **exactly zero**, the measured
> entropy-model term (0.3564) and the unsampled conformer-and-medium allowance
> (0.75) still combine to **0.8455 kcal/mol**"

then: "Those two numbers are themselves claims (`u-floor-without-method-kcal` =
0.845452, `u-floor-excess-over-tolerance-kcal` = 0.345452) standing on receipt
`d3bca43a`." Figure B4(b) is titled "the **two** non-electronic terms alone
exceed the tolerance". §7's caption: "the *measured* RRHO-versus-quasi-RRHO
sensitivity (0.356) and the never-sampled conformer-plus-neat-medium allowance
(0.750) combine to 0.845 kcal/mol".

Two terms do not give 0.845:

| combination | value |
|---|---|
| RSS(0.356444, 0.75) | **0.830393** |
| 0.356444 + 0.75 | 1.106444 |
| RSS(0.158862, 0.356444, 0.75) | **0.845452** ✓ |

The claim is a **three**-term quadrature — geometry convergence 0.158862, entropy
model 0.356444, conformer-plus-medium 0.750 — with only the functional term
(0.896554) zeroed. I confirmed this from the cycle-6 claim's own
`uncertainty_components` list and reproduced 0.845452117 to ten figures.
Figure B4's *bars* are all four and correct; only its title, the caption and the
prose are wrong.

The conclusion is untouched: 0.830 and 0.845 both exceed the 0.5 tolerance, and
the shape of the argument ("the floor survives the better method going to zero")
stands. But this is the number the settlement rests on, the report explicitly
offers it for checking, and as written it cannot be checked.

**Correction requested.** In §6.3, name all three terms in the floor (adding the
0.158862 geometry-convergence residual) and note that the agent's own quoted
sentence enumerates two of the three — quote it as it stands and correct it in
the report's voice immediately after. Retitle figure B4(b) to "the three
non-electronic terms alone exceed the tolerance" and repair the §7 caption the
same way. §8's "even if the electronic-structure term were driven to exactly
zero" needs no change.

### B3. SHOULD-FIX — the report never records that the final budget upgraded the DFT allowance from asserted to measured
§6.1 ("Two components `measured` with receipts (0.1589 geometry convergence,
0.3564 entropy model), two `asserted` and labelled as the session's judgement
(1.0 functional/basis, 0.75 conformer-plus-neat-medium)") and §0.4.

That is the **cycle-4** budget and correct for cycle 4. The delivered cycle-6
budget has three measured terms and no 1.0:

| term | magnitude | basis |
|---|---|---|
| geometry convergence | 0.158862 | measured |
| vibrational-entropy model | 0.356444 | measured |
| density functional (B3LYP→DSD-BLYP at def2-TZVP) | **0.896554** | **measured** |
| conformer + neat medium | 0.750 | asserted |
| total (RSS) | **1.232312** | |

I reproduced all four and the total. The host's own observation line on the
cycle-6 claim records only one model-authored literal
(`literal_value=0.75@lit-model-allowance`), where §6.1's quoted line —
correctly, for cycle 4 — records two. Figure B4(b) already shows 0.8966 labelled
"measured in cycle 5". The prose never says the upgrade happened, so a reader
following the report alone believes the final error bar rests half on assertion
when it rests three-quarters on measurement. This is the one place the report
undersells the work.

**Correction requested.** Mark §6.1's table and quote as cycle 4's, and add to
§6.3 the final four-term budget with three measured terms and the total
1.232312, noting that the measured functional shift replaced the asserted 1.0.

### B4. SHOULD-FIX — the pre-registered trigger for the refusal never fired
§6.2, §6.3 and §8 ("The refusal is narrow, verified and correctly aimed").

`ddg-functional-shift-tzvp-kcal` was declared at cycle 5 with this
`failure_update_rule`:

> "If the shift exceeds 1.0 kcal/mol the functional choice, not the basis or the
> geometry, is the dominant irreducible error: **I state that as the reason the
> task's 0.5 kcal/mol tolerance cannot be met in this envelope and refuse the
> tolerance on that measured receipt** rather than on an asserted allowance."

The shift came in at 0.896554 — inside the declared 0–1.0 band, scored `agreed`,
rule not triggered. The refusal was then reached by the separately constructed
floor argument at cycle 6. That is a legitimate and arguably stronger route, but
it is not the pre-registered one, and a report that (rightly) makes much of
pre-registered consequences firing should say when one did not.

**Correction requested.** One paragraph in §6.3 stating that the cycle-5
pre-registered refusal trigger did not fire (0.8966 < 1.0, graded `agreed`) and
that the refusal was reached instead by the cycle-6 floor construction.

### B5. SHOULD-FIX — the cost of the negative result is overstated by conflating it with two useful calculations
§0.2 ("Four of the eleven engine calls — 26 896 s, 79 % of all engine wall this
goal spent — bought a negative result. That is the largest single fact about this
case") and the figure B1 caption in §7 ("26 896 s of engine wall, 79 % of
everything this goal spent, on the wrong coordinate").

Cycle 1's four engine calls were two reactant optimisations and two scans. From
the per-node receipts:

| node | wall / s |
|---|---|
| `scan-esterc5-path` | 14 213.05 |
| `scan-esterc4-path` | 12 308.02 |
| `opt-benzyl-azide` | 266.98 |
| `opt-tfm-ynoate` | 107.59 |

The two reactant optimisations were *used* — they are the reference the absolute
barriers stand on, and the report itself reads `opt-tfm-ynoate` back in §9. The
wrong coordinate cost **two** of eleven calls and **26 521 s = 78.4 %**, not four
and 79 %.

**Correction requested.** In §0.2 and the figure B1 caption: "Two of the eleven
engine calls — 26 521 s, 78 % of all engine wall — bought a negative result",
with a clause noting the other two cycle-1 calls produced the reactant reference
both barriers use.

### B6. SHOULD-FIX — §1's tolerance argument is a non-sequitur as written
§1, closing paragraph.

The report sees the distinction and then argues past it:

> "Note what the tolerance *is*. It is not an accuracy target; it is the width of
> a decision. … That makes 'the calculation cannot resolve 0.5 kcal/mol' an
> answer to the question rather than an evasion of it"

If 0.5 kcal/mol is a decision width, then "cannot resolve 0.5" does not follow
from the task — it follows from the agent's additional choice to declare
`required_tolerance: 0.5, tolerance_origin: task`, which I confirmed in the
cycle-1 declaration. That choice is defensible (you need precision comparable to
a threshold to place a value against it) and it is recorded in `tolerance_basis`,
but it is an interpretation, and the settlement word
`unreachable_from_evidence` depends on it. Under the decision reading the case
answered the chemist outright: no level puts |ΔΔG‡| above 0.7 with an uncertainty
of 1.23, so by the chemist's own rule the thermal route is out. §8 makes that
argument well; §1 should not borrow the conclusion from the reading it has just
set aside.

**Correction requested.** Rewrite the paragraph to state both readings — decision
threshold vs required precision — say which one the agent declared and where, and
say that the chemist's decision is answered under either, while only the
precision reading produces `unreachable_from_evidence`.

### B7. SHOULD-FIX — the bare-id observable carries three inconsistent delivered values, two superseded, none retracted
§9 ("One nuance on 'the value was delivered'").

§9 correctly reports the cycle-6 miss and correctly says delivery under the exact
id happened at cycles 2, 4 and 5. What it does not say is what those three rows
contain. From `raw/po3-r19/.chemsmart-agent/workspace-record.jsonl`:

| cycle | `ddg-activation-regio-353k` | stated uncertainty | sufficiency | provenance |
|---|---|---|---|---|
| 2 | 0.614452 | 2.513852 | short | unconverged `ts-esterc4`; forming bonds mislabelled (§0.3) |
| 4 | 0.455590 | 1.309500 | short | converged pair, B3LYP/def2-SVP |
| 5 | 0.675214 | 1.179635 | short | DSD-BLYP/def2-TZVP; entropy term zeroed (§0.4) |

Three values under one id, spanning 48 % of the largest, all graded `short`, and
both of the superseded rows carry defects the report documents elsewhere. A
reader who greps the workspace record for the required observable — which §9
directs them to do — finds three mutually inconsistent numbers with no marker
saying which is current.

**Correction requested.** Give §9 the table above, note that all three are
`short`, and say which two are superseded and by what (the cycle-3 restart, and
the cycle-6 entropy re-measurement).

### B8. SHOULD-FIX — an unmarked elision inside a fenced "verbatim" engine block
§3.3.

The report shows:

```
orca: Use one of the keywords
orca: RIJCOSX: treat the HF exchange by chain-of-spheres
orca: RIJK : treat Coulomb and exchange both by RI
Error (ORCA_MAIN): ... aborting the run
```

The recorded `aborted_engine_lines` for each of the four nodes are:

```
Use one of the keywords
RIJDX : treat the HF exchange exactly (equals RIJONX)
RIJCOSX: treat the HF exchange by chain-of-spheres
RIJK : treat Coulomb and exchange both by RI
[file orca_main/main_input_check.cpp, line 3594]: Error (ORCA_MAIN): ... aborting the run
```

A line is dropped, a source-file prefix is dropped, and an `orca: ` prefix that
is not in the record is added — inside a block a reader will take as the
program's own words. The dropped line is not inert: it is a third legal option,
and the agent's own quoted reading ("B3LYP needs RIJCOSX (or RIJK)") also omits
it, which is worth noting in its own right.

**Correction requested.** Quote all five lines or mark the elision with an
ellipsis, drop the invented `orca: ` prefix, and add a clause noting that the
program named three options and the session's reading named two.

### B9. SHOULD-FIX — the convergence check's direction is the interesting part and is not stated
§5 ("That judgement was right in a checkable way: the restart moved the
ester-at-C4 353 K Gibbs energy by 0.158862 kcal/mol, which is **35 %** of the
difference being measured").

I reproduce 0.158862 kcal/mol and 34.87 %. The sign is missing, and it is what
makes the check matter: the restart *raised* G(ester-at-C4) by +0.158862 while
*lowering* its electronic energy by 0.073660 (the unconverged structure's
15.76 cm⁻¹ mode was inflating its RRHO entropy, exactly as the agent said). So
ΔΔG went from −0.6145 to −0.4556: standing by the unconverged structure would
have made the ester-at-C4 preference look **larger** — biased toward the isomer
the chemist said they need. A bias in that direction is the one worth naming.

**Correction requested.** State the sign and the consequence: "the restart raised
G(ester-at-C4) by 0.158862 kcal/mol while lowering its electronic energy by
0.074, moving ΔΔG‡ from −0.6145 to −0.4556 — the unconverged structure had
exaggerated the preference for the chemist's target isomer by 35 %."

### B10. SHOULD-FIX — figures B2 and B3 are not legible at print size
`figures/figB2-scan-structures.png`, `figures/figB3-saddles.png`.

B2: atom and distance labels collide with atoms and with each other in six of the
eight panels — path A points 5 and 7 have "C(ester)", "C(CF₃)" and the purple
distance overprinting one another ("C(CF3)" renders as "C(C…)"); path B point 1's
purple distance label is unreadable; path B point 7 piles "N1", "C(ester)",
"C(CF₃)" and "3.59" into one small region. The panels are also not on a common
camera or scale despite the caption saying they are "oriented on the alkyne
axis", so the caption's instruction — "Watch the azide in the top row" — cannot
be followed reliably: the molecule changes size and aspect between frames.
Roughly 30 % of the canvas below the second row is empty.

B3: the same collision problem in all three panels; in panel (b) "C(CF₃)" is
completely obscured by the F atoms. The lower third is empty whitespace.

**Correction requested.** Fix label placement (offset labels away from atoms, or
use leader lines); fix one camera and one scale per row in B2; reclaim the empty
canvas to enlarge the panels. Both figures carry the mechanical argument of the
case and need to survive a printed page.

### B11. SHOULD-FIX — "44 distinct" is a counting artefact
§8 ("the ledger holds 60 `approaches_recorded` entries (44 distinct)").

60 is right — I summed the six entries (9, 19, 10, 7, 7, 8). Distinct by exact
`approach` text is **51**. 44 is what you get by deduplicating on the first 60
characters, which merges the cycle-2 and cycle-3 restatements of the same
rejection (whose *reasons* differ and which the report elsewhere treats as
separate acts).

**Correction requested.** "60 entries, 51 distinct approach texts (44 distinct
opening clauses, because several rejections are restated with updated reasons)".

### B12. MINOR — a last-digit inconsistency with the report's own table
§6.3: "the electronic differences are 0.3257, 0.5600 and **1.4566** kcal/mol",
against the table's +1.4565. The value is 1.4565436. Correction: 1.4565.

### B13. MINOR — typography in the figures
Both `figB3` and `figB1` text blocks use "A" for Å and "cm^-1" for cm⁻¹. Correct
for print.

### B14. MINOR — figure B1's two energy panels use different y-ranges
Path A tops at 17.5 kcal/mol, path B above 20, so the two profiles cannot be
compared at a glance. The caption does not ask the reader to compare them, so
this is cosmetic — a shared axis would still be better.

### B15. MINOR — an out-of-plane distance whose definition is not given
§4.4 quotes "the ester carbonyl lies 0.22 Å out of the forming-ring plane". It
reproduces only under a particular convention: carbonyl carbon against the plane
of C4/C5/N3 gives 0.222 Å; against the best-fit plane of all five forming-ring
atoms it gives 0.112 Å. The report verifies the accompanying dihedrals (−4.1° /
168.6°) independently and correctly. Correction: state which plane, or mark the
number as the agent's with its convention unrecorded.

### B16. MINOR — the delegated approval is disclosed only in §10
Same as A15: §8's "one displayed human decision" should name
`claude-owner-delegated-reviewer` where the invariant is claimed.

---

## Things I re-derived and could **not** reproduce

Everything else I checked reproduced, most of it to machine precision. These did
not:

| # | Report's claim | Where | My value |
|---|---|---|---|
| 1 | `u-floor-without-method-kcal` = 0.845452 from "the entropy-model term (0.3564) and the … allowance (0.75)" | B §6.3, figB4(b) title, §7 caption | RSS(0.3564, 0.75) = **0.830393**; the claim's own three components RSS to 0.845452 (finding B2) |
| 2 | "60 `approaches_recorded` entries (**44 distinct**)" | B §8 | 60 total, **51** distinct by text (44 only if deduplicated on the first 60 characters) |
| 3 | "each row sums to 1.000 e within **2 × 10⁻⁶**" | A §7, figure A2 | max deviation **3.0 × 10⁻⁶** (B3LYP Mulliken, 0.9999969999999999) |
| 4 | "the electronic differences are 0.3257, 0.5600 and **1.4566**" | B §6.3 | **1.4565436** (the report's own table says 1.4565) |
| 5 | "the ester carbonyl lies 0.22 Å out of the forming-ring plane" | B §4.4 (agent quote) | **0.112 Å** against the 5-atom best-fit ring plane; 0.222 Å only against the C4/C5/N3 plane (finding B15) |
| 6 | "**Four** of the eleven engine calls — 26 896 s, 79 %" bought the negative result | B §0.2, §7 | **two** calls, **26 521 s**, **78.4 %** (finding B5) |
| 7 | five budget terms as "**1σ half-width**" | A §6, figA3(b) axis | the record calls them "the full spread" / "carried in full" / "the largest displacement"; the claim carries no `uncertainty_combination` at all (finding A2) |

Items 3, 4 and 5 are trivial; 1, 2, 6 and 7 are the substantive ones.

For the record, a partial list of what **did** reproduce, since a referee should
say what held: the entire ino3 uncertainty budget and its quadrature
(0.146959 / 0.078593 / 0.012906 / 0.008000 / 0.000983 → 0.167348, and the linear
sum 0.247442, all to ~1 × 10⁻¹⁵); the basis shift 0.0785927 from the four total
energies; the couple −0.472742 → −0.394149 V and the interval
[−0.561497, −0.226801]; the quartet gap 0.219078 eV and its 0.108760 eV T·S
component, from the host's own thermochemistry receipts rather than ORCA's
printed values (a distinction that matters — ORCA's printed Gibbs gives 0.2513 eV
and would have looked like a discrepancy); all twelve spin populations and three
⟨S²⟩ values in the figure-A2 table; the four ino3 structures' Ni–S, Ni–P, S–Ni–S
and P–Ni–P to three decimals; the three ferrocene constants and their 0.008 V gap
in `literature_constants.py`; both po3 scan tables point by point including the
retreating second forming bond and the C≡C series; the level ladder
(+0.3257/+0.5600/+1.4565 electronic; −0.4556/−0.2213/+0.6752 signed; spread
1.1308); the 353 K thermal difference −0.781329; both absolute barriers
(23.353043, 23.808633); every saddle forming bond, cross contact, regio margin
and imaginary mode; the qRRHO shift 0.356444 and the byte-identical cycle-5 bug;
the product-label transposition (independently, from the perceived graph); Guess
B's 1.770 Å N···C(CF₃) and 1.580 Å N···F clash; the 1.964 Å N···F in the rejected
dual-contact compose; the pre-saddles' exact 2.100 Å minimum interfragment
distance; the −4.1° / 168.6° ester–ring dihedrals; the byte-identical scan input
geometries differing only in the driven-coordinate line; `passed: 4, aborted: 4`
at cycle 1 and the `RI-MP2 needs an AuxC basis` abort at cycle 5; both envelopes
(43 200 s, 40 calls, 2 excursion calls, 12 revisions,
`claude-owner-delegated-reviewer`); and both settlements with their full
`decision_uncertainties`.

---

## What the reports do well

Specifically, and these are not small:

1. **The quoting is exact.** I checked every blockquote mechanically against the
   whole raw corpus on an alphanumeric-normalised basis: **31 of 31** in Case A
   and **50 of 50** in Case B are verbatim once the reports' own marked ellipses
   and `[ledger string cut]` markers are honoured. The only deviations are an
   HTML-subscript transliteration (`N<sub>benzyl</sub>` → `N_benzyl`) and two
   field-listing composites the reports present as field listings. Given how easy
   it is to strengthen an agent's hedge while paraphrasing, this discipline is
   the single thing that most earns a reader's trust — and it is exactly why
   findings B8 and A14 are worth fixing rather than shrugging at.
2. **Both figure scripts are genuinely reproducible.** Re-running each
   `make_figures.py` regenerated all seven PNGs **byte-identically**, and both
   print the tables they plot so a reader can diff numbers without reading
   matplotlib.
3. **Errors are led with, not buried.** Case B's §0.3 and §0.4 name two mistakes
   the run made and caught, and §0.4 records that the correction made the answer
   *worse* (uncertainty 1.1796 → 1.2323). Case A's §0.4 flags that cycle 2's
   claims never reached the workspace record. That is the right instinct and it
   is what makes findings A3 and A4 disappointing rather than fatal.
4. **Figure B1 is the best thing in either report.** Points numbered in sequence,
   the host-flagged grid-edge maximum in red with the anomaly's own signal id,
   the interior minimum in green, the second forming bond and the C≡C length on
   a twin axis below, and — decisively for the brief's criterion 6 — "no point on
   this surface was carried forward" written *on the figure*. It answers "which
   point was selected and why" with "none, and here is why", which is harder to
   draw than a selected point.
5. **Distinguishing host arithmetic from model judgement is done consistently.**
   Case A quotes the agent's own "That last sentence is my chemical judgement,
   not a computed quantity"; Case B labels the ruthenium caution as prior
   knowledge and says "I verified only that the label is there; I did not verify
   the chemistry, and neither did the case." Case A's §7 A2 is explicit that the
   functional dependence of the populations "is my reading of its own executed
   outputs, not a claim on its record."
6. **The calculation trajectory is followable in both** (criterion 3). Each cycle
   gets a section that says what it bought, why it followed the one before, and
   what it cost in engine calls. Case B's cycle 1 → 2 → 3 chain — scan the wrong
   coordinate, diagnose it from four distances and a bond length, replace it by
   cutting the product in half — reads as chemistry, not as a changelog.
7. **Case A §9 and Case B §9 exist at all.** A section headed "what I could not
   establish" is rare and both use it honestly (Case B's admission that it cannot
   break the regiochemical tie is the correct answer, not a hedge).

---

## Scientific concerns in the underlying work
*(kept separate from the writing; these are my judgements as a chemist, not
defects in the reports except where noted)*

**C1 — Case A: the whole answer rests on one functional's geometry, and it is the
outlier functional.** All three species were optimised at PBE0-D3BJ/def2-SVP, so
the adiabatic couple, the quartet gap and every thermochemical term inherit PBE0
geometry. The case's own three-functional table shows PBE0 is the most localising
of the three (Ni spin 0.810 vs B3LYP 0.670 and TPSSh 0.624), and the functional
uncertainty was measured *vertically* on those PBE0 geometries. Nothing tested
whether a less localising functional gives the same structure. This matters
physically because the doublet cation's geometry is not a mild perturbation: the
square plane buckles from S–Ni–S 179.9° / P–Ni–P 180.0° to **151.0° / 158.0°** on
one-electron oxidation (I reproduced both). Four-coordinate Ni(III) is genuinely
uncommon — most Ni(III) is five- or six-coordinate — and a 29° buckle is a large
second-order distortion for a non-degenerate d7 state. The Ni–P *lengthening*
(2.217 → 2.258 Å) is chemically sensible for a d7 centre with less density to
back-donate into the phosphine acceptor orbitals, so that part is fine. But the
distortion magnitude is exactly the sort of quantity that moves with exact
exchange, and it feeds the EPR question the chemist actually asked: the g-tensor
of a near-planar d7 and of a substantially flattened one are not the same
measurement. See finding A5 — the report should at least name this as an absent
term.

**C2 — Case A: the delivered spin population is the extreme of every level
computed in the workspace.** Delivered: Ni 0.809 (PBE0/def2-SVP Mulliken). Also
available in the same workspace: 0.670 (B3LYP/SVP), 0.624 (TPSSh/SVP), and
**0.708** (PBE0/def2-TZVP at the same geometry — I read it from
`nicat-sp-pme3`). So basis alone moves it 0.10 e and functional 0.19 e, and the
delivered number is the most metal-localised of the four. The report's figure-A2
conclusion gives the functional-only range ("62 % to 81 %") and puts the basis
effect in a closing parenthesis. The honest combined range for "how much of the
hole is on nickel" is ≈ 0.62–0.81 e, and the qualitative answer the chemist
needs — metal-centred, not a thiyl radical — is robust across all of it (no
sulfur exceeds 0.19 e anywhere, against the ~0.5 e a Ni(II)–S• would put on one
sulfur). **Requested:** state the combined basis-and-functional range where §6
states the spin answer, not only the functional range in §7.

**C3 — Case A: the reference-scale argument covers the solvent and not the
method.** The settlement argues the CPCM systematic is "common-mode to the
subtraction" because the same continuum model enters the ferrocene reference.
That is true of the solvation model and false of the electronic-structure method:
the registered 4.988 V is a G3(MP2)-RAD-Full-TZ/COSMO-RS number and the complex
is PBE0-D3BJ/CPCM, so nothing about the two legs' method error cancels. This is
the physical reason finding A1 matters rather than a bookkeeping quibble.

**C4 — Case B: the absolute barriers are quoted at an ideal-gas 1 atm standard
state.** ΔG‡ = 23.353 and 23.809 kcal/mol come from host thermochemistry whose
own recorded assumption is "ideal-gas translational standard state at 1 atm". For
a bimolecular association the 1 atm → 1 M correction at 353 K is
RT ln(RT/P°) ≈ **+2.4 kcal/mol**, so a solution-phase comparison would place
these near 25.8 kcal/mol. It cancels exactly in the ΔΔG‡ the chemist asked for,
and the agent labelled the absolutes "context rather than the deliverable" — but
§9's only caveat on them is a 4.43 cm⁻¹ mode in the free ynoate, which is an
order of magnitude smaller effect than the standard state. **Requested:** add the
standard-state clause to §9's caveat on the absolute barriers.

**C5 — Case B: the closing catalyst caution needs qualifying, and the report
correctly says it did not check it — so I did.** The agent's caution is that
"ruthenium gives 1,5-disubstituted triazoles and copper the 1,4 pattern, so
'straight to ruthenium' may deliver the isomer you do not want". Two problems on
this substrate. First, CuAAC requires a *terminal* alkyne; it does not operate on
a fully substituted alkyne like methyl 4,4,4-trifluorobut-2-ynoate, so the copper
half of the comparison is inapplicable rather than merely unhelpful. Second, for
an internal alkyne both ring carbons are substituted, so the 1,4/1,5 descriptor
does not pick out a product the way it does for a terminal alkyne: what RuAAC
controls here is *which* of the two substituents ends up adjacent to N1, and that
is set by the substituents' steric and electronic asymmetry, not by a fixed rule.
The chemist's instinct to reach for ruthenium is sound (RuAAC is the route that
works on internal alkynes); the reasoning offered for the caution is not. **This
is a defect in the underlying work, not in the report**, which flags the sentence
as unverified prior knowledge — but since the report chose to reproduce it, one
qualifying clause is owed.

**C6 — Case B: the cycle-1 mechanistic prediction is correctly reasoned and the
`indeterminate` grade is right.** For a conjugated ynoate the LUMO's larger
coefficient sits on the carbon distal from the carbonyl — here the CF₃-bearing
terminus — so nucleophilic attack by the terminal azide nitrogen there, placing
CF₃ at C4 and the ester at C5, is the standard conjugate-addition reading, and
the σ-induction-versus-π-resonance distinction the declaration draws is the right
one. It also contradicts the person who asked, was written before any
calculation, and carries a `method_resolution` of 1.0 kcal/mol that correctly
denied it credit for landing at 0.675. That one declaration is the genuinely
impressive act in either case, which is exactly why finding B1 matters: putting
it in a table beside four post-hoc restatements dilutes it.

**C7 — both cases: the levels of theory are defensible and neither report
oversells them.** B3LYP-D3BJ/def2-SVP geometries and frequencies for a
charge-transfer-flavoured asynchronous cycloaddition TS, and PBE0-D3BJ/def2-SVP
for a first-row d-manifold ionisation, are ordinary starting levels; both reports
say so and neither claims accuracy. Case A's own statement that the delivered
value "sits at the edge of its own functional spread, and the spread is one-sided"
is the correct thing to say about a number like this.

---

## One observation for the commissioning scientist (instrument, not manuscripts)

The `declared_after_evidence` flag behaves differently across the two windows.
In `po3-r19` (harness `e08c34fe`) the host stamps it `true` on every diagnostic
declared in cycles 3, 5 and 6, correctly scoping it to the goal rather than the
session — I verified the declarations were each the first typed act of their own
session, so session scoping would have missed them all. In `ino3-r12` (harness
`d692818c`, per `seals/SUFFICIENCY-3-INO3-R12.md`) the flag never appears, and
`u-reference-scale-pme3` — declared at cycle 3 over two cycles of existing
receipts — is not stamped. If the goal-scoping repair landed between the two
windows, then `ino3-r12`'s record is simply pre-repair, and any comparison of the
two cases' pre-registration discipline has to say so. That is the substance of
findings A8 and A9: without a harness pin in Case A's provenance table, a reader
comparing the two reports will attribute a harness difference to the agent.
