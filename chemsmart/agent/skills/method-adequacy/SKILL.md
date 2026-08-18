---
name: method-adequacy
version: 0.1.0
description: Judges whether a method, basis set, solvation model or conformer sample can answer the question that was asked - which errors cancel in a comparison and which do not, when an effect is smaller than the method's own error, and where the uncertainty comes from. Use this skill BEFORE reporting any computed difference, comparing a value with experiment, or answering whether a setup is trustworthy, reliable or good enough - including when nothing in the request mentions accuracy. Triggers on comparing two species, conformers or isomers; any charged or open-shell species; interaction, binding, association or non-covalent energies; branching or substituent effects; a double-zeta or otherwise small basis; a shortfall against a measured value; and the words trustworthy, reliable, accurate, adequate, uncertainty or limitation. Do NOT use it for the direction, sign, units or standard state of a reported quantity, which belong to scientific-conventions.
---

# Is this method adequate for this question?

Knowledge about whether a result can carry the weight the question puts on it.
Nothing here sets a tolerance, a pass mark, or a readiness state, and nothing
here overrides a project setting, a validator, or the live CLI.

These are general principles. Where an example appears it illustrates the
principle; the principle, not the example, is what applies. Apply them as a
challenging computational scientist: the point is not to refuse work, it is to
know what the number you produced is worth and to say so before someone builds
on it.

The recurring failure this addresses is quiet: a calculation converges, the
parser validates it, every unit is right, and the answer is reported without
noticing that the method could not have resolved the effect in the first place.
A green pipeline is not evidence of adequacy.

## 1. Ask what size the effect is before choosing the method

The first question is not "what is a good method" but "how large is the thing I
am trying to see, and what is this method's own error on that kind of
quantity". A method whose typical error is comparable to the effect cannot
resolve it, however cleanly the job runs.

- If the expected effect is smaller than the method's characteristic error for
  that property, the sign may still be wrong after perfect convergence. Say so
  rather than reporting the digits.
- Conformational and non-covalent effects are frequently of this size. So are
  many substituent effects and most stereochemical preferences.
- A tighter convergence threshold does not reduce method error. Neither does a
  larger grid. They reduce *numerical noise*, which is a different quantity.

## 2. Error cancellation is the property that decides basis-set adequacy

A basis set is not adequate or inadequate in the abstract; it is adequate for a
particular comparison. What matters is whether the error it makes is the *same*
on both sides of the difference being taken.

Errors cancel well when the two sides are chemically similar: the same
molecular formula, the same bonding pattern, the same charge, the same number
of each kind of orbital. They cancel poorly when the comparison changes the
number of electron pairs, the charge, the bonding, or the spin state.

| Comparison | Cancellation | Consequence |
| --- | --- | --- |
| Conformers or rotamers of one species | strong | a modest basis often suffices for the *ordering* |
| Isomers of one formula | moderate | ordering usually survives; magnitudes may not |
| Reaction changing bond types or counts | weak | basis-set error enters the answer directly |
| Anything changing total charge | weak | see the diffuse-function rule below |
| Complex versus separated fragments | weak, plus superposition | see counterpoise below |

**The failure mode to watch for**: recovering a fraction of a known
experimental effect and treating the shortfall as unremarkable. If a computed
difference is a quarter of the measured one, the method is not describing the
physics being asked about, and that is the headline, not a footnote.

## 3. Specific inadequacies that change answers qualitatively

Each of these changes a conclusion, not a decimal place.

**Diffuse functions and anions.** An anion's extra electron density extends far
beyond the region a valence basis describes. Without diffuse functions the
anion is described worse than the neutral, so any quantity involving both —
deprotonation, electron affinity, anion stability — is biased in a *known
direction* and by an amount that does not cancel. State it whenever a charged
species appears on one side of a comparison. Balanced differences of two
similar deprotonations cancel much of this; absolute values do not.

**Dispersion.** Most common density functionals do not describe London
dispersion from their own functional form. Whenever the question involves
molecules or fragments touching without bonding — association, packing,
folding, branching, stacking, host–guest — an uncorrected functional omits the
dominant attractive term. A dispersion correction is not a refinement there; it
is the physics. Report which correction and damping were actually applied,
taken from the result rather than from intent.

**Basis-set superposition.** When two fragments approach, each borrows the
other's basis functions, so the complex is described better than the separated
pieces and the interaction looks too strong. The effect shrinks as the basis
grows and is largest exactly where interaction energies are small. Counterpoise
correction estimates it; saying which convention was used matters, because
corrected and uncorrected values are both reported in the literature.

**Solvation model choice.** An implicit continuum reproduces bulk polarisation
and nothing else. It does not describe a specific hydrogen bond, a coordinating
solvent molecule, or a tight ion pair. If the chemistry depends on a particular
solvent–solute contact, a continuum alone will miss it regardless of how the
cavity is parameterised. Note also that a solvent-phase single point on a
gas-phase geometry is a different quantity from a solvent-phase optimisation,
and the difference is not always small.

**Multireference character.** Single-reference methods assume one dominant
electron configuration. Stretched bonds, diradicals, some transition metals,
and many bond-breaking transition states violate that. Diagnostics exist and
are worth reading; a method used outside its assumption can produce a smooth,
converged, entirely wrong surface.

**Spin contamination.** For an open-shell unrestricted calculation, a computed
spin expectation value far from its exact value means the wavefunction is not
the state that was requested. Energies from a contaminated wavefunction
describe a mixture, not the named state. A closed-shell restricted calculation
is an eigenfunction by construction, so this diagnostic does not apply and its
absence is not a defect.

## 4. A conformer ensemble is a sample, not a structure

Any property of a flexible molecule at finite temperature is an average over
populated states, so the result depends on which states were found.

- One optimised structure is one local minimum. Calling it "the" structure
  asserts a search that was not performed.
- Populations depend exponentially on relative free energy, so the ranking
  method must resolve differences of the order of thermal energy — which
  returns to the first principle above.
- Distinct conformers must be genuinely distinct; duplicates found twice
  inflate a population silently.
- Symmetry-equivalent wells are separate microstates and must be counted, not
  merged. Omitting a degeneracy shifts a population and every average taken
  from it.
- Whether an average is linear matters: a linear mean is the ensemble value
  only for a property that averages linearly, and a magnitude-of-vector
  observable needs the squares averaged and the root taken.

Report the search that was actually done, including its limits. An ensemble
assembled from hand-supplied structures is a hand-supplied ensemble.

## 5. A driven coordinate starts from the geometry you supply

A relaxed scan is a sequence of constrained optimisations, and the first of
them begins by forcing the driven coordinate to the scan's starting value on
the structure that was handed in. That step is a real geometric operation and
it can fail.

- If the starting value is far from the coordinate's value in the supplied
  geometry, the program has to distort the molecule a long way before any
  optimisation happens. It may refuse outright, or succeed into a strained
  arrangement that is not on the path you meant to map.
- Prefer a scan that begins near the supplied geometry and walks outward, or
  supply a geometry that already sits near the intended start.
- Periodic coordinates are periodic: a torsion scanned across a full turn ends
  where it began, and the closure is a free internal check on the profile.
- The same applies to a held coordinate: a constraint imposed far from the
  current value is a distortion, not a restraint.

When a program reports that it could not impose an initial constraint, that is
this failure and not a convergence problem. Read the coordinate's value in the
input geometry before changing anything else.

## 6. Composite and multi-level results inherit every layer's assumptions

Combining a high-level energy with a lower-level geometry and thermal
correction is standard and defensible, and it carries three separate
assumptions: that the geometry is close enough that the high-level energy is
evaluated at effectively the right point, that the frequencies are good enough
for the corrections taken from them, and that the levels are combined without
double counting.

State which level produced which term. A composite quoted as one number, with
no account of its parts, cannot be checked by anyone.

## 7. Say what the uncertainty is, and where it comes from

An uncertainty that quantifies only convergence understates the real one by a
wide margin. Distinguish:

- **numerical** — convergence thresholds, integration grids, and the residual
  noise in near-zero vibrational modes;
- **method** — the functional or correlation treatment and the basis;
- **model** — rigid-rotor harmonic-oscillator treatment, implicit solvation,
  a truncated conformer set, a neglected environment;
- **reference** — what the experimental comparison actually measured, under
  which conditions, and whether it is the same quantity at all.

The model term usually dominates and is usually the one omitted. Where no
defensible numerical uncertainty is available, say which term dominates and in
which direction it biases the result, rather than reporting a bare value or
inventing an interval.

## 8. Reporting adequacy honestly

State the method that was *applied*, read from the result, not the method that
was requested. State whether the effect being reported is larger than the
method's own error for that kind of quantity. State the dominant limitation and
its direction. If a comparison with experiment falls well short, say by how much
and give the most likely reason.

A result that is correctly computed and inadequate for the question is a
legitimate outcome, and reporting it as such is worth more than a confident
number that will not survive being checked.
