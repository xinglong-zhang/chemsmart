---
name: scientific-conventions
version: 0.3.0
description: How computational-chemistry quantities are conventionally defined and reported — direction of every difference quantity, adiabatic versus vertical geometry, which energy terms are included, when an established assignment may be stated, and thermochemistry standard states.
---

# Reporting conventions for computed quantities

Knowledge about how results are *expressed*, not about how accurate they must
be. Nothing here sets a tolerance, an error budget, or a readiness state.

These are general principles. Where an example appears it illustrates the
principle; the principle, not the example, is what applies.

Apply them as a challenging computational scientist. Treat an unexpected
number as a hypothesis-generating observation: verify geometry, electronic
state, energy definition, sign, unit, thermochemical convention and program
semantics before adding guidance. When a real calculation disproves or narrows
a rule, update the general rule and challenge it on a different chemical case.
Never learn a molecule-specific value or preferred DAG as a convention.

## 1. Every difference quantity needs a stated direction

A number that is a difference is meaningless until the reader knows which term
was subtracted. State the direction explicitly rather than relying on a
convention the reader may not share.

| Quantity class | Conventional direction | Sign that follows |
|---|---|---|
| Term value (`Te`, `T0`, excitation, ionization, electron affinity as `IE`-like) | upper/final state − ground state | non-negative by construction |
| Reaction energy / enthalpy / free energy | products − reactants | negative if exergonic |
| Activation barrier | transition state − reactants | positive |
| Interaction / binding energy | complex − separated fragments | negative if bound, but the opposite sign is also in common use — always say which |
| State-ordering gap (e.g. a singlet–triplet gap) | higher-lying state − ground state | non-negative when the ground state is correctly identified |

A term value is a spectroscopic band origin: it is measured upward from the
ground state and cannot be negative. A negative term value means the two states
were ordered the wrong way round, not that the quantity is negative.

## 2. Adiabatic versus vertical is a geometry convention

The distinction applies to any state-to-state quantity — ionization, electron
attachment, excitation, spin-state change — and it changes the shape of the
calculation, not only the number.

| | Geometry of each state | Shape of the calculation |
|---|---|---|
| **Adiabatic** | each state at **its own** relaxed geometry | relax both states, then compare |
| **Vertical** | both states at the **initial** geometry | relax the initial state only; evaluate the other state at that geometry |

A vertical quantity therefore contains **no optimization of the final state**.
Computing two independently relaxed geometries and calling the result vertical,
or using one geometry and calling the result adiabatic, reports something other
than what was asked. For the signed difference
``E_final - E_initial``, the vertical value is normally no smaller because the
final state is evaluated away from its own minimum. A positive electron
affinity is commonly defined with the opposite sign, so its numerical
inequality reverses; always state the direction before comparing values.

## 3. Say which energy terms are included

Distinguish, and name which one is being reported:

- **Electronic energy** — the converged SCF/post-SCF energy at a geometry.
- **Zero-point corrected (0 K)** — electronic plus harmonic zero-point energy.
  A molecular state needs its vibrational analysis, not only its electronic
  energy. A monatomic fragment has no vibrations or rotations and therefore
  has exactly zero vibrational ZPE; do not request a fictitious atomic
  optimization or frequency calculation to manufacture those modes.
- **Enthalpy / free energy at finite T** — electronic plus zero-point plus
  thermal corrections, with `G = H − TS` at the stated temperature.

Zero-point and thermal corrections belong to quantities defined at a relaxed
geometry. A vertical quantity is an energy difference at one fixed geometry and
carries no separate zero-point correction.

## 4. Established assignments may be stated as prior knowledge

Some facts — a ground-state spin multiplicity, a known point group, an
established conformer preference — are settled in the literature and do not
have to be computed before they can be stated. When one applies, state it,
mark it as prior knowledge rather than a computed result, and let the
calculation confirm or contradict it. Declining to state a settled fact is not
caution; it withholds information the requester asked for.

The reverse error matters as much: presenting an unsettled, substituent-
dependent, or method-sensitive assignment as settled. The test is whether the
assignment is established for **this** system, not for a superficially similar
one.

Two generative principles cover most ground-state spin assignments:

- Near-degenerate frontier orbitals follow **Hund's rule** — the high-spin
  configuration lies lowest, because exchange stabilisation outweighs the small
  orbital-energy gap.
- A substituent or ligand field that **splits those orbitals** far enough
  reverses the ordering in favour of the low-spin state. Strong π-donation into
  a formally empty frontier orbital, or a strong-field ligand set, does this.

Apply the principles to the system at hand rather than recalling a list.

## 5. Reference-method limitation for low-spin and stretched states

A single-determinant reference can be less balanced for low-spin,
near-degenerate, or stretched-bond states than for a corresponding high-spin
state. Open-shell singlets, diradicals, transition-metal spin states and
bond-breaking regions therefore require diagnostics rather than a universal
error direction. Depending on the system, reference, functional and amount of
exact exchange, a low-spin state may be biased upward or downward. Establish
the direction from spin diagnostics, reference stability, multireference
evidence or an applicable cited benchmark; do not infer it from spin alone.

## 6. Differences are more convergence-sensitive than absolutes

A difference of two total energies is more sensitive to convergence than either
energy alone, because the two errors cancel only when the states are converged
comparably. This applies to the SCF threshold, the integration grid, and the
geometry-optimization criteria alike. State the settings used rather than
leaving them at whatever default applies: the setting is part of what was
computed, and asymmetric settings between the two states invalidate the
cancellation the difference relies on.

## 7. Thermochemistry

- Partition functions are evaluated in the **rigid-rotor / harmonic-oscillator**
  approximation unless another treatment is named. Low-frequency modes are the
  approximation's weak point; if a quasi-harmonic treatment is used, name it.
- Rotational entropy requires the **rotational symmetry number** σ, taken from
  the molecular point group — not from the formula or the atom count. Omitting
  it inflates the entropy by `R ln σ`.

  | Point group | σ | | Point group | σ |
  |---|---|---|---|---|
  | `C1`, `Cs`, `Ci`, `C∞v` | 1 | | `D2h` | 4 |
  | `C2`, `C2v`, `C2h`, `D∞h` | 2 | | `D3h` | 6 |
  | `C3v` | 3 | | `D3`, `D3h` | 6 |
  | `Td` | 12 | | `Oh` | 24 |

- The modern standard state is **1 bar**. Older tables and some program defaults
  use 1 atm; the two differ by `R ln(1.01325)` in the standard entropy. When a
  request names a pressure, report the value at that pressure and say which one
  the calculation used.
- A reported free energy should be reconstructible from the electronic energy,
  the zero-point energy, the thermal corrections, and `G = H − TS` at the stated
  temperature.
- For a monatomic ideal-gas species confined to its ground electronic level,
  there is no rotational or vibrational contribution. Its molar enthalpy
  increment above the electronic energy is `5/2 RT`: `3/2 RT` translation plus
  `RT` from `pV`.
- Check a finite-temperature reaction enthalpy against its 0 K value by
  subtracting the per-species thermal increments with the same stoichiometric
  signs. Derive that difference from the planned quantities. A remembered
  literature number is an external comparison, not an internal consistency
  equation, unless its source and convention were supplied.
- When starting from a ZPE-corrected 0 K quantity, the increment to finite
  temperature is `H(T) - E_electronic - ZPE`. By contrast,
  `H(T) - E_electronic` already contains ZPE and belongs with an electronic
  energy difference. Do not add both corrections to the same reaction value.

## 8. A number without a unit and a convention is not a result

Report the unit explicitly and keep one unit per quantity within a document.
When converting, state the source unit. Energies computed in Hartree are
conventionally reported in kcal/mol or kJ/mol for chemical differences and in eV
for ionization, attachment, and excitation.
