---
name: typed-analysis-contract
version: 0.1.0
description: How the typed analysis layer expects a workflow to be expressed — what an identifier, a unit, a declared quantity kind and an evidence reference each mean, which parts the host derives for you, and how a rejected node is repaired without re-authoring the graph.
---

# Expressing an analysis so the host can carry the evidence

Knowledge about how the typed analysis layer wants intent *written*, not about
which chemistry to do. Nothing here chooses a method, a route, or a value, and
nothing here overrides the live schemas, which remain authoritative if this
document falls behind them.

The reason to read it: the analysis layer is where a plan turns into a number
with provenance, and its contract is about evidence rather than arithmetic.
Most rejections there are not chemistry errors. They are a mismatch between how
a quantity was described and what the host must be able to prove about it
afterwards.

## 1. Identifiers are labels you choose; the host resolves what they point at

A public identifier — a node id, an output id, an input id, a plan or workflow
id — is a name for one thing inside one request. It is lower case and starts
with a letter, and it exists so that a later reference is unambiguous.

Chemical notation is mixed case and these names are not. Fold the case rather
than dropping letters: a name for a change in Gibbs energy in kilojoules keeps
all of its letters when it is folded down. A compound name whose locants come
first needs a leading word, because a name may not begin with a digit.

These names are presentation. The host grades an identifier-independent
symbolic form, so a readable name costs nothing and buys a legible record.

## 2. A unit is a claim about dimension, and dimensionless is a unit

Every quantity carries a unit because the host checks dimensional consistency
before it will combine anything. Two consequences follow:

- A quantity that genuinely has no dimension — a count, a population or mole
  fraction, an oscillator strength, a degeneracy, a verdict — is dimensionless
  and says so with the dimensionless unit. It is not unitless-by-omission, and
  a word describing what is being counted is not a unit.
- A percentage is a presentation of a fraction, not a separate dimension.
  Choose the dimensionless form and present it however you like afterwards.

Declare the unit where the quantity is declared. A unit supplied only at
rendering time asks the host to check a dimension it was never given.

## 3. A declared quantity kind must be one the producing operation derives

An operation produces a fixed vocabulary of quantities. Naming a kind that
operation does not derive is refused, and the refusal lists what it does
derive — read that list rather than guessing a synonym. In particular, a
thermochemistry stage derives thermochemical quantities; a bare energy label is
not one of them.

## 4. Evidence references, and what the host derives for you

An analysis input points at a quantity inside a producing receipt. The host
must be able to say, afterwards, exactly which occurrence of which quantity
produced which number — that is the whole purpose of the typed layer.

This is why repetition matters. Comparing two structures repeats a source
quantity name by construction: several thermochemistry receipts each yield a
quantity under the same name, so that name alone no longer identifies one
occurrence.

**The host resolves this itself.** Each input already carries an id that is
unique within its request, so an occurrence is identified without you naming
it. An optional semantic role exists only to give an occurrence a more readable
label than its id. If you do supply roles, keep them distinct — two inputs
sharing one role would collapse onto a single slot and make the reference
ambiguous, which is refused.

The general rule this illustrates: **do not hand-author bookkeeping the host
can derive from what you have already told it.** Supply science; let the host
maintain identity.

## 5. Prefer a typed source quantity to a literal

When a number already exists as a quantity on a receipt, reference that
quantity rather than restating its value. A literal is for a constant that has
no typed source, and every literal is a value the host cannot trace. A
condition such as a temperature that a stage already carries is a typed
quantity, not a number to retype.

## 6. Say when something is not supported, rather than dropping it

A requested stage that cannot be materialised stays in the plan, marked as
blocked, with the reason. Deleting it is the cheapest way to make findings
clear and it silently changes the question that was asked. An honest plan and a
complete plan are the same plan, and a plan that says what it could not do is
more useful than one that quietly did less.

## 7. Repairing a rejected node

Rejections name the path, the value that was rejected, and the rule it broke;
several independent problems in one submission are reported together, so read
the whole rejection before revising. Where an inventory is offered — the
selectors a result actually resolves, the quantities a receipt actually
carries, the kinds an operation actually derives — that inventory is the
answer, not a hint.

Repair the named field. A rejection about how something is expressed is not a
reason to change the science, and changing method, structure, state, or the
requested observable in response to a naming rule replaces the question instead
of answering it.
