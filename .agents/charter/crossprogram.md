# Charter topic: crossprogram

> Evidence, not instruction: what this surface was qualified by and what it found, under the
> conditions stated. The live tree and ``chemsmart agent capabilities`` outrank it; doctrine is
> ``AGENTS.md``, ``CONDUCT.md`` and ``.agents/rsl/``.

This surface is qualified through one completed multi-program execution:
a single displayed approval covering seven nodes — one xTB optimisation
feeding three ORCA and three PySCF single points on that geometry at
three electronic states — all validated, with the geometry handoff
crossing programs twice and preserving atom identity and order. PySCF
answered through the merged plane on those live results, serving method,
basis, energies, orbital energies, the spin diagnostic and per-atom
Mulliken charges that close on each record's own formal charge; a
frequency selector is served for ``hess`` and refused for ``sp`` by the
job type detected from the stored spec. What the merge itself claimed
was no new execution surface: PySCF ``sp/opt/hess`` was already
release-qualified, and the response and correlated stages described
above came in a later round on top of it.

The first charged, open-shell xTB runs under the Agent found a defect
the merge had carried: the xTB result audit merged the project's
*resolved* settings over the bound identity, and a project that declares
no state resolves to charge 0, multiplicity 1, so three correct runs —
an anion, a neutral radical, a radical cation, each with the right
``--chrg`` and ``--uhf`` in its own receipt — were typed as state
mismatches while the neutral one passed. The bound identity is now the
only authority for charge and multiplicity in every program's audit; a
project field participates only when it is explicit.

A geometry may cross programs and a number may not follow it freely. The
optimised-geometry handoff is keyed on the producing program and refuses
any change of atom identity or order, so an xTB optimisation feeding an
ORCA or PySCF single point preserves that parent atom *i* is child atom
*i*. A typed value carries its unit and its dimension and not the method
that produced it, and the arithmetic checks only those, so a
tight-binding energy and a hybrid-DFT energy subtract without complaint.
That is displayed and never refused: a high-level single point on a
low-level geometry is an ordinary protocol and composite methods mix
levels deliberately, while two energies from one program at different
basis sets are equally unsubtractable and no program-identity check would
catch them. The displayed analysis chain names the level of theory behind
every input, resolved through the analysis chain rather than one hop, so
the reviewer decides whether a mixture is a method or a mistake.

Naming the level is necessary and it is not sufficient, and the release
says so rather than implying otherwise. In the qualifying run both
programs were asked for ``b3lyp`` with one basis and the chain displayed
two identical strings, while the total energies differed by 0.24 hartree
on every species: one program's ``B3LYP`` uses the VWN5 local
correlation and the other's uses VWN3, so five identical characters name
two functionals. The offset was nearly constant across the three charge
states and largely cancelled in the differences, which is why the derived
indices still agreed to about 0.1 eV while no total energy agreed at all.
A refusal would have compared the same two strings and been no wiser. The
host shows what the project asked each program for; whether two programs
mean the same thing by a keyword is a fact about the programs, and
checking it stays the scientist's.

A solvation term is comparable across ORCA and PySCF only where the model
is: both readers serve ``solvation_electrostatic_energy`` in hartree for
the same physical term, and the two programs' continua differ in cavity
construction — at one geometry and one dielectric PySCF's SMD and C-PCM
polarisation energies differ by 3.6 kcal/mol — so the model and the
solvent are read beside the number, never inferred from it. The frontier
pair of an unrestricted result is the extremum over both spin channels on
ORCA and on PySCF alike; Gaussian refuses an open-shell frontier value
rather than report one channel, and xTB serves the gap it prints from its
single orbital ladder.

An open-shell singlet is one typed request (`broken_symmetry: true`,
R10 Q18). Each program writes it as the mechanism measured to reach the
broken-symmetry solution on sites with one unpaired electron each:
Gaussian U + guess=mix; ORCA UHF with GuessMix
45; PySCF following its own RKS -> UKS instability, because PySCF's
documented `init_guess_breaksym` and an explicit HOMO/LUMO mix left H2 at
2.00 A and p-benzyne restricted (oracle O0, CUHK Slurm 2153330). The
request and its verification carry across programs; the particular SCF
solution does not. Each result's level states the reference that ran
(rks/uks) and whether symmetry broke (<S**2> against its target, with a
threshold), and a request that did not break raises
`spin.broken_symmetry_request_unbroken`. At the exactly degenerate
90-degree ethylene twist, Gaussian's U guess=mix reached a solution 26
kcal/mol above ORCA's with the same <S**2> ~ 1.00, while Gaussian's
stable=opt reached the lowest (O0b, 2153375). Equal <S**2> is therefore
not equal solutions, and no sensor yet compares them. On two S = 1
Ni(II) centres ORCA's GuessMix reached a state 33.40 mEh above the one a
flip of one centre from the quintet reaches (R11 truth-2, CUHK 2157086);
that flip is ORCA's typed `site_spin_flip: {atoms, final_ms}` on the
high-spin multiplicity, read at the Ms it converged to and measured
against S = |Ms| (R11 E2, CUHK 2157103); Gaussian and PySCF refuse it.
"broken" says the symmetry broke, not that an intended multi-site state
was reached.

A dispersion correction is one word only where each program means one
correction by it, and is written only where the program has parameters
for the pair (R10 Q20). Census D ran every functional + dispersion pair
the writers produce on the programs themselves (CUHK 2153534,
2153546): green previews died on the program's own lookup in Gaussian
(179, link 301), ORCA (96, after the SCF, past the input-check banner)
and PySCF (197, at run time). The writers now hold those measured tables
and refuse, with routes, exactly the pairs that died, writing every pair
that ran byte for byte; where programs share a pair, their dispersion
energies agree to 1e-9 Eh. The word `d3` was two corrections: ORCA's
bare D3 is D3(BJ), while the hub's Gaussian `gd3` is zero damping. ORCA
now refuses a bare `d3` and asks for `d3zero` or `d3bj`. ORCA's D2
applies a silent default C6 scaling (1.200) to B3LYP/G, the spelling the
hub writes for `b3lyp`, and that pair is refused. wB97X + D3(BJ) has no
same-physics translation: ORCA's and PySCF's are wB97X-V exchange and
correlation with D3(BJ), and Gaussian's wB97XD is its own functional.
PySCF still runs `wb97x` + `d3bj` with wB97X-V-fitted parameters on the
2008 functional, which its own table calls supported.
