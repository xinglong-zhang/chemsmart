---
id: matched-numerics-before-numerics
paths:
  - "chemsmart/jobs/settings.py"
  - "chemsmart/jobs/orca/settings.py"
  - "chemsmart/jobs/gaussian/settings.py"
  - "chemsmart/jobs/pyscf/settings.py"
  - "chemsmart/jobs/*/writer.py"
conditions: "ORCA 6.1.1, Gaussian 16 C.02, PySCF 2.14 on CUHK; closed- and open-shell molecules up to benzene; R10 episode Q2"
evidence:
  - "commit:d94706a2"
  - "commit:f84db238"
  - "note:CUHK Slurm 2149277 and 2149278 (121 matched runs), 2149487 (57 runs on the repaired tree), 2149622 (the D3BJ check)"
repeat_cost: "a real Hamiltonian difference (B3LYP's VWN form, 0.154 Eh in benzene) was handled as a convention and taught to the model for rounds instead of translated; the episode's own PBE0 and BP86 expectations failed where default numerics hid real differences"
falsifier: "a literal that agrees across programs to about 1e-6 Eh at matched tight numerics but whose relative energies at default settings differ by more than 0.05 kcal/mol"
home: prose
supersedes: []
earned: 2026-09-24
last_verified: "2026-09-28 @ 3dc11ef6"
---
Call a difference between programs "numerics" only after the same literal has been run at matched tight numerics in each program (no RI or density fitting, fine grids, tight SCF).
Agreement to about 1e-6 Eh means one Hamiltonian; the spread at default settings is 1e-4 to 3.5e-4 Eh, large enough to hide a real functional difference in totals and small enough to be mistaken for one.
A difference that survives matched numerics is a translation for the host to make, never a convention to teach the model.
Evidence: R10 Q2 oracles (CUHK 2149277/2149278, 2149487): B3LYP's VWN form and BP86's LDA part differed where default numerics had hidden them.
