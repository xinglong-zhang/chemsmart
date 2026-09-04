# Schema-gap proposal

The seed cases repeatedly contain scientifically meaningful settings that are not project-configurable today.

Highest-priority candidates:
1. pressure-coupling type (`pcoupltype`) for membrane semi-isotropic coupling;
2. mixed-force-field / parameter provenance;
3. heating or annealing schedules;
4. electrostatics/cutoff controls (`coulombtype`, `rcoulomb`, `rvdw`, `cutoff_scheme`, `nstlist`, `pbc`);
5. explicit box dimensions/vectors;
6. stage-specific restraints and LINCS details;
7. special-protocol settings such as NEMD driving pressure.

Rule: do not expose every possible `.mdp` keyword. Promote a setting when it recurs, changes scientific behavior, and users/agents should be able to choose it.
