# Worked example: C002 → T003

Source: Supporting Information for DOI `10.1039/D4TB01145G`.

## 1. Extract one stage, not one whole paper
The source gives a 100 ps NVT equilibration and a 50 ns NPT production. They become two records.

## 2. Normalize units for the production stage
- 2 fs = 0.002 ps
- 50 ns = 50,000 ps
- `nsteps = 50,000 / 0.002 = 25,000,000`

## 3. Map only supported fields
The production record maps:
- Amber14SB → `gromacs_settings.force_field`
- TIP3P → `gromacs_settings.water_model`
- 0.002 ps → `timestep`
- 300 K → `temperature`
- 1 bar → `pressure`
- v-rescale → `thermostat`
- Parrinello–Rahman → `barostat`
- hydrogen-bond constraints / LINCS → `constraints` / `constraint_algorithm`
- 25,000,000 → `nsteps`

## 4. Keep unsupported fields visible
PME and the 1.0 nm Coulomb/vdW cutoffs are recorded as schema gaps rather than silently dropped or falsely claimed to be represented in YAML.

See the `Field_Mapping` sheet and the matching `T003_...yaml`.
