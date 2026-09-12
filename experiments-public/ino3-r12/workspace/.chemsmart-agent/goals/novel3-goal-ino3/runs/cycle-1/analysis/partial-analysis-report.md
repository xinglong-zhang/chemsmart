# Host-validated structured analysis

Completion receipt: `a99408bec99461012253ff9816046fa6f570bd8c2466ae8f2da4cf2edf7c2950`
Toolchain plan: `e5335dbedc7774ce5dd07268af3350953d2f94228cb9329f461a54e72e1065db`

Partial analysis: 23 finding(s) over 30 analysis nodes; 2 claim record(s) rendered.

## Analysis nodes that did not execute

- ni0-opt-ph3 (calculation): did not validate
- ni0-opt-pme3 (calculation): did not validate
- ni0-sp-ph3 (calculation): did not validate
- ni0-sp-pme3 (calculation): did not validate
- nicat-opt-ph3 (calculation): did not validate
- nicat-sp-ph3 (calculation): did not validate
- niqrt-opt-ph3 (calculation): did not validate
- niqrt-opt-pme3 (calculation): did not validate
- ext-sp-ni0-ph3 (result_extraction): failed -- expected exactly one registered 'orca_output' result artifact for producer 'ni0-sp-ph3'; found 0
- ext-sp-ni0-pme3 (result_extraction): failed -- expected exactly one registered 'orca_output' result artifact for producer 'ni0-sp-pme3'; found 0
- ext-sp-nicat-ph3 (result_extraction): failed -- expected exactly one registered 'orca_output' result artifact for producer 'nicat-sp-ph3'; found 0
- thermo-ni0-ph3 (thermochemistry): failed -- !! ERROR: Detected imaginary frequencies in geometry optimization for /home/chemsmart/agent-campaigns/ax41-refine-100/novel-round-3/workspaces/ino3-nickel-thiolate-redox/nodes/ni0-opt-ph3/nickel-bis-thiolate-bis-phosphine_opt_opt.out. A valid optimized geometry should not contain imaginary frequencies. Please re-optimize the geometry to locate a true minimum.
- thermo-ni0-pme3 (thermochemistry): failed -- !! ERROR: Detected imaginary frequencies in geometry optimization for /home/chemsmart/agent-campaigns/ax41-refine-100/novel-round-3/workspaces/ino3-nickel-thiolate-redox/nodes/ni0-opt-pme3/geom-pme3v2-h18_opt_opt.out. A valid optimized geometry should not contain imaginary frequencies. Please re-optimize the geometry to locate a true minimum.
- thermo-nicat-ph3 (thermochemistry): failed -- !! ERROR: Detected imaginary frequencies in geometry optimization for /home/chemsmart/agent-campaigns/ax41-refine-100/novel-round-3/workspaces/ino3-nickel-thiolate-redox/nodes/nicat-opt-ph3/nickel-bis-thiolate-bis-phosphine_opt_opt.out. A valid optimized geometry should not contain imaginary frequencies. Please re-optimize the geometry to locate a true minimum.
- expr-redox-ph3 (quantity_expression): skipped -- upstream analysis did not execute: ext-sp-ni0-ph3, ext-sp-nicat-ph3, thermo-ni0-ph3, thermo-nicat-ph3
- expr-redox-pme3 (quantity_expression): skipped -- upstream analysis did not execute: ext-sp-ni0-pme3, thermo-ni0-pme3
- claims-redox (claim_rendering): skipped -- upstream analysis did not execute: expr-redox-ph3, expr-redox-pme3
- validate-freq-ni0-ph3 (scientific_validation): failed -- scientific validation input is not typed evidence from its planned producer
- validate-freq-ni0-pme3 (scientific_validation): failed -- scientific validation input is not typed evidence from its planned producer
- validate-freq-nicat-ph3 (scientific_validation): failed -- scientific validation input is not typed evidence from its planned producer
- validate-freq-nicat-pme3 (scientific_validation): failed -- scientific validation input is not typed evidence from its planned producer
- validate-state-nicat-ph3 (scientific_validation): failed -- scientific validation input is not typed evidence from its planned producer
- validate-state-nicat-pme3 (scientific_validation): failed -- scientific validation input is not typed evidence from its planned producer

## Quantities the results do not carry

| Quantity | Selector | Why the result does not carry it |
|---|---|---|
| hirshfeld_atomic_charges | hirshfeld_atomic_charges | orca result contains no 'hirshfeld_atomic_charges' value. This result resolves: alpha_homo, alpha_lumo, auxiliary_basis, auxiliary_basis_role, basis, beta_homo, beta_lumo, charge, connectivity, converged, dipole_moment, dipole_moment_magnitude, dispersion_energy, effective_multiplicity, energies, energy, entropy_times_temperature, functional, gap, gibbs_free_energy, homo, loewdin_atomic_charges, lumo, mulliken_atomic_charges, multiplicity, positions, reference_energy, scf_energy, solvation_cavity_surface_area, solvation_electrostatic_energy, solvation_model, solvent, spin_square, spin_square_deviation, spin_square_target, symbols, vibrational_frequencies, vibrational_mode_atom_participation, vibrational_mode_degeneracy_group |
| effective_multiplicity | effective_multiplicity | this result prints no <S^2> expectation value; a spin-restricted closed-shell calculation is an eigenfunction of S^2 by construction, so read 'multiplicity' for its electronic state. This result resolves: alpha_homo, alpha_lumo, auxiliary_basis, auxiliary_basis_role, basis, beta_homo, beta_lumo, charge, connectivity, converged, dipole_moment, dipole_moment_magnitude, dispersion_energy, energies, energy, entropy_times_temperature, functional, gap, gibbs_free_energy, homo, loewdin_atomic_charges, lumo, mulliken_atomic_charges, multiplicity, positions, reference_energy, scf_energy, solvation_cavity_surface_area, solvation_electrostatic_energy, solvation_model, solvent, spin_square_target, symbols, vibrational_frequencies, vibrational_mode_atom_participation, vibrational_mode_degeneracy_group |
| hirshfeld_atomic_charges | hirshfeld_atomic_charges | orca result contains no 'hirshfeld_atomic_charges' value. This result resolves: alpha_homo, alpha_lumo, auxiliary_basis, auxiliary_basis_role, basis, beta_homo, beta_lumo, charge, connectivity, converged, dipole_moment, dipole_moment_magnitude, dispersion_energy, energies, energy, entropy_times_temperature, functional, gap, gibbs_free_energy, homo, loewdin_atomic_charges, lumo, mulliken_atomic_charges, multiplicity, positions, reference_energy, scf_energy, solvation_cavity_surface_area, solvation_electrostatic_energy, solvation_model, solvent, spin_square_target, symbols, vibrational_frequencies, vibrational_mode_atom_participation, vibrational_mode_degeneracy_group |
| hirshfeld_atomic_charges | hirshfeld_atomic_charges | orca result contains no 'hirshfeld_atomic_charges' value. This result resolves: alpha_homo, alpha_lumo, auxiliary_basis, auxiliary_basis_role, basis, beta_homo, beta_lumo, charge, connectivity, converged, dipole_moment, dipole_moment_magnitude, dispersion_energy, effective_multiplicity, energies, energy, entropy_times_temperature, functional, gap, gibbs_free_energy, homo, loewdin_atomic_charges, lumo, mulliken_atomic_charges, multiplicity, positions, reference_energy, scf_energy, solvation_cavity_surface_area, solvation_electrostatic_energy, solvation_model, solvent, spin_square, spin_square_deviation, spin_square_target, symbols, vibrational_frequencies, vibrational_mode_atom_participation, vibrational_mode_degeneracy_group |
| effective_multiplicity | effective_multiplicity | this result prints no <S^2> expectation value; a spin-restricted closed-shell calculation is an eigenfunction of S^2 by construction, so read 'multiplicity' for its electronic state. This result resolves: alpha_homo, alpha_lumo, auxiliary_basis, auxiliary_basis_role, basis, beta_homo, beta_lumo, charge, connectivity, converged, dipole_moment, dipole_moment_magnitude, dispersion_energy, energies, energy, entropy_times_temperature, functional, gap, gibbs_free_energy, homo, loewdin_atomic_charges, lumo, mulliken_atomic_charges, multiplicity, positions, reference_energy, scf_energy, solvation_cavity_surface_area, solvation_electrostatic_energy, solvation_model, solvent, spin_square_target, symbols, vibrational_frequencies, vibrational_mode_atom_participation, vibrational_mode_degeneracy_group |
| hirshfeld_atomic_charges | hirshfeld_atomic_charges | orca result contains no 'hirshfeld_atomic_charges' value. This result resolves: alpha_homo, alpha_lumo, auxiliary_basis, auxiliary_basis_role, basis, beta_homo, beta_lumo, charge, connectivity, converged, dipole_moment, dipole_moment_magnitude, dispersion_energy, energies, energy, entropy_times_temperature, functional, gap, gibbs_free_energy, homo, loewdin_atomic_charges, lumo, mulliken_atomic_charges, multiplicity, positions, reference_energy, scf_energy, solvation_cavity_surface_area, solvation_electrostatic_energy, solvation_model, solvent, spin_square_target, symbols, vibrational_frequencies, vibrational_mode_atom_participation, vibrational_mode_degeneracy_group |

## Thermochemical conditions

| Stage | Temperature (K) | Standard state | Entropy model | Frequency scale |
|---|---:|---|---|---:|
| thermo-ni0-ph3 | `298.15` | 1 mol/L | rrho | `1` |
| thermo-ni0-pme3 | `298.15` | 1 mol/L | rrho | `1` |
| thermo-nicat-ph3 | `298.15` | 1 mol/L | rrho | `1` |
| thermo-nicat-pme3 | `298.15` | 1 mol/L | rrho | `1` |

## Literature constants

| Constant | Value | Unit | Family | Convention |
|---|---:|---|---|---|
| standard_hydrogen_electrode_absolute_potential_kelly2006 | `4.28` | V | tissandier1998_cluster_pair_proton_scale | absolute potential of the standard hydrogen electrode, 298.15 K, solutes at 1 mol/L and ideal gases at 1 bar, proton and electron on Boltzmann statistics; consistent with an intrinsic proton solvation free energy of -1105 kJ/mol, from which the surface potential is absent |

Claim record: `40e5598930bf504c794a3ecd4e564b8ec1c8b2fea335d5f33351363747705cbc`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| quartet-doublet-gap | `0.3546500680407081` | eV | `ba2563ac879b0fcd14d928180f9ef839d3fe312844bbe92177f8833619066090` |
| quartet-doublet-gap-ph3 | `0.35769508381908266` | eV | `ba2563ac879b0fcd14d928180f9ef839d3fe312844bbe92177f8833619066090` |

Claim record: `6b8667fa781c25002e403a6c52b5c1cea6b4738cbb42d8b2c45c79c2258cddd0`

## Host-rendered numerical claims

| Claim | Value | Unit | Source receipt |
|---|---:|---|---|
| d-ni-s1-cat-pme3 | `2.164832687808691` | angstrom | `c8942ed68cceb2ec44738f84e23ed86bae816ebd153e061d35beb8e21370d92d` |
| d-ni-s1-neut-pme3 | `2.218030592060668` | angstrom | `c8942ed68cceb2ec44738f84e23ed86bae816ebd153e061d35beb8e21370d92d` |
| d-ni-s2-cat-pme3 | `2.165138650427958` | angstrom | `c8942ed68cceb2ec44738f84e23ed86bae816ebd153e061d35beb8e21370d92d` |
| d-ni-s2-neut-pme3 | `2.2179360719414345` | angstrom | `c8942ed68cceb2ec44738f84e23ed86bae816ebd153e061d35beb8e21370d92d` |
| dq-ni-ph3 | `-0.034632` | e | `43af831331a5303f7cc98a6f67062ee017a470a63866f9193090eb0647f9bbb0` |
| dq-ni-pme3 | `-0.09814600000000001` | e | `f3f725b747d9ac88cd1fe08e4e7701ae6ffd3e6acd190da421b5bbd61bb4f18c` |
| dq-s1-ph3 | `-0.20242300000000002` | e | `43af831331a5303f7cc98a6f67062ee017a470a63866f9193090eb0647f9bbb0` |
| dq-s1-pme3 | `-0.17123499999999997` | e | `f3f725b747d9ac88cd1fe08e4e7701ae6ffd3e6acd190da421b5bbd61bb4f18c` |
| dq-s2-ph3 | `-0.202439` | e | `43af831331a5303f7cc98a6f67062ee017a470a63866f9193090eb0647f9bbb0` |
| dq-s2-pme3 | `-0.17088800000000004` | e | `f3f725b747d9ac88cd1fe08e4e7701ae6ffd3e6acd190da421b5bbd61bb4f18c` |
| spin-square-ph3 | `0.767524` | 1 | `6e86f3d35755b701b13ace8859c743a6690e2e4cd6fc12ae9defcaa01540055a` |
| spin-square-pme3 | `0.774059` | 1 | `17f584cdb3b46aa86912651ccaa5c15a3e887b0c6920c40206280e034517ea78` |
| spin-square-qrt-ph3 | `3.762007` | 1 | `179dbc8c698c3081b37dc62f0f4044d5393e4d0ac1abfc1f1d9852bf4fb348ae` |
| spin-square-qrt-pme3 | `3.764746` | 1 | `0bd0d63dcafffd6b6e20ef3183b460bcf96badfd7a4d67d4c7e9f1206c560ccd` |

## Surviving receipts (evidence, not claims)

| Quantity | Value | Unit | Evidence rung | Source receipt |
|---|---:|---|---|---|
| effective_multiplicity | `4.0073662173552345` | 1 | parsed | `0bd0d63dcafffd6b6e20ef3183b460bcf96badfd7a4d67d4c7e9f1206c560ccd` |
| energy | `-3304.657464221637` (= `-3304.657464221637` Eh as printed) | hartree | parsed | `0bd0d63dcafffd6b6e20ef3183b460bcf96badfd7a4d67d4c7e9f1206c560ccd` |
| multiplicity | `4` | 1 | parsed | `0bd0d63dcafffd6b6e20ef3183b460bcf96badfd7a4d67d4c7e9f1206c560ccd` |
| spin_square | `3.764746` | 1 | parsed | `0bd0d63dcafffd6b6e20ef3183b460bcf96badfd7a4d67d4c7e9f1206c560ccd` |
| effective_multiplicity | `4.005999001497629` | 1 | parsed | `179dbc8c698c3081b37dc62f0f4044d5393e4d0ac1abfc1f1d9852bf4fb348ae` |
| energy | `-3069.131212461491` (= `-3069.131212461491` Eh as printed) | hartree | parsed | `179dbc8c698c3081b37dc62f0f4044d5393e4d0ac1abfc1f1d9852bf4fb348ae` |
| multiplicity | `4` | 1 | parsed | `179dbc8c698c3081b37dc62f0f4044d5393e4d0ac1abfc1f1d9852bf4fb348ae` |
| spin_square | `3.762007` | 1 | parsed | `179dbc8c698c3081b37dc62f0f4044d5393e4d0ac1abfc1f1d9852bf4fb348ae` |
| charge | `1` | 1 | parsed | `17f584cdb3b46aa86912651ccaa5c15a3e887b0c6920c40206280e034517ea78` |
| effective_multiplicity | `2.023916006162311` | 1 | parsed | `17f584cdb3b46aa86912651ccaa5c15a3e887b0c6920c40206280e034517ea78` |
| energy | `-3304.670497371353` (= `-3304.670497371353` Eh as printed) | hartree | parsed | `17f584cdb3b46aa86912651ccaa5c15a3e887b0c6920c40206280e034517ea78` |
| mulliken_atomic_charges | `[0.134358,-0.269438,-0.172913,0.100104,0.10758,0.11547,-0.269742,-0.172983,0.100175,0.107565,0.115455,0.151968,0.151684,-0.209046,-0.193069,-0.178373,-0.208999,-0.193181,-0.1783,0.109485,0.116675,0.110099,0.108152,0.103742,0.113985,0.108698,0.103467,0.106585,0.109439,0.116681,0.110056,0.108152,0.103771,0.113949,0.108705,0.103448,0.106598]` | e | parsed | `17f584cdb3b46aa86912651ccaa5c15a3e887b0c6920c40206280e034517ea78` |
| multiplicity | `2` | 1 | parsed | `17f584cdb3b46aa86912651ccaa5c15a3e887b0c6920c40206280e034517ea78` |
| positions | `[[0.000409,1e-05,0.046015],[2.092439,0.135679,0.585927],[2.519631,-1.251124,1.671166],[2.650082,-2.196884,1.128419],[3.479899,-0.987284,2.137895],[1.76889,-1.374441,2.463334],[-2.091205,-0.135231,0.588865],[-2.51856,1.251132,1.674689],[-2.649594,2.19727,1.132648],[-3.478656,0.986554,2.141355],[-1.767795,1.374337,2.466855],[0.073177,2.21519,-0.384512],[-0.072988,-2.215642,-0.38462],[1.493808,2.577108,-1.462289],[-1.345348,2.782012,-1.368396],[0.264615,3.32968,1.031566],[-1.493196,-2.574493,-1.464026],[1.344968,-2.78436,-1.368281],[-0.266151,-3.330109,1.031161],[1.470361,3.643526,-1.733827],[1.42608,1.961716,-2.370629],[2.432087,2.350217,-0.938012],[-1.208969,3.844162,-1.622359],[-2.287384,2.647651,-0.823083],[-1.380992,2.18398,-2.290115],[0.414031,4.356403,0.665153],[1.144886,3.004376,1.604143],[-0.620807,3.294541,1.678605],[-1.472514,-3.641035,-1.735292],[-1.422174,-1.959577,-2.372446],[-2.431624,-2.344675,-0.941287],[1.206479,-3.846023,-1.623141],[2.287211,-2.652539,-0.822722],[1.382026,-2.185714,-2.289543],[-0.415409,-4.356917,0.664927],[-1.146676,-3.004532,1.603178],[0.618961,-3.294962,1.678678]]` (= `[[0.000409,1e-05,0.046015],[2.092439,0.135679,0.585927],[2.519631,-1.251124,1.671166],[2.650082,-2.196884,1.128419],[3.479899,-0.987284,2.137895],[1.76889,-1.374441,2.463334],[-2.091205,-0.135231,0.588865],[-2.51856,1.251132,1.674689],[-2.649594,2.19727,1.132648],[-3.478656,0.986554,2.141355],[-1.767795,1.374337,2.466855],[0.073177,2.21519,-0.384512],[-0.072988,-2.215642,-0.38462],[1.493808,2.577108,-1.462289],[-1.345348,2.782012,-1.368396],[0.264615,3.32968,1.031566],[-1.493196,-2.574493,-1.464026],[1.344968,-2.78436,-1.368281],[-0.266151,-3.330109,1.031161],[1.470361,3.643526,-1.733827],[1.42608,1.961716,-2.370629],[2.432087,2.350217,-0.938012],[-1.208969,3.844162,-1.622359],[-2.287384,2.647651,-0.823083],[-1.380992,2.18398,-2.290115],[0.414031,4.356403,0.665153],[1.144886,3.004376,1.604143],[-0.620807,3.294541,1.678605],[-1.472514,-3.641035,-1.735292],[-1.422174,-1.959577,-2.372446],[-2.431624,-2.344675,-0.941287],[1.206479,-3.846023,-1.623141],[2.287211,-2.652539,-0.822722],[1.382026,-2.185714,-2.289543],[-0.415409,-4.356917,0.664927],[-1.146676,-3.004532,1.603178],[0.618961,-3.294962,1.678678]]` Angstrom as printed) | angstrom | parsed | `17f584cdb3b46aa86912651ccaa5c15a3e887b0c6920c40206280e034517ea78` |
| spin_square | `0.774059` | 1 | parsed | `17f584cdb3b46aa86912651ccaa5c15a3e887b0c6920c40206280e034517ea78` |
| spin_square_deviation | `0.024059000000000053` | 1 | parsed | `17f584cdb3b46aa86912651ccaa5c15a3e887b0c6920c40206280e034517ea78` |
| vibrational_frequencies | `[28.3,47.79,52.71,82.29,93.84,111.57,114.71,136.76,137.16,143.58,149.46,158.14,165.2,169.83,171.5,180.31,186.87,196.19,198.15,199.5,205.47,208.19,217.41,218.93,251.65,268.24,272.72,281.08,291.71,308.02,336.35,375.01,406.7,685.29,685.33,734.63,735.77,754.43,754.91,769.92,770.39,793.11,793.21,843.67,847.79,847.97,850.22,934.17,940.44,942.08,946.21,947.78,954.66,954.72,955.92,956.29,962.02,1282.97,1283.0,1285.09,1285.8,1295.98,1299.24,1318.35,1321.37,1398.82,1398.89,1400.47,1400.48,1405.81,1405.89,1410.28,1410.34,1411.17,1411.37,1417.31,1417.71,1431.05,1431.19,1439.23,1439.42,3064.62,3064.64,3066.71,3066.86,3070.36,3070.41,3072.26,3072.31,3176.23,3176.47,3181.37,3181.4,3184.45,3184.46,3185.62,3185.68,3191.86,3192.18,3194.97,3195.01,3207.36,3207.67,3211.91,3211.94]` | cm^-1 | parsed | `17f584cdb3b46aa86912651ccaa5c15a3e887b0c6920c40206280e034517ea78` |
| charge | `0` | 1 | parsed | `418b88ca5a0aeb2a1b03e60807678289ee62014e4b47a3a464ed728bcd1fd2f9` |
| energy | `-3069.316727578493` (= `-3069.316727578493` Eh as printed) | hartree | parsed | `418b88ca5a0aeb2a1b03e60807678289ee62014e4b47a3a464ed728bcd1fd2f9` |
| mulliken_atomic_charges | `[-0.003298,-0.421204,-0.162329,0.064323,0.083345,0.083349,-0.421213,-0.162339,0.06431,0.083336,0.083358,0.153077,0.055231,0.072947,0.072949,0.153034,0.055208,0.072958,0.072959]` | e | parsed | `418b88ca5a0aeb2a1b03e60807678289ee62014e4b47a3a464ed728bcd1fd2f9` |
| multiplicity | `1` | 1 | parsed | `418b88ca5a0aeb2a1b03e60807678289ee62014e4b47a3a464ed728bcd1fd2f9` |
| positions | `[[3.8e-05,-1.6e-05,-0.0009],[2.204461,-0.088802,-0.000342],[2.997251,-1.731917,0.001088],[2.302548,-2.581776,0.001178],[3.634527,-1.827398,-0.890709],[3.63343,-1.826337,0.893782],[-2.204426,0.088893,-0.001218],[-2.997183,1.732039,0.001859],[-2.302471,2.581894,0.001432],[-3.63537,1.827788,-0.889258],[-3.632446,1.826182,0.895233],[0.260824,2.156237,-0.001726],[1.589892,2.646307,-0.002963],[-0.274996,2.896514,-1.083914],[-0.27325,2.897147,1.080893],[-0.260761,-2.156367,-0.001242],[-1.58979,-2.646558,-0.003306],[0.276002,-2.897241,-1.082529],[0.271719,-2.896588,1.082641]]` (= `[[3.8e-05,-1.6e-05,-0.0009],[2.204461,-0.088802,-0.000342],[2.997251,-1.731917,0.001088],[2.302548,-2.581776,0.001178],[3.634527,-1.827398,-0.890709],[3.63343,-1.826337,0.893782],[-2.204426,0.088893,-0.001218],[-2.997183,1.732039,0.001859],[-2.302471,2.581894,0.001432],[-3.63537,1.827788,-0.889258],[-3.632446,1.826182,0.895233],[0.260824,2.156237,-0.001726],[1.589892,2.646307,-0.002963],[-0.274996,2.896514,-1.083914],[-0.27325,2.897147,1.080893],[-0.260761,-2.156367,-0.001242],[-1.58979,-2.646558,-0.003306],[0.276002,-2.897241,-1.082529],[0.271719,-2.896588,1.082641]]` Angstrom as printed) | angstrom | parsed | `418b88ca5a0aeb2a1b03e60807678289ee62014e4b47a3a464ed728bcd1fd2f9` |
| vibrational_frequencies | `[-187.47,-186.51,-22.14,54.01,75.69,134.97,141.84,163.73,179.68,217.44,226.03,227.23,253.78,296.27,302.51,356.3,406.81,475.33,516.63,546.8,577.25,729.15,730.51,936.51,939.97,947.18,947.61,977.8,1011.12,1100.89,1101.98,1102.52,1104.36,1327.57,1329.83,1424.09,1424.29,1440.76,1440.94,2471.34,2473.16,2491.4,2491.73,2493.99,2494.14,3059.63,3059.92,3156.47,3156.52,3179.37,3179.64]` | cm^-1 | parsed | `418b88ca5a0aeb2a1b03e60807678289ee62014e4b47a3a464ed728bcd1fd2f9` |
| dq-ni-ph3 | `-0.034632` | e | derived | `43af831331a5303f7cc98a6f67062ee017a470a63866f9193090eb0647f9bbb0` |
| dq-s1-ph3 | `-0.20242300000000002` | e | derived | `43af831331a5303f7cc98a6f67062ee017a470a63866f9193090eb0647f9bbb0` |
| dq-s2-ph3 | `-0.202439` | e | derived | `43af831331a5303f7cc98a6f67062ee017a470a63866f9193090eb0647f9bbb0` |
| charge | `1` | 1 | parsed | `6e86f3d35755b701b13ace8859c743a6690e2e4cd6fc12ae9defcaa01540055a` |
| effective_multiplicity | `2.0174478927595625` | 1 | parsed | `6e86f3d35755b701b13ace8859c743a6690e2e4cd6fc12ae9defcaa01540055a` |
| energy | `-3069.144357513474` (= `-3069.144357513474` Eh as printed) | hartree | parsed | `6e86f3d35755b701b13ace8859c743a6690e2e4cd6fc12ae9defcaa01540055a` |
| mulliken_atomic_charges | `[0.031334,-0.218781,-0.164913,0.093189,0.117154,0.117153,-0.218774,-0.164926,0.093185,0.117144,0.117165,0.207872,0.101228,0.115727,0.115743,0.20782,0.101214,0.115714,0.115752]` | e | parsed | `6e86f3d35755b701b13ace8859c743a6690e2e4cd6fc12ae9defcaa01540055a` |
| multiplicity | `2` | 1 | parsed | `6e86f3d35755b701b13ace8859c743a6690e2e4cd6fc12ae9defcaa01540055a` |
| positions | `[[-2e-05,-3.5e-05,-0.001221],[2.147239,-0.0967,-0.000397],[2.955563,-1.720354,0.000542],[2.287286,-2.588652,0.00085],[3.591081,-1.770953,-0.894918],[3.590733,-1.770055,0.896301],[-2.147298,0.096784,-0.000771],[-2.955525,1.720498,0.001151],[-2.287216,2.588777,0.000943],[-3.591779,1.771185,-0.893779],[-3.589945,1.770159,0.897444],[0.273809,2.201244,-0.00107],[1.621082,2.626814,-0.001802],[-0.273798,2.888547,-1.106489],[-0.272553,2.888454,1.105025],[-0.273792,-2.201388,-0.00067],[-1.621044,-2.627036,-0.002058],[0.274487,-2.889343,-1.105346],[0.271689,-2.887945,1.106266]]` (= `[[-2e-05,-3.5e-05,-0.001221],[2.147239,-0.0967,-0.000397],[2.955563,-1.720354,0.000542],[2.287286,-2.588652,0.00085],[3.591081,-1.770953,-0.894918],[3.590733,-1.770055,0.896301],[-2.147298,0.096784,-0.000771],[-2.955525,1.720498,0.001151],[-2.287216,2.588777,0.000943],[-3.591779,1.771185,-0.893779],[-3.589945,1.770159,0.897444],[0.273809,2.201244,-0.00107],[1.621082,2.626814,-0.001802],[-0.273798,2.888547,-1.106489],[-0.272553,2.888454,1.105025],[-0.273792,-2.201388,-0.00067],[-1.621044,-2.627036,-0.002058],[0.274487,-2.889343,-1.105346],[0.271689,-2.887945,1.106266]]` Angstrom as printed) | angstrom | parsed | `6e86f3d35755b701b13ace8859c743a6690e2e4cd6fc12ae9defcaa01540055a` |
| spin_square | `0.767524` | 1 | parsed | `6e86f3d35755b701b13ace8859c743a6690e2e4cd6fc12ae9defcaa01540055a` |
| spin_square_deviation | `0.017523999999999984` | 1 | parsed | `6e86f3d35755b701b13ace8859c743a6690e2e4cd6fc12ae9defcaa01540055a` |
| vibrational_frequencies | `[-170.01,-168.44,24.97,77.75,98.07,134.84,145.97,166.42,187.26,240.6,242.14,252.6,259.7,292.87,315.88,383.11,423.03,488.1,511.71,557.52,578.5,732.8,734.86,926.4,929.05,940.64,955.87,956.15,975.46,1092.2,1093.82,1093.96,1095.78,1325.36,1329.38,1414.56,1414.58,1425.22,1425.44,2495.57,2495.76,2526.1,2526.22,2533.04,2533.1,3075.51,3075.54,3182.85,3182.92,3208.83,3209.3]` | cm^-1 | parsed | `6e86f3d35755b701b13ace8859c743a6690e2e4cd6fc12ae9defcaa01540055a` |
| charge | `1` | 1 | parsed | `7aac1ead740a7dabc4c237050bf3eb45ff5accd0022ec6fe2ad9f10f79fc6704` |
| energy | `-3305.741363794941` (= `-3305.741363794941` Eh as printed) | hartree | parsed | `7aac1ead740a7dabc4c237050bf3eb45ff5accd0022ec6fe2ad9f10f79fc6704` |
| multiplicity | `2` | 1 | parsed | `7aac1ead740a7dabc4c237050bf3eb45ff5accd0022ec6fe2ad9f10f79fc6704` |
| charge | `0` | 1 | parsed | `aba443d76804cbfabccfed5ebabc1d82c9f4c4b6822b642641caef04e1b9be6f` |
| energy | `-3304.819331805318` (= `-3304.819331805318` Eh as printed) | hartree | parsed | `aba443d76804cbfabccfed5ebabc1d82c9f4c4b6822b642641caef04e1b9be6f` |
| mulliken_atomic_charges | `[0.036212,-0.440673,-0.170251,0.068774,0.079324,0.079266,-0.44063,-0.170234,0.06876,0.079308,0.079294,0.149634,0.149677,-0.209493,-0.202255,-0.202297,-0.209488,-0.202234,-0.202314,0.086698,0.087384,0.087376,0.090015,0.095848,0.098269,0.090018,0.098271,0.095946,0.086704,0.087386,0.087378,0.089995,0.095844,0.098261,0.090029,0.098261,0.095939]` | e | parsed | `aba443d76804cbfabccfed5ebabc1d82c9f4c4b6822b642641caef04e1b9be6f` |
| multiplicity | `1` | 1 | parsed | `aba443d76804cbfabccfed5ebabc1d82c9f4c4b6822b642641caef04e1b9be6f` |
| positions | `[[-0.000118,-4.3e-05,-0.000419],[2.217682,-0.031957,0.00167],[3.201538,-1.571071,0.005445],[2.6254,-2.502711,0.007829],[3.846744,-1.581333,-0.886041],[3.846271,-1.576606,0.897311],[-2.217825,0.031808,0.000844],[-3.201438,1.571078,0.005852],[-2.625042,2.502549,0.007275],[-3.847668,1.581295,-0.88489],[-3.845084,1.576932,0.898497],[0.196074,2.205844,-0.001566],[-0.196236,-2.205855,-0.001517],[1.8337,3.019588,-0.00602],[-0.538666,2.97706,-1.489528],[-0.531189,2.977916,1.489723],[-1.833776,-3.019707,-0.006125],[0.538675,-2.976813,-1.489521],[0.531145,-2.978078,1.48964],[1.678793,4.10952,-0.007335],[2.400433,2.721294,-0.897684],[2.403925,2.72418,0.884373],[-0.641185,4.065135,-1.360324],[-1.506636,2.529354,-1.739354],[0.155852,2.77669,-2.319278],[-0.635928,4.065737,1.359976],[0.168199,2.779375,2.315776],[-1.497534,2.529743,1.744967],[-1.678867,-4.109636,-0.007495],[-2.400388,-2.721352,-0.897852],[-2.404128,-2.724365,0.884207],[0.640609,-4.064992,-1.360723],[1.506951,-2.52949,-1.738805],[-0.155406,-2.77577,-2.319478],[0.636376,-4.06582,1.359667],[-0.16843,-2.780102,2.315682],[1.497178,-2.529398,1.745219]]` (= `[[-0.000118,-4.3e-05,-0.000419],[2.217682,-0.031957,0.00167],[3.201538,-1.571071,0.005445],[2.6254,-2.502711,0.007829],[3.846744,-1.581333,-0.886041],[3.846271,-1.576606,0.897311],[-2.217825,0.031808,0.000844],[-3.201438,1.571078,0.005852],[-2.625042,2.502549,0.007275],[-3.847668,1.581295,-0.88489],[-3.845084,1.576932,0.898497],[0.196074,2.205844,-0.001566],[-0.196236,-2.205855,-0.001517],[1.8337,3.019588,-0.00602],[-0.538666,2.97706,-1.489528],[-0.531189,2.977916,1.489723],[-1.833776,-3.019707,-0.006125],[0.538675,-2.976813,-1.489521],[0.531145,-2.978078,1.48964],[1.678793,4.10952,-0.007335],[2.400433,2.721294,-0.897684],[2.403925,2.72418,0.884373],[-0.641185,4.065135,-1.360324],[-1.506636,2.529354,-1.739354],[0.155852,2.77669,-2.319278],[-0.635928,4.065737,1.359976],[0.168199,2.779375,2.315776],[-1.497534,2.529743,1.744967],[-1.678867,-4.109636,-0.007495],[-2.400388,-2.721352,-0.897852],[-2.404128,-2.724365,0.884207],[0.640609,-4.064992,-1.360723],[1.506951,-2.52949,-1.738805],[-0.155406,-2.77577,-2.319478],[0.636376,-4.06582,1.359667],[-0.16843,-2.780102,2.315682],[1.497178,-2.529398,1.745219]]` Angstrom as printed) | angstrom | parsed | `aba443d76804cbfabccfed5ebabc1d82c9f4c4b6822b642641caef04e1b9be6f` |
| vibrational_frequencies | `[-147.74,-146.39,-24.92,25.25,53.8,97.34,98.82,109.16,127.45,138.68,141.66,153.06,163.03,173.93,177.65,203.68,207.28,208.44,229.53,231.09,232.97,233.96,241.92,246.76,267.53,284.19,311.93,313.88,319.97,322.42,346.75,362.03,390.68,681.12,682.56,718.67,720.24,736.39,736.69,753.83,754.79,799.77,799.85,836.38,839.91,841.43,844.66,924.85,927.04,931.83,936.95,939.82,946.46,948.3,948.5,948.53,960.12,1278.48,1279.33,1279.43,1279.45,1294.28,1296.12,1313.67,1316.44,1398.69,1398.82,1406.17,1406.53,1406.55,1406.98,1410.24,1410.72,1411.04,1411.77,1427.73,1428.04,1442.17,1442.79,1445.82,1445.99,3058.02,3058.3,3066.04,3066.06,3066.44,3066.52,3067.24,3067.26,3148.2,3148.24,3171.04,3171.08,3171.48,3171.51,3181.42,3181.45,3199.24,3199.36,3201.27,3201.29,3217.55,3217.57,3218.43,3218.45]` | cm^-1 | parsed | `aba443d76804cbfabccfed5ebabc1d82c9f4c4b6822b642641caef04e1b9be6f` |
| dg-ph3-ev | `0.01314505198297411` (= `0.35769508381908266` eV as requested) | hartree | derived | `ba2563ac879b0fcd14d928180f9ef839d3fe312844bbe92177f8833619066090` |
| dg-pme3-ev | `0.013033149716193293` (= `0.3546500680407081` eV as requested) | hartree | derived | `ba2563ac879b0fcd14d928180f9ef839d3fe312844bbe92177f8833619066090` |
| d-ni-s1-neut-pme3 | `2.218030592060668` | angstrom | derived | `c8942ed68cceb2ec44738f84e23ed86bae816ebd153e061d35beb8e21370d92d` |
| d-ni-s2-neut-pme3 | `2.2179360719414345` | angstrom | derived | `c8942ed68cceb2ec44738f84e23ed86bae816ebd153e061d35beb8e21370d92d` |
| d-ni-s1-cat-pme3 | `2.164832687808691` | angstrom | derived | `c8942ed68cceb2ec44738f84e23ed86bae816ebd153e061d35beb8e21370d92d` |
| d-ni-s2-cat-pme3 | `2.165138650427958` | angstrom | derived | `c8942ed68cceb2ec44738f84e23ed86bae816ebd153e061d35beb8e21370d92d` |
| electronic_energy | `-3304.6704973713527` | hartree | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| zero_point_energy | `0.30833688870978565` | hartree | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| internal_energy | `-3304.3381460534756` | hartree | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| enthalpy | `-3304.337201868921` | hartree | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| entropy_times_temperature | `0.0750954092783307` | hartree | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| gibbs_free_energy | `-3304.4122972781993` | hartree | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| thermal_internal_energy_correction | `0.33235131787696354` | hartree | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| thermal_enthalpy_correction | `0.3332955024317726` | hartree | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| enthalpy_increment_above_zero_point | `0.024958613721986912` | hartree | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| thermal_gibbs_correction | `0.25820009315339565` | hartree | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| entropy | `0.0002518712368885819` | hartree K^-1 | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| temperature | `298.15` | K | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| pressure | `1.0` | atm | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| near_zero_mode_count | `0` | 1 | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| heat_capacity_cv | `0.00013254266460511536` | hartree K^-1 | derived | `e0b742911fff75742d4d74e79863d46f54a3f85aa05f1ff606a593dcd910cd84` |
| dq-ni-pme3 | `-0.09814600000000001` | e | derived | `f3f725b747d9ac88cd1fe08e4e7701ae6ffd3e6acd190da421b5bbd61bb4f18c` |
| dq-s1-pme3 | `-0.17123499999999997` | e | derived | `f3f725b747d9ac88cd1fe08e4e7701ae6ffd3e6acd190da421b5bbd61bb4f18c` |
| dq-s2-pme3 | `-0.17088800000000004` | e | derived | `f3f725b747d9ac88cd1fe08e4e7701ae6ffd3e6acd190da421b5bbd61bb4f18c` |

## Expected versus delivered

| Observable | Expected | Delivered | Unit | Agreement | Basis |
|---|---|---:|---|---|---|
| quartet-doublet-gap | positive 0.05..1.5 | `0.3546500680407081` | eV | agreed | Square-planar-derived Ni thiolate cations are normally low-spin doublets; a competitive quartet would require large exchange stabilization from metal-localized d electrons. |
| quartet-doublet-gap-ph3 | positive 0.05..1.5 | `0.35769508381908266` | eV | agreed | Same low-spin expectation as for the PMe3 model; PH3 is a worse donor, which should if anything stabilize the low-spin form further. |
| redox-potential-vs-fc | positive 0.0..0.7 | `` |  | not_comparable | The user's shelf oxidants (Fc+ 0.00 V, acetylferrocenium +0.27 V, tris(4-bromophenyl)aminium +0.70 V) were chosen to bracket the likely couple; Ni(II) bis(thiolate) bis(phosphine) complexes typically oxidize within this window in acetonitrile. |
| redox-potential-vs-fc-ph3 | positive 0.0..1.0 | `` |  | not_comparable | Weaker phosphine donors destabilize the oxidized (cationic) form, so the PH3 model should oxidize at higher (more positive) potential than the PMe3 model; the magnitude of the shift is being measured. |

An expectation is displayed, never scored: a diverging row settles nothing and means the chemistry disagreed with the reasoning, which is a result the reader owns.

Scientific decision: not recorded -- interpretation is a session act; this run executed extraction, thermochemistry, expressions, validation verdicts, and claim rendering only.

Recovery: the engine outputs and every receipt above are unchanged and remain readable. A later explicit analysis request -- a new session over this workspace planning an analysis-only toolchain on the registered results -- may re-extract, validate, and claim them without re-running any engine.

## Record delivery

| Record | Nodes (reached state) | Calculation | Analysis |
| --- | --- | --- | --- |
| 1 | ni0-opt-ph3 (engine_complete), ni0-sp-ph3 (not_attempted) | not_delivered | partial |
| 2 | nicat-opt-ph3 (engine_complete), nicat-sp-ph3 (not_attempted) | not_delivered | partial |
| 3 | niqrt-opt-ph3 (engine_complete) | not_delivered | executed |
| 4 | ni0-opt-pme3 (engine_complete), ni0-sp-pme3 (not_attempted) | not_delivered | partial |
| 5 | nicat-opt-pme3 (validated), nicat-sp-pme3 (validated) | validated | partial |
| 6 | niqrt-opt-pme3 (engine_complete) | not_delivered | executed |
| shared | claims-context (executed), claims-gaps (executed), claims-redox (skipped), expr-dists-pme3 (executed), expr-gaps (executed), expr-redox-ph3 (skipped), expr-redox-pme3 (skipped), expr-spin-context-ph3 (executed), expr-spin-context-pme3 (executed), spin-pop-analysis (blocked_unsupported) | - | partial |

A batch of N is N observations; each verdict above is one record's, and no aggregate quantity is rendered.
