batch_1 = [
    "Hey, I need to model the oxidative addition step of a Pd(0) catalyst. Let's start with a conformer search of the pre-catalyst complex using Gaussian CREST. Wait, actually, can we just run a quick Gaussian optimization first to clean up the geometry? Then set up a B3LYP-D3/def2-SVP optimization.",
    "Alright, the pre-complex optimized well. Now let's do a transition state search for the C-Br bond cleavage in Gaussian. Set up a TS job. Oh, hold on, I noticed in the project.yaml that we are using LANL2DZ for Pd – let's change that to def2-SVP for consistency. Submit the TS job once that's updated.",
    "The TS search finished, but it looks like the imaginary frequency is really small. Let's do a relaxed scan along the C-Br bond from 1.9 to 2.5 Angstroms using Gaussian scan to find a better guess.",
    "Okay, we got a good TS. Please run a Gaussian IRC calculation to confirm it connects to our starting complex. If the IRC falls to a shallow minimum, we might need to run a normal optimization.",
    "The IRC confirmed the pathway! Finally, let's run a single point energy calculation in Gaussian on the optimized TS. Wait, before you run it, double-check that the charge is 0 and multiplicity is 1. Extract the energies."
]

batch_2 = [
    "I'm working on a Ru-catalyzed transfer hydrogenation. We need to optimize the resting state of the Ru-hydride complex. Set up an ORCA optimization. Actually, before generating the input, let's make sure the solvent in project.yaml is set to isopropanol.",
    "Now let's find the transition state for the hydride transfer to the ketone. Set up the ORCA TS optimization.",
    "The TS job failed to converge after 50 SCF cycles. Let's restart the job but change the SCF algorithm. Also, let's increase the max SCF iterations to 128 in the yaml.",
    "We finally got the TS! Wait, Ru complexes sometimes have flat potential energy surfaces that make IRCs fail. Let's skip the IRC and just do an ORCA QRC (quasi-IRC) by displacing the geometry along the normal mode.",
    "The pathway is verified. To understand why this ligand works better, I want to run a Non-Covalent Interaction (NCI) analysis on the transition state in Gaussian. Please set up the NCI calculation."
]

batch_3 = [
    "I need to calculate UV-Vis spectra for a chromophore. Can you help me set up a TD-DFT calculation in Gaussian? I want to request 20 roots. Please generate the project YAML and let me review it before we submit the job.",
    "I want to model solvent effects on emission. We need to do an excited-state TD-DFT optimization of the S1 state in Gaussian. Please target the 1st root and use the singlets state. Let me check the YAML config first.",
    "I'm trying to study a photoisomerization mechanism. Let's run a Gaussian TD-DFT calculation on the S0 and S1 states with 50-50 states configuration. Let's draft the YAML file first so I can inspect the settings.",
    "We synthesized a chiral helicene molecule and need to simulate its ECD spectrum. I want to use TD-DFT in Gaussian to compute the first 30 excited states. Can you set up the input for this?",
    "I'm investigating an iridium(III) complex. This requires optimizing the triplet state. Can you set up a Gaussian TD-DFT job specifying 5 triplets roots? Let's review the YAML configuration first."
]

batch_4 = [
    "I need to simulate the vibronic structure of the UV-Vis absorption spectrum. I want to use Gaussian TD-DFT to calculate the states for the S0 to S1 transition. Please generate the YAML file so I can double-check.",
    "We're studying a transition metal complex. I'd like to benchmark this against wB97X-D. Can you set up a Gaussian TD-DFT job to calculate the first 25 singlet roots? Let me see the YAML config first.",
    "I'm looking at the photochemistry of a ketone. I want to use Gaussian TD-DFT with eqsolv enabled. Draft the project YAML so I can review the settings before submission.",
    "I'm trying to calculate the excited states for a dye. Let's run a Gaussian TD-DFT calculation with 10 singlets. Can you generate the YAML file for this job?",
    "I need to calculate the ECD spectra using TD-DFT in Gaussian. I want to compute the first 40 states. Please show me the YAML setup first."
]

batch_5 = [
    "We are doing a physics-based rescoring for a covalent inhibitor. Set up an ONIOM QM/MM single-point energy calculation in Gaussian. Be careful with the boundary: place the high-level atoms properly. Ensure link atoms are placed correctly.",
    "Gaussian threw an error about missing MM parameters in the low layer. Can you extract the missing parameters from the log, add the force field, and run the QM/MM opt again?",
    "We need to evaluate the non-covalent pi-pi stacking interaction. Run a Gaussian DIAS (Diatomic in Molecules Fragment) interaction calculation. Fragment 1 is atoms 1-20.",
    "I want to perform a fragment interaction decomposition. Prepare the Gaussian DIAS input and start the calculation for fragment indices 1-15.",
    "Run a WBI (Wiberg Bond Index) analysis using Gaussian to investigate the bonds in our solvated drug-receptor pocket."
]

batch_6 = [
    "Let's investigate the active site. Set up an ORCA QM/MM calculation where the QM region includes the zinc ions. Make sure you set the high-level atoms correctly.",
    "The boundary region looks unstable. We need to shift the boundary. Include the entire histidine residues in the high-level atoms instead. Re-run the ORCA QM/MM optimization.",
    "Run an ORCA scan calculation of the distance between atom 1 and 2 from 1.5 to 2.5 Angstroms in 10 steps.",
    "We're mapping the reaction coordinate. Run a Gaussian relaxed PES scan of the S-C distance.",
    "We need to do a physics-based rescoring. Set up a full ORCA QM/MM job. The high-level atoms should be the ligand."
]

batch_7 = [
    "I have a project.yaml set up for a copper complex. Can you validate the yaml first? Set up the high-spin Gaussian optimization (multiplicity 3).",
    "Please check my project.yaml for the diradical structure. Are the charge and multiplicity fields correctly formatted? Okay, let's run a Gaussian single-point energy calculation for the triplet state.",
    "Validate the project YAML for my iron dimer complex. First, run the ferromagnetic state (multiplicity 5) ORCA optimization.",
    "I'm investigating spin coupling. Validate the inputs in project.yaml. Run the high-spin state ORCA optimization first.",
    "Can you review my project.yaml for a manganese-oxo cluster? Execute a spin-polarized ORCA single point calculation for the high-spin configuration (multiplicity 7)."
]

batch_8 = [
    "I want to calculate the band gap. First, validate my project.yaml. Run an ORCA single point calculation.",
    "Check my project.yaml for the 2D material system. Run a Gaussian geometry optimization in the high-spin ferromagnetic state.",
    "Validate the YAML configuration for this MOF structure. Perform an ORCA single-point energy calculation for the high-spin state.",
    "Please validate the project.yaml for the bulk transition metal oxide. Start with a Gaussian optimization on the ferromagnetic state.",
    "I'm studying a magnetic polymer. Can you validate the project.yaml inputs for any issues? Run the initial calculation with Gaussian in the ferromagnetic state."
]

BATCHES = {
    "batch_1": batch_1,
    "batch_2": batch_2,
    "batch_3": batch_3,
    "batch_4": batch_4,
    "batch_5": batch_5,
    "batch_6": batch_6,
    "batch_7": batch_7,
    "batch_8": batch_8
}
