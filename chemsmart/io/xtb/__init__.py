class XTBRefs:
    XTB_GROUPS = [
        "chrg",  # Set the charge of the molecule
        "spin",  # Set Nalpha-Nbeta of the molecule
        "cma", # Shift molecule to center of mass
        "constrain", # Constrain the gradient by appling potentials
        "cube",
        "external",
        "fix", # Exact fixing
        "gbsa", # Generalized born (GB) model with solvent accessable surface area
        "gfn", # GFN Hamiltonian
        "hess",
        "metadyn",
        "md",
        "modef",
        "oniom",
        "opt",
        "path",
        "scan",
        "scc",
        "split",
        "stm",
        "symmetry",
        "thermo",
        "wall",
        "write"
    ]

    XTB_METHODS = ["GFN0-xTB", "GFN1-xTB", "GFN2-xTB", "GFN-FF"]

    XTB_OPT_LEVEL = [
        "crude",
        "sloppy",
        "loose",
        "lax",
        "normal",
        "tight",
        "vtight",
        "extreme"
    ]

    XTB_OPT_ENGINE = [
        "rf", # Approximate Normal Coordinate Rational Function optimizer
        "lbfgs", # L-BFGS Approximate Normal Coordinate optimizer
        "inertial" # Fast Inertial Relaxation Engine
    ]

    XTB_SOLVENT_MODELS = [
        "gbsa",  # Generalized born model with solvent accessable surface area contributions
        "alpb",  # analytical linearized Poisson-Boltzmann model
    ]