"""
Reference strings and citations for computational chemistry methods.

Contains formatted reference strings for various thermochemistry and
quantum chemistry methods used in ChemSmart. Includes citations for
QRRHO approximations, damping functions, and entropy calculations.
"""

# qrrho_header = f'"   " + "-" * 108 + "\n"'
#         log(
#             "   "
#             + " " * 32
#             + "Quasi-Rigid-Rotor-Harmonic-Oscillator Scheme"
#             + "\n"
#         )
#         log("   " + "-" * 108 + "\n")
#

qrrho_header = (
    "   "
    + "-" * 108
    + "\n"
    + "   "
    + " " * 32
    + "Quasi-Rigid-Rotor-Harmonic-Oscillator Scheme"
    + "\n"
    + "   "
    + "-" * 108
    + "\n"
)

head_gordon_damping_function_ref = (
    "   - Damping function: Chai and Head-Gordon\n"
    + "     REF: Chai, J.-D.; Head-Gordon, M. Phys. Chem. Chem. Phys. 2008, 10, 6615-6620\n\n"
)

grimme_quasi_rrho_entropy_ref = (
    "   - Grimme's quasi-RRHO entropy:\n"
    + "     REF: Grimme, S. J. Chem. Phys. 2011, 134, 064114\n\n"
)

truhlar_quasi_rrho_entropy_ref = (
    "   - Truhlar's quasi-RRHO entropy:\n"
    + "     REF: Ribeiro, R. F.; Marenich, A. V.; Cramer, C. J; Truhlar, D. G. J. Phys. Chem. B 2011, 115, 14556-14562\n\n"
)

head_gordon_quasi_rrho_enthalpy_ref = (
    "   - Head-Gordon's quasi-RRHO enthalpy:\n"
    + "     REF: Li, Y.; Gomes, J.; Sharada, S. M.; Bell, A. T.; Head-Gordon, M. J. Phys. Chem. C 2015, 119, 1840-1850\n\n"
)
