"""Generate isotope data for ChemSmart from NIST isotope information.

This script parses isotope data from NIST and generates a Python data file
containing isotopic masses, abundances, and calculated values for use in
quantum chemistry calculations.
"""

import re
from pathlib import Path

# Most stable mass numbers for radioactive elements without natural abundance
most_stable_mass_numbers = {
    43: 98,  # Tc-98
    61: 145,  # Pm-145
    84: 209,  # Po-209
    85: 210,  # At-210
    86: 222,  # Rn-222
    87: 223,  # Fr-223
    88: 226,  # Ra-226
    89: 227,  # Ac-227
    93: 237,  # Np-237
    94: 244,  # Pu-244
    95: 243,  # Am-243
    96: 247,  # Cm-247
    97: 247,  # Bk-247
    98: 251,  # Cf-251
    99: 252,  # Es-252
    100: 257,  # Fm-257
    101: 258,  # Md-258
    102: 259,  # No-259
    103: 262,  # Lr-262
    104: 267,  # Rf-267
    105: 268,  # Db-268
    106: 271,  # Sg-271
    107: 270,  # Bh-270
    108: 269,  # Hs-269
    109: 278,  # Mt-278
    110: 281,  # Ds-281
    111: 281,  # Rg-281
    112: 285,  # Cn-285
    113: 286,  # Nh-286
    114: 289,  # Fl-289
    115: 289,  # Mc-289
    116: 293,  # Lv-293
    117: 293,  # Ts-293
    118: 294,  # Og-294
}


def parse_isotope_file(filename):
    """Parse NIST isotope data file and extract isotopic information.
    
    This function processes the NIST isotope data file and extracts atomic
    numbers, mass numbers, relative atomic masses, and isotopic compositions
    for each isotope.
    
    Args:
        filename (str): Path to the NIST isotope data file.
        
    Returns:
        dict: Dictionary containing isotope data organized by atomic number
              and mass number, with calculated most abundant isotope and
              weighted atomic mass for each element.
    """
    isotopes = {}
    with open(filename, "r") as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if line.startswith("Atomic Number ="):
            atomic_number = int(line.split("=")[1].strip())
        elif line.startswith("Mass Number ="):
            mass_number = int(line.split("=")[1].strip())
        elif line.startswith("Relative Atomic Mass ="):
            mass = float(re.search(r"([0-9.]+)", line).group(1))
        elif line.startswith("Isotopic Composition ="):
            match = re.search(r"([0-9.]+)", line)
            abundance = float(match.group(1)) if match else 0.0
            if atomic_number not in isotopes:
                isotopes[atomic_number] = {}
            isotopes[atomic_number][mass_number] = {
                "mass": mass,
                "abundance": abundance,
            }

    # Calculate most abundant isotope and weighted atomic mass
    for atomic_number in isotopes:
        isotope_data = isotopes[atomic_number]
        all_isotopes = {}
        for mass_number in isotope_data:
            if isinstance(mass_number, int):
                all_isotopes[mass_number] = isotope_data[mass_number]

        # Find most abundant isotope by natural abundance
        if all_isotopes:
            total_abundance = sum(
                d["abundance"] for d in all_isotopes.values()
            )
            if total_abundance > 0.0:
                most_abundant_mass_number = max(
                    all_isotopes, key=lambda x: all_isotopes[x]["abundance"]
                )
            else:
                most_abundant_mass_number = most_stable_mass_numbers[
                    atomic_number
                ]
            most_abundant = all_isotopes[most_abundant_mass_number]
            isotope_data["most_abundant"] = {
                "mass_number": most_abundant_mass_number,
                "mass": most_abundant["mass"],
                "abundance": most_abundant["abundance"],
            }

            # Calculate natural abundance weighted atomic mass
            if total_abundance > 0:
                weighted_mass = (
                    sum(
                        d["mass"] * d["abundance"]
                        for d in all_isotopes.values()
                    )
                    / total_abundance
                )
            else:
                weighted_mass = isotope_data.get(
                    most_stable_mass_numbers[atomic_number]
                )["mass"]
            isotope_data["weighted_atomic_mass"] = weighted_mass

    return isotopes


if __name__ == "__main__":
    """Main execution for generating isotope data file."""
    input_file = Path(__file__).parent / "../utils/isotopes.txt"
    output_file = Path(__file__).parent / "../utils/isotopes_data.py"

    isotopes = parse_isotope_file(input_file)

    # Write formatted isotope data to Python file
    with open(output_file, "w") as f:
        f.write(
            '''"""Isotope data extracted from NIST public website.

    Source data has been compiled by NIST:

        https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses

    The atomic weights data were published in:

        J. Meija et al, Atomic weights of the elements 2013,
        Pure and Applied Chemistry 88, 265-291 (2016).
        https://doi.org/10.1515/pac-2015-0305
        http://www.ciaaw.org/atomic-weights.htm

    Isotopic compositions data were published in:

        Michael Berglund and Michael E. Wieser,
        Isotopic compositions of the elements 2009 (IUPAC Technical Report)
        Pure Appl. Chem., 2011, Vol. 83, No. 2, pp. 397-410
        https://doi.org/10.1351/PAC-REP-10-06-02

    The relative atomic masses of the isotopes data were published in:

        M. Wang, G. Audi, A.H. Wapstra, F.G. Kondev, M. MacCormick, X. Xu,
        and B. Pfeiffer, The AME2012 Atomic Mass Evaluation,
        Chinese Phys. C 36 1603
        https://doi.org/10.1088/1674-1137/36/12/003
        http://amdc.impcas.ac.cn/evaluation/data2012/ame.html
    """

isotopes = {\n'''
        )

        for Z in sorted(isotopes.keys()):
            f.write(f"    {Z}: {{\n")
            for A in sorted(k for k in isotopes[Z] if isinstance(k, int)):
                data = isotopes[Z][A]
                if "mass" in data and "abundance" in data:
                    f.write(
                        f"        {A}: {{\"mass\": {data['mass']}, \"abundance\": {data['abundance']}}},\n"
                    )
            if "most_abundant" in isotopes[Z]:
                ma = isotopes[Z]["most_abundant"]
                f.write(
                    f"        \"most_abundant\": {{\"mass_number\": {ma['mass_number']}, \"mass\": {ma['mass']}, \"abundance\": {ma['abundance']}}},\n"
                )
            if "weighted_atomic_mass" in isotopes[Z]:
                f.write(
                    f"        \"weighted_atomic_mass\": {isotopes[Z]['weighted_atomic_mass']},\n"
                )
            f.write("    },\n")
        f.write("}\n")
