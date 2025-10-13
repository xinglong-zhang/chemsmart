eV_pattern = r"([\d\.]+) eV"
nm_pattern = r"([\d\.]+) nm"
f_pattern = r"f=([\d\.]+)"
float_pattern = r"[-]?\d*\.\d+|\d+"
integer_pattern = r"^\+?\d+$"  # Non-negative integer (incl. 0)
                               # optional leading + allowed
float_pattern_with_exponential = r"^[+-]?(?:(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?|\d+(?:[eE][+-]?\d+))$"
raw_energy_value_pattern = r"(-\d+\.\d+)"

xyz_filename_pattern = r"([^\s\"']+\.xyz\b)"
# \b ensures that the match ends right after xyz
# and is not followed by something like: xyz1, xyzabc xyz_thing
# It will match if .xyz is followed by: a space, a quote, end of line, punctuation

normal_mode_pattern = r"\s*(\d+)\s+(\d+)((?:\s+[+-]?\d*\.\d+)+)\s*"
frozen_coordinates_pattern = (
    r"\s*([A-Z][a-z]?)\s+(-1|0)\s+(-?\d+\.\d*)\s+(-?\d+\.\d*)\s+(-?\d+\.\d*)"
)
scf_energy_pattern = r"SCF Done:\s+E\([^)]*\)\s*=\s*([-.\d]+)"
mp2_energy_pattern = r"EUMP2\s*=\s*(.*)"
oniom_energy_pattern = r"ONIOM:\s+extrapolated energy\s*=\s*(.*)"
xyz_energy_pattern = r"Energy\(Hartree\):\s*(-?\d+\.\d+)"

# standard coordinate pattern with (symbol x y z)
standard_coord_pattern = (
    r"\s*[A-Z][a-z]?\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*"
)

"""
An example of the relevant part of the output describing the structure is:
        | 20> * xyz 0 1
        | 21>   O   -0.00000000323406      0.00000000000000      0.08734060152197
        | 22>   H   -0.75520523910536      0.00000000000000     -0.50967029975151
"""
orca_input_coordinate_in_output = r"(?i)\|\s+(\d+)>\s+(\w+)\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)"

"""Given input example:
        Mayer bond orders larger than 0.100000
        B( 25-C , 27-C ) :   1.4263 B( 25-C , 29-H ) :   1.0007 B( 27-C , 58-C ) :   1.0067
        B(  0-O ,  1-H ) :   0.9959 B(  0-O ,  2-H ) :   0.9959.
"""

orca_nproc_used_line_pattern = (
    r"Program running with (\d+) parallel MPI-processes"
)

mayer_bond_order_segment_pattern = (
    r"B\(\s*(\d+)-([A-Z])\s*,\s*(\d+)-([A-Z])\s*\)\s*:\s*(\d+\.\d+)"
)

sdf_pattern = (
    r"\s*([\d\.-]+)\s+([\d\.-]+)\s+([\d\.-]+)\s+(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"
    r"\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*"
)

# regex pattern for constrained coordinates in orca input
# Examples of ORCA constraints patterns:
# Constraining bond distances : { B N1 N2 value C }
# Constraining bond angles : { A N1 N2 N1 value C }
# Constraining dihedral angles : { D N1 N2 N3 N4 value C }
# Constraining cartesian coordinates : { C N1 C }
# where the value is optional, by default it is the present value in the structure
constrained_bond_length_pattern_in_input = (
    r"\|\s*(\d+)>.*\{\s*B\s+(\d+)\s+(\d+)(?:\s+([\d.]+))?\s+C\s*\}"
)

constrained_bond_angles_pattern_in_input = (
    r"\|\s*(\d+)>.*\{\s*A\s+(\d+)\s+(\d+)\s+(\d+)(?:\s+([\d.]+))?\s+C\s*\}"
)

constrained_dihedrals_pattern_in_input = r"\|\s*(\d+)>.*\{\s*D\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)(?:\s+([\d.]+))?\s+C\s*\}"

# in orca output:
# Will constrain atom 2 coordinate 1
# Will constrain atom 2 coordinate 2
# Will constrain atom 2 coordinate 3
# etc.
orca_frozen_atoms_output_pattern = r"Will constrain atom \d+ coordinate \d"

# regex pattern for constrained coordinates in orca
# \s+[BAD] : B, A, or D with whitespace.
# example lines for match:
# 10. B(H   9,H   8)                  2.4714         0.000431 C
# 30. A(C   1,C   5,H   8)           69.0631         0.105398 C
# 49. D(H   9,C   4,C   3,C   2)   -180.0000         0.024745 C
orca_constrained_coordinates_pattern = r"^\d+\.\s+[BAD]\([A-Z][a-z]?\s+\d+,[A-Z][a-z]?\s+\d+(?:,[A-Z][a-z]?\s+\d+)?(?:,[A-Z][a-z]?\s+\d+)?\)\s+-?\d+\.\d{4}\s+\d+\.\d{6}\s+C$"

# filename pattern for orca output files

# filename matches with point pxx but not with fragment fx
orca_dias_filename_point_without_fragment = r".*_p(\d+)(?:_(?!f)(.+))?\.out"

# filename matches with point pxx and fragment f1
orca_dias_filename_point_with_fragment1 = r".*_p(\d+)_(f1)(?:_(.+)?)?\.out"

# filename matches with point pxx and fragment f2
orca_dias_filename_point_with_fragment2 = r".*_p(\d+)_(f2)(?:_(.+)?)?\.out"

# filename matches with reactant r1 or r2
orca_dias_filename_with_reactant = r".*_r([12])(?:_(.+)?)?\.out"


# filename pattern for gaussian output files

# filename matches with point pxx but not with fragment fx
# matches from the word "dias" onwards
gaussian_dias_filename_point_without_fragment = r"(?:.*dias_p(\d+)(?:_((?:(?!f\d).)+))?\.log)|(?:.*_p(\d+)(?:_((?:(?!f\d).)+))?\.log)"

# filename matches with point pxx and fragment f1
gaussian_dias_filename_point_with_fragment1 = r".*_p(\d+)_(f1)(?:_(.+)?)?\.log"

# filename matches with point pxx and fragment f2
gaussian_dias_filename_point_with_fragment2 = r".*_p(\d+)_(f2)(?:_(.+)?)?\.log"

# filename matches with reactant r1 or r2
gaussian_dias_filename_with_reactant = r".*_r([12])(?:_(.+)?)?\.log"


# Route string cleaning patterns for Gaussian link job settings

# Pattern to find optimization keywords from route string
# Matches: "opt", "opt=word", "opt=(parameters)", "opt = word", etc.
gaussian_opt_keywords_pattern = r"\bopt\s*(=\s*(\([^)]*\)|\w+))?\s*"

# Pattern to find frequency calculation keywords from route string
# Matches: "freq", "freq=numer", "freq = analytical", etc.
# The \b ensures that "freq" is matched as a whole word
# and not part of another word like "frequency"
# Although note that "frequency" is a valid Gaussian keyword,
# we'd want to avoid erroneous partial matches, eg.,
# freqency (spelling error eg)
gaussian_freq_keywords_pattern = r"\bfreq\b\s*(=\s*\w+)?\s*"

# Pattern to find multiple consecutive spaces in strings
multiple_spaces_pattern = r"\s+"


# PyMOL strings
pymol_isosurface_pattern = r"isosurface\s*=\s*[\d\.]+"
pymol_color_range_pattern = r"range\s*=\s*[\d\.]+"
