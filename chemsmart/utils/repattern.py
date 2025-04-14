eV_pattern = r"([\d\.]+) eV"
nm_pattern = r"([\d\.]+) nm"
f_pattern = r"f=([\d\.]+)"
float_pattern = r"[-]?\d*\.\d+|\d+"
normal_mode_pattern = r"\s*(\d+)\s+(\d+)((?:\s+[+-]?\d*\.\d+)+)\s*"
frozen_coordinates_pattern = (
    r"\s*([A-Z][a-z]?)\s+(-1|0)\s+(-?\d+\.\d*)\s+(-?\d+\.\d*)\s+(-?\d+\.\d*)"
)
scf_energy_pattern = r"SCF Done:\s+E\([^)]*\)\s*=\s*([-.\d]+)"
mp2_energy_pattern = r"EUMP2\s*=\s*(.*)"
oniom_energy_pattern = r"ONIOM:\s+extrapolated energy\s*=\s*(.*)"

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
# Examples of ORCA constraints patterns:
# Constraining bond distances : { B N1 N2 value C }
# Constraining bond angles : { A N1 N2 N1 value C }
# Constraining dihedral angles : { D N1 N2 N3 N4 value C }
# Constraining cartesian coordinates : { C N1 C }
# where the value is optional, by default it is the present value in the structure
constrained_bond_length_pattern = (
    r"\|\s*(\d+)>.*\{\s*B\s+(\d+)\s+(\d+)(?:\s+([\d.]+))?\s+C\s*\}"
)

constrained_bond_angles_pattern = (
    r"\|\s*(\d+)>.*\{\s*A\s+(\d+)\s+(\d+)\s+(\d+)(?:\s+([\d.]+))?\s+C\s*\}"
)

constrained_dihedrals_pattern = r"\|\s*(\d+)>.*\{\s*D\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)(?:\s+([\d.]+))?\s+C\s*\}"

orca_frozen_atoms_output_pattern = r"Will constrain atom \d+ coordinate \d"
