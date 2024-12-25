from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
from rdkit import RDLogger
import os
import csv
#import sys #for command line running of the script

#generate .sdf file to store steric information
def xyz_to_sdf(xyz_file, memo, num_atoms_per_conformer):
    file_name = os.path.splitext(os.path.basename(xyz_file))[0]
    sdf_file = f"{file_name}_{memo}.sdf"
    working_dir = os.path.dirname(os.path.abspath(xyz_file))
    os.chdir(working_dir)  # switch to the directory where the geometry is located
    xyz_file = f"{file_name}.xyz"

    with open(xyz_file, "r") as f:
        xyz_data = f.read()

    lines = xyz_data.strip().splitlines()
    num_conformers = len(lines) // (num_atoms_per_conformer + 2) #calculate number of conformers
    writer = Chem.SDWriter(sdf_file)
    RDLogger.DisableLog('rdApp.*') #mute rdKit warning messages

    for i in range(num_conformers):  #extract coordinates for each conformer, store in .sdf format
        start_idx = i * (num_atoms_per_conformer + 2) + 2
        end_idx = start_idx + num_atoms_per_conformer
        atoms = lines[start_idx:end_idx]
        mol = Chem.RWMol()
        conf_coords = []

        for atom_data in atoms:
            atom_info = atom_data.split()
            atom = Chem.Atom(atom_info[0])
            mol.AddAtom(atom)
            conf_coords.append((float(atom_info[1]), float(atom_info[2]), float(atom_info[3])))  # x, y, z

        conf = Chem.Conformer(len(conf_coords))
        for idx, (x, y, z) in enumerate(conf_coords):
            conf.SetAtomPosition(idx, Chem.rdGeometry.Point3D(x, y, z))
        mol.AddConformer(conf)

        try:
            Chem.SanitizeMol(mol, Chem.SANITIZE_ALL)
            mol.SetProp("_Name", f"Conformer_{i + 1}")
            writer.write(mol)
        except Exception as e:
            print(f"Failed to sanitize or write conformer {i + 1}: {e}")

    writer.close()
    print(f"generated temporary .sdf file from given geometry")
    return sdf_file


# Calculate 210 descriptors for conformers
def calculate_2d_descriptors(xyz_file, num_atoms_per_conformer):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    sdf_file=xyz_to_sdf(xyz_file, "2d", num_atoms_per_conformer)
    file_name = os.path.splitext(os.path.basename(xyz_file))[0]
    output_csv = f"{file_name}_2d.csv"
    supplier = Chem.SDMolSupplier(sdf_file)
    data = []

    for i, mol in enumerate(supplier):
        if mol is None:
            print(f"Failed to parse molecule {i + 1} in the SDF file.")
            continue

        descriptors = {desc[0]: desc[1](mol) for desc in Descriptors.descList}
        descriptors['Name'] = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Molecule_{i + 1}"
        data.append(descriptors)

    # write to csv
    df = pd.DataFrame(data)
    columns_order = ['Name'] + [col for col in df.columns if col != 'Name']
    df = df[columns_order]
    df.to_csv(output_csv, index=False)
    os.remove(sdf_file)
    print(f".sdf file removed.",f"Descriptors saved to {output_csv}")
    return output_csv

def create_file(file_name):
    # check if the file already exists
    if os.path.exists(file_name):
        # ask user if they want to overwrite
        decision = input(f"The file '{file_name}' already exists. Do you want to overwrite it? (yes/no): ").strip().lower()
        if decision not in ["yes", "y"]:
            print("Operation canceled. File was not overwritten.")
            return

    with open(file_name, "w") as file:    # create or overwrite the file
        file.write("")
    print(f"The file '{file_name}' has been created or overwritten.")

#find descriptors that differ between conformers
def find_steric_descriptors(csv_file, output_file):
    with open(csv_file, 'r') as file:
        reader = list(csv.reader(file))

    if len(reader) < 2:
        print("Please enter the correct file.")
    else:
        first_row = reader[0]
        second_row = reader[-2]
        last_row = reader[-1]
        different_columns = []
        # compare column values in the second-to-last and last rows
        for i, (val1, val2) in enumerate(zip(second_row, last_row)):
            if val1 != val2 and val1 is not None:
                different_columns.append(first_row[i])

        if different_columns and len(different_columns) !=1:
            create_file(output_file)
            for item in different_columns:

                if item != "Name":
                    with open(output_file, 'a') as out_file:
                        out_file.write(f"{item}\n")
            print(f"Values written to {output_file}")
        else:
            print("The descriptors values of conformers are the same")

# #Example
# num_atoms_per_conformer = 71
# # xyz_file = sys.argv[1].strip() #for command line running
# xyz_file = input("Please enter the path to the .xyz file: ").strip()
# # Calculate descriptors and save to CSV
# calculate_2d_descriptors(xyz_file, num_atoms_per_conformer)
## find_steric_descriptors(calculate_2d_descriptors(xyz_file, num_atoms_per_conformer), "steric_descriptors_2d.csv")
