from rdkit import Chem
from rdkit.Chem.Descriptors3D import CalcMolDescriptors3D
import pandas as pd
import os
from descriptors_2d import xyz_to_sdf, find_steric_descriptors
#import sys #for command line running of the script

#Calculate all 3D descriptors
def calculate_3d_descriptors(xyz_file, num_atoms_per_conformer):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    sdf_file = xyz_to_sdf(xyz_file, "3d", num_atoms_per_conformer)
    file_name = os.path.splitext(os.path.basename(xyz_file))[0]
    output_csv = f"{file_name}_3d.csv"
    supplier = Chem.SDMolSupplier(sdf_file)
    data = []

    for i, mol in enumerate(supplier):
        if mol is None:
            print(f"Failed to parse molecule {i + 1} in the SDF file.")
            continue

        descriptors = CalcMolDescriptors3D(mol)
        descriptors['Name'] = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Molecule_{i + 1}"
        data.append(descriptors)

    # write to csv
    df = pd.DataFrame(data)
    columns_order = ['Name'] + [col for col in df.columns if col != 'Name']
    df = df[columns_order]
    df.to_csv(output_csv, index=False)
    os.remove(sdf_file)
    print(f".sdf file removed.", f"Descriptors saved to {output_csv}")
    return output_csv


# Example
xyz_file = input("Please enter the path to the .xyz file: ").strip()
num_atoms_per_conformer = 71  # Replace with the number of atoms in your molecule
# Calculate descriptors and save to CSV
#calculate_3d_descriptors(xyz_file, num_atoms_per_conformer)
find_steric_descriptors(calculate_3d_descriptors(xyz_file, num_atoms_per_conformer), "steric_descriptors_3d.csv")
