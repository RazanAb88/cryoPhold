#!/usr/bin/env python3
import mdtraj as md
import argparse

def remove_hydrogens(input_pdb, output_pdb):
    # Load the PDB file
    traj = md.load(input_pdb)
    
    # Select non-hydrogen atoms
    non_hydrogen_indices = traj.top.select('not element H')

    # Create a new trajectory without hydrogens
    traj_no_h = traj.atom_slice(non_hydrogen_indices)

    # Save the new trajectory as PDB
    traj_no_h.save_pdb(output_pdb)
    print(f"Saved hydrogen-free PDB as: {output_pdb}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove hydrogens from a PDB file.")
    parser.add_argument("-i", "--input", required=True, help="Input PDB file")
    parser.add_argument("-o", "--output", required=True, help="Output PDB file without hydrogens")
    args = parser.parse_args()

    remove_hydrogens(args.input, args.output)
