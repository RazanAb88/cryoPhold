import mdtraj as md
import os
import glob
import argparse

def combine_pdb_to_xtc(input_pdbs, output_xtc, remove_hydrogens=False):
    """
    Combines multiple PDB files into a single XTC trajectory file using MDTraj.
    Optionally removes hydrogen atoms before combining.

    Parameters:
    input_pdbs (list): List of input PDB filenames.
    output_xtc (str): Output XTC filename.
    remove_hydrogens (bool): Whether to remove hydrogens before combining.
    """
    # Load all trajectories
    trajectories = []
    for pdb in input_pdbs:
        traj = md.load(pdb)
        print(f"Loaded {pdb}")

        # Remove hydrogen atoms if requested
        if remove_hydrogens:
            non_hydrogen_indices = traj.topology.select("not symbol H")
            traj = traj.atom_slice(non_hydrogen_indices)
            print(f"  Removed hydrogen atoms from {pdb}")

        trajectories.append(traj)
    
    # Combine trajectories
    combined_traj = md.join(trajectories)
    
    # Save as XTC
    combined_traj.save_xtc(output_xtc)

def main():
    parser = argparse.ArgumentParser(
        description="Combine all .pdb files in a specified folder into a single XTC trajectory."
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Path to the input folder containing PDB files."
    )
    parser.add_argument(
        "-o", "--output", type=str, default="combined.xtc",
        help="Output XTC filename (default: combined.xtc)."
    )
    parser.add_argument(
        "-r", "--remove-hydrogens", action="store_true",
        help="Remove hydrogen atoms from the PDB files before combining."
    )
    
    args = parser.parse_args()
    
    # Find all .pdb files in the provided folder and sort them.
    input_files = sorted(glob.glob(os.path.join(args.input, "*.pdb")))
    
    if not input_files:
        print(f"No .pdb files found in the folder: {args.input}")
        return
    
    print(f"Found {len(input_files)} PDB files. Combining them into {args.output}")
    
    try:
        combine_pdb_to_xtc(
            input_pdbs=input_files,
            output_xtc=args.output,
            remove_hydrogens=args.remove_hydrogens
        )
        print(f"\nSuccessfully created {args.output}")
    except Exception as e:
        print(f"Error occurred: {str(e)}")

if __name__ == "__main__":
    main()
