#!/usr/bin/env python3
"""
Streamlined Dihedral Featurization Tool

A simplified tool for computing dihedral angle features from MD trajectories.
"""

import numpy as np
import mdtraj as md
import os
import glob
import pandas as pd
import pickle
import argparse


class DihedralFeaturizer:
    """
    Compute dihedral angle features from MD trajectories.
    
    Parameters
    ----------
    types : list, default=['phi', 'psi']
        Dihedral angles to compute: 'phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4'
    sincos : bool, default=True
        Return sine and cosine of angles instead of raw angles
    """
    
    def __init__(self, types=['phi', 'psi'], sincos=True):
        if isinstance(types, str):
            types = [types]
        self.types = list(types)
        self.sincos = sincos
        
        known = {'phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4'}
        if not set(types).issubset(known):
            raise ValueError(f'angles must be subset of {known}, got {types}')

    def transform(self, traj):
        """
        Featurize trajectory into dihedral angle features.
        
        Parameters
        ----------
        traj : mdtraj.Trajectory
            Input trajectory
            
        Returns
        -------
        features : np.ndarray, shape=(n_frames, n_features)
            Dihedral angle features
        """
        features = []
        
        for angle_type in self.types:
            func = getattr(md, f'compute_{angle_type}')
            _, angles = func(traj)
            
            if self.sincos:
                features.extend([np.sin(angles), np.cos(angles)])
            else:
                features.append(angles)
        
        return np.hstack(features)

    def describe_features(self, traj):
        """
        Return a list of dictionaries describing the dihedral features.
        
        Returns
        -------
        feature_descs : list of dict
            Dictionary describing each feature with the following information
            about the atoms participating in each dihedral:
                - resnames: unique names of residues
                - atominds: the four atom indices
                - resseqs: unique residue sequence ids (not necessarily 0-indexed)
                - resids: unique residue ids (0-indexed)
                - featurizer: Dihedral
                - featuregroup: the type of dihedral angle
                - otherinfo: sin/cos/raw transformation applied
        """
        feature_descs = []
        top = traj.topology
        
        for dihed_type in self.types:
            # Get the compute function for this dihedral type
            func = getattr(md, f'compute_{dihed_type}')
            # Get atom indices participating in each dihedral
            aind_tuples, _ = func(traj)
            
            # Extract residue information for each dihedral
            for ainds in aind_tuples:
                resid = set(top.atom(ai).residue.index for ai in ainds)
                resids = list(resid)
                reseq = set(top.atom(ai).residue.resSeq for ai in ainds)
                resseqs = list(reseq)
                resname = set(top.atom(ai).residue.name for ai in ainds)
                resnames = list(resname)
                
                base_info = {
                    'resnames': resnames,
                    'atominds': ainds,
                    'resseqs': resseqs,
                    'resids': resids,
                    'featurizer': 'Dihedral',
                    'featuregroup': dihed_type
                }
                
                if self.sincos:
                    # Add sin feature
                    sin_info = base_info.copy()
                    sin_info['otherinfo'] = 'sin'
                    feature_descs.append(sin_info)
                    
                    # Add cos feature
                    cos_info = base_info.copy()
                    cos_info['otherinfo'] = 'cos'
                    feature_descs.append(cos_info)
                else:
                    base_info['otherinfo'] = 'raw'
                    feature_descs.append(base_info)
        
        return feature_descs


def featurize_pdb(pdb_path, angles=['phi', 'psi'], sincos=True, output_dir=None):
    """
    Featurize a single PDB file.
    
    Parameters
    ----------
    pdb_path : str
        Path to PDB file
    angles : list, default=['phi', 'psi']
        Dihedral angles to compute
    sincos : bool, default=True
        Use sin/cos transformation
    output_dir : str, optional
        Output directory (default: based on PDB filename)
    
    Returns
    -------
    features : np.ndarray
        Computed features
    descriptions : list
        Feature descriptions
    """
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
    
    print(f"Loading: {os.path.basename(pdb_path)}")
    traj = md.load(pdb_path)
    print(f"Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")
    
    # Featurize
    featurizer = DihedralFeaturizer(types=angles, sincos=sincos)
    features = featurizer.transform(traj)
    descriptions = featurizer.describe_features(traj)
    
    print(f"Features shape: {features.shape}")
    
    # Save results
    if output_dir is None:
        base_name = os.path.splitext(os.path.basename(pdb_path))[0]
        output_dir = f"{base_name}_features"
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save as pickle and numpy
    np.save(os.path.join(output_dir, 'features.npy'), features)
    with open(os.path.join(output_dir, 'features.pkl'), 'wb') as f:
        pickle.dump(features, f)
    
    # Save descriptions as DataFrame pickle (matching original format)
    desc_df = pd.DataFrame(descriptions)
    with open(os.path.join(output_dir, 'feature_descriptions.pkl'), 'wb') as f:
        pickle.dump(desc_df, f)
    
    # Also save as CSV for readability
    desc_df.to_csv(os.path.join(output_dir, 'feature_descriptions.csv'), index=False)
    
    print(f"Results saved to: {output_dir}")
    print("Files created:")
    print(f"  - features.pkl (Python pickle format)")
    print(f"  - features.npy (NumPy array format)")
    print(f"  - feature_descriptions.pkl (DataFrame pickle)")
    print(f"  - feature_descriptions.csv (CSV format)")
    
    return features, descriptions


def featurize_directory(directory, angles=['phi', 'psi'], sincos=True, 
                       stride=1, output_dir=None):
    """
    Featurize all .xtc files in a directory using available .pdb files as topology.
    
    Parameters
    ----------
    directory : str
        Directory containing .xtc and .pdb files
    angles : list, default=['phi', 'psi'] 
        Dihedral angles to compute
    sincos : bool, default=True
        Use sin/cos transformation
    stride : int, default=1
        Load every stride-th frame
    output_dir : str, optional
        Output directory (default: based on input directory)
        
    Returns
    -------
    all_features : list
        List of feature arrays for each trajectory
    descriptions : list
        Feature descriptions
    """
    if not os.path.isdir(directory):
        raise ValueError(f"Directory not found: {directory}")
    
    # Find files
    xtc_files = glob.glob(os.path.join(directory, '*.xtc'))
    pdb_files = glob.glob(os.path.join(directory, '*.pdb'))
    
    if not xtc_files:
        raise ValueError(f"No .xtc files found in {directory}")
    if not pdb_files:
        raise ValueError(f"No .pdb files found in {directory}")
    
    print(f"Found {len(xtc_files)} .xtc and {len(pdb_files)} .pdb files")
    
    featurizer = DihedralFeaturizer(types=angles, sincos=sincos)
    all_features = []
    descriptions = None
    
    for xtc_file in xtc_files:
        print(f"Processing: {os.path.basename(xtc_file)}")
        
        # Find matching topology
        base_name = os.path.splitext(os.path.basename(xtc_file))[0]
        specific_pdb = os.path.join(directory, f"{base_name}.pdb")
        
        if os.path.exists(specific_pdb):
            pdb_file = specific_pdb
        else:
            # Try to find compatible PDB by atom count
            pdb_file = None
            try:
                temp_traj = md.load_frame(xtc_file, 0, top=pdb_files[0])
                xtc_atoms = temp_traj.n_atoms
                
                for candidate in pdb_files:
                    if md.load(candidate).n_atoms == xtc_atoms:
                        pdb_file = candidate
                        break
            except:
                pass
        
        if pdb_file is None:
            print(f"  No compatible topology found, skipping")
            continue
            
        try:
            # Load and featurize
            traj = md.load(xtc_file, top=pdb_file, stride=stride)
            features = featurizer.transform(traj)
            all_features.append(features)
            
            # Get descriptions from first successful trajectory
            if descriptions is None:
                descriptions = featurizer.describe_features(traj)
            
            print(f"  Loaded {traj.n_frames} frames → {features.shape}")
            
        except Exception as e:
            print(f"  Error: {e}")
            continue
    
    if not all_features:
        raise ValueError("No trajectories were successfully processed")
    
    # Save results
    if output_dir is None:
        output_dir = f"{os.path.basename(directory)}_features"
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save individual and concatenated features
    for i, features in enumerate(all_features):
        np.save(os.path.join(output_dir, f'traj_{i}_features.npy'), features)
    
    concatenated = np.vstack(all_features)
    np.save(os.path.join(output_dir, 'all_features.npy'), concatenated)
    
    with open(os.path.join(output_dir, 'all_features.pkl'), 'wb') as f:
        pickle.dump(all_features, f)
    
    # Save descriptions as DataFrame pickle (matching original format)  
    if descriptions:
        desc_df = pd.DataFrame(descriptions)
        with open(os.path.join(output_dir, 'feature_descriptions.pkl'), 'wb') as f:
            pickle.dump(desc_df, f)
        
        # Also save as CSV for readability
        desc_df.to_csv(os.path.join(output_dir, 'feature_descriptions.csv'), index=False)
    
    print(f"Processed {len(all_features)} trajectories")
    print(f"Total frames: {concatenated.shape[0]}, Features: {concatenated.shape[1]}")
    print(f"Results saved to: {output_dir}")
    print("Files created:")
    print(f"  - all_features.pkl (list of trajectory features)")
    print(f"  - all_features.npy (concatenated features)")
    print(f"  - feature_descriptions.pkl (DataFrame pickle)")
    print(f"  - feature_descriptions.csv (CSV format)")
    print(f"  - traj_*_features.npy (individual trajectory features)")
    
    return all_features, descriptions


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description="Compute dihedral angle features from MD trajectories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s protein.pdb                           # Featurize single PDB
  %(prog)s --directory ./trajectories            # Featurize directory
  %(prog)s protein.pdb --angles phi psi chi1     # Custom angles
  %(prog)s protein.pdb --sincos false            # Raw angles (no sin/cos)
        """
    )
    
    parser.add_argument('pdb_file', nargs='?', help='Single PDB file to process')
    parser.add_argument('--directory', '-d', help='Directory with .xtc/.pdb files')
    parser.add_argument('--angles', '-a', nargs='+', default=['phi', 'psi'],
                       choices=['phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4'],
                       help='Dihedral angles to compute (default: phi psi)')
    parser.add_argument('--sincos', type=lambda x: x.lower() == 'true', default=True,
                       help='Use sin/cos transformation (default: true)')
    parser.add_argument('--stride', '-s', type=int, default=1,
                       help='Load every Nth frame (default: 1)')
    parser.add_argument('--output', '-o', help='Output directory')
    
    args = parser.parse_args()
    
    # Validate input
    if not args.pdb_file and not args.directory:
        parser.error("Must specify either a PDB file or --directory")
    if args.pdb_file and args.directory:
        parser.error("Cannot specify both PDB file and --directory")
    
    try:
        if args.pdb_file:
            featurize_pdb(args.pdb_file, args.angles, args.sincos, args.output)
        else:
            featurize_directory(args.directory, args.angles, args.sincos, 
                              args.stride, args.output)
                              
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
