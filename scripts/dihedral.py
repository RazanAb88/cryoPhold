#!/usr/bin/env python3
"""
Trajectory Loading and Dihedral Featurization - Standalone Implementation

This script demonstrates how to load a PDB file and perform dihedral angle featurization.
"""

import numpy as np
import mdtraj as md
import os
import glob
import pandas as pd
import pickle
import itertools


class DihedralFeaturizer:
    """
    Standalone featurizer based on dihedral angles.

    This featurizer transforms MD trajectories into vector datasets by
    representing each frame with dihedral angles (backbone or side-chain).

    Parameters
    ----------
    types : list
        One or more of ['phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4']
    sincos : bool
        Instead of outputting the angle, return the sine and cosine of the
        angle as separate features.
    """

    def __init__(self, types=['phi', 'psi'], sincos=True):
        if isinstance(types, str):
            types = [types]
        self.types = list(types)  # force a copy
        self.sincos = sincos

        known = {'phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4'}
        if not set(types).issubset(known):
            raise ValueError('angles must be a subset of %s. you supplied %s' %
                             (str(known), str(types)))

    def partial_transform(self, traj):
        """
        Featurize an MD trajectory into a vector space via calculation
        of dihedral (torsion) angles

        Parameters
        ----------
        traj : mdtraj.Trajectory
            A molecular dynamics trajectory to featurize.

        Returns
        -------
        features : np.ndarray, dtype=float, shape=(n_samples, n_features)
            A featurized trajectory is a 2D array of shape
            `(length_of_trajectory x n_features)` where each `features[i]`
            vector is computed by applying the featurization function
            to the `i`th snapshot of the input trajectory.
        """
        x = []
        for a in self.types:
            func = getattr(md, 'compute_%s' % a)
            _, y = func(traj)

            if self.sincos:
                x.extend([np.sin(y), np.cos(y)])
            else:
                x.append(y)

        return np.hstack(x)

    def transform(self, traj_list):
        """
        Featurize several trajectories.

        Parameters
        ----------
        traj_list : list(mdtraj.Trajectory)
            Trajectories to be featurized.

        Returns
        -------
        features : list(np.ndarray), length = len(traj_list)
            The featurized trajectories.  features[i] is the featurized
            version of traj_list[i] and has shape (n_samples_i, n_features)
        """
        return [self.partial_transform(traj) for traj in traj_list]

    def describe_features(self, traj):
        """
        Return a list of dictionaries describing the dihedral features.

        Parameters
        ----------
        traj : mdtraj.Trajectory
            The trajectory to describe

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
                - featuregroup: the type of dihedral angle and whether sin or cos has been applied.
        """
        feature_descs = []
        
        for dihed_type in self.types:
            # Get the compute function for this dihedral type
            func = getattr(md, 'compute_%s' % dihed_type)
            # Get atom indices participating in each dihedral
            aind_tuples, _ = func(traj)
            
            top = traj.topology
            
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
                    sin_info['featuregroup'] = dihed_type
                    sin_info['otherinfo'] = 'sin'
                    feature_descs.append(sin_info)
                    
                    # Add cos feature
                    cos_info = base_info.copy()
                    cos_info['featuregroup'] = dihed_type
                    cos_info['otherinfo'] = 'cos'
                    feature_descs.append(cos_info)
                else:
                    base_info['otherinfo'] = 'raw'
                    feature_descs.append(base_info)

        return feature_descs


class TrajectoryDihedralAnalyzer:
    """
    A class to handle trajectory loading and dihedral featurization
    """
    
    def __init__(self, topology_file=None):
        """
        Initialize the analyzer
        
        Parameters
        ----------
        topology_file : str, optional
            Path to topology file (PDB, PSF, etc.)
        """
        self.topology_file = topology_file
        self.trajectories = []
        self.featurizer = None
        self.features = None
        self.feature_descriptions = None
        
    def load_trajectories(self, trajectory_files, topology=None, stride=1, chunk=1000):
        """
        Load multiple trajectory files
        
        Parameters
        ----------
        trajectory_files : list of str
            List of trajectory file paths
        topology : str, optional
            Topology file path (overrides class topology)
        stride : int, default=1
            Load every stride-th frame
        chunk : int, default=1000
            Chunk size for memory-efficient loading
        
        Returns
        -------
        trajectories : list of mdtraj.Trajectory
            Loaded trajectories
        """
        if topology is None:
            topology = self.topology_file
            
        self.trajectories = []
        
        print(f"Loading {len(trajectory_files)} trajectory files...")
        
        for i, traj_file in enumerate(trajectory_files):
            print(f"  Loading trajectory {i+1}/{len(trajectory_files)}: {os.path.basename(traj_file)}")
            
            try:
                if traj_file.endswith('.h5'):
                    # HDF5 files already contain topology
                    traj = md.load(traj_file, stride=stride)
                else:
                    # Other formats need topology
                    if topology is None:
                        raise ValueError(f"Topology file required for {traj_file}")
                    traj = md.load(traj_file, top=topology, stride=stride)
                
                self.trajectories.append(traj)
                print(f"    Loaded {traj.n_frames} frames, {traj.n_atoms} atoms")
                
            except Exception as e:
                print(f"    Error loading {traj_file}: {e}")
                continue
        
        print(f"Successfully loaded {len(self.trajectories)} trajectories")
        return self.trajectories
    
    def setup_dihedral_featurizer(self, types=['phi', 'psi'], sincos=True):
        """
        Setup dihedral featurizer
        
        Parameters
        ----------
        types : list of str, default=['phi', 'psi']
            Types of dihedral angles to compute
            Options: 'phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4'
        sincos : bool, default=True
            Whether to compute sin/cos of angles instead of raw angles
        """
        self.featurizer = DihedralFeaturizer(types=types, sincos=sincos)
        print(f"Dihedral featurizer setup with types: {types}, sincos: {sincos}")
        
    def featurize_trajectories(self):
        """
        Apply dihedral featurization to loaded trajectories
        
        Returns
        -------
        features : list of np.ndarray
            Featurized trajectories
        """
        if self.featurizer is None:
            raise ValueError("Featurizer not setup. Call setup_dihedral_featurizer() first.")
        
        if not self.trajectories:
            raise ValueError("No trajectories loaded. Call load_trajectories() first.")
        
        print("Featurizing trajectories...")
        self.features = []
        
        for i, traj in enumerate(self.trajectories):
            print(f"  Featurizing trajectory {i+1}/{len(self.trajectories)}")
            features = self.featurizer.partial_transform(traj)
            self.features.append(features)
            print(f"    Shape: {features.shape}")
        
        # Get feature descriptions from first trajectory
        if self.trajectories:
            self.feature_descriptions = self.featurizer.describe_features(self.trajectories[0])
            print(f"Total features: {len(self.feature_descriptions)}")
        
        return self.features
    
    def get_feature_info(self):
        """
        Get information about the computed features
        
        Returns
        -------
        df : pd.DataFrame
            DataFrame with feature information
        """
        if self.feature_descriptions is None:
            raise ValueError("Features not computed yet. Call featurize_trajectories() first.")
        
        df = pd.DataFrame(self.feature_descriptions)
        return df
    
    def plot_time_series(self, traj_idx=0, feature_indices=None, max_features=6, figsize=(15, 8)):
        """
        Plot time series of dihedral features
        
        Parameters
        ----------
        traj_idx : int, default=0
            Index of trajectory to plot
        feature_indices : list, optional
            Specific feature indices to plot
        max_features : int, default=6
            Maximum number of features to plot if feature_indices not given
        figsize : tuple, default=(15, 8)
            Figure size
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Matplotlib not available for plotting")
            return
            
        if self.features is None:
            raise ValueError("Features not computed yet. Call featurize_trajectories() first.")
        
        features = self.features[traj_idx]
        
        if feature_indices is None:
            feature_indices = list(range(min(features.shape[1], max_features)))
        
        fig, axes = plt.subplots(len(feature_indices), 1, figsize=figsize, sharex=True)
        if len(feature_indices) == 1:
            axes = [axes]
        
        time = np.arange(features.shape[0])
        
        for i, feat_idx in enumerate(feature_indices):
            axes[i].plot(time, features[:, feat_idx])
            axes[i].set_ylabel(f'Feature {feat_idx}')
            axes[i].grid(True, alpha=0.3)
        
        axes[-1].set_xlabel('Frame')
        plt.suptitle(f'Dihedral Feature Time Series (Trajectory {traj_idx})')
        plt.tight_layout()
        plt.show()
    
    def save_features(self, output_dir='featurized_data'):
        """
        Save featurized data to pickle files
        
        Parameters
        ----------
        output_dir : str, default='featurized_data'
            Directory to save the data
        """
        if self.features is None:
            raise ValueError("Features not computed yet. Call featurize_trajectories() first.")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Save individual trajectory features as pickle files
        for i, features in enumerate(self.features):
            filename = os.path.join(output_dir, f'features_traj_{i}.pkl')
            with open(filename, 'wb') as f:
                pickle.dump(features, f)
            print(f"Saved trajectory {i} features to {filename}")
        
        # Save concatenated features as pickle
        all_features = np.vstack(self.features)
        concat_filename = os.path.join(output_dir, 'features_all.pkl')
        with open(concat_filename, 'wb') as f:
            pickle.dump(all_features, f)
        print(f"Saved concatenated features to {concat_filename}")
        
        # Save feature descriptions as DataFrame pickle
        if self.feature_descriptions:
            desc_df = pd.DataFrame(self.feature_descriptions)
            desc_filename = os.path.join(output_dir, 'feature_descriptions.pkl')
            with open(desc_filename, 'wb') as f:
                pickle.dump(desc_df, f)
            print(f"Saved feature descriptions as DataFrame to {desc_filename}")
            
            # Also save as CSV for readability
            csv_filename = os.path.join(output_dir, 'feature_descriptions.csv')
            desc_df.to_csv(csv_filename, index=False)
            print(f"Saved feature descriptions to {csv_filename}")
        
        # Save complete analyzer state
        analyzer_data = {
            'features': self.features,
            'feature_descriptions': self.feature_descriptions,
            'featurizer_types': self.featurizer.types if self.featurizer else None,
            'featurizer_sincos': self.featurizer.sincos if self.featurizer else None,
            'n_trajectories': len(self.trajectories),
            'trajectory_shapes': [traj.xyz.shape for traj in self.trajectories]
        }
        
        state_filename = os.path.join(output_dir, 'analyzer_state.pkl')
        with open(state_filename, 'wb') as f:
            pickle.dump(analyzer_data, f)
        print(f"Saved complete analyzer state to {state_filename}")
        
        print(f"All data saved to {output_dir}/")
    
    def load_features(self, filepath):
        """
        Load features from a pickle file
        
        Parameters
        ----------
        filepath : str
            Path to the pickle file containing features
            
        Returns
        -------
        features : np.ndarray
            Loaded features
        """
        with open(filepath, 'rb') as f:
            features = pickle.load(f)
        print(f"Loaded features from {filepath}, shape: {features.shape}")
        return features


def featurize_all(filenames, featurizer, topology, chunk=1000, stride=1):
    """
    Load and featurize many trajectory files.

    Parameters
    ----------
    filenames : list of strings
        List of paths to MD trajectory files
    featurizer : DihedralFeaturizer
        The featurizer to be invoked on each trajectory trajectory as
        it is loaded
    topology : str, Topology, Trajectory
        Topology or path to a topology file, used to load trajectories with
        MDTraj
    chunk : {int, None}
        If chunk is an int, load the trajectories up in chunks using
        md.iterload for better memory efficiency (less trajectory data needs
        to be in memory at once)
    stride : int, default=1
        Only read every stride-th frame.

    Returns
    -------
    data : np.ndarray, shape=(total_length_of_all_trajectories, n_features)
    indices : np.ndarray, shape=(total_length_of_all_trajectories)
    fns : np.ndarray shape=(total_length_of_all_trajectories)
        These three arrays all share the same indexing, such that data[i] is
        the featurized version of indices[i]-th frame in the MD trajectory
        with filename fns[i].
    """
    data = []
    indices = []
    fns = []

    for file in filenames:
        kwargs = {} if file.endswith('.h5') else {'top': topology}
        count = 0
        for t in md.iterload(file, chunk=chunk, stride=stride, **kwargs):
            x = featurizer.partial_transform(t)
            n_frames = len(x)

            data.append(x)
            indices.append(count + (stride * np.arange(n_frames)))
            fns.extend([file] * n_frames)
            count += (stride * n_frames)
    
    if len(data) == 0:
        raise ValueError("No data loaded!")

    return np.concatenate(data), np.concatenate(indices), np.array(fns)


def featurize_best_frames_pdb(pdb_path, output_dir='dihedral_features', 
                             angles=['phi', 'psi'], sincos=True):
    """
    Perform dihedral featurization on best_frames.pdb file
    
    Parameters
    ----------
    pdb_path : str
        Path to the best_frames.pdb file or directory containing it
    output_dir : str, default='dihedral_features'
        Directory to save the featurized data
    angles : list of str, default=['phi', 'psi']
        Types of dihedral angles to compute
        Options: 'phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4'
    sincos : bool, default=True
        Whether to compute sin/cos of angles instead of raw angles
        
    Returns
    -------
    features : np.ndarray
        Featurized dihedral angles
    feature_descriptions : list
        Descriptions of each feature
    """
    # Handle both file path and directory path
    if os.path.isdir(pdb_path):
        pdb_file = os.path.join(pdb_path, 'best_frames.pdb')
    else:
        pdb_file = pdb_path
    
    # Check if file exists
    if not os.path.exists(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    print(f"Loading PDB file: {pdb_file}")
    
    try:
        # Load the PDB file (no separate topology needed for PDB)
        trajectory = md.load(pdb_file)
        print(f"Loaded {trajectory.n_frames} frames with {trajectory.n_atoms} atoms")
        
        # Setup featurizer with specified dihedral angles
        featurizer = DihedralFeaturizer(types=angles, sincos=sincos)
        print(f"Dihedral featurizer setup with types: {angles}, sincos: {sincos}")
        
        # Featurize the trajectory
        print("Computing dihedral features...")
        features = featurizer.partial_transform(trajectory)
        print(f"Featurized trajectory shape: {features.shape}")
        
        # Get feature descriptions
        feature_descriptions = featurizer.describe_features(trajectory)
        print(f"Number of features: {len(feature_descriptions)}")
        
        # Save results
        os.makedirs(output_dir, exist_ok=True)
        
        # Save features as pickle (as list containing single numpy array for the trajectory)
        features_list = [features]  # Wrap the features array in a list
        features_file = os.path.join(output_dir, 'best_frames_features.pkl')
        with open(features_file, 'wb') as f:
            pickle.dump(features_list, f)
        print(f"Saved features as list containing trajectory array to {features_file}")
        
        # Save feature descriptions as DataFrame pickle
        desc_df = pd.DataFrame(feature_descriptions)
        desc_file = os.path.join(output_dir, 'feature_descriptions.pkl')
        with open(desc_file, 'wb') as f:
            pickle.dump(desc_df, f)
        print(f"Saved feature descriptions as DataFrame to {desc_file}")
        
        # Save as CSV for readability
        csv_file = os.path.join(output_dir, 'feature_descriptions.csv')
        desc_df.to_csv(csv_file, index=False)
        print(f"Saved feature descriptions to {csv_file}")
        
        # Save as numpy array for easy loading
        npy_file = os.path.join(output_dir, 'best_frames_features.npy')
        np.save(npy_file, features)
        print(f"Saved features as numpy array to {npy_file}")
        
        # Print summary
        print(f"\nFeaturization Summary:")
        print(f"Input file: {pdb_file}")
        print(f"Number of frames: {trajectory.n_frames}")
        print(f"Number of features: {features.shape[1]}")
        print(f"Feature range: [{np.min(features):.3f}, {np.max(features):.3f}]")
        print(f"Output directory: {output_dir}")
        
        return features, feature_descriptions
        
    except Exception as e:
        print(f"Error during featurization: {e}")
        raise


def example_usage():
    """
    Example usage of the TrajectoryDihedralAnalyzer
    """
    # Initialize analyzer
    analyzer = TrajectoryDihedralAnalyzer(topology_file='topology.pdb')
    
    # Example trajectory files (replace with your actual files)
    trajectory_files = [
        'trajectory1.xtc',
        'trajectory2.xtc',
        'trajectory3.dcd'
    ]
    
    try:
        # Load trajectories
        analyzer.load_trajectories(trajectory_files, stride=10)  # Load every 10th frame
        
        # Setup dihedral featurizer
        analyzer.setup_dihedral_featurizer(types=['phi', 'psi', 'chi1'], sincos=True)
        
        # Featurize trajectories
        features = analyzer.featurize_trajectories()
        
        # Get feature information
        feature_info = analyzer.get_feature_info()
        print("\nFeature Information:")
        print(feature_info.head())
        
        # Plot time series (optional - requires matplotlib)
        try:
            analyzer.plot_time_series(traj_idx=0, max_features=4)
        except ImportError:
            print("Matplotlib not available - skipping plots")
        
        # Save results
        analyzer.save_features('dihedral_features')
        
        # Print summary statistics
        all_features = np.vstack(features)
        print(f"\nSummary Statistics:")
        print(f"Total frames: {all_features.shape[0]}")
        print(f"Total features: {all_features.shape[1]}")
        print(f"Feature range: [{np.min(all_features):.3f}, {np.max(all_features):.3f}]")
        print(f"Feature means: {np.mean(all_features, axis=0)[:5]}...")  # First 5 features
        
    except Exception as e:
        print(f"Error in example usage: {e}")
        print("Make sure to provide valid trajectory and topology files!")


def batch_featurize_directory(traj_dir, topology_file, output_dir='batch_features'):
    """
    Batch featurize all trajectories in a directory
    
    Parameters
    ----------
    traj_dir : str
        Directory containing trajectory files
    topology_file : str
        Path to topology file
    output_dir : str, default='batch_features'
        Output directory for features
    """
    # Find all trajectory files
    extensions = ['*.xtc', '*.dcd', '*.trr', '*.h5']
    traj_files = []
    
    for ext in extensions:
        traj_files.extend(glob.glob(os.path.join(traj_dir, ext)))
    
    if not traj_files:
        print(f"No trajectory files found in {traj_dir}")
        return
    
    print(f"Found {len(traj_files)} trajectory files")
    
    # Use the standalone featurizer
    featurizer = DihedralFeaturizer(types=['phi', 'psi'], sincos=True)
    
    try:
        # Featurize all trajectories at once
        features, indices, filenames = featurize_all(
            traj_files, 
            featurizer, 
            topology_file,
            chunk=1000,
            stride=1
        )
        
        print(f"Featurized {len(traj_files)} trajectories")
        print(f"Total features shape: {features.shape}")
        
        # Save results
        os.makedirs(output_dir, exist_ok=True)
        
        # Save as pickle files
        features_file = os.path.join(output_dir, 'features.pkl')
        indices_file = os.path.join(output_dir, 'indices.pkl')
        filenames_file = os.path.join(output_dir, 'filenames.pkl')
        
        with open(features_file, 'wb') as f:
            pickle.dump(features, f)
        with open(indices_file, 'wb') as f:
            pickle.dump(indices, f)
        with open(filenames_file, 'wb') as f:
            pickle.dump(filenames, f)
        
        print(f"Results saved to {output_dir}/ as pickle files")
        
        return features, indices, filenames
        
    except Exception as e:
        print(f"Error in batch featurization: {e}")
        return None, None, None


if __name__ == "__main__":
    # Run featurization on best_frames.pdb
    print("Standalone Dihedral Featurization")
    print("=" * 40)
    
    # Parse command line arguments
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Featurize molecular dynamics structures using dihedral angles',
        epilog='''
Examples:
  %(prog)s best_frames.pdb
    Basic usage with default settings (phi, psi angles with sin/cos)
  
  %(prog)s --angles phi psi chi1 chi2
    Compute backbone phi/psi and first two sidechain chi angles
  
  %(prog)s --sincos false --output raw_angles
    Use raw angles instead of sin/cos transformation
  
  %(prog)s /path/to/protein.pdb --angles phi psi omega chi1 chi2 chi3 chi4
    Compute all available dihedral angles from a custom path
  
  %(prog)s --angles chi1 chi2 chi3 chi4 --output sidechain_features
    Focus only on sidechain conformations

Available dihedral angles:
  phi    : Backbone phi (Ï†) dihedral angle (C-N-CA-C)
  psi    : Backbone psi (Ïˆ) dihedral angle (N-CA-C-N) 
  omega  : Backbone omega (Ï‰) dihedral angle (CA-C-N-CA)
  chi1   : First sidechain (Ï‡1) dihedral angle
  chi2   : Second sidechain (Ï‡2) dihedral angle
  chi3   : Third sidechain (Ï‡3) dihedral angle
  chi4   : Fourth sidechain (Ï‡4) dihedral angle

Output files:
  best_frames_features.pkl     : Features as pickle (for Python)
  best_frames_features.npy     : Features as NumPy array
  feature_descriptions.csv     : Human-readable feature descriptions
  feature_descriptions.pkl     : Feature descriptions as pickle

Notes:
  - Sin/cos transformation is recommended for machine learning (periodic angles)
  - Raw angles are in radians and range from -Ï€ to Ï€
  - PDB files contain topology, so no separate topology file is needed
  - Chi angles are only computed for residues that have them
        ''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('pdb_path', nargs='?', default='best_frames.pdb',
                        help='Path to PDB file or directory containing best_frames.pdb (default: %(default)s)')
    
    parser.add_argument('--angles', nargs='+', default=['phi', 'psi'],
                        choices=['phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4'],
                        metavar='ANGLE',
                        help='Dihedral angles to compute (default: %(default)s). '
                             'Choices: {%(choices)s}')
    
    parser.add_argument('--sincos', type=str, default='true',
                        choices=['true', 'false'],
                        metavar='BOOL',
                        help='Compute sin/cos of angles instead of raw angles (default: %(default)s). '
                             'Sin/cos is recommended for machine learning. Choices: {%(choices)s}')
    
    parser.add_argument('--output', default='dihedral_features',
                        metavar='DIR',
                        help='Output directory for featurized data (default: %(default)s)')
    
    parser.add_argument('--version', action='version', version='%(prog)s 1.0',
                        help='Show version information and exit')
    
    # Handle case where script is run without argparse (for backwards compatibility)
    try:
        args = parser.parse_args()
        pdb_path = args.pdb_path
        angles = args.angles
        sincos = args.sincos.lower() == 'true'
        output_dir = args.output
    except SystemExit:
        # Fallback to simple command line parsing if argparse fails
        if len(sys.argv) > 1:
            pdb_path = sys.argv[1]
        else:
            pdb_path = 'best_frames.pdb'
        angles = ['phi', 'psi']
        sincos = True
        output_dir = 'dihedral_features'
    
    print(f"Configuration:")
    print(f"  PDB path: {pdb_path}")
    print(f"  Angles: {angles}")
    print(f"  Sin/Cos: {sincos}")
    print(f"  Output: {output_dir}")
    print()
    
    try:
        features, descriptions = featurize_best_frames_pdb(
            pdb_path, 
            output_dir=output_dir,
            angles=angles, 
            sincos=sincos
        )
        print("Featurization completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        print("\nUsage:")
        print("  python script.py [pdb_path] [--angles angle1 angle2 ...] [--sincos true/false] [--output dir]")
        print("\nExamples:")
        print("  python script.py best_frames.pdb")
        print("  python script.py --angles phi psi chi1 chi2")
        print("  python script.py --sincos false")
        print("  python script.py /path/to/dir/ --angles phi psi omega --output my_features")
        print("  python script.py best_frames.pdb --angles phi psi chi1 chi2 chi3 chi4 --sincos true")
        print("\nAvailable angles: phi, psi, omega, chi1, chi2, chi3, chi4")
