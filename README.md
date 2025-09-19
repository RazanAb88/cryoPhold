# CryoPhold: CryoEM meets AlphaFold and molecular simulation to reveal protein dynamics

CryoPhold is a computational pipeline that performs Bayesian reweighting of AlphaFold predicted structural ensembles against experimental cryo-EM data. The pipeline determines optimal weights for each AlphaFold structure and identifies the minimal subset that best reproduces experimental density maps, enabling the elucidation of protein dynamics through the integration of AI-predicted structures and experimental data.

## 🧬 Scientific Background

CryoPhold implements a maximum entropy approach (BioEN algorithm) to:

1. **Generate simulated density maps** from each AlphaFold structure
2. **Optimize structural weights** to match experimental cryo-EM data
3. **Posterior structural ensemble** while preserving fit quality and conformational heterogeneity 
  
## 🚀 Key Features

- **Bayesian Ensemble Reweighting**: Optimally combines AlphaFold structures using experimental cryo-EM constraints
- **Structure Selection**: Identifies minimal structural subsets that maintain fit quality
- **Comprehensive Analysis**: Provides correlation analysis, Fourier Shell Correlation, RMSF, and PCA
- **Automated Pipeline**: End-to-end workflow from AlphaFold generation to final analysis

## 📦 Installation

```bash
git clone git@github.com:strauchlab/cryoPhold.git
cd cryoModule
chmod +x setup.sh
./setup.sh
conda activate cryophold
```

### Test Installation
```bash
python ./cryoModule/cryoPhold.py -h
```

**Boom! You should be ready to go!** 🎉

## 🔧 Usage Workflow

### Step 1: Generate AlphaFold Structural Ensemble

Check if `colabfold` is installed in the `cryophold` enviornment using:

```bash
colabfold_batch -h
```

Generate a conformational ensemble using ColabFold with dropout to sample structural diversity:

```bash
colabfold_batch \
  --num-recycle 3 \
  --recycle-early-stop-tolerance 0.5 \
  --num-ensemble 5 \
  --model-type auto \
  --templates \
  --use-gpu-relax \
  --amber \
  --max-seq 8 \
  --max-extra-seq 16 \
  --num-seeds 16 \
  --use-dropout \
  --relax-max-iterations 200 \
  ./input/ \
  ./results/
```

**Quick start**: Use the provided script `./scripts/colabfold.sh`

### Step 2: Align Ensemble to Cryo-EM Map

Align the AlphaFold generated conformational ensemble to your experimental cryo-EM map using Situs:

- Use the **Situs standalone tool**: https://situs.biomachina.org
- Run the alignment script: `./scripts/align.sh`
- **Important**: Update the script with your specific:
  - Map resolution
  - Path to cryo-EM map  
  - Number of processors

**Alternative**: If you have a corresponding PDB, align the AlphaFold ensemble using PyMOL or UCSF Chimera and save as `combined.pdb`. Transform the `combined.pdb` to a .xtc format using `gmx`, `mdtraj` or `mdanalysis`

### Step 3: Align and Combine Structures

Process the aligned structures:

```bash
# Refine aligned structures
./scripts/refine.sh

# Visually inspect alignment quality with cryo-EM map
# Combine all aligned PDBs into trajectory format
python ./scripts/combine-pdb.py --remove-hydrogens  # if needed
```

This generates `combined.xtc` for the next step.

We have already provided `AF2` and `Boltz` generated structural ensemble as examples.

### Step 4: Bayesian Reweighting 

Run the main CryoPhold pipeline:

```bash
python cryoPhold.py \
  --path 'path_to_directory_containing_combined.xtc' \
  --threshold 2.95 \
  --resolution 3.05 \
  --mask sim
```

#### Parameters:
- `--path`: Directory containing AlphaFold output (`combined.xtc`, `prot-masses.pdb`, `reference.map`)
- `--threshold`: Cryo-EM map density threshold (default: 2.95)  
- `--resolution`: Map resolution in Å (default: 3.05)
- `--mask`: Masking strategy - `'exp'` (experimental) or `'sim'` (simulated, default)

## 📊 Output Files

### Primary Results
| File | Description |
|------|-------------|
| `map_posterior.mrc` | Bayesian reweighted density map |
| `map_prior.mrc` | Uniform weighted (prior) map |
| `weights.dat` | Optimal weights for each AlphaFold structure |
| `best.pdb` | Best single structure with RMSF as B-factors |
| `statistics.dat` | Comprehensive analysis statistics |

### Analysis Plots
| File | Description |
|------|-------------|
| `lcurve.svg` | L-curve showing entropy vs χ² trade-off |
| `weights.svg` | Structure weight distribution |
| `FSC.svg` | Fourier Shell Correlation analysis |
| `RMSF.svg` | Root Mean Square Fluctuation comparison |
| `PCA.svg` | Principal Component Analysis |

### Iterative Refinement (`iter/` directory)
The pipeline automatically identifies the minimal structural subset:
- `iter/best_frames.pdb` - Final optimized structural ensemble
- `iter/map_posterior_iter.mrc` - Refined posterior map
- `iter/statistics_iter.dat` - Iteration-specific statistics
- Updated plots for each refinement iteration

## Featurize CryoPhold generated PDB and ML Module
```bash
python ./scripts/dihedral.py -h
```

Use: 
```bash
python dihedral.py iter/best_frames.pdb --angles phi psi chi1
```

This will create a folder `dihedral_features`

Go to `mlmodule` and install the software train the features using ML and it will create a free energy surface, residue hotspots and `plumed.dat` file. You can use the `plumed.dat` file to run Enhanced sampling calculations

## Simulation Module

Now you want to use the `best_frames.pdb` generated by `cryoPhold` to launch new simulations.

### Step 1: Split the PDBs

```bash
python ./scripts/split_pdb.py -h
usage: split_pdb.py [-h] -i INPUT [-o OUTDIR]

Split a multi‐MODEL PDB into one file per MODEL block.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to input multi‐model PDB
  -o OUTDIR, --outdir OUTDIR
                        Directory to write individual PDBs
```
### Step 2: Prepare the PDBs to launch simulations

Install `Gromacs`, `AmberTools` and `acpype`. Go to `simmodule` and run the `process.sh` script to prepare the PDB files.

**You are ready to launch simulations from cryoPhold ensemble!** 🎉
**Try launching multiple independent metadynamics simulations from each structures using `plumed.dat` file!** 🎉

### Step 3: Featurize the MD data

```bash
python ./scripts/dihedral.py -h
```

Use: 
```bash
python dihedral.py --directory ./trajectories  --angles phi psi chi1
```

### Step 4: Train ML models on the MD data

Create a new `conda` enviornment using `MDML package (https://github.com/svats73/mdml)` or use the `mlmodule` to train the featurized data using ML model.

**Boom! Build MSM, capture Hotspots, perform Clustering and find allosteric pockets!** 🎉

## 🤝 Contributing

We welcome contributions! Please feel free to submit issues and pull requests.

## 📄 License

Please refer to the repository license for usage terms.

## 💜 Builders

Soumendranath Bhakat, Shray Vats, Andreas Mardt and Eva M. Strauch

## 🙏 Citation

```bash
@article {Bhakat2025.09.12.675912,
	author = {Bhakat, Soumendranath and Vats, Shray and Mardt, Andreas and Strauch, Eva M.},
	title = {CryoPhold: CryoEM meets AlphaFold and molecular simulation to reveal protein dynamics},
	elocation-id = {2025.09.12.675912},
	year = {2025},
	doi = {10.1101/2025.09.12.675912},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2025/09/13/2025.09.12.675912},
	eprint = {https://www.biorxiv.org/content/early/2025/09/13/2025.09.12.675912.full.pdf},
	journal = {bioRxiv}
}
```

**Bridge the gap between AI-predicted structures and cryo-EM? Start your CryoPhold journey today!** ✨
