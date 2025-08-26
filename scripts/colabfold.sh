#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=6
#SBATCH --time=0-10:00:00
#SBATCH --mem=30G
#SBATCH --gres=gpu:1
#SBATCH --output=slurm_%j.out
#SBATCH --job-name=Colabfold

# Load required modules or activate your environment
# module load colabfold  (if available as a module)
# source activate my_env  (if using conda)

# Run the colabfold_batch command with specified parameters
# Generate multiple conformations by enabling the following
## --max-msa 16:32 or 8:16\
## --num-seeds 16 \
## --use-dropout True \

export LD_LIBRARY_PATH=/home/bhakat/localcolabfold/colabfold-conda/lib:$LD_LIBRARY_PATH

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
