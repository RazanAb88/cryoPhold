#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# 1. Install the conda environment from env.yaml
echo "Creating conda environment from env.yml..."
conda env create -f env.yml

# 2. Activate the conda environment
echo "Activating conda environment..."
source activate cryophold  # or: conda activate your_env_name

# 3. Clone the denss repository into the current directory
echo "Cloning denss repository..."
git clone https://github.com/tdgrant1/denss.git

# 4. Change to the denss directory
cd denss

# 5. Checkout the icosahedral_sym branch
echo "Switching to icosahedral_sym branch..."
git checkout icosahedral_sym

# 6. Install the denss package
echo "Installing denss package..."
pip install .

echo "Setup complete!"
