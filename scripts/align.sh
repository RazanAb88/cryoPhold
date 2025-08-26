#!/bin/bash
#SBATCH --mem=60g
#SBATCH -p long
#SBATCH -c 12
#SBATCH --job-name=Fitting-glp1

# Create the Map directory if it doesn't exist
#mkdir -p Map

# Main processing loop
for i in {001..080}; do
    # Create the protein directory if it doesn't exist
    mkdir -p Prot${i}
    
    # Check if the directory was created successfully
    if [ ! -d "Prot${i}" ]; then
        echo "Error: Failed to create directory Prot${i}"
        continue
    fi
    
    # Find the corresponding AlphaFold2 PDB file using pattern matching
    source_file=$(find . -maxdepth 1 -name "prot_unrelaxed_rank_${i}_alphafold2_ptm_model_3_seed_*.pdb" | head -1)
    
    # Check if the source PDB file exists
    if [ -z "$source_file" ] || [ ! -f "$source_file" ]; then
        echo "Warning: Source file prot_unrelaxed_rank_${i}_alphafold2_ptm_model_3_seed_*.pdb not found"
        continue
    fi
    
    echo "Found source file: $source_file"
    
    # Copy the corresponding PDB file to the directory
    cp "$source_file" Prot${i}/prot.pdb
    
    # Verify the copy was successful
    if [ ! -f "Prot${i}/prot.pdb" ]; then
        echo "Error: Failed to copy $source_file to Prot${i}/prot.pdb"
        continue
    fi
    
    # Change to the protein directory
    cd Prot${i}
    
    # Check if the required reference map file exists
    if [ ! -f "../Map/reference.map" ]; then
        echo "Warning: Reference map file not found for Prot${i}"
        cd ..
        continue
    fi
    
    # Run the colores command
    # Change the resolution (-res) and number of processors according to your need
    # Give a correct path to the reference.map 
    echo "Processing Prot${i}..."
    ~/Situs/Situs_3.2/bin/colores ../Map/reference.map prot.pdb -res 3.8 -nprocs 12
    
    # Check if the command completed successfully
    if [ $? -eq 0 ]; then
        echo "Successfully processed Prot${i}"
    else
        echo "Error: colores command failed for Prot${i}"
    fi
    
    # Return to the parent directory
    cd ..
done

echo "Processing complete."
