#!/bin/bash

# Process all PDB files in the current directory
for pdb_file in *.pdb; do
    # Skip if no PDB files found (when glob doesn't match)
    [[ ! -f "$pdb_file" ]] && continue
    
    # Extract base name without extension (e.g., prot1.pdb -> prot1)
    base_name=$(basename "$pdb_file" .pdb)
    
    # Create directory named after the PDB file
    mkdir "Model-${base_name}"
    cd "Model-${base_name}"
    
    # Copy the specific PDB file as prot.pdb
    cp "../${pdb_file}" prot.pdb
    cp ../tleap.all .
    
    echo "Processing ${pdb_file} in Model-${base_name}..."
    
    tleap -s -f tleap.all
    
    # Energy minimization with sander
    sander -O -i ../inputs/full-min.in -o full-min.out -p com_solvated.top -c com_solvated.crd -ref com_solvated.crd -r Full_Mini.rst
    
    # Convert to PDB format
    ambpdb -p com_solvated.top -c Full_Mini.rst > Full_Mini.pdb
    
    # Convert to GROMACS format
    echo 0 | gmx trjconv -s Full_Mini.pdb -f Full_Mini.pdb -o Full-min.gro
    
    # Convert AMBER topology to GROMACS
    acpype -p com_solvated.top -x com_solvated.crd -b gmx
    
    # GROMACS processing
    cd gmx.amb2gmx
    cp ../../insert-posre.sh .
    
    # Create index and restraints
    echo q | gmx make_ndx -f gmx_GMX.gro
    echo 4 | gmx genrestr -f gmx_GMX.gro -fc 500 500 500 -n index.ndx
    
    # Modify topology
    bash insert-posre.sh gmx_GMX.top
    
    # Prepare and run energy minimization
    gmx grompp -f ../../inputs/em.mdp -c ../Full-min.gro -p gmx_GMX.top -o em.tpr -maxwarn 1
    gmx mdrun -ntmpi 1 -v -deffnm em
    
    # Clean up and final processing
    rm -f em-use.gro
    echo 0 | gmx trjconv -s ../com_solvated.pdb -f em.gro -o em-use.gro -ur compact -pbc nojump
    
    # Return to parent directory
    cd ../..
    
    echo "Completed processing ${pdb_file}"
done

echo "All PDB files processed successfully!"
