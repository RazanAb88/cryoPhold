#!/bin/bash

# Script to insert POSRES block above moleculetype section in gmx_GMX.top

file="gmx_GMX.top"

# Check if file exists
if [ ! -f "$file" ]; then
    echo "Error: File $file not found!"
    exit 1
fi

# Create backup
cp "$file" "${file}.backup"
echo "Created backup: ${file}.backup"

# POSRES block to insert
posres_block="#ifdef POSRES
#include \"posre.itp\"
#endif"

# Find all [ moleculetype ] lines and insert before the second-to-last one
# First, find line numbers of all [ moleculetype ] occurrences
moleculetype_lines=$(grep -n '^\[ moleculetype \]' "$file" | cut -d: -f1)
line_count=$(echo "$moleculetype_lines" | wc -l)

if [ "$line_count" -lt 2 ]; then
    echo "Error: Need at least 2 [ moleculetype ] sections to insert before second-to-last"
    exit 1
fi

# Get the second-to-last line number
second_to_last_line=$(echo "$moleculetype_lines" | tail -2 | head -1)

echo "Found $line_count [ moleculetype ] sections"
echo "Inserting POSRES block before line $second_to_last_line"

# Use awk to insert the POSRES block at the specific line
awk -v line="$second_to_last_line" '
NR == line {
    print "#ifdef POSRES"
    print "#include \"posre.itp\""
    print "#endif"
    print ""
}
{print}
' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"

echo "POSRES block inserted successfully in $file"
echo "Original file backed up as ${file}.backup"
