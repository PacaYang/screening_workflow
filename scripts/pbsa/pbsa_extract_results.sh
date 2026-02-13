#!/bin/bash

pbsa_dir=$1
outdir=$2

# Output file
output="${outdir}/tmp.csv"

# Write header
echo "Folder,Delta_TOTAL,SD(Prop.),SD,SEM(Prop.),SEM" > "$output"

# Loop over each file
find $pbsa_dir/. -type f -name "FINAL_RESULTS_MMPBSA.dat" | while read -r file; do
    # Extract folder name (you can adjust this if you want just basename)
    folder=$(dirname "$file")
    # Extract the ΔTOTAL line in the Delta section only
    line=$(awk '/Delta \(Complex - Receptor - Ligand\)/, /Using Interaction Entropy Approximation/' "$file" | grep 'ΔTOTAL')

    # Extract values from the line (assuming fixed column widths)
    if [[ -n "$line" ]]; then
        # Remove leading/trailing spaces and extract fields
        values=$(echo "$line" | awk '{print $2","$3","$4","$5","$6}')
        echo "$folder,$values" >> "$output"
    fi
done
