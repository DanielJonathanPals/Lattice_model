#!/bin/bash

echo "Reading in the parameters..."
echo ""

# Initialize an array to store the right-hand sides
rhs_array=()

# Read the file line by line
while IFS= read -r line || [ -n "$line" ]; do
    echo "$line"
    # Use awk to extract the right-hand side of each equation
    rhs=$(echo "$line" | awk -F'=' '{print $2}' | tr -d '[:space:]')
    
    # Add the right-hand side to the array
    rhs_array+=($rhs)
done < "parameters.txt"

echo ""
echo "Executing julia script..."
echo ""

printf "%s\n" "${rhs_array[@]}" | xargs julia src/Run_simulation.jl
