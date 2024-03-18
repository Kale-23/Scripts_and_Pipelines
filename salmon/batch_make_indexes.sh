#!/bin/bash

# Take directory name as an argument
directory="$1"

# Check if the argument is provided
if [ -z "$directory" ]; then
	echo "Please provide a directory path as an argument."
	exit 1
fi

# Check if the directory exists
if [ ! -d "$directory" ]; then
	echo "Directory does not exist."
	exit 1
fi

# Loop through each file in the directory
for file in "$directory"/*; do
	if [ -f "$file" ]; then
		# Extract the filename without the extension
		filename=$(basename -- "$file")
		name_no_extension=${filename%%.*}

		# Run Salmon index maker command using the filename without extension
		sbatch ./salmon_make_index_new.sh $file "${name_no_extension}_index"
	fi
done
