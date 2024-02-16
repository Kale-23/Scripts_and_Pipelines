#! /bin/bash

# check if a directory is provided as an argument
if [ -z "$1" ]; then
  echo "Please provide a directory as an argument."
  exit 1
fi

# Use the find command to search for FASTA files
# Starting from the directory specified by $1
# and including all subdirectories
find "$1" -type f -name "*.fasta" -o -name "*.fa" -o -name "*.fna" | while read -r file
do
  # Copy each found file to the current directory
  cp "$file" .
done

