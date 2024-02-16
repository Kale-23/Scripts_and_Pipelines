#! /bin/bash

head_dir=$1

if [[ $# -lt 1 ]]
then
	echo "Enter directory of species directories to rename all fasta filenames as species.TRINITY.fasta"
	exit 1
fi

for sub_dir in "$head_dir"/*
do
	dir_name=$(basename "$sub_dir")
	new_file="$dir_name".TRINITY.fasta
	
	cat "$sub_dir"/*.fasta > "$new_file"	
done
