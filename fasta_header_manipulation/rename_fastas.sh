#! /bin/bash

head_dir=$1

counter=1

for sub_dir in "$head_dir"/*
do
	for fasta in "$sub_dir"/*
	do
		if [[ -f $fasta ]]
		then

			dir_name=$(basename "$sub_dir")
			file_name=$(basename "$fasta")
			new_file_name="${dir_name}_${counter}.Trinity.fasta"
			counter=$(($counter + 1))
			new=$(dirname $fasta)/$new_file_name
			mv $fasta $new 
		fi		
	done
done
