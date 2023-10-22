#! /bin/bash

# one argument is provided that is used as target directory
if [ "$#" -ne 1 ]; then
    echo "provide target directory as argument"
    exit 1
fi

# loops through each fasta file in directory
for fasta_file in "$1"/*.fasta
do
	if [[ -f "$fasta_file" ]]
        then
            base=$(basename $fasta_file)
            species="${base%%.*}"
            # changes headers to ">species_sequence#"
            echo "working on $fasta_file"
            awk "/^>/{printf(\">"$species"_%d\\n\"), ++i; next}{print}" $fasta_file > tmp
            mv tmp $fasta_file
        fi
done
