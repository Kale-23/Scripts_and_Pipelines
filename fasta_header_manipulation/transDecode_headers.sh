#! /bin/bash

# one argument is provided that is used as target directory
if [ "$#" -ne 1 ]; then
    echo "provide target directory as argument"
    exit 1
fi

for file in $1/*pep
do
	if [[ -f "$file" ]]
        then
            base=$(basename $file)
            file_start="${base%%.*}"
            echo "working on $file"
		
            awk -F "::" '/^>/{printf(">%s\n", $2); next}{print}' $file > "$file_start".pep.fasta
            #mv tmp $fasta_file
        fi
done
