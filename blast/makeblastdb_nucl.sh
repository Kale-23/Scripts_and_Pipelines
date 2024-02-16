#!/bin/bash

if [[ $# -lt 1 ]]; then
	echo "enter directory of fastas to run makeblastdb on all, creates new directory \"blastdb\" and moves outputs there"
	exit 1
fi

for i in $1/*.fasta; do
	echo "formating for BLAST  $i ..."
	makeblastdb -in $i -parse_seqids -dbtype nucl
done

mkdir blastdb

# moves everything to blastdb
mv $1/* blastdb

# copies .fasta files back to original directory
cp blastdb/*.fasta $1/
