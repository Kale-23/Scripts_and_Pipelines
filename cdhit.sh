#! /bin/bash

if [[ $# -lt 1 ]]
then
	echo "Enter directory of fastas to run CD-HIT on all fastas within"
	exit 1
fi

for fasta in $1/*.fasta
do
	if [[ -f $fasta ]]
	then
		cd-hit -i $fasta -o ${fasta%.TRINITY.fasta}.TRINITY.hit.fasta -c 0.98 -n 5 -T 0 -M 120000
	fi
done
