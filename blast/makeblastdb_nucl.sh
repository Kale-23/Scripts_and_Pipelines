#!/bin/bash

if [[ $# -lt 1 ]]; then
	echo "enter directory of fastas to run makeblastdb on all"
	exit 1
fi

for i in $1/*.fa; do
	if [[ -f $i ]]; then
		echo "formating for BLAST  $i ..."
		makeblastdb -in $i -parse_seqids -dbtype nucl
	fi
done
