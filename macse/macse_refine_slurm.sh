#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --job-name="macse refine"
#SBATCH --output=mace_refine-%j.log
#SBATCH --exclude=node117,node118

filename="${1%.*}"

java -jar ~/progs/macse/macse_v2.07.jar -prog refineAlignment -align $1 -out_NT "$filename"_aligned_NT.fasta -out_AA "$filename"_aligned_AA.fasta
