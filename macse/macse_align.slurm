#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --job-name="macse align"
#SBATCH --output=mace_align-%j.log
#SBATCH --exclude=node117,node118
#
filename="${1%.*}"

java -jar ~/scripts/macse/macse_v2.07.jar -prog alignSequences -seq $1 -out_NT "$filename"_aligned_NT.fasta -out_AA "$filename"_aligned_AA.fasta
