#! /bin/bash

fasta=$1

base_fasta=$(basename $fasta)

java -jar macse_v2.07.jar -prog exportAlignment -align $1 -charForRemainingFS - -out_NT "$base_fasta".noFS_NT.fasta -out_AA "$base_fasta".noFA_AA.fasta


