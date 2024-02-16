#! /bin/bash

fasta=$1

file=${1%.*}

java -jar ~/scripts/macse/macse_v2.07.jar -prog exportAlignment -align $1 -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS - -out_NT "$file"_NT_noFS.fasta -out_AA "$file"_AA_noFS.fasta -out_stat_per_seq "$file"_noFS_stats.csv

#echo "sep=;" > ./temp
#cat "$file"_noFS_stats.csv >> ./temp
#mv ./temp "$file"_noFS_stats.csv
