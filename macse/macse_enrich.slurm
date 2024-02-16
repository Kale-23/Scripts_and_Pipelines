#! /bin/bash

#SBATCH --ntasks=1

## default 24 threads.
##SBATCH --cpus-per-task=2

## More than 128G
##SBATCH --mem=512000

#SBATCH --job-name="macse enrich"
#SBATCH --output=macse-enrich-%j.log
#BATCH --exclude=node117,node118

## repeats script name and arguments
echo "$0 $@"

## reads out exact macse arguments each time it is run
grep -x -A 1 "##search" $0 | tail -n 1
##search
java -jar ~/scripts/macse/macse_v2.07.jar -prog enrichAlignment -align $1 -seq $2 ##-maxSTOP_inSeq 2 -maxFS_inSeq 5 ##-maxINS_inSe:wq 20 ##-maxDEL_inSeq 40
