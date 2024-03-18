#! /bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name="salmon_index"
#SBATCH --output=salmon-%J.log
#SBATCH --exclude=node117,node118
module purge
module load anaconda/colsa
conda activate salmon
#module load linuxbrew/colsa
salmon index -p 24 -t $1 -i $2
