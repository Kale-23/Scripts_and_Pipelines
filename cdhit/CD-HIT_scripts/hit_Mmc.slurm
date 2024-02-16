#! /bin/bash

#SBATCH --job-name="cdhit"
#SBATCH --output="hit_Mmc.log"

##SBATCH --cpus-per-task=2
##SBATCH --mem=128000
#SBATCH --ntasks=1

module purge
module load linuxbrew/colsa

cd-hit -i ./M_mcmillanae.TRINITY.fasta -o ./M_mcmillanae.TRINITY.hit.fasta -c 0.98 -n 5 -T 0 -M 120000
