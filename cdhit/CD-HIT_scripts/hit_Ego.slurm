#! /bin/bash

#SBATCH --job-name="cdhit"
#SBATCH --output="hit_Ego.log"

##SBATCH --cpus-per-task=2
##SBATCH --mem=128000
#SBATCH --ntasks=1

module purge
module load linuxbrew/colsa

cd-hit -i ./temp_fastas/E_goslinei.TRINITY.fasta -o ./temp_fastas/E_goslinei.TRINITY.hit.fasta -c 0.98 -n 5 -T 0 -M 120000
