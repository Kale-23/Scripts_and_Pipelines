#! /bin/bash

#SBATCH --job-name="TransDecode"
#SBATCH --output="TransDecode_individuals_1.log"

##SBATCH --cpus-per-task=2
##SBATCH --mem=128000
#SBATCH --ntasks=1

module purge
module load linuxbrew/colsa

for file in "$1"/*
do
	if [[ -f $file ]]
	then
		TransDecoder.LongOrfs -t $file
		TransDecoder.Predict -t $file
	fi
done 
