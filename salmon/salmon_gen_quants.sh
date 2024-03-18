#! /bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name="salmon_quant"
#SBATCH --output=salmon-%J.log
#SBATCH --exclude=node117,node118

show_usage() {
	echo "Usage: $0 <index_directory> <ind_raw_reads_directory> <quant_output_directory>"
}

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
	show_usage
	exit 0
fi

if [ "$#" -ne 3 ]; then
	echo "incorrect number of arguments"
	show_usage
	exit 1
fi

module purge
module load anaconda/colsa # testing
conda activate salmon      # testing
#module load linuxbrew/colsa # testing conda instead of old salmon in this module

for dir in $2/*; do
	if [[ -d $dir ]]; then
		sample="$(basename $dir)"
		echo "processing sample $sample"

		salmon quant -p 24 \
			--seqBias \
			--gcBias \
			-i "$1" \
			-l A \
			-1 "$dir"/*_1.fq.gz \
			-2 "$dir"/*_2.fq.gz \
			-o "$(basename $3)"/"$(basename $2)"/"$sample"

	fi
done
