#!/bin/bash
#SBATCH --ntasks=1

#SBATCH --job-name="align_pipeline"
#SBATCH --output=%j_align.log
#BATCH --exclude=node117,node118

usage() {
	echo "Usage: $0 QUERY_SEQ OUTPUT_DIR BLAST_DB FASTA_DIR ALIGNMENT_TYPE"
	echo "  - QUERY_SEQ: Path to the query sequence file."
	echo "  - OUTPUT_DIR: Directory to store the outputs."
	echo "  - BLAST_DB: Path to the BLAST database."
	echo "  - FASTA_DIR: Directory containing FASTA files"
	echo "  - ALIGNMENT_TYPE: (single/multi) either combine fastas"
	echo "                    to create single alignment, or seperate by individual"
	echo "                    Default=multi"
}

if [ "$#" -ne 5 ]; then
	usage
	exit 1
fi

QUERY_SEQ="$1"
OUTPUT_DIR="$2"
BLAST_DB="$3"
FASTA_DIR="$4"

if [[ $5 -ne "single" ]]; then
	ALIGNMENT_TYPE="multi"
else
	ALIGNMENT_TYPE="$5"
fi

echo "$(date): begin running script:"
echo "$0 $@"

# blasting query
BLAST_DIR="$OUTPUT_DIR"/blast_output
mkdir -p "$BLAST_DIR"
for file in "$BLAST_DB"/*.fasta "$BLASTDB"/*.fa; do
	if [[ -f $file ]]; then
		BASE_FILE="$(basename $file)"
		echo "Blasting "$BASE_FILE""
		blastn -query "$QUERY_SEQ" \
			-db $file \
			-out "$BLAST_DIR"/"${BASE_FILE%%.*}"_ref_blastout \
			-outfmt 6 \
			-max_target_seqs 10 \
			-evalue 0.00001 \
			-num_threads 24
	fi
done

# pulling headers out of blast output files
HIT_DIR="$OUTPUT_DIR"/hit_files
mkdir -p "$HIT_DIR"
for file in "$BLAST_DIR"/*; do
	if [[ -f $file ]]; then
		BASE_FILE="$(basename $file)"
		echo "pulling hits from $BASE_FILE and deleting duplicates"
		cut -f2 "$file" | sort | uniq >"$HIT_DIR"/${BASE_FILE%%.*}_hits
	fi
done

# pulling fasta sequences from header files
RECOVERED_FASTAS="$OUTPUT_DIR"/recovered_fastas
mkdir -p "$RECOVERED_FASTAS"
for file in "$HIT_DIR"/*; do
	if [[ -f $file ]]; then
		BASE_FILE="$(basename $file)"
		echo "pull FASTA seqs from  $BASE_FILE"
		REF_FILE=$(ls "$FASTA_DIR"/"${BASE_FILE%_ref*}"*)
		echo $REF_FILE
		~/progs/sequence_extraction/selectSeqs.pl -f "$file" $REF_FILE >"$RECOVERED_FASTAS"/${BASE_FILE%%_ref_blastout_hits*}.fasta
	fi
done

# combine all sequences into one output fasta
FULL_FASA="$OUTPUT_DIR"/full_fasta
mkdir -p $FULL_FASA
cat $RECOVERED_FASTAS/*.fasta >$FULL_FASA/all_recovered.fasta

# create alignment output locations
ALIGNMENTS="$OUTPUT_DIR"/alignments
mkdir -p "$ALIGNMENTS"/NT
mkdir -p "$ALIGNMENTS"/AA
mkdir -p "$ALIGNMENTS"/stats

if [[ $ALIGNMENT_TYPE == "single" ]]; then
	#SINGLE ALIGNMENT
	java -jar ~/progs/macse/macse_v2.07.jar -prog alignSequences \
		-seq "$FULL_FASA"/all_recovered.fasta \
		-out_NT "$ALIGNMENTS"/NT/${BASE_FILE%%.*}_NT.fasta \
		-out_AA "$ALIGNMENTS"/AA/${BASE_FILE%%.*}_AA.FASTA
else
	#INDIVIDUAL SPECIES ALIGNMENTS
	for file in "$RECOVERED_FASTAS"/*; do
		BASE_FILE="$(basename $file)"
		echo "running macse command:"
		grep -x -A 9 "##MACSE_COMMAND" $0 | tail -n 1
		##MACSE_COMMAND
		java -jar ~/progs/macse/macse_v2.07.jar -prog enrichAlignment \
			-align "$QUERY_SEQ" -seq "$file" \
			-maxSTOP_inSeq 2 \
			-maxFS_inSeq 5 \
			-out_NT "$ALIGNMENTS"/NT/${BASE_FILE%%.*}_NT.fasta \
			-out_AA "$ALIGNMENTS"/AA/${BASE_FILE%%.*}_AA.fasta \
			-out_tested_seq_info "$ALIGNMENTS"/stats/${BASE_FILE%%.*}_output_stats.csv #\
		#-maxINS_inSeq 20 \
		#-maxDEL_inSeq 40

		#	java -jar ~/progs/macse/macse_v2.07.jar -prog alignSequences \
		#		-seq "$file" \
		#		-out_NT "$ALIGNMENTS"/NT/${BASE_FILE%%.*}_NT.fasta \
		#		-out_AA "$ALIGNMENTS"/AA/${BASE_FILE%%.*}_AA.fasta
	done
fi
echo "$(date): finished script"
