#! /bin/bash
JUST_STATS=$1

for i in ~/Desktop/newer/*; do
	if [[ -d $i ]]; then
		name=$(basename "$i" | cut -d'_' -f1,2)
		ident=$(basename "$i" | cut -d'_' -f3)
		echo
		echo "### Running R Script on "$name"_"$ident" ###"
		echo
		#mkdir -p "$i"/plots
		#mkdir -p "$i"/data
		#Rscript ~/progs/complete_analysis.r $i "$name".TRINITY.hit.fasta.transdecoder_fixed_headers.pep 9 $JUST_STATS
		echo
		echo "### Running Python Script on "$name"_"$ident" ###"
		echo
		#~/progs/select_seq_single_list.py -f "$i"/"$name"_missing_genes.tsv -s "$i"/"$name".TRINITY.hit.fasta -o "$i"/"$name"_"$ident"_missing_genes.fasta
		echo
		echo "### Running orfipy on "$name"_"$ident" ###"
		echo
		#orfipy "$i"/"$name"_"$ident"_missing_genes.fasta --outdir "$i"/orfipy_out --procs 5 --bed "$name"_"$ident"_orfs.bed --pep "$name"_"$ident"_orfs.pep.fasta --longest
		echo
		echo "### Running blast on "$name"_"$ident" missing labeled seqs ###"
		echo
		mkdir -p "$i"/blastout
		source /opt/homebrew/anaconda3/bin/activate blast_osx-64
		blastp -query "$i"/orfipy_out/"$name"_"$ident"_orfs.pep.fasta -db ~/Desktop/Orthofinder/swissprot_db/uniprotkb_reviewed_true_AND_reviewed_tr_2024_02_27.fasta -out "$i"/blastout/"$name"_"$ident"_blastout -outfmt 6 -evalue 0.00001 -num_threads 10
		source /opt/homebrew/anaconda3/bin/deactivate
	fi
done
