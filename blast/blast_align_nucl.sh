#! /bin/bash

fasta_query=$1

for fasta in ../../blastdb/*.fasta
do
    echo "BLASTING  $fasta ..."
    blastn -query $fasta_query -db $fasta -out ${fasta%.}_ref_blastout -outfmt 6 -max_target_seqs 10 -evalue 0.00001 -num_threads 24
done

mkdir blastout

mv ../../blastdb/*blastout ./blastout/

#2 pull accs from blastout

for i in ./blastout/*blastout
do
  echo "getting hit accessions  $i ..."
  cut -f2 $i > ${i%.}_hits
done

mkdir hit1_accessions

mv ./blastout/*hits ./hit1_accessions

#3 pull seqs from blastout

for hitfile in ./hit1_accessions/*_ref_blastout_hits
do
    perl -p -e 's/\>//g' $hitfile > temp
    filenamefull=${hitfile##*/}
    filename="${filenamefull%_ref_blastout_hits*}"
    echo "pull FASTA seqs from  $filename ..."
    selectSeqs.pl -f temp ../../fastas/$filename >> $filename
    rm temp
done

mkdir hits_fasta

mv *.fasta ./hits_fasta

cat hits_fasta/*fasta > ed1.fa

cd-hit -i ed1.fa -o ed2.fa -c 1.00 -n 5    #cdhit --> looks for and removes redundant seq - c is threshold 

cat ed2.fa  $1 > ed3.fa  #adding in bait file 

#mafft --auto ed3.fa > ed4.fa    #mafft aligning my total file that just when thru cdhit

#mkdir rax     #the $1/rax will make a nested directory within $1 directory.... name it rax because you are going to run raxml with this alignment 

#mv *.fa rax

#rm *.clstr


