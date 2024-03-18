#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

# Function to translate peptide sequences to nucleotide sequences
def pep_to_nucl(input_file, output_file):
    with open(input_file, "r") as f:
        peptide_records = list(SeqIO.parse(f, "fasta"))

        nucleotide_records = []
        for record in peptide_records:
            peptide_seq = record.seq
            dna_seq = peptide_seq.back_transcribe()
            print(dna_seq)
            nucleotide_record = SeqRecord(dna_seq, id=record.id, description=record.description)
            nucleotide_records.append(nucleotide_record)

    # Writing nucleotide sequences to a new FASTA file
    with open(output_file, "w") as f:
        SeqIO.write(nucleotide_records, f, "fasta")

parser = argparse.ArgumentParser(description="turns peptide fasta to nucleotide fasta")

parser.add_argument("-i", required=True, help="input peptide fasta")
parser.add_argument("-o", required=True, help="name of output nucleotide fasta")

args = parser.parse_args()
# Call the function
pep_to_nucl(args.i, args.o)
