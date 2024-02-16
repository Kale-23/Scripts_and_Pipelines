#! /usr/bin/env python3

import os
import argparse
from Bio import SeqIO

# taking arguments
parser = argparse.ArgumentParser(description="select sequences from blast output")

parser.add_argument("-f", required=True, type=str, help="file of fasta headers")
parser.add_argument("-s", "-seqs", required=True, type=str, help="directory of fasta files to search headers for")
parser.add_argument("-o", "-outfile", required=True, type=str, help="location of output fasta")

args = parser.parse_args()

### getting headers from blast output ###

seq_headers_to_find = []

try:
    with open(args.f, "r") as fasta_in:
        for line in fasta_in:
            line_stripped = line.strip()
            if line_stripped.startswith(">"):
                line_stripped = line.lstrip(">").strip()
            seq_headers_to_find.append(line_stripped)
except IOError as e:
    print(f"IOError {e}")

### finding headers in fastas + output ###

try:
    with open(args.o, "w") as fasta_out:
        for entry in os.scandir(args.s):
            if not entry.is_file():
                continue
            for record in SeqIO.parse(entry.path, "fasta"):
                header = record.id
                if header in seq_headers_to_find:
                    SeqIO.write(record, fasta_out, "fasta")
                    seq_headers_to_find = [x for x in seq_headers_to_find if x != header]
                    
except IOError as e:
    print(f"IOError {e}")

# print if couldnt find sequences in args.s
if len(seq_headers_to_find) > 0:
    print("could not find sequences for {0}".format(", ".join(seq_headers_to_find)))
else:
    print("all sequences found and added")
