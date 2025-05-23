#! /usr/bin/env python3
import argparse
from Bio import SeqIO
import os

# taking arguments
parser = argparse.ArgumentParser(description="select sequences from blast output")

parser.add_argument("-f", required=True, type=str,
                    help="file of fasta headers")
parser.add_argument("-s", "-seqs", required=True, type=str,
                    help="single fasta file to search headers for")
parser.add_argument("-o", "-outfile", required=True, type=str,
                    help="location of output fasta")

args = parser.parse_args()

# getting headers from blast output #
seq_headers_to_find = {}

try:
    with open(args.f, "r") as fasta_in:
        for line in fasta_in:
            if line.startswith("names"):
                continue
            line_stripped = line.strip()
            if line_stripped.startswith(">"):
                line_stripped = line.lstrip(">").strip()
            split_list = line_stripped.split(",")
            if split_list[1].upper().startswith("MUC"):
                print(split_list)
                seq_headers_to_find[split_list[0]] = split_list[1]
except IOError as e:
    print(f"IOError {e}")

# finding headers in fastas + output #
name = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(args.f))))
try:
    with open(args.o, "w") as fasta_out:
        for record in SeqIO.parse(args.s, "fasta"):
            header = record.id
            if header in seq_headers_to_find.keys():
                record.id = seq_headers_to_find[header] + "_" + name + "_" + record.id
                SeqIO.write(record, fasta_out, "fasta")
                seq_headers_to_find.pop(header)
except IOError as e:
    print(f"IOError {e}")

# print if couldnt find sequences in args.s
if len(seq_headers_to_find) > 0:
    print("could not find sequences for:\n{0}".format("\n".join(seq_headers_to_find)))
else:
    print("all sequences found and added")
