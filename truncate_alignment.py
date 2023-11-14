#! /usr/bin/env python3

import Bio import SeqIO
import argparse

# taking arguments
parser = argparse.ArgumentParser(description="truncates alignments at position specified")

parser.add_argument("-i", required=True, type=str, help="input alignment")
parser.add_argument("-o", required=True, type=str, help="name of output alignment")
parser.add_argument("-p", required=True, type=int, help="position of truncation (remove >= value")

args = parser.parse_args()

# open input and output files
try:
	with open(args.i, "r") as input_handle, open(args.o, "w") as output_handle:
		records = list(SeqIO.parse(input_handle, "fasta"))
		
		# replace record with truncated record and output each to output
		for record in records:
		    truncated_sequence = record.seq[:args.p]
		    truncated_record = record
		    truncated_record.seq = truncated_sequence
		    SeqIO.write(truncated_record, output_handle, "fasta")
except IOError as e:
	print(f"IOError {e}")

print(f"sequences in {args.i} truncated at position {args.p}, new alignment at {args.o}")
