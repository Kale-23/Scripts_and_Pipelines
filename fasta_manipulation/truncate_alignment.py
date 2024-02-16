#! /usr/bin/env python3

from Bio import SeqIO
import argparse

# taking arguments
parser = argparse.ArgumentParser(description="trims alignments from the front or back at the position specified")

parser.add_argument("-i", required=True, type=str, help="input alignment")
parser.add_argument("-o", required=True, type=str, help="name of output alignment")
parser.add_argument("-l", required=True, type=str, choices=["front", "back"], help="location of trimming: front|back")
parser.add_argument("-p", required=True, type=int, help="position of truncation (front: remove <= position| back: remove >= position)")

args = parser.parse_args()


def is_only_dashes(string):
    #determines if sequence is completely removed after trimming

    for char in string:
        if char != "-":
            return False
    return True

# open input and output files
try:
    with open(args.i, "r") as input_handle, open(args.o, "w") as output_handle:
        records = list(SeqIO.parse(input_handle, "fasta"))
        
        # replace record with truncated record and output each to output
        for record in records:
            if args.l == "back":
                truncated_sequence = record.seq[:args.p - 1]
            if args.l == "front":
                truncated_sequence = record.seq[args.p:]
            if is_only_dashes(truncated_sequence):
                print(f"full sequence of record {record.id} removed, removing from alignment")
                continue
            truncated_record = record
            truncated_record.seq = truncated_sequence
            SeqIO.write(truncated_record, output_handle, "fasta")
except IOError as e:
    print(f"IOError {e}")

print(f"sequences in {args.i} truncated at position {args.p}, new alignment at {args.o}")
