#! /usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import argparse
import os

# taking arguments
parser = argparse.ArgumentParser(description="trims alignments from the front or back at the position specified")

parser.add_argument("-i", required=True, type=str, help="input getorf fasta output")
parser.add_argument("-o", required=True, type=str, help="name of output fasta")
parser.add_argument("-s", "-seqs", required=True, type=str, help="directory of fasta files to search headers for")

args = parser.parse_args()

### Pull headers
seq_headers_to_find = {}
removed_headers = ()
pattern = r'\[(\d+) - (\d+)\]'
try:
    with open(args.i, "r") as fasta_in:
        for line in fasta_in:
            line_stripped = line.strip()
            if line_stripped.startswith(">"):
                line_stripped = line.lstrip(">").strip()
                split_header = line_stripped.split(" ")[0] 
                matches = re.findall(pattern, line_stripped) # should only find 1 match
                #print(matches)
                reverse = True if "REVERSE SENSE" in line_stripped else False
                if int(matches[0][0]) < int(matches[0][1]):
                    start = matches[0][0]
                    end = matches[0][1]
                else:
                    start = matches[0][1]
                    end = matches[0][0]
                #print(f"start: {start}, end: {end}, reverse: {reverse}, matches: {matches}")
                seq_headers_to_find[split_header] = (int(start), int(end), reverse)
except IOError as e:
    print(f"IOError {e}")
#print(seq_headers_to_find)

### finding headers in fastas + output ###
removed_headers = []
try:
    with open(args.o, "w") as fasta_out:
        for entry in os.scandir(args.s):
            if not entry.is_file():
                continue
            for record in SeqIO.parse(entry.path, "fasta"):
                header = record.id
                for seq_header in seq_headers_to_find:
                    cut_spot = seq_header.rfind("_")
                    if cut_spot >= 0:
                        seq_header_reduced = seq_header[:cut_spot]
                    if header == seq_header_reduced:
                        new_seq = str(record.seq)[:seq_headers_to_find[seq_header][1] + 1] if not seq_headers_to_find[seq_header][2] else str(record.seq)[seq_headers_to_find[seq_header][0]:seq_headers_to_find[seq_header][1]:-1]
                        print(f"adding {seq_header}\n{new_seq}")
                        SeqIO.write(SeqRecord(Seq(new_seq), id=record.id, description=""), fasta_out, "fasta")
                        removed_headers.append(seq_header)                    
except IOError as e:
    print(f"IOError {e}")

for header in removed_headers:
    seq_headers_to_find.pop(removed_headers)

# print if couldnt find sequences in args.s
if len(seq_headers_to_find) > 0:
    print("could not find sequences for:\n{0}".format("\n".join(seq_headers_to_find)))
else:
    print("all sequences found and added")
# open input and output files
