#! /usr/bin/env python3

import os
import argparse
from Bio import SeqIO

# taking arguments
parser = argparse.ArgumentParser(description="takes transdecoded sequences and outputs original headers with gene number")

parser.add_argument("-s", "-seqs", required=True, type=str, help="location of transdecoder .pep files for input")
parser.add_argument("-o", "-output", required=True, type=str, help="location of output")

args = parser.parse_args()

# make output directory and format output path correctly
os.makedirs(args.o, exist_ok=True)
if not args.o.endswith("/"):
	args.o = args.o + "/"

# for every .pep file in input directory, for every record, fix the name and output to new file
for entry in os.scandir(args.s):
	if not entry.is_file():
		continue
	if entry.name.endswith(".pep"):
		print(f"Working on {entry.name}")

		# seperate file name with extension and add "_fixed_headers" between, and creates new file with that name for output
		last_dot_index = entry.name.rfind('.')
		outfile = entry.name[:last_dot_index] + "_fixed_headers"  + entry.name[last_dot_index:]
		with open(f"{args.o}{outfile}", "w") as out_handle:

			# for every record take the 1/2 data locations (original header/gene number) and save as new headers in output file
			for index, record in enumerate(SeqIO.parse(entry.path, "fasta"), start=1):
				fixed_header = ("_").join(record.id.split("::")[1:3])
				record.id = fixed_header
				record.description = ""
				SeqIO.write(record, out_handle, "fasta")
