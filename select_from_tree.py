#! /usr/bin/env python3

from typing import Tuple

from Bio import Phylo
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import re
import argparse

def seq_sim(seq1: str, seq2: str) -> float:
	# jaccardian similarity to compare two strings
	s1 = set(seq1)
	s2 = set(seq2)
	intersection = len(s1.intersection(s2))
	union = len(s1.union(s2))
	return intersection / union 

def seq_sim_score(seq1: str, seq2: str, sim_score: float) -> Tuple[bool, float]:
	# returns whether two strings are within threshold using jaccardian similarity
	score = seq_sim(seq1, seq2)
	return (score >= sim_score, score)

# taking arguments
parser = argparse.ArgumentParser(description="use tree tip names to select sequences from fasta file")

parser.add_argument("-tree", "-t", required=True, type=str, help="input tree file")
parser.add_argument("-filetype", "-ft", type=str, default="newick", help="tree file type, defaults to newick")
parser.add_argument("-seqs", "-s", required=True, type=str, help="input fasta file to select sequences from")
parser.add_argument("-outfile", "-o", required=True, type=str, help="name of output fasta")
parser.add_argument("-p", type=float, help="enables partial matches and sets threshold for jaccardian similarity between headers and tips")

args = parser.parse_args()

# reading in tree and selecting tip names
tree = Phylo.read(args.tree, args.filetype) 
seq_names = set(tip.name for tip in tree.get_terminals())

# reading in fasta file and searching for matching headers
selected_seqs = []
for record in SeqIO.parse(args.seqs, "fasta"):
	header = record.id
	
	# looking for full matches between tree tips and headers
	if header in seq_names:
		selected_seqs.append(SeqRecord.SeqRecord(id=header, seq=record.seq))
		seq_names.remove(header)
		continue

	# if user specifies partial matches, uses jaccardian similarity with threshold args.p to determine which to add. Adds highest jaccardian score tip if multiple
	if args.p is not None:
		highest_score = None
		for seq in seq_names:
			higher, score = seq_sim_score(header, seq, args.p)
			if higher: 
				print(f"partial match found: {header}, {seq}, {score}")
	
				if highest_score is None:
					highest_score = (score, seq) 
					continue
				if highest_score[0] < score:
					highest_score = (score, seq) 
					continue
		# if at least 1 sequence beat threshold, add header and remove seqence from the search
		if highest_score is not None:
			print(f"using {header}, {seq}, as match")
			selected_seqs.append(SeqRecord.SeqRecord(id=header, seq=record.seq, description=""))
			seq_names.remove(highest_score[1])
			continue

# lists unmatched tips
if len(seq_names) > 0:
	print("no sequences found for tips: {0}".format(", ".join(seq_names)))
else:
	print("all sequences found!")

# write output file		
SeqIO.write(selected_seqs, args.outfile, "fasta")


