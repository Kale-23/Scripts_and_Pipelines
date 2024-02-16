#! /usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description = "removes ! from alignment files")

parser.add_argument("infile", type=str, help="input file")
parser.add_argument("outfile", type=str, help="output file")
parser.add_argument("replace", type=str, help="character to replace with \"\"")

args = parser.parse_args()

with open(args.infile, "r") as infile, open(args.outfile, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            # Write sequence headers as they are
            outfile.write(line)
        else:
            # Remove "!" characters from sequence lines
            sequence = line.replace("*", "")   #args.replace, "")
            outfile.write(sequence)
