#!/usr/bin/env python3
import argparse

def remove_gaps_and_exclamations(input_file, output_file):
    sequences = {}
    current_sequence = None
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_sequence = line[1:]
                sequences[current_sequence] = "" 
            else:
                if current_sequence is not None:
                    sequences[current_sequence] += line.replace("-", "").replace("!", "")

    cleaned_fasta = ""
    for header, sequence in sequences.items():
        cleaned_fasta += f">{header}\n{sequence}\n"

    with open(output_file, "w") as f:
        f.write(cleaned_fasta)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove '-' and '!' characters from sequences in a FASTA file.")
    parser.add_argument("-i", required=True, help="Path to the input FASTA file")
    parser.add_argument("-o", required=True, help="Path to the output FASTA file")
    args = parser.parse_args()

    remove_gaps_and_exclamations(args.i, args.o)

