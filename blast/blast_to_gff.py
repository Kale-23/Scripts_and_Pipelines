#!/usr/bin/env python

import argparse

def blast_to_gff(blast_output_file, gff_output_file):
    with open(blast_output_file, 'r') as blast_file, open(gff_output_file, 'w') as gff_file:
        for line in blast_file:
            # Parse each line of the BLAST output
            fields = line.strip().split('\t')
            query_id = fields[0]
            subject_id = fields[1]
            alignment_length = int(fields[3])
            start = int(fields[8])
            end = int(fields[9])
            e_value = float(fields[10])
            bit_score = float(fields[11])
            strand = '+' if int(fields[8]) < int(fields[9]) else '-'

            # Write GFF entry
            gff_file.write('\t'.join([
                subject_id,  # Sequence ID
                'BLAST',  # Source
                'alignment',  # Feature type
                str(start),  # Start position
                str(end),  # End position
                str(bit_score),  # Score
                strand,  # Strand
                '.',  # Phase
                f'ID={query_id};Name={query_id};E-value={e_value};Alignment_length={alignment_length}'  # Attributes
            ]) + '\n')

def main():
    parser = argparse.ArgumentParser(description="Convert BLAST output to GFF format")
    parser.add_argument("blast_output", help="Path to the BLAST output file")
    parser.add_argument("gff_output", help="Path to the output GFF file")

    args = parser.parse_args()

    blast_to_gff(args.blast_output, args.gff_output)

if __name__ == "__main__":
    main()
