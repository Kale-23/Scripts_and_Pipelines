#!/bin/bash

java -jar ~/scripts/macse_v2.07.jar -prog alignTwoProfiles -p1 $1 -p2 $2 -out_NT combined_alignment_NT.fasta -out_AA combined_alignment_AA.fasta 
