#! /bin/bash

# takes salmon quant file output and configures to work with RNAseq_plots
setup_library.py -d .

# runs RNAseq_plots with default values
Rscript ~/progs/RNAseq_plots.R ./libraries.tsv -all skin slime

Rscript ~/progs/complete_analysis.r . "$1".TRINITY.hit.fasta.transdecoder_fixed_headers.pep 9
