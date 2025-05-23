#! /bin/bash 
for file in ~/Desktop/Hagfish_stuff/Salmon_only_annotated/*; do 
    thing=$(basename $file)
    output=~/Desktop/Hagfish_stuff/Salmon_only_annotated/${thing}
    #./rnadeap.R \
    #    -i $output \
    #    -f skin,slime -e \
    #    ~/Desktop/Hagfish_stuff/EnTAP/${thing}/final_results/entap_results.tsv \
    #    -x 10 
    mkdir -p ${output}/fastas
    python3 ./select_from_rnadeap.py \
        -f ${output}/RNAdeap_output_skin-slime/list/${thing}_vonWillebrand_slime.txt \
        -s ~/Desktop/Hagfish_stuff/final_fastas/${thing}_annotated.fasta \
        -o ${output}/fastas/${thing}_vonWillebrand_slime.fasta 
done
