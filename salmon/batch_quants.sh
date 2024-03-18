#! /bin/bash

for dir in $1/*; do
	name=$(basename ${dir%_*})
	sbatch ./salmon_gen_quants.sh indexes/"$name"_index $dir quants/
done
