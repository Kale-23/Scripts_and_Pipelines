#! /bin/bash

# check arguments
if [ "$#" -ne 2 ]; then
    echo "first arg: source directory, second arg: target directory"
    exit 1
fi

# setting source and target directories
source_dir="$1"
target_dir="$2"

# for every directory in source directory, make new directory in target, and place all fasta files within
for dir in "$source_dir"/*
do
    if [[ -d $dir ]]
    then
        new_dir="$target_dir/$(basename $dir)"

        mkdir -p "$new_dir"
        
        cp "$dir"/*.fasta "$new_dir"

    fi
done
