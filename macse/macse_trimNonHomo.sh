#!/bin/bash

file=${1%.*}

java -jar ~/scripts/macse_v2.07.jar -prog trimNonHomologousFragments -seq $1 -out_trim_info "$file"_stats.csv -min_trim_in 40 -min_trim_ext 20
