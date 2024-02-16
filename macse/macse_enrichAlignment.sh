#!/bin/bash

#$file=${1%.*}

java -jar ~/scripts/macse_v2.07.jar -prog enrichAlignment -align $1 -seq $2 
