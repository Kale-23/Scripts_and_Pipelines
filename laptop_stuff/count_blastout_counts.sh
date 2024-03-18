#!/bin/bash

cat $1 | cut -f1 | sort | uniq -c | awk '{$1=$1};1' | cut -d" " -f1 | sort -n | uniq -c
