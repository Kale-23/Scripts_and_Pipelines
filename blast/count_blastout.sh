#! /bin/bash

IFS=$'\n'

for item in blastout/*
do
 	counter=0
 	if [[ -e $item ]]
 	then
 		for line in $(cat $item)
 		do
 			counter=$(($counter + 1))
 		done
	file=$(basename $item)
	echo ${file%%.*}: $counter
 	fi
done	
