#!/bin/bash

#Script that will automate frag gene scan over the different metagenomes at different time points. It takes 3 arguments to run
#1st argument is the input directory for the genomes or contigs depending
#2nd argument is the number of threads
#3rd argument is the output directory where it outputs the genes
#example run: sh runFGS.sh -i megahit -t 9 -o FGS_contigs
#at the end of processing each sample this script also calls another python script called prepend_sample2FGS_headers.py
#this script will modify the fasta headers of the FGS gene predictions and will prepend the name/id of the sample that it comes from
#this will later be usefull to keep track of the predicted genes whenever I merge all the predicted genes from the different samples

while getopts i:t:o:e: option
do
case "${option}"
in
    i) in_dir=${OPTARG};;
    t) n_threads=${OPTARG};;
    o) out_dir=${OPTARG};;
    e) extension=${OPTARG};;
esac
done

if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

for file in $(ls $in_dir*$extension)
do
        echo $file
	echo $extension
	#fname=$(echo $file | cut -d'/' -f2)
	fname=$(echo $file | rev | cut -d'/' -f1 | rev)
	echo $fname
	echo $out_dir$fname".FGS"
	run_FragGeneScan.pl -genome=$file -out=$out_dir$fname".FGS" -complete=1 -train=complete -thread=$n_threads
done

