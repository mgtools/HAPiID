#!/bin/bash
#
#Script that will call and run hmmscan, part of HMMER tools, over a scpecied directory, typically where the FGS genes/proteins are predicted.
#The script takes four command line arguments.
#1st argument is the input directory, which is the location for protein sequences for genes in fasta format, each sample distinguished by sample name/id
#2nd argument is the number of threads, specifying the number of CPUs hmmscan would utilize.
#3rd argument is the output directory, hmmscan returns two outputs, first is a tabular format with columns specifying the accession numbers and for PFAM hmm and the ID for the fasta header of the sequence with the most siginifcant hit. The text output contains more detailed information about all the hits mapping to a profile from the list of protein sequences.
#4th argument specifies the location for the marker gene profiles. i.e. its the database which HMMSCAN searches these protein sequences against.
##Note the e-value is hardcoded here at 1e-10 feel free to change that or even make it a command line parameter.
#5th argumment is the extension of the files, i.e. tells the script what files to search for ending with what extensions.

while getopts i:t:o:m:e: option
do
    case "${option}"
	in
	    i) in_dir=${OPTARG};;
	    t) n_threads=${OPTARG};;
	    o) out_dir=${OPTARG};;
	    m) marker_profiles=${OPTARG};;
	    e) extension=${OPTARG};;
    esac
done

if [ ! -d $out_dir ]; then
    mkdir $out_dir
fi

files=( $in_dir*$extension)
n="${#files[@]}"

for y in $(eval echo "{0..$n}");
do
   printf "%s\0"  "${files[$y]}"; done| xargs -0 -I @@ -P $n_threads sh runHMMSCANcommand.sh -i @@ -o $out_dir -m $marker_profiles
