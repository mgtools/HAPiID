#!/bin/bash

while getopts i:o:m: option
do
    case "${option}"
	in
	    i) in_f=${OPTARG};;
	    o) out_dir=${OPTARG};;
	    m) marker_profiles=${OPTARG};;
    esac
done

fname=$(echo $in_f | rev | cut -d'/' -f1 | rev)
hmmscan -o $out_dir$fname".txt" --tblout $out_dir$fname".hmm" -E 1e-10 $marker_profiles $in_f
