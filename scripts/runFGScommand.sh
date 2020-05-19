#!/bin/bash


while getopts i:o: option
do
case "${option}"
in
    i) file=${OPTARG};;
    o) out_dir=${OPTARG};;
esac
done

fname=$(echo $file | rev | cut -d'/' -f1 | rev)
echo $fname

run_FragGeneScan.pl -genome=$file -out=$out_dir$fname".FGS" -complete=1 -train=complete