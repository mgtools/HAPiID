#!/bin/bash
#
#script where I will issue the search6Frames.sh pipeline in parallel
#


while getopts i:o:e:c:f:n:h: option
do
    case "${option}"
	in
	    i) in_dir=${OPTARG};;
	    o) mzid_out_dir=${OPTARG};;
	    e) MSGFplus_db=${OPTARG};;
	    c) n_threads=${OPTARG};;
	    f) in_f=${OPTARG};;
	    n) n_topGenomes=${OPTARG};;
	    h) HAPiID_n_topGenomes=${OPTARG};;
    esac
done


IFS=$'\n' read -d '' -r -a files < $in_f

l="${#files[@]}"

echo $l
l="$(($l-1))"

for y in $(eval echo "{0..$l}");
do
    printf "%s\0"  "${files[$y]}"; done| xargs -0 -I @@ -P $n_threads sh search6Frames.sh -i $in_dir@@ -o $mzid_out_dir -e $MSGFplus_db -t 5 -e $MSGFplus_db -n $n_topGenomes -h $HAPiID_n_topGenomes

