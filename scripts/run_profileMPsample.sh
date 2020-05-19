#!/bin/bash                                                                                                                                                                      
#                                                                                                                                                                                
#script where I will issue the HAPiID pipeline in parallel
#
#


while getopts i:o:d:e:c:f:p:x: option
do
    case "${option}"
        in
            i) in_dir=${OPTARG};; #the directory where the mgf files are present
            o) mzid_out_dir=${OPTARG};; #directory to dump the MSGF+ search results
            e) MSGFplus_db=${OPTARG};; #directory of where to save the protein database for search
	    d) ribP_elonF_db_file=${OPTARG};; #directory of the ribosomal proteins and elongation factors database
            c) n_threads=${OPTARG};; #number of threads to perform parallel runs
            f) in_f=${OPTARG};;  #the file containing the list of mgf files to run the pipeline on
	    p) percent_spectra=${OPTARG};; #percentage of spectra covered during profiling stage
	    x) protein_extension=${OPTARG};; #extensions of the protein sequences

    esac
done


IFS=$'\n' read -d '' -r -a files < $in_f

l="${#files[@]}"

echo $l
l="$(($l-1))"

for y in $(eval echo "{0..$l}");
do
    printf "%s\0"  "${files[$y]}"; done| xargs -0 -I @@ -P $n_threads sh profileMPsample.sh -i $in_dir@@ -o $mzid_out_dir -d $ribP_elonF_db_file -t 5 -e $MSGFplus_db -p $percent_spectra -x .faa
