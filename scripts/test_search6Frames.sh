#!/bin/bash
#
#script where I will first get the six frame translations of the contigs for the n most abundant genomes ranked by HAPiID
#and then create a protein database using these 6 frame translations for the top n genomes and search for peptides
#within this newly created translated contig database
#

START=$(date +%s)

MSGFPlus='../MSGF+/MSGFPlus.jar'
mod_f='../MSGF+/Mods_normal.txt'
pyscripts='../MSGF+/'
#directory of all the proteomes needed to construct databases                                                                                                                    
genomes_dir='../data/genomes/'
#directory of the ribosomal proteins and elongation factor proteins for all genomes
ribP_elonF_proteins_dir='../data/ribP_elonF_proteins/'

ribP_elonF_db_file='../data/ribP_elonF_MSGF_search/ribP_elonF_all_cdHit_100_proteinSeqs.fasta'
ribP_elonF_clstr_file='../data/ribP_elonF_all_cdHit_100_proteinSeqs.fasta.clstr'


while getopts i:o:e:t:n:h: option
do
case "${option}"
in
    i) mgf_file=${OPTARG};;
    o) mzid_out_dir=${OPTARG};;
    e) MSGFplus_db=${OPTARG};;
    t) n_threads=${OPTARG};;
    n) n_topGenomes=${OPTARG};;
    h) HAPiID_n_topGenomes=${OPTARG};;

esac
done

mgf_fname="$(echo $mgf_file | awk -F/ '{print $NF}')"
mgf_fname="$(echo $mgf_fname | rev | cut -d'.' -f2- | rev)"


genomes_rank_f=$mzid_out_dir$mgf_fname"_ribP_elonF/"$mgf_fname".GenomeCoveringAllHEGpeptides.txt"

echo $mgf_file
echo $mgf_fname
echo $genomes_rank_f
echo $n_threads
echo $n_topGenomes
echo $HAPiID_n_topGenomes


END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Exec time $DIFF seconds"
