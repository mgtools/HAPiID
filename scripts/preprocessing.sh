#!/bin/bash
#
#if your starting with contigs and dont have a database yet there are some preprcessing steps that needs to be taken.
#First predict all the genes using frag gene scan
#script used to do that:
#- runFGS_parallel.sh
# sh runFGS_parallel.sh -i dereplicated_genomes/ -t 20 -o dereplicated_genomes_FGS/ -e .fa
#second, use the protein sequences from the protein coding genes and then scan them using hmmscan against the pfamDB and elonF DB that I manually curated and created.
#script used to do that:
#- runHMMSCAN_parallel.sh
# sh runHMMSCAN_parallel.sh -i dereplicated_genomes_FGS/ -t 20 -o dereplicated_genomes_FGS_HMMSCAN/ -m ribP_elonF_pfam_db_refined_manually/ribP_elonF_profiles_refined_manually.hmm -e .faa

#third you need to extract all the protein sequences for ribosomal proteins and elongation factors for all these genomes, to do  that I used the script
# - extract_ribProteins_elonFacror_genes.py
# - which requires the script hmmscan_out_reader.py
# python3 extract_ribProteins_elonFacror_genes.py dereplicated_genomes_FGS_HMMSCAN/ dereplicated_genomes_FGS/ dereplicated_genomes_ribP_elonF/

#fourth step is to CD-hit these sequences with 100% sequence identity to remove all the redundancy while doing a sample profiling step.
#before we issue the CD-hit command we first concatinated all the sequences together using the linux cat command:
#then I issued a cd-hit command over the concatinated sequence
#- cd-hit -i dereplicated_genomes_ribP_elonF_all.fasta -o dereplicated_genomes_ribP_elonF_all_cdhit_100.fasta -c 1.00 -n 5 -M 30000 -d 0 -T 8

#both the combined protein fasta file and the cd-hit .clstr output will be used by the profiling database.

#The last step is to create a dictionary mapping between the proteins to genome names IDs without the extensions:
#to do that I created the script proteins2genomes.py
#- python3  proteins2genomes.py all_proteomes/
#The following steps after this will be just to modify the pipeline slightly to be able to incorporate a general version of the task using the information just processed in these steps
START=$(date +%s)
while getopts i:t:e: option
do
case "${option}"
in
    i) genomes_dir=${OPTARG};;
    t) n_cores=${OPTARG};;
    e) contigs_extension=${OPTARG};;
esac
done


if [[ $genomes_dir = *"/" ]]; then
    FGS_out_dir="$(echo $genomes_dir | rev | cut -d'/' -f2- | rev).FGS/"
    HMM_out_dir="$(echo $FGS_out_dir | rev | cut -d'/' -f2- | rev).HMMScan/"
    echo $FGS_out_dir
    echo $HMM_out_dir
else
    FGS_out_dir=$genomes_dir".FGS/"
    HMM_out_dir="$(echo $FGS_out_dir | rev | cut -d'/' -f2- | rev).HMMScan/"
    echo $FGS_out_dir
    echo $HMM_out_dir
fi

if [ ! -d $FGS_out_dir ]; then
    mkdir $FGS_out_dir
fi

if [ ! -d $HMM_out_dir ]; then
    mkdir $HMM_out_dir
fi
#predict protein coding genes using frag gene scan
sh runFGS_parallel.sh -i $genomes_dir -t $n_cores -o $FGS_out_dir -e $contigs_extension

if [ ! -d '../data/' ]; then
    mkdir '../data/'
fi

if [ ! -d '../data/proteomes/' ]; then
    mkdir '../data/proteomes/'
    mkdir '../data/ribP_elonF_proteins/'
    mkdir '../data/ribP_elonF_MSGF_search/'
fi

#move all the proteomes to a proteomes folder under the data folder
mv $FGS_out_dir*.faa ../data/proteomes/
#scan all the proteomes to extract ribosomal and elongation factor proteins
sh runHMMSCAN_parallel.sh -i ../data/proteomes/ -t $n_cores -o $HMM_out_dir -m ../data/ribP_elonF_pfam_db_refined_manually/ribP_elonF_profiles_refined_manually.hmm -e .faa
#from the HMMScan results extract the protein sequences for the ribosomal proteins and elongation factor proteins
python3 extract_ribProteins_elonFacror_genes.py $HMM_out_dir ../data/proteomes/ ../data/ribP_elonF_proteins/
#combine all the ribosomal proteins and elongation factor proteins into one big multifasta file
cat ../data/ribP_elonF_proteins/*.faa > ../data/ribP_elonF_all.fasta
#we need to dereplicate the ribosomal proteins and elongation factors by issuing a CD-hit command at 100% sequence identity to avoid redundant search
cd-hit -i ../data/ribP_elonF_all.fasta -o ../data/ribP_elonF_all_cdHit_100_proteinSeqs.fasta -c 1.00 -n 5 -M 30000 -d 0 -T $n_cores
mv ../data/ribP_elonF_all_cdHit_100_proteinSeqs.fasta ../data/ribP_elonF_MSGF_search/
#create a dictionary mapping between the protein sequence IDs to their corresponding genome ID names without extensions
python3 proteins2genomes.py ../data/proteomes/ ../data/


END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Exec time $DIFF seconds"
