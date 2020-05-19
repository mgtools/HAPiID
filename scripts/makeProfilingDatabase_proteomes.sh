#!/bin/bash
#
#script that will be run before running the pipeline that will generate all 
#the necesarry files for the pipeline to run
#

START=$(date +%s)


n_cores=20
#directory of all the proteomes needed to construct databases                                                        
proteomes_dir='../data/proteomes/'

while getopts i:t:e: option
do
case "${option}"
in
    i) proteomes_dir=${OPTARG};;
    t) n_cores=${OPTARG};;
    e) proteomes_extension=${OPTARG};;
esac
done
 
#directory of the ribosomal proteins and elongation factor proteins for all genomes
ribP_elonF_proteins_dir='../data/ribP_elonF_proteins/'

data_dir='../data/'
cd_hit_out_f=$data_dir'ribP_elonF_all_cdHit_100_proteinSeqs.fasta'
ribP_elonF_db_file='../data/ribP_elonF_MSGF_search/ribP_elonF_all_cdHit_100_proteinSeqs.fasta'
ribP_elonF_clstr_file=$cd_hit_out_f".clstr"

ribP_elonF_proteins2genome_fname='../data/ribP_elonF_proteins2genome_dic.json'
proteins2genomes_fname='../data/proteins2genome_dic.json'
proteomes_suffix='.fasta.FGS.faa'
ribP_elonF_suffix='_ribP_elonF_proteins.faa'
ribP_elonF_proteins_dir="../data/ribP_elonF_proteins/"

if [[ $proteomes_dir = *"/" ]]; then
    HMM_out_dir="$(echo $proteomes_dir | rev | cut -d'/' -f2- | rev).HMMScan/"
    echo $HMM_out_dir
else
    HMM_out_dir="$(echo $proteomes_dir | rev | cut -d'/' -f2- | rev).HMMScan/"
    echo $HMM_out_dir
fi

if [ -f "$ribP_elonF_clstr_file" ]; then
    rm -r $ribP_elonF_clstr_file
fi

if [ -f "$proteins2genomes_fname" ]; then
    rm -r $proteins2genomes_fname
fi

if [ -f "$ribP_elonF_proteins2genome_fname" ]; then
    rm -r $ribP_elonF_proteins2genome_fname
fi

if [ -f "$data_dir""ribP_elonF_all.fasta" ]; then
    rm -r $data_dir"ribP_elonF_all.fasta"
fi

if [ -d "$ribP_elonF_proteins_dir" ]; then
    rm -rf $ribP_elonF_proteins_dir
fi

if [ -d "$HMM_out_dir" ]; then
    rm -rf $HMM_out_dir
fi

if [ ! -d "$ribP_elonF_proteins_dir" ]; then
    mkdir $ribP_elonF_proteins_dir
fi

if [ -d "$data_dir""ribP_elonF_MSGF_search/" ] ; then
    rm -r "$data_dir""ribP_elonF_MSGF_search/"    
fi
mkdir $data_dir"ribP_elonF_MSGF_search/"

#unify proteomes extension
python3 unifyGenomeExtensions.py $proteomes_dir $proteomes_extension

#scan all the proteomes to extract ribosomal and elongation factor proteins                                 
sh runHMMSCAN_parallel.sh -i ../data/proteomes/ -t $n_cores -o $HMM_out_dir -m ../data/ribP_elonF_pfam_db_refined_manually/ribP_elonF_profiles_refined_manually.hmm -e .faa
#from the HMMScan results extract the protein sequences for the ribosomal proteins and elongation factor proteins  
python3 extract_ribProteins_elonFacror_genes.py $HMM_out_dir $proteomes_dir $ribP_elonF_proteins_dir $proteomes_extension

echo "creating proteome 2 genome map file"
python3 proteins2genomes.py $proteomes_dir $proteins2genomes_fname $proteomes_extension

echo "creating ribosomal and elongation factor proteome 2 genome map file"
python3 proteins2genomes.py $ribP_elonF_proteins_dir $ribP_elonF_proteins2genome_fname $ribP_elonF_suffix

#clean protein sequences
python3 cleanSequences.py $ribP_elonF_proteins_dir

#append a new line after each file
for file in $ribP_elonF_proteins_dir*$ribP_elonF_suffix
do
    echo "" >> "$file"
done

cat $ribP_elonF_proteins_dir*$ribP_elonF_suffix > $data_dir"ribP_elonF_all.fasta"

cd-hit -i $data_dir"ribP_elonF_all.fasta" -o $data_dir"ribP_elonF_all_cdHit_100_proteinSeqs.fasta" -c 1.0 -n 5 -T 30 -d 200 -M 50000

mv $data_dir"ribP_elonF_all_cdHit_100_proteinSeqs.fasta" $data_dir"ribP_elonF_MSGF_search/"

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Exec time $DIFF seconds"
