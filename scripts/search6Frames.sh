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


while getopts i:o:t:e:n:h: option
do
case "${option}"
in
    i) mgf_file=${OPTARG};;
    o) mzid_out_dir=${OPTARG};;
    t) n_cores=${OPTARG};;
    e) MSGFplus_db=${OPTARG};;
    n) n_topGenomes=${OPTARG};;
    h) HAPiID_n_topGenomes=${OPTARG};;
esac
done

if [ ! -d $mzid_out_dir ]; then
    mkdir $mzid_out_dir
fi

if [ ! -d $MSGFplus_db ]; then
    mkdir $MSGFplus_db
fi

mgf_fname="$(echo $mgf_file | awk -F/ '{print $NF}')"
mgf_fname="$(echo $mgf_fname | rev | cut -d'.' -f2- | rev)"

echo $mgf_fname

six_frame_db_dir=$MSGFplus_db$mgf_fname"_top"$n_topGenomes"_6frameTranslation/"

if [ ! -d $six_frame_db_dir ]; then
    mkdir $six_frame_db_dir
fi

MSGF_out_dir=$mzid_out_dir$mgf_fname"_top"$n_topGenomes"_6frameTranslation_search/"

if [ ! -d $MSGF_out_dir ]; then
    mkdir $MSGF_out_dir
fi

peptideMatches_out_dir=$mzid_out_dir"peptideMatches/"

if [ ! -d $peptideMatches_out_dir ]; then
    mkdir $peptideMatches_out_dir
fi

expressedORFs_out_dir=$mzid_out_dir"expressedORFs/"

if [ ! -d $expressedORFs_out_dir ]; then
    mkdir $expressedORFs_out_dir
fi

genomes_rank_f=$mzid_out_dir$mgf_fname"_ribP_elonF/"$mgf_fname".GenomeCoveringAllHEGpeptides.txt"


#create the six frame translation protein datbasae
python3 make_6_frame_db.py $mgf_fname $genomes_dir $genomes_rank_f $n_topGenomes $six_frame_db_dir

java -Xmx32000M -jar $MSGFPlus -s $mgf_file  -o $MSGF_out_dir$mgf_fname".mzid" -d $six_frame_db_dir"top_"$n_topGenomes"_most_abundant_6FrameTranslated.fasta" -inst 1 -t 15ppm -ti -1,2 -mod $mod_f -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1 -thread $n_cores

java -Xmx32000M -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -i $MSGF_out_dir$mgf_fname".mzid" -showDecoy 1

python3 parseFDR_o.py $MSGF_out_dir$mgf_fname".tsv" 0.01

python3 get_uniquePeptides_6frames.py $mzid_out_dir $mgf_fname $n_topGenomes

python3 getUnmappedPeptides.py $mzid_out_dir $mgf_fname "top"$HAPiID_n_topGenomes "top"$n_topGenomes $n_cores

python3 get_expressed_ORFs.py $mzid_out_dir $MSGFplus_db $mgf_fname "top"$HAPiID_n_topGenomes "top"$n_topGenomes


END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Exec time $DIFF seconds"
