#!/bin/bash
#
#script in which I will put the different steps involved in metaprotemic profiling together. This first should involve
#starting with an .MGF sample and then searching the MGF peptides against the databse that contains the HEG genes
#i.e. ribosomal proteins and elongation factors. Then based on this it will call the greedy approach to basically
#keep the minimum amount of genomes needed to explain the peptides. And then rerun MSGF+ search over the expanded
#selected genomes.
#The steps involved in this pipeline, is first search the given mas spec data against the ribosomal proteins and elongation factor
#database coming from all the bins.
#then use the greedy approach to keep enough bins to explain all the peptides mapped to these HEG.
#select the top N bins covering most of the peptides and then expand their genomes by creating a new database
#search the mass spec data against the new database, and then report the number of peptides mapped.
#
#
#
START=$(date +%s)

MSGFPlus='../MSGF+/MSGFPlus.jar'
mod_f='../MSGF+/Mods_normal.txt'
n_cores=20
pyscripts='../MSGF+/'
#directory of all the proteomes needed to construct databases
proteomes_dir='../data/proteomes/'
#directory of the ribosomal proteins and elongation factor proteins for all genomes

ribP_elonF_proteins_dir='../data/ribP_elonF_proteins/'

ribP_elonF_db_file='../data/ribP_elonF_MSGF_search/ribP_elonF_all_cdHit_100_proteinSeqs.fasta'
ribP_elonF_clstr_file='../data/ribP_elonF_all_cdHit_100_proteinSeqs.fasta.clstr'
#MSGFplus_db='/data/mstambou/human_gut_microbiota/MSGFplus_db/'

proteins2genomes_fname='../data/proteins2genome_dic.json'

while getopts i:o:d:t:e:n: option
do
case "${option}"
in
    i) mgf_file=${OPTARG};;
    o) mzid_out_dir=${OPTARG};;
    d) ribP_elonF_db_file=${OPTARG};;
    t) n_cores=${OPTARG};;
    e) MSGFplus_db=${OPTARG};;
    n) n_topGenomes=${OPTARG};;
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

ribP_elonF_dir=$mzid_out_dir$mgf_fname"_ribP_elonF/"
extended_db_dir=$mzid_out_dir$mgf_fname"_extended_db_search/"

if [ ! -d $ribP_elonF_dir ]; then
    mkdir $ribP_elonF_dir
fi

if [ ! -d $extended_db_dir ]; then
    mkdir $extended_db_dir
fi

if [ ! -d $mzid_out_dir'data/'$mgf_file'/' ]; then
    mkdir -p $mzid_out_dir'data/'$mgf_fname'/'
fi

if [ ! -d $mzid_out_dir'unique_peptides/' ]; then
    mkdir -p $mzid_out_dir'unique_peptides/'
fi


#initial screening against the HEGs (should be run once only if using the same MGF file in the subsequent runs)

if [ ! -f $ribP_elonF_dir$mgf_fname'.GenomeCoveringAllHEGpeptides.txt' ]; then    

    java -Xmx32000M -jar $MSGFPlus -s $mgf_file  -o $ribP_elonF_dir$mgf_fname'.mzid' -d $ribP_elonF_db_file -inst 1 -t 15ppm -ti -1,2 -mod $mod_f -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1 -thread $n_cores

    java -Xmx32000M -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -i $ribP_elonF_dir$mgf_fname'.mzid' -showDecoy 1

    python3 parseFDR_o.py $ribP_elonF_dir$mgf_fname'.tsv' 0.01

    python3 get_genome2peptideMappintgsMatrix.py $ribP_elonF_dir$mgf_fname'.tsv.0.01.tsv' $ribP_elonF_clstr_file $proteins2genomes_fname $mzid_out_dir'data/'$mgf_fname'/'

    python3 coverAllPeptide_greedy.py $mzid_out_dir'data/'$mgf_fname'/'$mgf_fname'_genome2peptide_dic.json' $ribP_elonF_dir$mgf_fname'.GenomeCoveringAllHEGpeptides.txt'

fi

#search the MGF files with the expanded genomes using the top n genomes
python3 extract_representative_genomes_abundanceBased.py $ribP_elonF_dir$mgf_fname'.GenomeCoveringAllHEGpeptides.txt' $n_topGenomes $MSGFplus_db$mgf_fname'/' $proteomes_dir
#augment the ribosomal proteins and elongation factor proteins to the database
python3 extract_representative_ribP_elonF_proteins.py $ribP_elonF_dir$mgf_fname'.GenomeCoveringAllHEGpeptides.txt' $n_topGenomes $MSGFplus_db$mgf_fname'/' $ribP_elonF_proteins_dir

#delete the exisiting database index files if exists under same name, MSGF+ doesnt overwrite
if [ -f $MSGFplus_db$mgf_fname'/top_'$n_topGenomes'_most_abundant.revCat.canno' ]; then
    rm -r $MSGFplus_db$mgf_fname'/top_'$n_topGenomes'_most_abundant.revCat.'*
fi

java -Xmx32000M -jar $MSGFPlus -s $mgf_file -o $extended_db_dir$mgf_fname'_top'$n_topGenomes'CoveringMostPeptides.mzid' -d $MSGFplus_db$mgf_fname'/top_'$n_topGenomes'_most_abundant.fasta' -inst 1 -t 15ppm -ti -1,2 -mod $mod_f -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1 -thread $n_cores

java -Xmx32000M -cp $MSGFPlus edu.ucsd.msjava.ui.MzIDToTsv -i $extended_db_dir$mgf_fname'_top'$n_topGenomes'CoveringMostPeptides.mzid' -showDecoy 1

python2 parseFDR_o.py $extended_db_dir$mgf_fname'_top'$n_topGenomes'CoveringMostPeptides.tsv' 0.01

python3 get_uniquePeptides.py $mzid_out_dir $mgf_fname $n_topGenomes

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Exec time $DIFF seconds"
