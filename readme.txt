read me file for reference based spectral search pipeline

the folder under metaproteomics should contain 4 folders:

1) data
2) MSGF+
3) MSGF+_latest
4) pyscript

the folder pyscript contains all the necesarry scripts for the pipeline. The 
script that puts together the steps for the pipeline is called 'profileMPsample.sh', which is a bash script and can be run by the sh profileMPsample.sh command

This script requires 6 command line arguments to run, which are the following:

-i mgf_file : which is the mgf spectral file as input that we are going to search against
-o mzid_out_dir : is the output directory where the search results by MSGF+ will be stored
-d ribP_elonF_db_file : is the database created using the ribosomal and elongation factor proteins which is found under /Metaproteomics/data/ribP_elonF_MSGF_search/ribP_elonF_all_cdHit_100_proteinSeqs.fasta
-t n_cores : specifies the number of threads for the script to use
-e MSGFplus_db : specifies the output folder where the new database (using the expanded genomes) will be created
-n n_topGenomes : the top n genomes to be selected to do the spectral search over using the expanded genomes

example to run this script:

sh profileMPsample.sh -i path/to/mgf_file -o directory/to/store/output/ -d ../data/ribP_elonF_MSGF_search/ribP_elonF_all_cdHit_100_proteinSeqs.fasta -t 20 -e dir/to/store/exapanded/genomes/database/ -n 5

PEASE make sure that when moving this pipeline to another folder or to another server, all you have to do is to copy paste the
main folder containing all the four folders inside them to another server. 
The main folder being Metaproteomics/ folder, without changing the internal folder organizations.
