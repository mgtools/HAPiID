# HAPiID
HAPiID (pronunciation: Happy ID). High Abundance Protein guided Metaproteomics identification


An overview of the pipeline is summarized in this figure:

<p align="center">
  <img src="pipeline_schema.png" width="40%"/>
 </p>
 
 
 We develop a new approach that uses two steps to optimize  the  use  of  reference  genomes  as  the  universal  reference  for  human  gut metaproteomics  identification.   The  first  step  is  to  use  only  the  high  abundance  pro-teins (HAPs) (i.e., ribosomal proteins and elongation factors) for metaproteomic MS/MSdatabase search and derive the taxonomic composition of the microbiome based on theidentification results.  The second step is to expand the search database by including all proteins from identified abundant species.  We call our approach HAPiID (HAPs guidedmetaproteomics IDentification). While HAPiID was originally designed for human gut metaproteomics, we expand out pipeline and include scripts allowing users to precompile their own custom made protein database and perform spectral search using that.
 
 * [Using HAPiID with the precompiled human gut bacteria database please follow the steps here:](#using-hapiid-with-the-precompiled-human-gut-database)
* [Using HAPiID with user defined database starting from contigs please follow the steps here:](#using-hapiid-with-user-defined-database-starting-from-contigs)

# Using HAPiID with the precompiled human gut database

HAPiID is a reference based peptide identification pipeline. It uses 3,300 reference and metagenome binned genomes in order to profile and construct protein database for spectral search. To run this pipeline in one command we made the script [profileMPsample.sh](scripts/profileMPsample.sh).
This script takes 6 command line arguments:
```
-i specifies the mgf formatted spectra file as an input to the script

-o specifies the output directory for the script to store the mzid file after processing the mgf files using MSGF+

-d specifies the sequence file (in fasta format) for the initial profiling database (i.e. the ribosomal and elongation factors)

-t specifies the number of threads used by the script to run

-e specifies the directory for the new database files to be stored (the ones with the expanded genomes)

-n specifies the top N genomes, i.e. the first n most abundant species to be used to expand their genomes in the final database.
```


exampe to run script:
```
sh profileMPsample.sh -i /mgf/file -o /dir/to/output/mzid/fies/ -d /ribosomal and elongation factor/fasta/sequence/file -t number of cores -e /output/directory/to/store/database/files/ -n top_N_genomes
```

The pipeline has two main steps, the "profiling" step and "targetted search step". The program MSGF+ is used for peptide spectral matching, which could be found [here](MSGF+/).

The initial "profiling" step starts by issuing the an MSGF+ search over the database made from ribosomal proteins and elongation factors. This is done by the [MSGFPlus.jar](MSGF+/MSGFPlus.jar) program.

After this step the script [get_genome2peptideMappintgsMatrix.py](scripts/get_genome2peptideMappintgsMatrix.py) is used to map all the peptides identified through the ribosomal and elongation factor proteins, to their corresponding genomes.

The script [coverAllPeptide_greedy.py](scripts/coverAllPeptide_greedy.py) is then used to list genomes in decreasing order of their contribution needed to explain all the peptides identified in the "profiling" step.

The top N genomes are then extracted and their genomes expanded to create the final database for spectral identification using the scripts [extract_representative_genomes_abundanceBased.py](scripts/extract_representative_genomes_abundanceBased.py) & [extract_representative_ribP_elonF_proteins.py](scripts/extract_representative_ribP_elonF_proteins.py)


A final round of spectral search is performed using [MSGFPlus.jar](MSGF+/MSGFPlus.jar) over the final (expanded) database. 
Unique identified peptides reported by the script [get_uniquePeptides.py](scripts/get_uniquePeptides.py)


# Using HAPiID with user defined database starting from contigs

We have included helper scripts that will allow the user to compile their own database. The user has to provide their own contings (and/or complete reference genomes) in fasta format, and the scripts will extract the necesarry information from these contigs in order to compile a protein database and run HAPiID. To run HAPiID over a custom database, the user has to precompile their own protein database first. We have created the script [preprocessing.sh](scripts/preprocessing.sh) that precompiles a protein database in one command. The script takes three command line arguments:
```
-i specifies the directory where the contigs (in fasta format) are located

-t specifies the number of threads used by the script to run

-e specifies the extension of the contigs (i.e. ".fasta", ".fa" etc.) it is important that the user renames all the contigs to have the same extension.
```
example to run script

```
sh preprocessing.sh -i /dir/to/contigs/ -t n_threads -e .fasta
```

We start by predicting the protein coding genes using the script [runFGS_parallel.sh](scripts/runFGS_parallel.sh), which runs frag gene scan.

After this step we extract all the ribosomal and elongation factor protein sequences from each genome/bin, by first scanning all the predicted genes from the previous step against a precompiled database containing profiles for [ribosomal and elongation factor proteins](data/ribP_elonF_pfam_db_refined_manually/), using HMMSCAN through the script [runHMMSCAN_parallel.sh](scripts/runHMMSCAN_parallel.sh).

The final step is to de-duplicate these ribosomal and elongation factor genes using CD-hit at 100% sequence similarity and create a dictionary mapping between the protein sequence IDs to their bin IDs using the script [proteins2genomes.py](scripts/proteins2genomes.py).

After compiling a custom protein database the users can search for peptides using the script [profileMPsample_custom.sh](scripts/profileMPsample_custom.sh), which follows similar steps as the script [profileMPsample.sh](scripts/profileMPsample.sh) defined in the previous section, by using the same input arguments.

exampe to run script:
```
sh profileMPsample_custom.sh -i /mgf/file -o /dir/to/output/mzid/fies/ -d /ribosomal and elongation factor/fasta/sequence/file -t number of cores -e /output/directory/to/store/database/files/ -n top_N_genomes
```

