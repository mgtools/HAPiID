# HAPiID
HAPiID (pronunciation: Happy ID). High Abundance Protein guided Metaproteomics identification


An overview of the pipeline is summarized in this figure:

<p align="center">
  <img src="pipeline_schema.png" width="40%"/>
 </p>

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



