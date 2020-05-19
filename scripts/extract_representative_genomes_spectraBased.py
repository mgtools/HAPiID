# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 17:02:08 2019

@author: mstambou

Script where it takes two things as inputs. The first input file is a TSV file which is resulted 
from an abundance qunatification results, which quantifies genomes based on the abundances of marker 
genes, the second input should be the cutoff number to include genomes. i.e. Top 50 most abundant genomes
for instance
"""

import os
import pandas as pd
from shutil import copyfile
import shutil
import fileinput
import glob
import sys


if len(sys.argv) != 6:
    print('please enter 4 command line arguments, i.e. example to run\n python3 search_abundance_results/t1d/131211-28623-12-ATH-F01-03.tsv.0.01_sortedAbundance.tsv 5 MSGFplus_db/131211-28623-12-ATH-F01-03/ proteomes_dir/ protein_extension')


representative_p_list_f  = sys.argv[1]
percent_spectra = str(sys.argv[2])
out_dir = str(sys.argv[3])
representative_genomes_dir = out_dir+'covering_'+percent_spectra+'_percent_ribP_elonF_spectra/'

proteomes_dir = str(sys.argv[4])

prot_ext = str(sys.argv[5])

if os.path.isdir(representative_genomes_dir):
    shutil.rmtree(representative_genomes_dir)
os.makedirs(representative_genomes_dir)

representative_p_list_df = pd.read_csv(representative_p_list_f, sep = '\t')

genome_list = list()
for i, row in representative_p_list_df.iterrows():
    genome_list.append(row['genome'])
    if row['nSpectra_%'] > float(percent_spectra):
        break

proteomes = [item for item in os.listdir(proteomes_dir) if item.endswith(prot_ext)]

for genome in genome_list:
    copyfile(proteomes_dir + genome + prot_ext, representative_genomes_dir + genome + prot_ext)
