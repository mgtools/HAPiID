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


if len(sys.argv) != 5:
    print('please enter 4 command line arguments, i.e. example to run\n python3 search_abundance_results/t1d/131211-28623-12-ATH-F01-03.tsv.0.01_sortedAbundance.tsv 5 MSGFplus_db/131211-28623-12-ATH-F01-03/ proteomes_dir/')


representative_p_list_f  = sys.argv[1]
top_n_genomes = str(sys.argv[2])
out_dir = str(sys.argv[3])
representative_genomes_dir = out_dir+'top_'+top_n_genomes+'_most_abundant/'

proteomes_dir = str(sys.argv[4])

if os.path.isdir(representative_genomes_dir):
    shutil.rmtree(representative_genomes_dir)
os.makedirs(representative_genomes_dir)

representative_p_list_df = pd.read_csv(representative_p_list_f, sep = '\t')

genome_list = list(representative_p_list_df['genome'])[:int(top_n_genomes)]

proteomes = [item for item in os.listdir(proteomes_dir) if item.endswith('faa')]

for genome in genome_list:
    if 'UMGS' in genome:
        copyfile(proteomes_dir+genome+'.faa', representative_genomes_dir+genome+'.fa.FGS.faa')
    elif 'GCF_' in genome:
        fname = [item for item in os.listdir(proteomes_dir) if item.startswith(genome) and item.endswith('faa')]
        if len(fname) > 1:
            print('more than one: '+genome)
        else:
            copyfile(proteomes_dir+fname[0], representative_genomes_dir+fname[0])
    else:
        copyfile(proteomes_dir+genome+'.fa.FGS.faa', representative_genomes_dir+genome+'.fa.FGS.faa')


