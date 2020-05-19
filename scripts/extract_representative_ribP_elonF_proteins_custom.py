# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 17:02:08 2019

@author: mstambou
"""

import os
import pandas as pd
from shutil import copyfile
import shutil
import fileinput
import glob
import sys

if len(sys.argv) != 5:
    print('please enter 4 command line arguments, i.e. example to run\n python3 extract_representative_ribP_elonF_proteins.py search_abundance_results/t1d/131211-28623-12-ATH-F01-03.tsv.0.01_sortedAbundance.tsv 5 MSGFplus_db/131211-28623-12-ATH-F01-03/ ribP_elonF_dir/')


representative_p_list_f  = sys.argv[1]
top_n_genomes = str(sys.argv[2])
out_dir = str(sys.argv[3])

representative_p_list_df = pd.read_csv(representative_p_list_f, sep = '\t')

representative_ribP_elonF_dir = out_dir+'top_'+top_n_genomes+'_most_abundant/'

ribP_elonF_dir = str(sys.argv[4])


remaining_HEG_list = list(representative_p_list_df['genome'])[int(top_n_genomes):]
print(ribP_elonF_dir)
print(representative_ribP_elonF_dir)
for genome in remaining_HEG_list:
    copyfile(ribP_elonF_dir+genome+'_marker_gene.faa', representative_ribP_elonF_dir + genome+'_marker_gene.faa')

file_list = glob.glob(representative_ribP_elonF_dir+"*")
with open('/'.join(representative_ribP_elonF_dir.split('/')[:-1])+'.fasta', 'w') as file:
    print('IM CONCATINATING THE FASTA FILES HERE!!!!')
    input_lines = fileinput.input(file_list)
    file.writelines(input_lines)
