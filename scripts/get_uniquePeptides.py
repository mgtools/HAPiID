# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:08:54 2019

@author: mstambou

script that will go over two FDR filtered TSVs, the first one being the TSV
for the HEG, and the second one being the TSV for the top N expanded genomes.

"""

import pandas as pd
import sys
import os

if len(sys.argv) != 4:
    print('please enter three command line argument to run this script, example to run python3 get_uniquePeptides.py /dir/to/MSGF+/output/ mgf_fname top_n')

else:
    out_dir = sys.argv[1]
    mgf_fname = sys.argv[2]
    top_n = sys.argv[3]

ribP_elonF_dir = out_dir+mgf_fname+'_ribP_elonF/'
extended_db_dir = out_dir+mgf_fname+'_extended_db_search/'

ribP_elonF_file = [item for item in os.listdir(ribP_elonF_dir) if item.endswith('.tsv.0.01.tsv')][0]
extended_db_file = [item for item in os.listdir(extended_db_dir) if item.endswith(top_n+'CoveringMostPeptides.tsv.0.01.tsv')][0]

ribP_elonF_df = pd.read_csv(ribP_elonF_dir + ribP_elonF_file, sep = '\t', header = None)
extended_db_df = pd.read_csv(extended_db_dir + extended_db_file, sep = '\t', header = None)

ribP_elonF_peptides = [''.join(filter(str.isalpha, item)) for item in  list(ribP_elonF_df[9])]
extended_db_peptides = [''.join(filter(str.isalpha, item)) for item in list(extended_db_df[9])]

ribP_elonF_peptides = list(set(ribP_elonF_peptides))
extended_db_peptides = list(set(extended_db_peptides))

all_peptides = list()
all_peptides.extend(ribP_elonF_peptides)
all_peptides.extend(extended_db_peptides)

all_peptides = list(set(all_peptides))

with open(out_dir + 'unique_peptides/'+mgf_fname+'_ribP_elonF_peptides.fasta', 'w') as out_f:
    for i, peptide in enumerate(ribP_elonF_peptides):
        out_f.write('>peptide'+str(i+1)+'\n')
        out_f.write(peptide+'\n')

with open(out_dir + 'unique_peptides/'+mgf_fname+'top'+str(top_n)+'_extended.fasta', 'w') as out_f:
    for i, peptide in enumerate(extended_db_peptides):
        out_f.write('>peptide'+str(i+1)+'\n')
        out_f.write(peptide+'\n')

with open(out_dir + 'unique_peptides/'+mgf_fname+'_all_peptides_top'+str(top_n)+'.fasta', 'w') as out_f:
    for i, peptide in enumerate(all_peptides):
        out_f.write('>peptide'+str(i+1)+'\n')
        out_f.write(peptide+'\n')

print('number of ribP elonF peptides identified: ', len(ribP_elonF_peptides))
print('number of peptides identified from extended genomes: ', len(extended_db_peptides))
print('number of peptides identified from combining both: ', len(all_peptides))
