# -*- coding: utf-8 -*-                                                                                                                                                          
"""
script that will get the unique peptides from the MSGF+ search over all the six frame translation databsae
"""

import pandas as pd
import sys
import os


out_dir = sys.argv[1]
mgf_fname = sys.argv[2]
top_n = sys.argv[3]

msgf_out_dir = out_dir + mgf_fname+"_top"+str(top_n)+"_6frameTranslation_search/"

out_file = [item for item in os.listdir(msgf_out_dir) if item.endswith('.tsv.0.01.tsv')][0]

out_df = pd.read_csv(msgf_out_dir + out_file, sep = '\t', header = None)

peptides = [''.join(filter(str.isalpha, item)) for item in  list(out_df[9])]

peptides = list(set(peptides))

print(len(peptides), ' unique peptides identified')

with open(out_dir + 'unique_peptides/' + mgf_fname + '_top'+str(top_n)+'_6frameTranslationPeptides.fasta', 'w') as out_f:
    for i, peptide in enumerate(peptides):
        out_f.write('>peptide'+str(i+1)+'\n')
        out_f.write(peptide+'\n')
