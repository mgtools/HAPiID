# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 16:23:59 2019

@author: mstambou

script that will go over the FDR filtered file at 10% and will extract a list
of proteins in that file that have at least one spectral match and will
create a new fasta file that will include only these proteins with at least
one hit.
"""

import os
import pandas as pd
from Bio import SeqIO
import sys

if len(sys.argv) != 6:
    print('please enter 5 command line arguments, i.e. example to run\n python3 extract_proteins_with_hits.py output/dir/for/MSGF+_searchResults/ mgf_fname top_n dir/to/msgfplus_db/ FDR')

else:
    out_dir = sys.argv[1]
    mgf_fname = sys.argv[2]
    top_n = str(sys.argv[3])
    msgfplus_db = sys.argv[4]
    FDR = sys.argv[5]
    
    tsv_fname = [item for item in os.listdir(out_dir+'/'+mgf_fname+'_extended_db_search/') if item.endswith('top'+top_n+'CoveringMostPeptides.tsv.'+FDR+'.tsv')][0]
    print(tsv_fname)
    tsv_df = pd.read_csv(out_dir+'/'+mgf_fname+'_extended_db_search/'+tsv_fname, sep = '\t')
    
    protein_cols = list(tsv_df['Protein'])
    
    protein_ids = list()
    for protein_col in protein_cols:
        protein_col = protein_col.split(';')
        for protein in protein_col:
            protein_ids.append((protein.split('(pre=')[0]).replace('XXX_', ''))
            
    unique_protein_ids = list(set(protein_ids))
    
    seqs = SeqIO.to_dict(SeqIO.parse(msgfplus_db+mgf_fname+'/top_'+top_n+'_most_abundant.fasta', 'fasta'))
    with open(msgfplus_db+mgf_fname+'/top_'+top_n+'_proteins_with_hits_FDR_'+str(FDR)+'.fasta', 'w') as out_f:
        for seqID in unique_protein_ids:
            seq = seqs[seqID]
            out_f.write('>'+seq.description+'\n')
            out_f.write(str(seq.seq)+'\n')

    print('creating new database with '+str(len(unique_protein_ids)) + ' protein sequences from original '+str(len(seqs)) + ' sequences\n')
