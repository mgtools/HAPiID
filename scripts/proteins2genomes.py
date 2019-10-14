# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 11:54:42 2018

@author: mstambou

script that will go over all the protein sequences in folder and create a dictionary mapping between the protein sequence ID to its corresponding genome ID without the extension
"""

import os
import json
from Bio import SeqIO
import sys

proteome_dir  = sys.argv[1]
out_dir = sys.argv[2]

proteins2genomes_dic = dict()
counter = 0
for file in [item for item in os.listdir(proteome_dir) if item.endswith('.faa')]:
    f_id = file.split('.fa.FGS.faa')[0]
    seqs = list(SeqIO.parse(proteome_dir + file, 'fasta'))
    new_dic = {str(k.id):f_id for k in seqs}
    counter += len(seqs)
    proteins2genomes_dic.update(new_dic)

with open(str(out_dir)+'proteins2genome_dic.json', 'w') as out_f:
    json.dump(proteins2genomes_dic, out_f)
    
    
print('total number of proteins: ', str(counter))
