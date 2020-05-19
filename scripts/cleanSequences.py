#!/bin/python3
"""
script that will remove premature stop codon from protein sequences that are to be used 
o be clustered by CD-hit
"""

import sys
import os
from Bio import SeqIO

proteins_dir = sys.argv[1]

for file in os.listdir(proteins_dir):
    seqs = list(SeqIO.parse(proteins_dir + file, 'fasta'))
    with open(proteins_dir + file, 'w') as out_f:
        for seq in seqs:
            out_f.write('>'+str(seq.description)+'\n')
            out_f.write(str(seq.seq).replace('*', '')+'\n')
