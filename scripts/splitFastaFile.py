"""
script that will split one big multi fasta sequence flile into n small files so that I can multiprocess 
them in parallel.
"""

from Bio import SeqIO
import os
import sys
import numpy as np
import sys

#in_f = '/data/mstambou/proteome_landscapes/small_proteins/small_proteins.faa'
#out_dir = '/data/mstambou/proteome_landscapes/small_proteins/small_proteins_split/'
#n_partitions = 20

in_f = sys.argv[1]
out_dir = sys.argv[2]
n_partitions = int(sys.argv[3])

f_name = in_f.rsplit('/', 1)[-1].split('.')[0]

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
    
protein_seqs = list(SeqIO.parse(in_f, 'fasta'))

indexes = np.arange(0, len(protein_seqs), len(protein_seqs)//n_partitions)
indexes[-1] = len(protein_seqs)

for i in range(len(indexes)-1):
    _from = indexes[i]
    to = indexes[i+1]
    with open(out_dir + f_name + '_partition'+str(i)+'.fasta', 'w') as out_f:
        out_f.write(''.join(['>'+seq.id+'\n'+str(seq.seq)+'\n' for seq in protein_seqs[_from: to]   ]))
