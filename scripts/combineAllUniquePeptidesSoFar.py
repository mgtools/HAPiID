"""
script that will go over all the identified unmapped peptudes for each sample
and will create one fasta file combining all the unique peptide sequences from
al the metaproteome samples searched so far
"""

from Bio import SeqIO
import os
import sys

in_dir = sys.argv[1]
suffix = sys.argv[2]
out_fname = sys.argv[3]

files = [item for item in os.listdir(in_dir) if item.endswith(suffix)]


all_unique_peptides = list()

for _file in files:
    peptides = [str(item.seq) for item in list(SeqIO.parse(in_dir + _file, 'fasta'))]
    if not all_unique_peptides:
        all_unique_peptides.extend(peptides)
    else:
        new_peps = set(peptides).difference(set(all_unique_peptides))
        all_unique_peptides.extend(list(new_peps))
all_unique_peptides = list(set(all_unique_peptides))

with open(out_fname, 'w') as out_f:
    for i, peptide in enumerate(all_unique_peptides):
        out_f.write('>peptide'+str(i)+'\n')
        out_f.write(peptide+'\n')
    
