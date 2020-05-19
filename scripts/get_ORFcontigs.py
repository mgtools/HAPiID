"""
Script that will extract the DNA contigs of the ORFs with protein evidence and put them on one file
These DNA contigs will lated be used by the gene predictors to predict all the proteins in them
and thenmake a protein blast database so that I can blast the protential ORFs againast them
"""

from Bio import SeqIO
import json 
import os
import sys

seqs_f = '/data/mstambou/proteome_landscapes/HAPiID_results/verifyORFs/allExpressedORFsSoFar_c1.0.fasta'
contig2genome_dic_f = '/data/mstambou/proteome_landscapes/HAPiID/data/contigs2genomes.json'
genomes_dir = '/data/mstambou/proteome_landscapes/HAPiID/data/genomes/'
out_dir = '/data/mstambou/proteome_landscapes/HAPiID_results/verifyORFs/'

seqs_f, out_dir = sys.argv[1], sys.argv[2]

with open(contig2genome_dic_f, 'r') as in_f:
    contig2genome_dic = json.load(in_f)

orf_seqs = SeqIO.to_dict(SeqIO.parse(seqs_f, 'fasta'))
seqIDs = list(orf_seqs.keys())
        
contig_names = [item.split('|')[1].rsplit('_',1)[0] for item in seqIDs]
contig_names = list(set(contig_names))
genomes = [contig2genome_dic[item] for item in contig_names]

with open(out_dir + 'contigsWithPeptides.fasta', 'w') as out_f:
    for genome, contig in zip(genomes, contig_names):
        genome_seqs = SeqIO.to_dict(SeqIO.parse(genomes_dir + genome+'.fasta', 'fasta'))
        contig_seq = genome_seqs[contig]
        out_f.write('>'+contig+'\n')
        out_f.write(str(contig_seq.seq)+'\n')
