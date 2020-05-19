# -*- coding: utf-8 -*-
"""
script where I will use to extract the 6 frame translations of the contigs for the top N species 
suggested by HAPiID and then will form a DB using these 6 fram translations of the participating
bins and will look for peptides and to see if we miss or not miss any gene predictions that are 
actually expressed. Possibly by identifying novel peptides.
"""

from Bio import SeqIO
import os
import sys
import pandas as pd
from Bio import SeqUtils
import fileinput
import glob

#sample_name = 'HM403'
#genomes_dir = '/data/mstambou/proteome_landscapes/HAPiID/data/genomes/'
#genome_ranks_f = '/data/mstambou/proteome_landscapes/HAPiID_results/HM403_ribP_elonF/HM403.GenomeCoveringAllHEGpeptides.txt'
#top_n = 10
#db_out_dir = '/data/mstambou/proteome_landscapes/HAPiID_db/'+sample_name+ 'top'+str(top_n)+'_6frameTranslation/'

sample_name = sys.argv[1]
genomes_dir = sys.argv[2]
genome_ranks_f = sys.argv[3]
top_n = int(sys.argv[4])
db_out_dir = sys.argv[5]

if not os.path.isdir(db_out_dir):
    os.mkdir(db_out_dir)

translated_genomes_dir = db_out_dir + 'top_'+str(top_n)+'_most_abundant_6FrameTranslated/'

if not os.path.isdir(translated_genomes_dir):
    os.mkdir(translated_genomes_dir)

genome_ranks_df = pd.read_csv(genome_ranks_f, sep = '\t')

top_n_genomes = list(genome_ranks_df['genome'])[:top_n]

gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

basepairs = {'A':'T', 'C':'G', 'G':'C', 'T':'A'} #defined for complimentarity

def translate_frameshifted( sequence ):
      translate = ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])
      #// to return int after division
      #return 'X' if key value pair not found
      return translate

def reverse_complement( sequence ):
      reversed_sequence = (sequence[::-1])
      # This ::-1 syntax in python means reverse the string.
      rc = ''.join([basepairs.get(reversed_sequence[i], 'X') for i in range(len(sequence))])
      return rc
    
def six_frame_translation(seq):
    f1 = translate_frameshifted(seq[0:]) #first frame
    f2 = translate_frameshifted(seq[1:]) #second frame
    f3 = translate_frameshifted(seq[2:]) #third frame
    f4 = translate_frameshifted(reverse_complement(seq)) #first frame negative strand
    f5 = translate_frameshifted(reverse_complement(seq[:len(seq)-1])) #second frame negative strand
    f6 = translate_frameshifted(reverse_complement(seq[:len(seq)-2])) #third frame negative strand
    
    return f1, f2, f3, f4, f5, f6


for genome in top_n_genomes:
    genome_seqs = list(SeqIO.parse(genomes_dir+genome+'.fasta', 'fasta'))
    with open(translated_genomes_dir + genome + '_6frameTranslated.fasta', 'w') as out_f:
        for seq in genome_seqs:
            _id = seq.id
            sequence = str(seq.seq)
            six_frame_translations = six_frame_translation(sequence)
            out_f.write('>'+_id+'_f1\n')
            out_f.write(six_frame_translations[0]+'\n')
            out_f.write('>'+_id+'_f2\n')
            out_f.write(six_frame_translations[1]+'\n')
            out_f.write('>'+_id+'_f3\n')
            out_f.write(six_frame_translations[2]+'\n')
            out_f.write('>'+_id+'_f4\n')
            out_f.write(six_frame_translations[3]+'\n')
            out_f.write('>'+_id+'_f5\n')
            out_f.write(six_frame_translations[4]+'\n')
            out_f.write('>'+_id+'_f6\n')
            out_f.write(six_frame_translations[5]+'\n')

file_list = glob.glob(translated_genomes_dir + '*')
with open('/'.join(translated_genomes_dir.split('/')[:-1]) + '.fasta', 'w') as file:
    input_lines = fileinput.input(file_list)
    file.writelines(input_lines)
