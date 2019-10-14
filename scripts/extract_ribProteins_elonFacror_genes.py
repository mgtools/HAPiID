# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 18:48:06 2019

@author: mstambou
"""

import hmmscan_out_reader as hmm_r
import os
from Bio import SeqIO
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
import sys

if len(sys.argv) != 4:
    print('please enter 3 command line arguments: dir/to/hmmscan/ dir/to/fgs/out/ out/dir/for/results/')
else:
    hmm_scan_dir = sys.argv[1]
    fgs_out_dir = sys.argv[2]
    ribP_elonF_prot_out_dir = sys.argv[3]

plot = True

hmm_scan_files = [item for item in os.listdir(hmm_scan_dir) if item.endswith('.hmm')]

if os.path.isdir(ribP_elonF_prot_out_dir) == False:
    os.mkdir(ribP_elonF_prot_out_dir)

n_ribosomals, n_elongations = [], []
for hmm_scan_file in hmm_scan_files:
    n_ribosomal, n_elongation = 0, 0
    prefix = hmm_scan_file.split('.fa.FGS')[0]
    ntd_seq = fgs_out_dir+prefix+'.fa.FGS.faa'
    target_seqs_description_dic = dict()
    hmm_scan_df = hmm_r.pfam_2_query_2_description(hmm_scan_dir+hmm_scan_file)
    target_seq_ids = list(set(hmm_scan_df['query']))
    print('number of ribosomal proteins and elongation factors for genome: ', prefix,len(target_seq_ids))
    target_seqs_description_dic = {target_seq_id:', '.join(list(hmm_scan_df[hmm_scan_df['query'] == target_seq_id]['description'])) for target_seq_id in target_seq_ids}
    for item in target_seqs_description_dic.values():
        if 'elong' in item.lower():
            n_elongation += 1
        elif 'ribosom' in item.lower():
            n_ribosomal += 1
    all_seqs = list(SeqIO.parse(ntd_seq, 'fasta'))
    target_seqs = [seqs for seqs in all_seqs if seqs.id in target_seq_ids]
    with open(ribP_elonF_prot_out_dir+prefix+'_ribP_elonF.faa', 'w') as out_f:
        SeqIO.write(target_seqs, out_f, 'fasta')
    n_ribosomals.append(n_ribosomal)
    n_elongations.append(n_elongation)
    
frequencies_df = pd.DataFrame(columns = ['ribosomal_proteins', 'elongation_factors'] )
frequencies_df['ribosomal_proteins'] = n_ribosomals
frequencies_df['elongation_factors'] = n_elongations

if plot == True:
    plt.figure(figsize=(12,8))
    frequencies_df.boxplot()
    plt.title('Distribuiton of frequencies for ribosomal proteins, and elongation factors \nin '+str(len(hmm_scan_files))+' bacterial genomes from HBC dataset', fontsize=20)
    plt.ylabel('number of genes / genomes', fontsize = 18)
    plt.savefig('ribProteins_elongationFactors_counts.pdf')
