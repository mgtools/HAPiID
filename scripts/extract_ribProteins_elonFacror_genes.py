# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 18:48:06 2019

@author: mstambou
"""

import hmmscan_out_reader as hmm_r
import os
from Bio import SeqIO
import pandas as pd
import sys

if len(sys.argv) != 5:
    print('please enter 4 command line arguments: dir/to/hmmscan/ dir/to/fgs/out/ out/dir/for/results/ suffix')
else:
    hmm_scan_dir = sys.argv[1]
    fgs_out_dir = sys.argv[2]
    ribP_elonF_prot_out_dir = sys.argv[3]
    suffix = sys.argv[4]

hmm_scan_files = [item for item in os.listdir(hmm_scan_dir) if item.endswith('.hmm')]

if os.path.isdir(ribP_elonF_prot_out_dir) == False:
    os.mkdir(ribP_elonF_prot_out_dir)

for hmm_scan_file in hmm_scan_files:
    n_ribosomal, n_elongation = 0, 0
    prefix = hmm_scan_file.split(suffix)[0]
    ntd_seq = fgs_out_dir+prefix+suffix
    target_seqs_description_dic = dict()
    hmm_scan_df = hmm_r.pfam_2_query_2_description(hmm_scan_dir+hmm_scan_file)
    target_seq_ids = list(set(hmm_scan_df['query']))
    print('number of ribosomal proteins and elongation factors for genome: ', prefix,len(target_seq_ids))
    target_seqs_description_dic = {target_seq_id:', '.join(list(hmm_scan_df[hmm_scan_df['query'] == target_seq_id]['description'])) for target_seq_id in target_seq_ids}

    all_seqs = list(SeqIO.parse(ntd_seq, 'fasta'))
    target_seqs = [seqs for seqs in all_seqs if seqs.id in target_seq_ids]
    with open(ribP_elonF_prot_out_dir+prefix+'_ribP_elonF_proteins.faa', 'w') as out_f:
        SeqIO.write(target_seqs, out_f, 'fasta')

    
