"""
script that will take the table produced by the blast out processing script and 
will go through the hits that have 100% PID with a gene from reference sequences
And will instead get the genes annotaed from the reference sequences instead of 
the ones predicted using the three predictors here (FGS, prodigal, genemark)
"""
import pandas as pd
import json
from Bio import SeqIO
import os
import sys

#blastSummary_f = '/data/mstambou/proteome_landscapes/HAPiID_results/verifyORFs/blast_out/blastHitsSummary.tsv'
#contig2genome_dic_f = '/data/mstambou/proteome_landscapes/HAPiID/data/contigs2genomes.json'
#proteomes_dir = '/data/mstambou/proteome_landscapes/HAPiID/data/proteomes/'

#out_dir = '/data/mstambou/proteome_landscapes/HAPiID_results/verifyORFs/refseq_out/'

blastSummary_f, contig2genome_dic_f, proteomes_dir, out_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

contigs_seqs = SeqIO.to_dict(SeqIO.parse(out_dir + 'allExpressedORFsSoFar_c1.0.fasta', 'fasta'))

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

with open(contig2genome_dic_f) as in_f:
    contig2genome_dic = json.load(in_f)

blastSummary_df = pd.read_csv(blastSummary_f, sep = '\t', header = [0, 1])
hits100 = blastSummary_df[blastSummary_df['FGS']['% identity'] == 100]

refseq_out_dir = out_dir+'refseq_out/'

with open(refseq_out_dir + 'contigsWithPeptides.fasta.refseq.faa', 'w') as out_f, open(out_dir+'ORFsWithPeptidesWith100%hit.fasta', 'w') as out_f2:
    for i, row in hits100.iterrows():
        contig_id = row['FGS']['subject acc.ver'].rsplit('_', 3)[0]
        genome_id = contig2genome_dic[contig_id]
        if genome_id.startswith('GCF_'):
            proteome = list(SeqIO.parse(proteomes_dir + genome_id+'.fasta.FGS.faa', 'fasta'))
            out_f.write('\n'.join(['>'+seq.id+'\n'+str(seq.seq).replace('*', '') for seq in proteome]))
            orf_id = row['ORFs']['orf_id']
            orf_seq = contigs_seqs[orf_id]
            out_f2.write('>' + str(orf_seq.id)+ '\n')
            out_f2.write(str(orf_seq.seq)+'\n')

