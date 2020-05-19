"""
parse BlastP table output and for each query report results from the three blast sources
coming from the three different predictors
"""

import pandas as pd
import os
import sys

pd.set_option('display.max_columns', 5400)
pd.set_option('display.max_rows', 5400)

#FGS_blast_out_f = '/data/mstambou/proteome_landscapes/HAPiID_results/verifyORFs/blast_out/allExpressedORFsSoFar_c1.0_FGS_blastout.txt'
#prodigal_blast_out_f = '/data/mstambou/proteome_landscapes/HAPiID_results/verifyORFs/blast_out/allExpressedORFsSoFar_c1.0_prodigal_blastout.txt'
#genemark_blast_out_f = '/data/mstambou/proteome_landscapes/HAPiID_results/verifyORFs/blast_out/allExpressedORFsSoFar_c1.0_genemark_blastout.txt'

orf2peptides_dir = '/data/mstambou/proteome_landscapes/HAPiID_results/expressedORFs/'

FGS_blast_out_f = sys.argv[1]
prodigal_blast_out_f = sys.argv[2]
genemark_blast_out_f = sys.argv[3]
orf2peptides_dir = sys.argv[4]

orf2peptides_tables = [item for item in os.listdir(orf2peptides_dir) if item.endswith('.tsv')]

allORF2peptides_df = pd.DataFrame(columns = ['ORF_id', 'peptide_seqs', 'from-to'])
for file in orf2peptides_tables:
    df = pd.read_csv(orf2peptides_dir + file, sep = '\t')
    allORF2peptides_df = pd.concat([allORF2peptides_df, df], ignore_index = True)
    
allORF2peptides_dic = dict()

for i, row in allORF2peptides_df.iterrows():
    orf = row['ORF_id']
    peptide = row['peptide_seqs']
    coord = row['from-to']
    if orf not in allORF2peptides_dic:
        allORF2peptides_dic[orf] = dict()
        allORF2peptides_dic[orf]['peptide'] = list()
        allORF2peptides_dic[orf]['coord'] = list()
    allORF2peptides_dic[orf]['peptide'].append(peptide)
    allORF2peptides_dic[orf]['coord'].append(coord)
    
allORF2peptides_df = pd.DataFrame(columns = ['orf_id', 'peptides', 'coords'])
for i, orf in enumerate(allORF2peptides_dic):
    allORF2peptides_df.loc[i] = [orf, '|'.join(allORF2peptides_dic[orf]['peptide']), '|'.join(allORF2peptides_dic[orf]['coord'])]
    
allORF2peptides_df = allORF2peptides_df.sort_values(by=['orf_id']).reset_index(drop=True)

allORF_ids = list(allORF2peptides_dic.keys())

out_dir = genemark_blast_out_f.rsplit('/',1)[0]+'/'
blast_cols = ['query acc.ver', 'subject acc.ver', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score']

FGS_blast_out_df = pd.read_csv(FGS_blast_out_f, sep = '\t', comment = '#', header = None)
prodigal_blast_out_df = pd.read_csv(prodigal_blast_out_f, sep = '\t', comment = '#', header = None)
genemark_blast_out_df = pd.read_csv(genemark_blast_out_f, sep = '\t', comment = '#', header = None)

FGS_blast_out_df.columns = blast_cols
prodigal_blast_out_df.columns = blast_cols
genemark_blast_out_df.columns = blast_cols

FGS_best_blast_out_df = FGS_blast_out_df.drop_duplicates(subset = ['query acc.ver'], keep = 'first', inplace = False)
prodigal_best_blast_out_df = prodigal_blast_out_df.drop_duplicates(subset = ['query acc.ver'], keep = 'first', inplace = False)
genemark_best_blast_out_df = genemark_blast_out_df.drop_duplicates(subset = ['query acc.ver'], keep = 'first', inplace = False)

#append rows for missing ORFS(if orf did not receive any signifcant hit)
missing_vals = ['NA']*11
cols = list(FGS_best_blast_out_df.columns)

FGS_missing = set(allORF_ids).difference(set(list(FGS_best_blast_out_df['query acc.ver']) ))
FGS_missing_df = pd.DataFrame(  [[orf_id] + missing_vals for orf_id in FGS_missing], columns =  cols)
FGS_best_blast_out_df = FGS_best_blast_out_df.append(FGS_missing_df, ignore_index = True).reset_index(drop = True)

prodigal_missing = set(allORF_ids).difference(set(list(prodigal_best_blast_out_df['query acc.ver']) ))
prodigal_missing_df = pd.DataFrame(  [[orf_id] + missing_vals for orf_id in prodigal_missing], columns =  cols)
prodigal_best_blast_out_df = prodigal_best_blast_out_df.append(prodigal_missing_df, ignore_index = True).reset_index(drop = True)

genemark_missing = set(allORF_ids).difference(set(list(genemark_best_blast_out_df['query acc.ver']) ))
genemark_missing_df = pd.DataFrame(  [[orf_id] + missing_vals for orf_id in genemark_missing], columns =  cols)
genemark_best_blast_out_df = genemark_best_blast_out_df.append(genemark_missing_df, ignore_index = True).reset_index(drop = True)

FGS_best_blast_out_df = FGS_best_blast_out_df.sort_values(by=['query acc.ver']).reset_index(drop=True)
FGS_best_blast_out_df.drop(['query acc.ver', 'mismatches', 'gap opens',  'evalue', 'bit score'], axis = 1, inplace = True)
prodigal_best_blast_out_df = prodigal_best_blast_out_df.sort_values(by=['query acc.ver']).reset_index(drop=True)
prodigal_best_blast_out_df.drop(['query acc.ver', 'mismatches', 'gap opens',  'evalue', 'bit score'], axis = 1, inplace = True)
genemark_best_blast_out_df = genemark_best_blast_out_df.sort_values(by=['query acc.ver']).reset_index(drop=True)
genemark_best_blast_out_df.drop(['query acc.ver', 'mismatches', 'gap opens',  'evalue', 'bit score'], axis = 1, inplace = True)


d = {}
d['ORFs'] = allORF2peptides_df
d['FGS'] = FGS_best_blast_out_df
d['prodigal'] = prodigal_best_blast_out_df
d['genemark'] = genemark_best_blast_out_df

concatinated_df = pd.concat(d, axis =1)
concatinated_df.to_csv(out_dir + 'blastHitsSummary.tsv', sep = '\t', index = None)
