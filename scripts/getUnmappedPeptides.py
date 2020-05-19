"""
script where I will filter the peptides identified by the 6 frame translation method further.
This is basically done by looking into those peptides and searching them against the proteome
of the genomes that were used to construct the 6 frame translated database
"""

from Bio import SeqIO
import sys
import multiprocessing as mp
import os
import time
import pandas as pd
import re


#arg1: the output directory for the MSGF+ search results; mzid out dir
#arg2: name of the sample, should be the same as the names of the subfolders in the arg1 folder
#arg3: top N species used in HAPiID original search
#arg4: top N species used for the six frame translation, database construction and search
#arg5: the number of threads used for direct peptide matching between peptides and genes.

if len(sys.argv) != 6:
    print('please enter 5 command line arguments')
else:

    results_out, sample, HAPiID_topN, sixFrame_topN, n_threads = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]

ranking_f = results_out + sample + '_ribP_elonF/'+sample+'.GenomeCoveringAllHEGpeptides.txt'
unique_peptides = results_out+'unique_peptides/'
proteomes_dir = '../data/proteomes/'
topN = int(re.sub("[^0-9]", "", sixFrame_topN))


out_dir = results_out+'peptideMatches/'+sample + '_' + HAPiID_topN + '_' + sixFrame_topN + '/'


ranking_df = pd.read_csv(ranking_f, sep = '\t')

topN_genomes = list(ranking_df['genome'])[:topN]

all_prot_files = [file+'.fasta.FGS.faa' for file in topN_genomes]


print('processing sample ', sample)

compareWith_f = [file for file in [item for item in os.listdir(unique_peptides) if item.endswith('_6frameTranslationPeptides.fasta')] if sample in file and sixFrame_topN in file][0]
twoStepAppraoch_f = [file for file in [item for item in os.listdir(unique_peptides) if item.endswith('_all_peptides_'+str(HAPiID_topN)+'.fasta')] if sample in file][0]

with open(unique_peptides + compareWith_f, 'r') as six_frame_in_f:
    new2oldPeptide_dic = {str(seq.seq).replace('I', 'L'):str(seq.seq) for seq in list(SeqIO.parse(six_frame_in_f, 'fasta'))}

with open(unique_peptides + compareWith_f, 'r') as compWith_f, open(unique_peptides + twoStepAppraoch_f, 'r') as twoStep_f:
    compareWith_seqs = set([str(seq.seq).replace('I', 'L') for seq in list(SeqIO.parse(compWith_f, 'fasta'))]) #6 frame translated peptides
    twoStep_seqs = set([str(seq.seq).replace('I', 'L') for seq in list(SeqIO.parse(twoStep_f, 'fasta'))]) #FGS protein peptides

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

peptide_list = compareWith_seqs.difference(twoStep_seqs) #Peptides identified using 6 frame translation as reference DB but not using FGS protein DB
pep_dic = {'peptide'+str(i):peptide_seq for i, peptide_seq in enumerate(peptide_list)}
peptide_seq_list = [(item[0], item[1])  for item in pep_dic.items()]


new_peptide_list = list()

print('converting peptides I 2 L')
for peptide in peptide_seq_list:
    new_peptide = peptide[1].replace('I', 'L')
    new_tuple = (peptide[0], new_peptide)
    new_peptide_list.append(new_tuple)


print('searching peptides against the proteins')

def peptide2protMatch(prot_f, out_dir = out_dir, prot_seqs_dir = proteomes_dir):
    start_time = time.time()
    #print('processing '+prot_f+' ...')
    with open(out_dir+prot_f.rsplit('.', 1)[0]+'_peptideMatches.txt', 'w') as out_f:
        for i, peptide in enumerate(new_peptide_list):
            for seq in list(SeqIO.parse(prot_seqs_dir + prot_f, 'fasta')):
                prot_seq = str(seq.seq).replace('I', 'L')
                prot_id = str(seq.id)
                if peptide[1] in prot_seq:
                    out_f.write(str(peptide[0])+'\t'+str(new2oldPeptide_dic[peptide[1]])+'\t'+str(prot_id)+'\n')
   # print('processed', prot_f, 'in: ', time.time() - start_time, 'seconds')
pool = mp.Pool(int(n_threads))

zip([*pool.map(peptide2protMatch, all_prot_files)])
if out_dir.endswith('/'):
    out_2_dir = out_dir.rsplit('/', 2)[0]+'/'
    prefix = out_dir.split('/')[-2]
else:
    out_2_dir = out_dir.rsplit('/', 1)[0]+'/'
    prefix = out_dir.split('/')[-1]

os.system("cat "+out_dir+"* > " + out_2_dir + prefix + "_peptide2allProteinMatches.txt")

with open(out_2_dir + prefix + '_peptide2allProteinMatches.txt', 'r') as in_f:
    fline = "peptide_ID\tpeptide\tprotein_match_ID\n"
    oline = in_f.readlines()
    oline.insert(0, fline)

with open(out_2_dir + prefix +'_peptide2allProteinMatches.txt', 'w') as out_f:
    out_f.writelines(oline)
    
mapped_ids = list()
match_df = pd.read_csv(out_2_dir + prefix +'_peptide2allProteinMatches.txt', sep = '\t')
mapped_ids = list(set(match_df['peptide_ID']))

unmapped_ids = set(list(pep_dic.keys())).difference(set( mapped_ids))

with open(out_2_dir + prefix + '_unmappedPeptides.fasta', 'w') as out_f:
    for _id in unmapped_ids:
        out_f.write('>'+str(_id)+'\n')
        out_f.write(new2oldPeptide_dic[pep_dic[_id]]+'\n')
print('number of unmapped peptides in', sample, len(unmapped_ids))    
