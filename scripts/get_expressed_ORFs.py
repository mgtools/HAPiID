"""
script where, for each contig predicted 6 frame contig translated I will get all the peptides
that are not mapped to any of the FGS predicted genes which is previously calculated using the script
getUnmappedPeptides.py, then I will map these peptides on the respective contigs that they were found 
during  MSGF+ search and then I will find the closest ORF encapsulating them, and filter out 
the proper ORFs, i.e. the ones starting with start codon and ending with stop codon and no 
stop codon in between, Thus the resulting list of ORFs are hence expressed with proteomics evidence.
This script will output all these potential ORFs to a file under the folder expressedORFs
"""

from Bio import SeqIO
import pandas as pd
import copy
import re
import sys

results_out = '../../HAPiID_results/'
db_out = '../../HAPiID_db/'
sample = 'HM403'
HAPiID_topN = 'top50'
sixFrame_topN = 'top10'

if len(sys.argv) != 6:
    print('please enter 5 command line arguments')
else:
    results_out, db_out, sample, HAPiID_topN, sixFrame_topN = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]


tsv_f = results_out + sample +'_' + sixFrame_topN + '_6frameTranslation_search/' + sample+'.tsv.0.01.tsv'
tsv_df = pd.read_csv(tsv_f, sep = '\t')


topn = int(re.sub("[^0-9]", "", sixFrame_topN))
db_f = db_out + sample + '_' + sixFrame_topN + '_6frameTranslation/top_'+str(topn)+'_most_abundant_6FrameTranslated.fasta'
db_seqs = SeqIO.to_dict(SeqIO.parse(db_f, 'fasta')) #protein sequences used to construct the 6 frame translatede protein databases

#convert the identified peptides to AA sequences from the tsv file
peptides_col = [''.join(filter(str.isalpha, item)) for item in list(tsv_df['Peptide'])]
tsv_peptide_df = copy.deepcopy(tsv_df)
tsv_peptide_df['Peptide'] = [item.replace('I', 'L') for item in peptides_col] #convert all the Isoleucines to Leucines in the peptide identification table

unmapped_peptides_seq_f = results_out + 'peptideMatches/' + sample + '_' + HAPiID_topN + '_' + sixFrame_topN + '_unmappedPeptides.fasta'

new2oldPeptide_dic = {str(item.seq).replace('I', 'L'):str(item.seq) for item in list(SeqIO.parse(unmapped_peptides_seq_f, 'fasta'))} #map the new peptides to the original sequences

unmapped_peptides_seqs = [str(item.seq).replace('I', 'L') for item in list(SeqIO.parse(unmapped_peptides_seq_f, 'fasta'))] #list of unmapped peptides created from the previous script 
                                                                                                                           #(not mapped to FGS predicted proteins), also convert I to L 
                                                                                                                           #to match the peptides in the TSV file

#this line will only keep the rows where unmapped peptides are mapped to real contigs (non decoy sequences in the database)
#all isoleucinces are converted to leucines on both ends
tsv_unmapped_peptides_df = tsv_peptide_df[tsv_peptide_df['Peptide'].isin(unmapped_peptides_seqs) &  ~tsv_peptide_df['Protein'].str.startswith('XXX_')]

contig2UnmappedPeptides_dic = dict()
for i, row in tsv_unmapped_peptides_df.iterrows():
    contigID = row['Protein'].split('(pre')[0]
    if contigID not in contig2UnmappedPeptides_dic:
        contig2UnmappedPeptides_dic[contigID] = list()
    peptide = row['Peptide']
    if peptide not in contig2UnmappedPeptides_dic[contigID]:
        contig2UnmappedPeptides_dic[contigID].append(peptide)

orfs = list()
new2oldContig_dic = dict()
new2oldORF_dic = dict()
for contig in contig2UnmappedPeptides_dic:
    original_contig_seq = str(db_seqs[contig].seq)
    contig_seq = str(db_seqs[contig].seq).replace('I', 'L') #temporarily change all the Isoleucines to leucines for peptide lookup
    new2oldContig_dic[contig_seq] = original_contig_seq
    contig_peptides = contig2UnmappedPeptides_dic[contig]
    contig_id = contig
    for peptide in contig_peptides:
        coordinates = [(i.start(), i.end()) for i in re.finditer(peptide, contig_seq)]
        for coord in coordinates:
            start, end  = coord[0], coord[1]
            nearset_methionine = [i.start() for i in re.finditer('M', contig_seq[: start])]
            if nearset_methionine:
                nearset_methionine = nearset_methionine[-1]
            else:
                nearset_methionine = 0
            nearest_stop = [i.start() for i in re.finditer('\*', contig_seq[end:])]
            if nearest_stop:
                nearest_stop = nearest_stop[0]
            else:
                nearest_stop = end #the end coordinate for the peptide
            
            orf = contig_seq[nearset_methionine: end + nearest_stop+1]#make sure you get the ORF after the stop AA coordinate
            original_orf = original_contig_seq[nearset_methionine: end + nearest_stop+1]
            new2oldORF_dic[orf] = original_orf
            orfs.append((contig_id, orf, peptide))
            

out_fname = results_out + 'expressedORFs/' + sample +'_' +HAPiID_topN + '_' + sixFrame_topN +'_translatedORFs.fasta'
out_fname_2 = results_out + 'expressedORFs/' + sample +'_' +HAPiID_topN + '_' + sixFrame_topN +'_ORF_peptides.tsv'
            
orf_counter = 0
filtered_orfs = list()

with open(out_fname, 'w') as out_f, open(out_fname_2, 'w') as out_f2:
    out_f2.write('ORF_id\tpeptide_seqs\tfrom-to\n')
    orf2peptide2coords_dic = dict()
    for orf_tuple in orfs:
        orf = orf_tuple[1]
        peptide = orf_tuple[2]
        contig_id = orf_tuple[0]
        if orf.startswith('M') and orf.endswith('*'): #make sure starts with M and ends with *
            if '*' not in orf[:-1]: #make sure there's not premature stop in between                
                filtered_orfs.append(orf)
                if orf not in orf2peptide2coords_dic:
                    orf2peptide2coords_dic[orf] = dict()
                    orf_id = sample+'|'+contig_id+'|orf'+ str(orf_counter+1)
                    orf2peptide2coords_dic[orf]['id'] = orf_id
                    orf2peptide2coords_dic[orf]['coord'] = list()
                    orf2peptide2coords_dic[orf]['peptide'] = list()
                    orf_counter += 1
                    
                coordinates = [(i.start(), i.end()) for i in re.finditer(peptide, orf)]
                if len(coordinates) == 1:
                    start = coordinates[0][0]
                    stop = coordinates[0][1]
                    orf2peptide2coords_dic[orf]['coord'].append(str(start)+'-'+str(stop))
                    orf2peptide2coords_dic[orf]['peptide'].append(peptide)
                else:
                    coord_lst = list()
                    peptide_lst = list()
                    for coord in coordinates:
                        coord_lst.append(str(coord[0])+'-'+str(coord[1]))
                        peptide_lst.append(peptide)
                    orf2peptide2coords_dic[orf]['coord'].extend(coord_lst)
                    orf2peptide2coords_dic[orf]['peptide'].extend(peptide_lst)

    for orf in orf2peptide2coords_dic:
        out_f.write('>'+orf2peptide2coords_dic[orf]['id']+'\n')
        out_f.write(new2oldORF_dic[orf]+'\n')
        peptides, coords = orf2peptide2coords_dic[orf]['peptide'], orf2peptide2coords_dic[orf]['coord']
        for peptide, coord in zip(peptides, coords):
            out_f2.write(orf2peptide2coords_dic[orf]['id']+'\t'+str(new2oldPeptide_dic[peptide])+'\t'+coord+'\n')
                        

print(str(orf_counter)+ ' expressed ORFs identified')
