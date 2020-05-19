"""
script that will perform an exact string matching between the peptides
and protein sequences from the genomes and will rerturn a table
reporting these matches.

this script requires three command line arguments
the first argument being the file that has the peptides to be searched in fasta format
these peptides could be obtained by any peptide searching methods such as MSGF+ or some other way
second command line argument is the file containing all the proteins coming from all the genomes
I pre-compiled a file by predicting from all the contigs the protein coding genes and then translating these
genes into their relative amino acid sequences. If you have not done so you need to precomple such a file
before runing this script
the third command line argument is the directory for the output file:
"""

from Bio import SeqIO
import sys
import multiprocessing as mp
import os
import time

if len(sys.argv) != 5:
    print('please enter 4 command line arguments to run this script, example to run script i.e. python3 peptides2sequenceStringMatching.py dir/2/peptides/file dir/to/allProteinsSeqs/ dir/to/output/directory n_threads')

else:
    peptides_seq_f = sys.argv[1]
    prot_seqs_dir = sys.argv[2]
    out_dir = sys.argv[3]    
    n_threads = int(sys.argv[4])
    
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    
    peptide_seqs = SeqIO.to_dict(SeqIO.parse(peptides_seq_f, 'fasta'))
    print('reading all proteins into memory')
    all_prot_files = [file for file in os.listdir(prot_seqs_dir) if file.endswith('.fasta.FGS.faa')]
    
    peptide_seq_list = [(peptide_id ,str(peptide_seqs[peptide_id].seq)) for peptide_id in peptide_seqs]
    
    new2oldPeptide_dic = dict()
    new_peptide_list = list()
    
    print('converting peptides I 2 L')
    for peptide in peptide_seq_list:
        new_peptide = peptide[1].replace('I', 'L')
        new_tuple = (peptide[0], new_peptide)
        new2oldPeptide_dic[new_peptide] = peptide[1]
        new_peptide_list.append(new_tuple)
    
    print(len(new_peptide_list))
    print(len(all_prot_files))
    print('searching peptides against the proteins')    
    
    def peptide2protMatch(prot_f, out_dir = out_dir, prot_seqs_dir = prot_seqs_dir):
        if os.path.exists(out_dir+prot_f.rsplit('.', 1)[0]+'_peptideMatches.txt') == False: #check first whether or not this genome is processed
            start_time = time.time()
            print('processing '+prot_f+' ...')
            with open(out_dir+prot_f.rsplit('.', 1)[0]+'_peptideMatches.txt', 'w') as out_f:
                for i, peptide in enumerate(new_peptide_list):
                    for seq in list(SeqIO.parse(prot_seqs_dir + prot_f, 'fasta')):
                        prot_seq = str(seq.seq).replace('I', 'L')
                        prot_id = str(seq.id)
                        if peptide[1] in prot_seq:
                            out_f.write(str(peptide[0])+'\t'+str(new2oldPeptide_dic[peptide[1]])+'\t'+str(prot_id)+'\n')
            print('processed', prot_f, 'in: ', time.time() - start_time, 'seconds')

    pool = mp.Pool(int(n_threads))
    
    zip([*pool.map(peptide2protMatch, all_prot_files)])
    if out_dir.endswith('/'):
        out_2_dir = out_dir.rsplit('/', 2)[0]+'/'
    else:
        out_2_dir = out_dir.rsplit('/', 1)[0]+'/'

    os.system("cat "+out_dir+"* > " + out_2_dir + "peptide2allProteinMatches.txt")
    
    with open(out_2_dir + 'peptide2allProteinMatches.txt', 'r') as in_f:
        fline = "peptide_ID\tpeptide\tprotein_match_ID\n"
        oline = in_f.readlines()
        oline.insert(0, fline)

    with open(out_2_dir + 'peptide2allProteinMatches.txt', 'w') as out_f:
        out_f.writelines(oline)
