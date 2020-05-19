"""
script that will perform an exact string matching between the peptides
and protein sequences from the genomes and will rerturn a table
reporting these matches.

this script requires three command line arguments
the first argument being the file that has the peptides to be searched in fasta format
these peptides could be obtained by any peptide searching methods such as MSGF+ or some other way

second command line argument is the file containing all the proteins coming from all the genomes 
to parallelize this process then each proteome is kept in a separate physical file that way many 
prpteomes are searched against the peptides in parallel, depending on how many cores are utilized.
I pre-compiled a file by predicting from all the contigs the protein coding genes and then translating these
genes into their relative amino acid sequences. If you have not done so you need to precomple such a file
before runing this script
the third command line argument is the directory for the output file:
the fourth command line argument is the number of threads used to parellize this prorcess
the fifth command line argument is specifying the suffix where this 
"""

from Bio import SeqIO
import sys
import multiprocessing as mp
import os
import time

if len(sys.argv) != 6:
    print('please enter 5 command line arguments to run this script, example to run script i.e. python3 peptides2sequenceStringMatching.py dir/2/peptides/file dir/to/allProteinsSeqs/ dir/to/output/directory n_threads suffix for the files to search into')

else:
    peptides_seq_f = sys.argv[1]
    prot_seqs_dir = sys.argv[2]
    out_dir = sys.argv[3]    
    n_threads = int(sys.argv[4])
    suffix = sys.argv[5]
    
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    print('reading all peptide sequences into memory')
    peptide_seqs = SeqIO.to_dict(SeqIO.parse(peptides_seq_f, 'fasta'))
    print('reading all proteins into memory')
    all_prot_files = [file for file in os.listdir(prot_seqs_dir) if file.endswith(suffix)]
    
    peptide_seq_list = [(peptide_id ,str(peptide_seqs[peptide_id].seq)) for peptide_id in peptide_seqs]
    
    new2oldPeptide_dic = dict()
    new_peptide_list = list()
    
    print('converting peptides I 2 L')
    for peptide in peptide_seq_list:
        new_peptide = peptide[1].replace('I', 'L')
        new_tuple = (peptide[0], new_peptide)
        new2oldPeptide_dic[new_peptide] = peptide[1]
        new_peptide_list.append(new_tuple)
    
    print('searching', len(new_peptide_list), 'peptides')
    print('against', len(all_prot_files), 'proteomes')
    print('searching peptides against the proteins')    
    
    def peptide2protMatch(prot_f, out_dir = out_dir, prot_seqs_dir = prot_seqs_dir):
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
