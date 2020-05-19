"""
script that will parse the gene mark output and will split the file into 2 other files
one containing the predicted proteins and another containing the respective genes
"""

import sys

#genemark_f = '/data/mstambou/proteome_landscapes/HAPiID_results/verifyORFs/genemark_out/contigsWithPeptides.fasta.genemark'

genemark_f = sys.argv[1]

f_name = genemark_f.rsplit('/', 1)[-1]
proteins_f = f_name+'.faa'
genes_f = f_name+'.ffn'
out_dir = genemark_f.rsplit('/', 1)[0]+'/'

protein_flag = False
gene_flag = False

with open(genemark_f, 'r') as in_f, open(out_dir + proteins_f, 'w') as proteins_out, open(out_dir + genes_f, 'w') as genes_out:
    for line in in_f:
        line = line.strip()
        if line.startswith('Predicted proteins:'):            
            protein_flag = True
            continue
        if line.startswith('Nucleotide sequence of predicted genes:'):
            gene_flag = True
            protein_flag = False
            continue
        if protein_flag:
            if line:
                proteins_out.write(line+'\n')
        if gene_flag:
            if line:
                genes_out.write(line + '\n')
        
