# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 18:57:49 2019

@author: mstambou

script in which I will quantify the abundance of genomes composed of only marker proteins in this case.
i.e. the genome length here is the sum of the lengths for those marker genes coming from that bin.
Hence the program needs first genome2marker genes fasta lists. 
The quantification in this case is done following the qin et al paper for the liver cirrhosis.
In the case of proteomics you treat the spec IDs as the reads, i.e. if the same peptide is mentioned
twice, i.e. in two different rows but has different spectrum IDs then it should be counted twice 
in the quantification step.
if the same peptide but different spec IDs is mapped twice to the same protein then this should be 
treated as unique maps

in the end once the script is done calculating the abundances it will generate a CSV file that has 
unique abundances and multi abundances and total abundances (unique + multi) as columns. in PMKB
it will also produce a boxplot summarizing the distribution of these maps in unique and multiple reads.

This script also takes the cd-hit output and uses it as a mapping between clsuters and members
and expands on the hits for the cases where genes are redundant.
"""


def getGenomeLength(genome2ribPelonF_dic, prot_seqs):
    """
    function that returns the length of a genome based on the proteins that are mapped to it
    (i..e. it could be just by using some marker genes/proteins)
    """
    genome2proteome_lengths_dic = dict()
    for genome in genome2ribPelonF_dic:
        genome_proteins = genome2ribPelonF_dic[genome]
        total_length = sum([len(str(prot_seqs[prot].seq)) for prot in genome_proteins])
        genome2proteome_lengths_dic[genome] = total_length
    return genome2proteome_lengths_dic

def get_uniquelyMappedPeptideSpectra(peptideSpectrum_list, peptideSpectrum2allGenomesWithUniqueSpectra_dic):
    """
    function that goes over all the pepetide Spectra mapped to genomes in the library and returns 
    the ones that are only uniquely maped to one proteins and no where else.
    """
    uniquely_mapped_peptideSpectra = []
    for peptideSpectrum in peptideSpectrum_list:
        if len(peptideSpectrum2allGenomesWithUniqueSpectra_dic[peptideSpectrum]) == 1:
            uniquely_mapped_peptideSpectra.append(peptideSpectrum)
    return uniquely_mapped_peptideSpectra

def get_multiMappedPeptideSpectra(peptideSpectrum_list, peptideSpectrum2allGenomesWithUniqueSpectra_dic):
    """
    function that goes over all the pepetide Spectra mapped to genomes in the library and returns 
    the ones that are mapped to more than one protein
    """
    multi_mapped_peptideSpectra = []
    for peptideSpectrum in peptideSpectrum_list:
        if len(peptideSpectrum2allGenomesWithUniqueSpectra_dic[peptideSpectrum]) > 1:
            multi_mapped_peptideSpectra.append(peptideSpectrum)
    return multi_mapped_peptideSpectra

def getGenomeUniquePeptideSpectrumCounts(genomesWithUniqueMappedPeptideSpectra, unique_proteins_expanded, unique_mapped_proteins, protein2peptideSpectrum_dic, peptideSpectrum2protein_dic, peptideSpectrum2allGenomesWithUniqueSpectra_dic):
    """
    Function that for each genome returns the number of peptide Spectra that are uniquely mapped to them and nothing else.
    Function that for each genome get its proteome, then get proteins that have shown abundance in this sample
    i.e. proteins that have at least one peptide spectrum being mapped to them, then for each of thes genomes
    get the uniquely mapped peptide spectra to them, and for each genome return the counts, i.e. the number
    of unique peptide spectra that are mapped to them.
    this function will also return genome to all its peptide spectra maps (both unqiue and multimaps)
    this function will also return peptide spectra to genomes mapping (i.e. a dictionary listing all the genomes that a certain peptide is mapped to
    will be one genome if that peptide is uniquely mapped)
    """
    genome2uniquePeptideSpectrum_counts_dic = dict()
    genome2allPeptideSpectra_dic = dict()    
    for genome in genomesWithUniqueMappedPeptideSpectra:
        genome_prots = genome2ribPelonF_dic[genome]
        #proteins in the genome that have peptides matched to them in this sample
        genome_abundant_prots = list(set(genome_prots).intersection(set(unique_proteins_expanded)))        
        genome_mapped_peptideSpectra = list(set([peptideSpectrum for prot in genome_abundant_prots for peptideSpectrum in protein2peptideSpectrum_dic[prot]]))
        genome2allPeptideSpectra_dic[genome] = genome_mapped_peptideSpectra
        uniquely_mapped_genome_peptideSpectra = get_uniquelyMappedPeptideSpectra(genome_mapped_peptideSpectra, peptideSpectrum2allGenomesWithUniqueSpectra_dic)
        genome2uniquePeptideSpectrum_counts_dic[genome] = len(uniquely_mapped_genome_peptideSpectra)        
                        
    return genome2uniquePeptideSpectrum_counts_dic, genome2allPeptideSpectra_dic
    
def get_uniqueAbundance_RPKM(genome2uniquePeptideSpectrum_counts_dic, genome2proteome_lengths_dic, total_n_peptideSpectrum_mapped_from_library):
    """
    function that calculates genome abundances using unique maps only. it reports the Reads Per Million Killo-Base (RPKM) values
    hence it normalizes for both sequencing depth and genome lengths, when it reports abundances. For sequencing depth I divide 
    the values by the total number of peptide Spectra in my library that are being mapped to the genomes.
    """
    genome_uniqueAbundanceRPKM_dic = dict()
    for genome in genome2uniquePeptideSpectrum_counts_dic:
        n_unique_peptideSpectra = genome2uniquePeptideSpectrum_counts_dic[genome]
        genome_length = genome2proteome_lengths_dic[genome]
        unique_RPKM = ((n_unique_peptideSpectra)*(10**3)*(10**6))/(genome_length*total_n_peptideSpectrum_mapped_from_library)
        genome_uniqueAbundanceRPKM_dic[genome] = unique_RPKM
    return genome_uniqueAbundanceRPKM_dic         
    
def get_speciesSpecificCoefficient(genome, multi_peptideSpectrum_genomes, genome_uniqueAbundanceRPKM_dic):
    unique_abundances = list()
    for multi_genome in multi_peptideSpectrum_genomes:
        unique_abundances.append(genome_uniqueAbundanceRPKM_dic[multi_genome])
    coefficient = genome_uniqueAbundanceRPKM_dic[genome]/float(sum(unique_abundances))
    return coefficient

def get_multiAbundance_RPKM(genomesWithUniqueMappedPeptideSpectra, genome2allPeptideSpectra_dic, peptideSpectrum2allGenomesWithUniqueSpectra_dic, genome_uniqueAbundanceRPKM_dic, genome2proteome_lengths_dic, total_n_peptideSpectrum_mapped_from_library):
    genome_multiAbundanceRPKM_dic = dict()
    for genome in genomesWithUniqueMappedPeptideSpectra:        
        genome_all_peptideSpectra = genome2allPeptideSpectra_dic[genome]
        genome_multi_peptideSpectra = get_multiMappedPeptideSpectra(genome_all_peptideSpectra, peptideSpectrum2allGenomesWithUniqueSpectra_dic)
        coefficient_lst = list()
        for multi_peptideSpectrum in genome_multi_peptideSpectra:            
            multi_peptideSpectrum_genomes = peptideSpectrum2allGenomesWithUniqueSpectra_dic[multi_peptideSpectrum]          
            coefficient = get_speciesSpecificCoefficient(genome, multi_peptideSpectrum_genomes, genome_uniqueAbundanceRPKM_dic)
            coefficient_lst.append(coefficient)    
        genome_length = genome2proteome_lengths_dic[genome]
        n_multi_peptideSpectra = sum(coefficient_lst)
        multi_RPKM = ((n_multi_peptideSpectra)*(10**3)*(10**6))/(genome_length*total_n_peptideSpectrum_mapped_from_library)
        genome_multiAbundanceRPKM_dic[genome] = multi_RPKM
    
    return genome_multiAbundanceRPKM_dic      

def get_log10(val):
    if val == 0:
        return 0
    else:
        return math.log10(val)

#######################main###############################
import pickle
from Bio import SeqIO
import pandas as pd
import utils
import matplotlib
import copy
import math
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import json

if len(sys.argv) != 6:
    print('please enter 5 command line arguments, i.e. example to run\n python3 quantify_peptide_mapping.py data/genome2ribPelonF_dic.pickle ribP_elonF_db/ribP_elonF_all_proteinSeqs.faa 131211-28623-ATH-F_MSGF+/131211-28623-12-ATH-F01-03.tsv.0.01.tsv data/ribP_elonF_all_cdHit_100.fasta.clstr data/umgs_hgg_proteins2genomes_dic.json' )
else:
    #precalculated dictionary containing genome/bin-ID to it's gene IDs mappings
    #genome2ribPelonF_dic = pickle.load(open('genome2ribPelonF_dic.pickle', 'rb'))
    genome2ribPelonF_dic = pickle.load(open(sys.argv[1], 'rb'))
    #file containing full length protein sequences for ribosomal proteins and elongation factors
    #protein_seqs_f = 'ribP_elonF_db/ribP_elonF_all_proteinSeqs.faa'
    protein_seqs_f = sys.argv[2]
    #The actual TSV file for the FDR filtered peptide identification
    #all_tsv_f = '131211-28623-ATH-F_MSGF+/131211-28623-12-ATH-F01-03.tsv.0.01.tsv'
    all_tsv_f = sys.argv[3]
    #cdhit output that maps cluster representatives to cluster members.
    #cdhit_out_f = '131211-28623-ATH-F_MSGF+/ribP_elonF_all_cdHit_100.fasta.clstr'
    cdhit_out_f = sys.argv[4]
    #sample_name = (all_tsv_f.rsplit('/', 1)[-1]).split('.tsv')[0]
    with open(sys.argv[5], 'r') as in_f:
        umgs_hgg_proteins2genomes_dic = json.load(in_f)
    #sample_name = (all_tsv_f.rsplit('/', 1)[-1]).split('.tsv')[0]
    
    abundance_out_f = '.'.join(all_tsv_f.split('.')[:-1])+'_sortedAbundance.tsv'
    all_tsv_df = pd.read_csv(all_tsv_f, sep = '\t')
    #this step makes sure if there were multiple rows of the same spectrum ID then keep the first one (i.e. the best one)
    all_tsv_df = all_tsv_df.drop_duplicates(subset = 'SpecID', keep = 'first')
    
    prot_seqs = SeqIO.to_dict(SeqIO.parse(protein_seqs_f, 'fasta'))
    
    unique_peptides, unique_proteins, protein2peptideSpectrum_dic, peptideSpectrum2protein_dic = utils.get_peptideSpectra_proteins_mappings(all_tsv_df)
    
    protein2peptideSpectrum_expanded_dic, peptideSpectrum2protein_expanded_dic = utils.expand_cdhit_clustered_prots(cdhit_out_f, protein2peptideSpectrum_dic, peptideSpectrum2protein_dic)
    
    unique_proteins_expanded = list(protein2peptideSpectrum_expanded_dic.keys())
    
    unique_mapped_proteins = utils.get_proteinsWithUniquePeptideSpectrum(peptideSpectrum2protein_expanded_dic, umgs_hgg_proteins2genomes_dic)
    genomesWithUniqueMappedPeptideSpectra = utils.get_genomes_with_unique_peptides(peptideSpectrum2protein_expanded_dic, umgs_hgg_proteins2genomes_dic)
    
    peptideSpectrum2allGenomesWithUniqueSpectra_dic = utils.get_peptideSpectrum2genomeWithUniqueSpectra_maps(peptideSpectrum2protein_expanded_dic, genomesWithUniqueMappedPeptideSpectra, umgs_hgg_proteins2genomes_dic)
    
    genome2uniquePeptideSpectrum_counts_dic, genome2allPeptideSpectra_dic = getGenomeUniquePeptideSpectrumCounts(genomesWithUniqueMappedPeptideSpectra, unique_proteins_expanded, unique_mapped_proteins, protein2peptideSpectrum_expanded_dic, peptideSpectrum2protein_expanded_dic, peptideSpectrum2allGenomesWithUniqueSpectra_dic)
    genome2proteome_lengths_dic = getGenomeLength(genome2ribPelonF_dic, prot_seqs)
    
    total_n_peptideSpectrum_mapped_from_library = len(set(all_tsv_df['SpecID']))
    
    genome_uniqueAbundanceRPKM_dic = get_uniqueAbundance_RPKM(genome2uniquePeptideSpectrum_counts_dic, genome2proteome_lengths_dic, total_n_peptideSpectrum_mapped_from_library)
    genome_multiAbundanceRPKM_dic = get_multiAbundance_RPKM(genomesWithUniqueMappedPeptideSpectra, genome2allPeptideSpectra_dic, peptideSpectrum2allGenomesWithUniqueSpectra_dic, genome_uniqueAbundanceRPKM_dic, genome2proteome_lengths_dic, total_n_peptideSpectrum_mapped_from_library)
    
    abundance_RPKM_df = pd.DataFrame(columns = ['genome', 'unique_abundance_RPKM', 'multi_abundance_RPKM', 'total_abundance_RPKM'])
    for i, genome in enumerate(genome_uniqueAbundanceRPKM_dic):
        abundance_RPKM_df.loc[i] = [genome, genome_uniqueAbundanceRPKM_dic[genome], genome_multiAbundanceRPKM_dic[genome], genome_uniqueAbundanceRPKM_dic[genome]+ genome_multiAbundanceRPKM_dic[genome]]
        
    sorted_abundance_RPKM_df = abundance_RPKM_df.sort_values(by = ['total_abundance_RPKM'], ascending = False)
    sorted_abundance_RPKM_df.to_csv(abundance_out_f, sep = '\t', index = False)
    
    genomeUniquePeptideCounts_df = pd.DataFrame(columns = ['genome', 'uniquePeptideCounts'])
    for i, genome in enumerate(genome2uniquePeptideSpectrum_counts_dic):
        genomeUniquePeptideCounts_df.loc[i] = [genome, genome2uniquePeptideSpectrum_counts_dic[genome]]
        
    genomeUniquePeptideCounts_df = genomeUniquePeptideCounts_df.sort_values(by = ['uniquePeptideCounts'], ascending = False)
    #genomeUniquePeptideCounts_df.to_csv(all_tsv_f.split('/')[0]+'/'+'.'.join(all_tsv_f.split('/')[-1].split('.')[:-1])+'_genomesWithUniquePeptides.txt', sep = '\t', index = False)
    genomeUniquePeptideCounts_df.to_csv('.'.join(all_tsv_f.split('.')[:-1]) +'_genomesWithUniquePeptides.txt', sep = '\t', index = False)
    #abundance plots
    log_abundance_RPKM_df = copy.deepcopy(sorted_abundance_RPKM_df)
    log_abundance_RPKM_df['unique_abundance_RPKM'] = list(map(get_log10, list(sorted_abundance_RPKM_df['unique_abundance_RPKM'])))
    log_abundance_RPKM_df['multi_abundance_RPKM'] = list(map(get_log10, list(sorted_abundance_RPKM_df['multi_abundance_RPKM'])))
    log_abundance_RPKM_df['total_abundance_RPKM'] = list(map(get_log10, list(sorted_abundance_RPKM_df['total_abundance_RPKM'])))
    
    plt.figure(figsize=(12,8))
    
    boxplot = log_abundance_RPKM_df.boxplot(column = ['unique_abundance_RPKM', 'multi_abundance_RPKM'])
    plt.ylabel('abundance distribution in PMKB', fontsize = 16)
    plt.title('distribution of abundances across genomes having at least one unqiue peptide mapped to them', fontsize = 18)
    plt.savefig('.'.join(all_tsv_f.split('.')[:-1])+'_unique_multi_boxplots_abundance.pdf', bbox_inches = 'tight')
    plt.figure(figsize=(12,8))
    x_range = np.arange(0, len(log_abundance_RPKM_df))
    plt.plot(x_range, log_abundance_RPKM_df['total_abundance_RPKM'], 'r')
    plt.xlabel('genomes/assembled bins with at least one unique peptide mapped to them', fontsize = 16)
    plt.ylabel('total abundance (PMKB, log10 scale)', fontsize = 16)
    plt.title('total abundance for genomes/bins that have at least one peptide mapped back to them (in PMKB, log10 scale)', fontsize = 18)
    plt.savefig('.'.join(all_tsv_f.split('.')[:-1])+'_total_abundance.pdf', bbox_inches = 'tight')