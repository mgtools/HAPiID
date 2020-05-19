# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 11:30:27 2019

@author: mstambou

script containing usefull functions that will be used in other scritps to perform quantification
"""


import pandas as pd
import cdhit_out_reader as cdhit
import copy
import json 

def readProteins2genomes_dic(proteins2genomes_dic_dir):
    with open(proteins2genomes_dic_dir, 'r') as in_f:
        umgs_hgg_proteins2genomes_dic = json.load(in_f)
        return umgs_hgg_proteins2genomes_dic
    
def extract_genome_from_prot(prot, umgs_hgg_proteins2genomes_dic):      
    print(prot)
    if prot in umgs_hgg_proteins2genomes_dic:
        genome = umgs_hgg_proteins2genomes_dic[prot]    
    elif ':' in prot:
        genome = prot.split(':')[0]        
    elif 'UMGS' in prot:
        genome = prot.split('_')[0]
        
    return genome               

def get_peptideSpectra_proteins_mappings(all_tsv_df):
    """
    function that takes the initial tab separated file for peptide identification and 
    returns all the unique proteins, all the unique peptides, protein2peptideSpectrm mappings
    and Spectra2proteins mappings. In this case I'm using spectra as a proxy perform quantificatin
    i.e. even if the same peptide is matched to the same protein twice but if they are coming from
    different spectrum then you treat these as two different matches and count them twice.
    """
    unique_peptides = list(set(all_tsv_df['Peptide']))
    peptideSpectrum2protein_dic = dict()
    protein2peptideSpectrum_dic = dict()
    proteins_n_mapped = 0
    unique_proteins = []
    for i, row in all_tsv_df.iterrows():
        item = row['Protein']
        peptideSpectrum = row['SpecID']
        for protein in item.split(';'):
            proteins_n_mapped += 1
            protein = protein.replace('XXX_', '').replace('XXX', '')
            if '(pre' in protein:
                protein = protein.split('(pre')[0]
            unique_proteins.append(protein)
            if protein in protein2peptideSpectrum_dic:
                protein2peptideSpectrum_dic[protein].append(peptideSpectrum)
            else:
                protein2peptideSpectrum_dic[protein] = [peptideSpectrum]
            if peptideSpectrum in peptideSpectrum2protein_dic:
                peptideSpectrum2protein_dic[peptideSpectrum].append(protein)
            else:
                peptideSpectrum2protein_dic[peptideSpectrum] = [protein]
    
    unique_proteins = list(set(unique_proteins))
    protein2peptideSpectrum_dic = {k:list(set(v)) for k,v in protein2peptideSpectrum_dic.items()}
    peptideSpectrum2protein_dic = {k:list(set(v)) for k,v in peptideSpectrum2protein_dic.items()}
    return unique_peptides, unique_proteins, protein2peptideSpectrum_dic, peptideSpectrum2protein_dic



def expand_cdhit_clustered_prots(cdhit_out_f, protein2peptideSpectrum_dic, peptideSpectrum2protein_dic):
    """
    function that expands on the list of the representative proteins that have peptide spectrum hits
    for each representative protein with peptide Spectrum hits, it's information i.e. number of peptide Spectra 
    mapped to the proteins and then peptide spectra 2 protein mappings are updated by adding all the 
    proteins within the cd-hit cluster.
    """
    cdhit_out = cdhit.cdhit2dic(cdhit_out_f)
    cdhit_2repr = cdhit.cdhit2representative2dic(cdhit_out_f)
    cdhit_repr2members = {cdhit_2repr[k]:cdhit_out[k] for k in cdhit_2repr}
    
    protein2peptideSpectrum_expanded_dic = copy.deepcopy(protein2peptideSpectrum_dic)
    peptideSpectrum2protein_expanded_dic = copy.deepcopy(peptideSpectrum2protein_dic)
    
    for protein in protein2peptideSpectrum_dic:
        protein_cluster = cdhit_repr2members[protein]
        if len(protein_cluster) > 1:        
            protein_cluster.remove(protein)
            for member in protein_cluster:        #here I am updating the old dictionaries and expanding on them    
                member_peptideSpectra = protein2peptideSpectrum_dic[protein]
                protein2peptideSpectrum_expanded_dic[member] = member_peptideSpectra
                for peptideSpectrum in member_peptideSpectra:
                    peptideSpectrum2protein_expanded_dic[peptideSpectrum].append(member)
                
    protein2peptideSpectrum_expanded_dic = {k:list(set(v)) for k,v in protein2peptideSpectrum_expanded_dic.items()}
    peptideSpectrum2protein_expanded_dic = {k:list(set(v)) for k,v in peptideSpectrum2protein_expanded_dic.items()}
    
    return protein2peptideSpectrum_expanded_dic, peptideSpectrum2protein_expanded_dic

def get_proteinsWithUniquePeptideSpectrum(peptideSpectrum2protein_expanded_dic, umgs_hgg_proteins2genomes_dic):
    """
    function returning proteins that have at least one peptide Spectrum that's unique to them
    """
    unique_mapped_proteins = list()
    for peptideSpectrum in peptideSpectrum2protein_expanded_dic:
        if len(peptideSpectrum2protein_expanded_dic[peptideSpectrum]) == 1:
            unique_mapped_proteins.append(peptideSpectrum2protein_expanded_dic[peptideSpectrum][0])
    return list(set(unique_mapped_proteins))

def get_genomesWithUniqueMappedProts(unique_mapped_proteins, umgs_hgg_proteins2genomes_dic):
    """
    returns a the list of genome names containing at least one protein having unique peptides
    """
    return list(set([extract_genome_from_prot(prot, umgs_hgg_proteins2genomes_dic) for prot in unique_mapped_proteins]))
        
def get_genomes_with_unique_peptides(peptideSpectrum2protein_expanded_dic, umgs_hgg_proteins2genomes_dic):
    """
    function that goes over every peptide and uses peptide to protein dictionaries
    for each peptide to protein maping list it sees whether or not if the genomes that all these proteins are coming 
    from is the same genome (i.e. coming from one genome or more) different genomes, and accordingly returns a list 
    of genomes that have at least one unique peptide.
    """
    genomesWithUniqueMappedPeptideSpectrum = list()
    uniqueSpectrum2genome = dict()
    uniqueGenome2Spectrum = dict()
    for peptideSpectrum in  peptideSpectrum2protein_expanded_dic:
        peptideSpectrumProteins = peptideSpectrum2protein_expanded_dic[peptideSpectrum]
        
        genomes = [extract_genome_from_prot(prot, umgs_hgg_proteins2genomes_dic) for prot in peptideSpectrumProteins]
        if len(set(genomes)) == 1:
            genomesWithUniqueMappedPeptideSpectrum.append(genomes[0])
            uniqueSpectrum2genome[peptideSpectrum] = genomes[0]
            uniqueGenome2Spectrum[genomes[0]] = peptideSpectrum
    genomesWithUniqueMappedPeptideSpectrum = list(set(genomesWithUniqueMappedPeptideSpectrum))
    return genomesWithUniqueMappedPeptideSpectrum

def get_peptideSpectrum2genomeWithUniqueSpectra_maps(peptideSpectrum2protein_expanded_dic, genomesWithUniqueMappedPeptideSpectra, umgs_hgg_proteins2genomes_dic):
    """
    function that returns a mapping between peptides and repsective genomes that they get mapped to.
    This is done through peptide to protein mappings
    """
    genomesWithUniqueMappedPeptideSpectra = set(genomesWithUniqueMappedPeptideSpectra)
    peptideSpectrum2genome_expanded_dic = dict()
    for peptideSpectrum in peptideSpectrum2protein_expanded_dic:
        peptideSpectrum_proteins = peptideSpectrum2protein_expanded_dic[peptideSpectrum]
        peptideSpectrum_genomes = [extract_genome_from_prot(peptideSpectrum_protein, umgs_hgg_proteins2genomes_dic) for peptideSpectrum_protein in peptideSpectrum_proteins]
        peptideSpectrum_genomesWithUniqueSpectra = list(set(peptideSpectrum_genomes).intersection(genomesWithUniqueMappedPeptideSpectra))
        peptideSpectrum2genome_expanded_dic[peptideSpectrum] = list(set(peptideSpectrum_genomesWithUniqueSpectra))
    return peptideSpectrum2genome_expanded_dic

def get_peptide_proteins_mappings(all_tsv_df):
    """
    function that takes the initial tab separated file for peptide identification and 
    returns all the unique proteins, all the unique peptides, protein2peptide mappings
    and peptide2proteins mappings. 
    """
    unique_peptides = list(set(all_tsv_df['Peptide']))
    peptide2protein_dic = dict()
    protein2peptide_dic = dict()
    proteins_n_mapped = 0
    unique_proteins = []
    for i, row in all_tsv_df.iterrows():
        item = row['Protein']
        peptide = row['Peptide']
        for protein in item.split(';'):
            proteins_n_mapped += 1
            protein = protein.replace('XXX_', '').replace('XXX', '')
            if '(pre' in protein:
                protein = protein.split('(pre')[0]
            unique_proteins.append(protein)
            if protein in protein2peptide_dic:
                protein2peptide_dic[protein].append(peptide)
            else:
                protein2peptide_dic[protein] = [peptide]
            if peptide in peptide2protein_dic:
                peptide2protein_dic[peptide].append(protein)
            else:
                peptide2protein_dic[peptide] = [protein]
    
    unique_proteins = list(set(unique_proteins))
    protein2peptide_dic = {k:list(set(v)) for k,v in protein2peptide_dic.items()}
    peptide2protein_dic = {k:list(set(v)) for k,v in peptide2protein_dic.items()}
    return unique_peptides, unique_proteins, protein2peptide_dic, peptide2protein_dic


def get_spectrum_proteins_mappings(all_tsv_df):
    """                                                                                                                                                                         
    function that takes the initial tab separated file for peptide identification and                                                                                           
    returns all the unique proteins, all the unique spectra, protein2spectra mappings                                                                                          
    and spectra2proteins mappings.                                                                                                                                              
    """
    unique_spectra = list(set(all_tsv_df['SpecID']))
    spectrum2protein_dic = dict()
    protein2spectrum_dic = dict()
    proteins_n_mapped = 0
    unique_proteins = []
    for i, row in all_tsv_df.iterrows():
        item = row['Protein']
        spectrum = row['SpecID']
        for protein in item.split(';'):
            proteins_n_mapped += 1
            protein = protein.replace('XXX_', '').replace('XXX', '')
            if '(pre' in protein:
                protein = protein.split('(pre')[0]
            unique_proteins.append(protein)
            if protein in protein2spectrum_dic:
                protein2spectrum_dic[protein].append(spectrum)
            else:
                protein2spectrum_dic[protein] = [spectrum]
            if spectrum in spectrum2protein_dic:
                spectrum2protein_dic[spectrum].append(protein)
            else:
                spectrum2protein_dic[spectrum] = [protein]

    unique_proteins = list(set(unique_proteins))
    protein2spectrum_dic = {k:list(set(v)) for k,v in protein2spectrum_dic.items()}
    spectrum2protein_dic = {k:list(set(v)) for k,v in spectrum2protein_dic.items()}
    return unique_spectra, unique_proteins, protein2spectrum_dic, spectrum2protein_dic
