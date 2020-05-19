# -*- coding: utf-8 -*-
"""
Created on Thu May 16 15:52:55 2019

@author: mstambou

eventually this script will produce a matrix to be an input to the integer linear programming 
it will get the tsv after FDR filtering and then extract all the proteins to peptide mappings
expand on them using the cd-hit mappings and then get the expanded protein to peptide mappings
and eventually combine the proteins into genomes and get genome to peptide mappings
this also produces genome to peptide mapping dictionaries
"""

import utils
import pandas as pd
import sys
import numpy as np
import json
import os



if len(sys.argv) != 5:
    print('please enter 5 command line arguments to run this script, i.e. example to run\n python3 get_basicPeptideSpectra2Portein2genome_stats.py 131211-28623-ATH-F_MSGF+/131211-28623-12-ATH-F01-03.tsv.0.01.tsv 131211-28623-ATH-F_MSGF+/ribP_elonF_all_cdHit_100.fasta.clstr umgs_hgg_proteins2genomes_dic.json dir/to/write/outputs')

else:
    all_tsv_f = sys.argv[1]
    #all_tsv_f = '131211-28623-ATH-F_MSGF+/131211-28623-12-ATH-F01-03.tsv.0.01.tsv'
    sample_name = (all_tsv_f.rsplit('/', 1)[-1]).split('.')[0]
    
    cdhit_out_f = sys.argv[2]
    #cdhit_out_f = '131211-28623-ATH-F_MSGF+/ribP_elonF_all_cdHit_100.fasta.clstr'
    
    proteins2genomes_dic_dir = sys.argv[3]
    #proteins2genomes_dic_dir = 'umgs_hgg_proteins2genomes_dic.json'
    
    out_dir = sys.argv[4]
    
    def getGenome2Spectrum_mat(protein2spectrum_expanded_dic, umgs_hgg_proteins2genomes_dic, out_dir, tsv_in = all_tsv_f):
        """
        function that takes in the expanded dictionary mapping from proteins to spectra and 
        also a dictionary mapping of all the protein IDs to their respective genome/bin IDs
        and returns a binary matrix whose rows are all the spectra that are indentified from
        the input sample and columns are all the genomes that have at least one of these spectra
        the spectra/genome intersection in the matrix will be a 1 denoting that that spectra was 
        associated with that genome otherwise it will be a 0
        """
        genome2spectrum_dic = dict()
        all_spectra = list()
        for protein in protein2spectrum_expanded_dic:
            genome = utils.extract_genome_from_prot(protein, umgs_hgg_proteins2genomes_dic)
            spectra = protein2spectrum_expanded_dic[protein]
            all_spectra.extend(spectra)
            if genome in genome2spectrum_dic:
                genome2spectrum_dic[genome].extend(spectra)
            elif genome not in genome2spectrum_dic:
                genome2spectrum_dic[genome] = spectra
                
        genome2spectrum_dic = {k:list(set(v)) for k,v in genome2spectrum_dic.items()}
        all_spectra = list(set(all_spectra))
        spectrum2idx_dic, idx2spectrum_dic = dict(), dict()
        genome2idx_dic, idx2genome_dic = dict(), dict()
        
        spectrum2genomeMat = np.ndarray((len(all_spectra), len(genome2spectrum_dic)))
        
        for i, genome in enumerate(genome2spectrum_dic):
            genome2idx_dic[genome] = i
            idx2genome_dic[i] = genome
        
        for i, spectrum in enumerate(all_spectra):
            spectrum2idx_dic[spectrum] = i
            idx2spectrum_dic[i] = spectrum
        
        for i in range(len(all_spectra)): #loop over all spectra
            for j in  range(len(genome2spectrum_dic)): #loop over all genomes
                if idx2spectrum_dic[i] in genome2spectrum_dic[idx2genome_dic[j]]:
                    spectrum2genomeMat[i][j] = 1
                    
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
                    
        out_f = (tsv_in.rsplit('/', 1)[-1]).rsplit('.')[0]+'_spectrum2genome_mat'
        
        np.save(out_dir+out_f, spectrum2genomeMat)
        
        with open(out_dir+sample_name+'_genome2matIdx_spectra.json', 'w') as out_f:
            json.dump(genome2idx_dic, out_f)
        with open(out_dir+sample_name+'_matIdx2genome_spectra.json', 'w') as out_f:
            json.dump(idx2genome_dic, out_f)
        with open(out_dir+sample_name+'_spectrum2matIdx_dic.json', 'w') as out_f:
            json.dump(spectrum2idx_dic, out_f)
        with open(out_dir+sample_name+'_matIdx2spectrum_dic.json', 'w') as out_f:
            json.dump(idx2spectrum_dic, out_f)
        with open(out_dir+sample_name+'_genome2spectrum_dic.json', 'w') as out_f:
            json.dump(genome2spectrum_dic, out_f)
        
    umgs_hgg_proteins2genomes_dic = utils.readProteins2genomes_dic(proteins2genomes_dic_dir)
    
    
    all_tsv_df = pd.read_csv(all_tsv_f, sep = '\t')
    #this step makes sure if there were multiple rows of the same spectrum ID then keep the first one (i.e. the best one)
    all_tsv_df = all_tsv_df.drop_duplicates(subset = 'SpecID', keep = 'first')
    
    unique_spectra, unique_proteins, protein2spectrum_dic, spectrum2protein_dic = utils.get_spectrum_proteins_mappings(all_tsv_df)
    #in expanding the cd-hit clusters I can reuse the same function in the utlis, it's basically doing the same thing.
    protein2spectrum_expanded_dic, spectrum2protein_expanded_dic = utils.expand_cdhit_clustered_prots(cdhit_out_f, protein2spectrum_dic, spectrum2protein_dic)
    
    getGenome2Spectrum_mat(protein2spectrum_expanded_dic, umgs_hgg_proteins2genomes_dic, out_dir)
