# -*- coding: utf-8 -*-
"""
Created on Mon May 20 13:50:51 2019

@author: mstambou

script in which I will keep genomes that are able to explain all the peptides,
that were first recovered from searching using the highly expressed genes.
Here I will take the greedy approach that basically will first get the genome
that has the most peptides mapped back to it, and then remove it, then get the 
genome that covers the next most peptides after removing peptides that were 
previously covered.
"""

import json
import sys
import copy
import pandas as pd

def getNextBestGenome(genome2remainingSpectrum_dic):
    """
    function that will take a dictionary mappings of genome2remaining peptides and will
    retunr the genome ID of the next genome that covers the biggest number of peptides 
    remaining.
    """
    genome2nSpectrum_tuple = [(k,len(v)) for k,v in genome2remainingSpectrum_dic.items()]
    genome2nSpectrum_tuple = sorted(genome2nSpectrum_tuple, key=lambda x:x[1], reverse = True)
    bestGenome = genome2nSpectrum_tuple[0]
    return bestGenome[0]

def updateGenome2remainingSpectrum_dic(genome2remainingSpectrum_dic, remaining_spectra, covered_spectra):
    """
    function that updates the genome to remaining spectra mappings everytime a 
    genome gets selected by the greedy approach
    """
    remaining_spectra = list(set(remaining_spectra).difference(covered_spectra))
    for genome in genome2remainingSpectrum_dic:
        genome2remainingSpectrum_dic[genome] = set(genome2remainingSpectrum_dic[genome]).intersection(set(remaining_spectra))
    return genome2remainingSpectrum_dic, remaining_spectra

if len(sys.argv) != 3:
    print('please enter 2 command line argument to run this script, i.e. example to run\n python3 coverAllPeptide_greedy.py genome/2peptide/dic_f output/file/name')

else:
    genome2spectrum_dic_f = sys.argv[1]
    out_f = sys.argv[2]
    
    with open(genome2spectrum_dic_f, 'r') as in_f:
        genome2spectrum_dic = json.load(in_f)
        
    all_spectra = list()
    for spectra in genome2spectrum_dic.values():
        all_spectra.extend(spectra)
        
    all_spectra = list(set(all_spectra))
    remaining_spectra = copy.deepcopy(all_spectra)
    genome2remainingSpectrum_dic = copy.deepcopy(genome2spectrum_dic)
    
    selected_genomes = list()
    covered_spectra_soFar = list()
    with open(out_f, 'w') as _out_f:
        _out_f.write('genome\tnSpectraCovered\n')
        while(remaining_spectra):
            nextBestGenome = getNextBestGenome(genome2remainingSpectrum_dic)
            selected_genomes.append(nextBestGenome)
            covered_spectra = genome2remainingSpectrum_dic[nextBestGenome]
            covered_spectra_soFar.extend(covered_spectra)
            covered_spectra_soFar = list(set(covered_spectra_soFar))
            genome2remainingSpectrum_dic, remaining_spectra = updateGenome2remainingSpectrum_dic(genome2remainingSpectrum_dic, remaining_spectra, covered_spectra)
            _out_f.write(nextBestGenome+'\t'+str(len(covered_spectra_soFar))+'\n')
    
    covered_spectra_df = pd.read_csv(out_f, sep = '\t')
    covered_spectra_df['nSpectra_%'] = [(item*100)/list(covered_spectra_df['nSpectraCovered'])[-1] for item in covered_spectra_df['nSpectraCovered']]
    covered_spectra_df.to_csv(out_f, sep = '\t', index = False)
