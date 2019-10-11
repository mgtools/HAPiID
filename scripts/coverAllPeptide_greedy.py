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

def getNextBestGenome(genome2remainingPeptide_dic):
    """
    function that will take a dictionary mappings of genome2remaining peptides and will
    retunr the genome ID of the next genome that covers the biggest number of peptides 
    remaining.
    """
    genome2nPeptide_tuple = [(k,len(v)) for k,v in genome2remainingPeptide_dic.items()]
    genome2nPeptide_tuple = sorted(genome2nPeptide_tuple, key=lambda x:x[1], reverse = True)
    bestGenome = genome2nPeptide_tuple[0]
    return bestGenome[0]

def updateGenome2remainingPeptide_dic(genome2remainingPeptide_dic, remaining_peptides, covered_peptides):
    """
    function that updates the genome to remaining peptides mappings everytime a 
    genome gets selected by the greedy approach
    """
    remaining_peptides = list(set(remaining_peptides).difference(covered_peptides))
    for genome in genome2remainingPeptide_dic:
        genome2remainingPeptide_dic[genome] = set(genome2remainingPeptide_dic[genome]).intersection(set(remaining_peptides))
    return genome2remainingPeptide_dic, remaining_peptides

if len(sys.argv) != 3:
    print('please enter 2 command line argument to run this script, i.e. example to run\n python3 coverAllPeptide_greedy.py genome/2peptide/dic_f output/file/name')

else:
    genome2peptide_dic_f = sys.argv[1]
    out_f = sys.argv[2]
    
    with open(genome2peptide_dic_f, 'r') as in_f:
        genome2peptide_dic = json.load(in_f)
        
    all_peptides = list()
    for peptides in genome2peptide_dic.values():
        all_peptides.extend(peptides)
        
    all_peptides = list(set(all_peptides))
    remaining_peptides = copy.deepcopy(all_peptides)
    genome2remainingPeptide_dic = copy.deepcopy(genome2peptide_dic)
    
    selected_genomes = list()
    covered_peptides_soFar = list()
    with open(out_f, 'w') as _out_f:
        _out_f.write('genome\tnPeptidesCovered\n')
        while(remaining_peptides):
            nextBestGenome = getNextBestGenome(genome2remainingPeptide_dic)
            selected_genomes.append(nextBestGenome)
            covered_peptides = genome2remainingPeptide_dic[nextBestGenome]
            covered_peptides_soFar.extend(covered_peptides)
            covered_peptides_soFar = list(set(covered_peptides_soFar))
            genome2remainingPeptide_dic, remaining_peptides = updateGenome2remainingPeptide_dic(genome2remainingPeptide_dic, remaining_peptides, covered_peptides)
            _out_f.write(nextBestGenome+'\t'+str(len(covered_peptides_soFar))+'\n')
    
       