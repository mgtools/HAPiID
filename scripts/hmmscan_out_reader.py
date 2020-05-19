# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 15:45:59 2018

@author: mstambou

helper script that will be imported by other scripts in order to facilitate reading the hmmscan file.
because the HMMSCAN file did not have a good structure and was not well structured I am jsut extracting
the accession ID for the PFAM HMM model and the query fasta ID for the sequence that is mapped to it
from each line. some of these fasta headers were bigger than the actual space alotted by the HMMSCAN
hence it was going over the other columns if I used predifined spaces as the delimiter as defined in 
the 3rd line of the HMMSCAN output file. I therefore only extract a data fram that maps PFAM accession
to query fasta ID. This dataframe could be used as a mapper to map the sequences to PFAM domains.
"""

import pandas as pd
import re

def pfam_2_query(hmmscan_out_f):    
    with open(hmmscan_out_f) as in_f:
        line_count = 0
        accession_query_df = pd.DataFrame(columns = ['accession', 'query'])
        for line in in_f:
            if not line.startswith('#'):
               line = list(filter(None, line.rstrip('\n').split(' ')))
               accession_query_df.loc[line_count] = [line[1], line[2]]
               line_count += 1    
    return accession_query_df

def pfam_2_query_2_description(hmmscan_out_f):    
    with open(hmmscan_out_f) as in_f:
        line_count = 0
        accession_query_df = pd.DataFrame(columns = ['accession', 'query', 'description'])
        for line in in_f:
            if not line.startswith('#'):
               line = list(filter(None, line.rstrip('\n').split(' ')))               
               accession_query_df.loc[line_count] = [line[1], line[2], ' '.join(line[18:])]
               line_count += 1    
    return accession_query_df

def pfam_2_query_old(hmmscan_out_f):    
    #hmmscan_out_f = 'extract_marker_genes/HMMSCAN_out/D01.hmm'
    with open(hmmscan_out_f) as in_f:
        lines = in_f.readlines()
        header = lines[2]
        del lines
    header = header.split('   ')
    header_1 = header[0]    
    space_indices_1 = [match.start() for match in re.finditer(re.escape(' '), header_1)]
    space_indices_1.insert(0,0)
    space_indices_1.append(space_indices_1[-1]+6)
    hmmscan_colspecs_1 = [[space_indices_1[i] ,space_indices_1[i+1]] for i in range(len(space_indices_1)-1)]
    hmmscan_colspecs_1.insert(2, [hmmscan_colspecs_1.pop(2)[0], hmmscan_colspecs_1.pop(2)[1]])
    accession_query_name_cols = hmmscan_colspecs_1[1:3]    
    accession_query_df = pd.read_fwf(hmmscan_out_f, header = None, colspecs = accession_query_name_cols, comment = '#')
    accession_query_df = accession_query_df[3:-10]
    accession_query_df.columns = ['accession', 'query']
    def get_query_only(s):
        return s.split(' ')[0]
    
    accession_query_df['query'] = list(map(get_query_only, accession_query_df['query']))
    
    return accession_query_df

def accession_2_targetName(hmmscan_out_f):
    
    #hmmscan_out_f = 'extract_marker_genes/HMMSCAN_out/D01.hmm'
    with open(hmmscan_out_f) as in_f:
        lines = in_f.readlines()
        header = lines[2]
        del lines
    header = header.split('   ')
    header_1 = header[0]    
    space_indices_1 = [match.start() for match in re.finditer(re.escape(' '), header_1)]
    space_indices_1.insert(0,0)
    space_indices_1.append(space_indices_1[-1]+6)
    hmmscan_colspecs_1 = [[space_indices_1[i] ,space_indices_1[i+1]] for i in range(len(space_indices_1)-1)]
    #hmmscan_colspecs_1.insert(2, [hmmscan_colspecs_1.pop(2)[0], hmmscan_colspecs_1.pop(2)[1]])
    accession_query_name_cols = hmmscan_colspecs_1[0:2]    
    targetName_2_accession_query_df = pd.read_fwf(hmmscan_out_f, header = None, colspecs = accession_query_name_cols, comment = '#')
    targetName_2_accession_query_df = targetName_2_accession_query_df [3:-10]
    targetName_2_accession_query_df.columns = ['target name', 'accession']
    targetName_2_accession_query_df = targetName_2_accession_query_df .drop_duplicates()
    accession_2_targetName_dic = dict(zip(targetName_2_accession_query_df['accession'], targetName_2_accession_query_df['target name']))
    

    return accession_2_targetName_dic
