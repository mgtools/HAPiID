# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 11:54:42 2018

@author: mstambou

helper script that will later be imported in other script to read and parse CDhit outputs

this script should return a dictionary of clusters to sequence ID mappings
"""


def cdhit2dic(cdhit_out_dir):
    with open(cdhit_out_dir, 'r') as cdhit_f:
        cdhit_dic = dict()
        for line in cdhit_f:
            if line.startswith('>'):
                k = '_'.join(line[1:].split(' ')).strip('\n')
                cdhit_dic[k] = []
            else:
                line = (line.split('>')[1]).split('...')[0]
                cdhit_dic[k].append(line)
        return cdhit_dic            

def cdhit2representative2dic(cdhit_out_dir):
    with open(cdhit_out_dir, 'r') as cdhit_f:
        cdhit_dic = dict()
        for line in cdhit_f:
            line = line.strip('\n')
            if line.startswith('>'):
                k = '_'.join(line[1:].split(' ')).strip('\n')                
            elif line.endswith('*'):
                line = (line.split('>')[1]).split('...')[0]
                cdhit_dic[k] = line
        return cdhit_dic        