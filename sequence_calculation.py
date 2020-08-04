#!/usr/bin/python

import sys
import numpy as np
from expression_split import *

def sequence_calculation(total_no_seq, desired_seq_percent,depth,seed,gff_df):
    n_seq = total_no_seq*desired_seq_percent/100
    n_seq_per_chr = expression_split(n_seq,24,'uniform',seed,depth)
    no_mir_chr = []
    chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
               'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
               'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    actual_no_mir_chr = [len(list(gff_df[gff_df['chr'].str.contains(chr_name+'$')].index)) for chr_name in chr_list]
        
    for v1,v2 in zip(n_seq_per_chr,actual_no_mir_chr):
        no_mir = len(expression_split(v1,v2,'poisson',seed,depth))
        if no_mir >= 1:
            no_mir_chr.append(no_mir)
        else:
            print('Error : Minimum number of miRs per chromosome should be %d. With your provided information it is %d'%(1, no_mir))
            print('Hint: Try reducing minimum depth so that minimum number of miRs per chromosome becomes more than 1.')
            sys.exit(1)
    return no_mir_chr,n_seq_per_chr