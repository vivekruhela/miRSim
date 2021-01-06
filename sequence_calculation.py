#!/usr/bin/python

import sys
import numpy as np
import pandas as pd
from expression_split import *
from sort_by_chromosome import *

def sequence_calculation(total_no_seq, desired_seq_percent,depth,seed,gff_df):
    
#     chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
#                'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
#                'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY','chrMT']
    chr_list = sort_by_chromosome(list(gff_df['chr'].unique()))
    if not desired_seq_percent == 0:
        n_seq = total_no_seq*desired_seq_percent/100
        chr_count = pd.DataFrame(pd.Series(gff_df['chr']).value_counts())
        chr_count = pd.DataFrame(chr_count, index=chr_list)
        n_seq_per_chr = [n_seq/len(chr_list) for i in list(chr_count['chr']/gff_df.shape[0])]
        n_seq_per_chr = list(map(int,n_seq_per_chr))
        n_seq_per_chr[-1] = int(n_seq - sum(n_seq_per_chr[:-1]))

        no_mir_chr = []    
        actual_no_mir_chr = [len(list(gff_df[gff_df['chr'].str.contains(chr_name+'$')].index)) for chr_name in chr_list]
        for v1,v2 in zip(n_seq_per_chr,actual_no_mir_chr):
            if v1 > depth:
                no_mir = len(expression_split(v1,v2,'poisson',seed,depth))
                if no_mir >= 1:
                    no_mir_chr.append(no_mir)
                else:
                    print('Error : Minimum number of miRs per chromosome should be %d. With your provided information it is %d'%(1, no_mir))
                    print('Hint: Try reducing minimum depth so that minimum number of miRs per chromosome becomes more than 1.')
                    sys.exit(1)
            else:
                no_mir = 1
                no_mir_chr.append(no_mir)
    else:
        no_mir_chr = []
        n_seq_per_chr = []
    return no_mir_chr,n_seq_per_chr