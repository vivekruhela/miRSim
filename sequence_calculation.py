#!/usr/bin/python

import sys
import numpy as np
from expression_split import *

def sequence_calculation(total_no_seq, desired_seq_percent,depth,seed):
    n_seq = total_no_seq*desired_seq_percent/100
    n_seq_per_chr = expression_split(n_seq,24,'uniform',seed,depth)
    no_mir_chr = [] 
    for v in n_seq_per_chr:
        no_mir = int(np.ceil(v/depth))
        if no_mir >= 1:
            no_mir_chr.append(no_mir)
        else:
            print('Error : Minimum number of miRs per chromosome should be %d. With your provided information it is %d'%(1, no_mir))
            print('Hint: Try reducing minimum depth so that minimum number of miRs per chromosome becomes more than 1.')
            sys.exit(1)
    return no_mir_chr