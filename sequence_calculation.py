#!/usr/bin/python

import sys

def sequence_calculation(total_no_seq, desired_seq_percent,depth,seq_error,gff_df):
    n_seq = total_no_seq*desired_seq_percent/100
    n_seq_per_chr = int(np.floor(n_seq/23))
    n_seq_Y_chr = n_seq - 23*n_seq_per_chr  # remaining sequence goes to last chromosome
    no_mir_chr = int(np.ceil(n_seq_per_chr/depth))
    print('no_mir_chr',no_mir_chr)
    no_mir_chr_Y = int(np.ceil(n_seq_Y_chr/depth))
    print('no_mir_chr_Y',no_mir_chr_Y)
    if seq_error == 'None':
        n_min_mir = 4
    else:
        n_min_mir = 2
    n_min_mir_Y = 2
    
    if no_mir_chr < n_min_mir : 
        print('Error : Minimum number of miRs per chromosome should be %d. With your provided information it is %d'%(n_min_mir, no_mir_chr))
        print('Hint: Try reducing minimum depth so that minimum number of miRs per chromosome becomes more than 4.')
        sys.exit(1)
    if gff_df.shape[0] == 0:
        warnings.warn("Warning : There are no unique RNAs left in gff file. All RNAs have been considered in earlier cases. So no new sequences can be generated for %s. Try higher value of Depth parameter" %seq_error)

    
    if no_mir_chr_Y < n_min_mir_Y : 
#         print('Error : Minimum number of miRs in chromosome Y should be %d. With your provided information it is %d'%(n_min_mir_Y, no_mir_chr_Y))
#         print('Hint: Try reducing minimum depth so that minimum number of miRs per chromosome becomes more than 4.')
        no_mir_chr_Y = n_min_mir_Y
#         sys.exit(1)
    return no_mir_chr, no_mir_chr_Y