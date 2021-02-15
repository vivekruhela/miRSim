#!/usr/bin/python

import random
from alter_nt import *

def sequence_alteration(mir_seq,site,no_mismatch_seed,no_mismatch_xseed, seed):
    if isinstance(seed, int):
        random.seed(seed)
    if site == 'seed':
        seq_seed = mir_seq[1:7]
        seq_seed = list(seq_seed)
        idx = list(range(1,7))
        for i in range(no_mismatch_seed):            
            x1 = random.randint(0, len(idx)-1)
            seq_seed[x1] = alter_nt(seq_seed,x1)
            idx.pop(x1)
            new_seed = ""
            new_seed = new_seed.join(seq_seed)
        mir_seq_new = mir_seq[:1] + new_seed + mir_seq[7:]
        
    elif site == 'xseed':
        Xseq_seed = mir_seq[7:]
        Xseq_seed = list(Xseq_seed)
        idx = list(range(0,len(Xseq_seed)))
        for  j in range(no_mismatch_xseed):            
            x3 = random.randint(0, len(idx)-1)
            Xseq_seed[x3] = alter_nt(Xseq_seed,x3)
            idx.pop(x3)
            new_xseed = ""
            new_xseed = new_xseed.join(Xseq_seed)
        mir_seq_new = mir_seq[:7] + new_xseed
        
    elif site == 'both':
        no_mismatch_both = no_mismatch_seed + no_mismatch_xseed
        for k in range(no_mismatch_both):
            mir_seq1 = sequence_alteration(mir_seq,'seed',no_mismatch_seed , no_mismatch_xseed, seed)
            mir_seq_new = sequence_alteration(mir_seq1,'xseed',no_mismatch_seed , no_mismatch_xseed, seed)
    
    return mir_seq_new