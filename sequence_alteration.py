#!/usr/bin/python

import random
from alter_nt import *

def sequence_alteration(mir_seq,site,seed):
    random.seed(seed)
    if site == 'seed':
        seq_seed = mir_seq[1:7]
        seq_seed = list(seq_seed)
        idx = list(range(1,7))
        x1 = random.randint(0, len(idx)-1)
        seq_seed[x1] = alter_nt(seq_seed,x1)
        idx.pop(x1)
        x2 = random.randint(0, len(idx)-1)
        seq_seed[x2] = alter_nt(seq_seed,x2)
        new_seed = ""
        new_seed = new_seed.join(seq_seed)
        mir_seq_new = mir_seq[:1] + new_seed + mir_seq[7:]
    elif site == 'xseed':
        Xseq_seed = mir_seq[7:]
        Xseq_seed = list(Xseq_seed)
        idx = list(range(0,len(Xseq_seed)))
        x3 = random.randint(0, len(idx)-1)
        Xseq_seed[x3] = alter_nt(Xseq_seed,x3)
        idx.pop(x3)
        x4 = random.randint(0, len(idx)-1)
        Xseq_seed[x4] = alter_nt(Xseq_seed,x4)
        new_xseed = ""
        new_xseed = new_xseed.join(Xseq_seed)
        mir_seq_new = mir_seq[:7] + new_xseed
    elif site == 'both':
        mir_seq1 = sequence_alteration(mir_seq,'seed',seed)
        mir_seq_new = sequence_alteration(mir_seq1,'xseed',seed)
    
    return mir_seq_new