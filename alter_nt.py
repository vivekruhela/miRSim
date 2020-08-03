#!/usr/bin/python

def alter_nt(seq_seed,position):
    if seq_seed[position] == 'A':
        seq_seed[position] = 'C'
    elif seq_seed[position] == 'T':
        seq_seed[position] = 'G'
    elif seq_seed[position] == 'G':
        seq_seed[position] = 'A'
    else:
        seq_seed[position] = 'T'
        
    return seq_seed[position]