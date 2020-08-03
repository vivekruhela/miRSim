#!/usr/bin/python

def cigar_generation(actual_seq,mod_seq):
    seq1 = list(actual_seq)
    seq2 = list(mod_seq)
    cigar = []
    if len(seq1) == len(seq2):
        for idx in range(len(seq1)):
            if seq1[idx] == seq2[idx]:
                cigar.append('-')
            else:
                cigar.append(seq2[idx])
    else:
        cigar = ['-' for l in range(len(seq1))]
        
    return ''.join(cigar).upper()