#!/usr/bin/python

def check_sequence(rna_dict, seq):
    all_sequences = [v for _,v in rna_dict.items()]
    if seq.upper() in all_sequences:
        unique_flag = False
        return unique_flag
    else:
        unique_flag = True
        return unique_flag 