#!/usr/bin/python

import random
import os
from expression_split import *
from cigar_generation import *
from mir_location import *
from sequence_alteration import *

def generate_sequence(fasta_seq, gff_df, rna_dict, no_mir_chr, n_seq_per_chr, depth, seq_error, out, out_file,write_mode,repeat,distribution,seed):
    
    if isinstance(seed, int):
        random.seed(seed)
    if write_mode == 'write':
        rna_ground_truth = open(os.path.join(out,out_file),'w')
        header = 'RNA_ID\tImpure_Region\tCigar_String\tref_Sequence\tsynthetic_sequence\tchr\tchr_start\tchr_end\tExpression_count\n'
        rna_ground_truth.write(header)
    elif write_mode == 'append':
        rna_ground_truth = open(os.path.join(out,out_file),'a')
    
    chr_list = list(gff_df['chr'].unique())
    if no_mir_chr and n_seq_per_chr:
        for i in range(len(chr_list)):
            chr_name = chr_list[i] + '$'
            mir_complete_list = list(gff_df[gff_df['chr'].str.contains(chr_name)].index)   
            if mir_complete_list:
                total_exp = n_seq_per_chr[i]
                if not no_mir_chr[i] > len(mir_complete_list):                    
                    mir_idx = [random.randint(0,len(mir_complete_list)-1) for v in range(int(no_mir_chr[i]))]
                    mir_list = [mir_complete_list[v] for v in mir_idx]
                else:                
                    if repeat:
                        complete_flag = True
                        while complete_flag == True:
                            mir_complete_list += mir_complete_list
                            if no_mir_chr[i] < len(mir_complete_list):
                                complete_flag = False
                        mir_idx = [random.randint(0,len(mir_complete_list)-1) for v in range(int(no_mir_chr[i]))]
                        mir_list = [mir_complete_list[v] for v in mir_idx]
                    else:                    
                        mir_idx = [random.randint(0,len(mir_complete_list)-1) for v in range(int(mir_complete_list))]
                        mir_list = [mir_complete_list[v] for v in mir_idx]

                complete_flag = True
                while complete_flag:
                    try:
                        expression_counts = expression_split(total_exp,len(mir_list),distribution,seed,depth)
                        if expression_counts:
                            complete_flag = False
                    except:
                        print('The number of RNAs and total available depth in %s  is %d and %d too less with error %s' %(chr_name,len(mir_list),total_exp,seq_error))
                        expression_counts = expression_split(total_exp,len(mir_list),distribution,seed,int(depth*0.5))
                        if expression_counts:
                            complete_flag = False

                for mir,exp in zip(mir_list,expression_counts):
                    if seq_error == 'None':
                        mir_seq_new = rna_dict[mir]
                        mir_cigar = cigar_generation(rna_dict[mir],mir_seq_new)
                    elif seq_error == 'Seed_region':
                        mir_seq_new = sequence_alteration(rna_dict[mir],'seed',seed)
                        mir_cigar = cigar_generation(rna_dict[mir],mir_seq_new)
                    elif seq_error == 'Outside_Seed_region':
                        mir_seq_new = sequence_alteration(rna_dict[mir],'xseed',seed)
                        mir_cigar = cigar_generation(rna_dict[mir],mir_seq_new)
                    elif seq_error == 'Both_region':
                        mir_seq_new = sequence_alteration(rna_dict[mir],'both',seed)
                        mir_cigar = cigar_generation(rna_dict[mir],mir_seq_new)

                    loc = mir_location(gff_df,mir,seed)
                    mir_depth = int(exp)
                    line = ''
                    line += mir + '\t' + seq_error + '\t' + mir_cigar + '\t' + rna_dict[mir] + '\t' + mir_seq_new + '\t' + loc[0] + '\t' + str(loc[1]) + '\t' + str(loc[2]) + '\t' + str(mir_depth) + '\n'
                    rna_ground_truth.write(line)            
                    for dep in range(mir_depth):
                        fasta_seq.append('>' + mir)
                        fasta_seq.append(mir_seq_new)

    return fasta_seq