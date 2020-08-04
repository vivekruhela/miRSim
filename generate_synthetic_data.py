#!/usr/bin/python

import time
import threaded
from gff import *
from generate_sequence import *
from write_fastq import *
from sequence_calculation import *


def generate_synthetic_data(fa_file,total_no_seq,pure_seq_percent,seed_error_percent,xseed_error_percent,
                      both_error_percent,depth,ascii_base,out,out_file_name,ground_truth_filename,
                      gff_file,mir_type,out_file_type,repeat,distribution,seed,adaptor,parallel_thread):
    start = time.time()
    file_read = open(fa_file).read().split('\n')
    file_read = file_read[:-1]
    rna_dict = {}
    for line in file_read:
        if ">" in line:
            line_no = file_read.index(line)
            rna_dict[line[1:]] = file_read[line_no+1]

    fasta_seq = []
    # for pure sequences only
    gff_df = gff(gff_file,mir_type)
    no_mir_chr = sequence_calculation(total_no_seq, pure_seq_percent, depth,seed)    
    fasta_seq = generate_sequence(fasta_seq,gff_df, rna_dict, no_mir_chr, depth, 'None', out, ground_truth_filename,'write',repeat,distribution,seed)
    print('Pure Sequence are generated.')

    # for seed region error sequences only
    if not repeat:
        gff_df = update_gff(gff_df,out,ground_truth_filename)  # Update gff dataframe so that there is no repitition.
    no_mir_chr = sequence_calculation(total_no_seq, seed_error_percent, depth,seed)    
    fasta_seq = generate_sequence(fasta_seq,gff_df, rna_dict, no_mir_chr, depth, 'Seed_region', out, ground_truth_filename,'append',repeat,distribution,seed)
    print('Sequence with error in seed region are generated.')
    
    # for outside seed region error sequences only
    if not repeat:
        gff_df = update_gff(gff_df,out,ground_truth_filename)
    no_mir_chr = sequence_calculation(total_no_seq, xseed_error_percent, depth,seed)    
    fasta_seq = generate_sequence(fasta_seq,gff_df, rna_dict, no_mir_chr, depth, 'Outside_Seed_region', out, ground_truth_filename,'append',repeat,distribution,seed)
    print('Sequence with error in xseed region are generated.')
    
    # for both seed region error and outside seed region error sequences
    if not repeat:
        gff_df = update_gff(gff_df,out,ground_truth_filename)
    no_mir_chr = sequence_calculation(total_no_seq, both_error_percent, depth,seed)    
    fasta_seq = generate_sequence(fasta_seq,gff_df, rna_dict, no_mir_chr, depth, 'Both_region', out, ground_truth_filename,'append',repeat,distribution,seed)
    print('Sequence with error in both seed and xseed region are generated.')

    write_fastq(fasta_seq, adaptor, ascii_base, out, out_file_name,out_file_type,seed,parallel_thread)
    
#     merge similar synthetic sequences
    ground_truth_df = pd.read_csv(os.path.join(out,ground_truth_filename),sep='\t')
    ground_truth_df = ground_truth_df.groupby(['RNA_ID','Impure_Region','Cigar_String','ref_Sequence','synthetic_sequence','chr','chr_start','chr_end']).sum().reset_index()
    ground_truth_df = ground_truth_df.sort_values('Impure_Region')
    ground_truth_df.to_csv(os.path.join(out,ground_truth_filename),encoding='utf-8',index=False)
    
    # deleting temporary folder
    cmd = 'rm -rf '+ os.path.join(out,'temp','..?* .[!.]*') + ' && rm -rf '+ os.path.join(out,'temp') + ' 2> /dev/null'
    os.system(cmd)
    end = time.time()
    time_taken = '%.2f' % ((end-start)/60)
    print(f"Runtime of the program is {time_taken} minutes")
        
        
        