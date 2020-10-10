#!/usr/bin/python

import os
from write_small_fastq_chunks import *
import numpy as np

def execute_parallel_thread_for_file_write(fasta_seq_chunk,out_file,counter,out,out_file_type,quality_char,adaptor,seed):
    max_no_line = len(fasta_seq_chunk)
    if max_no_line > 10000:
        chunk_size = 10000
        no_of_splitted_file = int(np.ceil(max_no_line/chunk_size))
    else:
        chunk_size = max_no_line
        no_of_splitted_file = 1
    for j in range(no_of_splitted_file):
        if out_file_type == 'fastq':
            tmp_file_name = out_file + '_' + str(j) + '.fastq'
        else:
            tmp_file_name = out_file + '_' + str(j) + '.fasta'
            
        out_file1 = open(os.path.join(out,'temp',tmp_file_name),'w')
        fasta_seq_small_chunk = fasta_seq_chunk[j*chunk_size:(j+1)*chunk_size]
        counter1 = counter + j*len(fasta_seq_small_chunk)
        write_small_fastq_chunks(fasta_seq_small_chunk,out_file1,counter1,out_file_type,quality_char,adaptor,seed)