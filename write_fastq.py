#!/usr/bin/python

import threading
import gzip
import numpy as np
import os
import time
from execute_parallel_thread_for_file_write import *

def write_fastq(fasta_seq, adaptor, ascii_base, out, out_file_name,out_file_type,seed,parallel_thread):
                
    # Not considering phred score less then 20 so that mean quality score is always greater than 20 for smooth bioinformatics analysis
    if ascii_base == 33:
        quality_char = ['5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']
    elif ascii_base == 64:
        quality_char = ['T','U','V','W','X','Y','Z', '^',' ', '-', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
    else:
        quality_char = []
    
    max_no_line = len(fasta_seq)
    no_of_splitted_file = parallel_thread
    chunk_size = int(np.ceil(max_no_line/no_of_splitted_file))
    print('Writing %d temporary files with chunk size of %d'%(no_of_splitted_file,chunk_size))
    if 'temp' in list(os.listdir(out)):
        cmd = 'rm -rf '+ os.path.join(out,'temp','..?* .[!.]*') + ' && rm -rf '+ os.path.join(out,'temp') + ' 2> /dev/null'
        os.system(cmd)
#         shutil.rmtree(os.path.join(out,'temp'))
        os.mkdir(os.path.join(out,'temp'))
    else:
        os.mkdir(os.path.join(out,'temp'))
        
    for j in range(no_of_splitted_file):
        tmp_file_name = 'tmp_'+str(j)            
        out_file = tmp_file_name 
        fasta_seq_chunk = fasta_seq[j*chunk_size:(j+1)*chunk_size]
        if j >=1:
            counter = j*len(fasta_seq_chunk)
            print('------------------------------')
            print('writing %d chunk'%j)
            threading.Thread(target=execute_parallel_thread_for_file_write,args=(fasta_seq_chunk,out_file,counter,out,out_file_type,quality_char,adaptor,seed,)).start()
        else:
            pass

    counter = 0
    tmp_file_name = 'tmp_'+str(0)
    out_file = tmp_file_name 
    fasta_seq_chunk = fasta_seq[0*chunk_size:(0+1)*chunk_size]
    print('------------------------------')
    print('writing %d chunk'%(j+1))
    execute_parallel_thread_for_file_write(fasta_seq_chunk,out_file,counter,out,out_file_type,quality_char,adaptor,seed)
    time.sleep(20)
    print('-----------------------------------------------------------------------')
    print('All chunks are written. Now merging all the temporary files and compressing it.')
    if out_file_type == 'fastq':
        cmd = 'pigz '+ os.path.join(out,'temp') +'/tmp_*.fastq -c > ' + os.path.join(out,out_file_name)
    else:
        cmd = 'pigz '+ os.path.join(out,'temp') + '/tmp_*.fasta -c > ' + os.path.join(out,out_file_name)
    os.system(cmd)
    