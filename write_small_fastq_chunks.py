#!/usr/bin/python

import random

def write_small_fastq_chunks(fasta_seq_chunk,out_file,counter,out_file_type,quality_char,adaptor,seed):

    for line in fasta_seq_chunk:
        line1 = ""
        line_no = fasta_seq_chunk.index(line)
        if ">" in line:
            if out_file_type == 'fastq':
                line1 += "@" + line[1:] + '-' + str(counter) + '\n'
                seq = fasta_seq_chunk[line_no+1] + adaptor # "TGGAATTCTCGGGTGCCAAGG"
                line1 += fasta_seq_chunk[line_no+1] + adaptor + '\n' # "TGGAATTCTCGGGTGCCAAGG" + '\n'
                line1 += '+' + '\n'
                quality_string = ''.join(char*random.randint(0,15) for char in quality_char)
                quality_string = list(quality_string)
                quality_string = random.sample(quality_string,len(quality_string))
                quality_string = ''.join(c for c in quality_string)
                quality_flag = False
                while quality_flag:
                    if len(quality_string) < len(seq):
                        quality_string += quality_string
                        quality_string = list(quality_string)
                        quality_string = random.sample(quality_string,len(quality_string))
                        quality_string = ''.join(c for c in quality_string)
                    else:
                        quality_flag = True
                        break
                line1 += quality_string[:len(seq)] + '\n'
                out_file.write(line1)
                counter += 1
            elif out_file_type == 'fasta':
                line1 += line + '-' + str(counter) + '\n'
                line1 += fasta_seq_chunk[line_no+1] + adaptor + '\n'
                out_file.write(line1)
                counter += 1

#     out_file.close()