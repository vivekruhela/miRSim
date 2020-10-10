#!/usr/bin/python

import random

def write_small_fastq_chunks(fasta_seq_chunk,out_file,counter,out_file_type,quality_char,adaptor,seed):

    primer_nt = ['A','T','G','C']
    for line in fasta_seq_chunk:
        line1 = ""
        line_no = fasta_seq_chunk.index(line)
        if ">" in line:
            if out_file_type == 'fastq':
                line1 += "@" + line[1:] + '-' + str(counter) + '\n'
                try:
                    seq = fasta_seq_chunk[line_no+1] + adaptor # "TGGAATTCTCGGGTGCCAAGG"
                except:
                    seq = fasta_seq_chunk[line_no] + adaptor # "TGGAATTCTCGGGTGCCAAGG"
                primer_len = 75 - len(seq)
                primer_string = ''.join(char*random.randint(0,100) for char in primer_nt)
                primer_string = list(primer_string)
                primer_string = random.sample(primer_string,len(primer_string))
                primer_string = ''.join(c for c in primer_string)[:primer_len]
                seq += primer_string                
                try:
                    line1 += fasta_seq_chunk[line_no+1] + adaptor + primer_string + '\n' # "TGGAATTCTCGGGTGCCAAGG" + '\n'
                except:
                    line1 += fasta_seq_chunk[line_no] + adaptor + primer_string + '\n' # "TGGAATTCTCGGGTGCCAAGG" + '\n'
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
                try:
                    line1 += fasta_seq_chunk[line_no+1] + adaptor + '\n'
                except:
                    line1 += fasta_seq_chunk[line_no] + adaptor + '\n'
                out_file.write(line1)
                counter += 1

#     out_file.close()