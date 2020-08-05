#!/usr/bin/python

import pandas as pd

def gff(gff_file,rna_type):
    gff_file = open(gff_file).read().split('\n')
    df = pd.DataFrame()
    if rna_type == 'miRNA':
        gff_file = gff_file[13:-1]
        chr_loc = [c.split('\t')[0] for c in gff_file]
        mirna_type = [c.split('\t')[2] for c in gff_file]
        chr_start = [c.split('\t')[3] for c in gff_file]
        chr_end = [c.split('\t')[4] for c in gff_file]
        strand = [c.split('\t')[6] for c in gff_file]
        mir_id = [c.split('\t')[8].split(';')[-2].split('=')[-1] for c in gff_file]        
#         df = pd.DataFrame()
        df['mir_id'] = mir_id
        df['chr'] = chr_loc
        df['chr_start'] = chr_start
        df['chr_end'] = chr_end
        df['strand'] = strand
        df['mirna_type'] = mirna_type
        df = df.set_index('mir_id')
        df = df[~df['mirna_type'].str.contains('miRNA_primary_transcript')]
        df = df.drop('mirna_type',axis=1)
    else:
        gff_file = gff_file[1:-1]
        rna_id = [c.split('\t')[-1] for c in gff_file]
        chr_loc = [c.split('\t')[0] for c in gff_file]
        chr_start = [c.split('\t')[3] for c in gff_file]
        chr_end = [c.split('\t')[4] for c in gff_file]
        strand = [c.split('\t')[6] for c in gff_file]
        df['rna_id'] = rna_id
        df['chr'] = chr_loc
        df['chr_start'] = chr_start
        df['chr_end'] = chr_end
        df['strand'] = strand
        df = df.set_index('rna_id')
    
    return df
    
    