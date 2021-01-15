#!/usr/bin/python

import argparse
from collections import defaultdict
import time
import warnings
import os
import sys
import random
from generate_synthetic_data import * 

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i","--input",help="Reference fasta file input.",type=str)
parser.add_argument("-t","--total_seq",help="Total number of sequence to be generated.",nargs="?", default=50000,type=int)
parser.add_argument("-st","--std_seq",help="Fraction of stadnard sequence out of total generated sequence.",type=int)
parser.add_argument("-nr","--non_rna",help="Fraction of Non RNA sequence out of total generated sequence.",type=int)
parser.add_argument("-s","--seed_error_seq",help="Fraction of sequence having impurity in seed region out of total generated sequence.",type=int)
parser.add_argument("-x","--xseed_error_seq",help="Fraction of sequence having impurity in xseed region (extra region outside seed region) out of total generated sequence.",type=int)
parser.add_argument("-b","--both_seed_xseed_error_seq",help="Fraction of sequence having impurity in both seed and xseed region outo of total generated sequence.",type=int)
parser.add_argument("-d","--min_depth",help="Minimum depth of sequence to be generated.", default=5,type=int)
parser.add_argument("-e","--encoding_quality",help="Quality score encoding for fastq file (33/64 for fastq, 0 for fasta).", default=33,type=int)
parser.add_argument("-se","--seed",help="Seed (random/fixed-prided by user).",type=int)
parser.add_argument("-o","--out_path",help="Path of saving output (fastq/fasta) file.")
parser.add_argument("-n","--out_file_name",help="Name of output sequence file (fastq/fasta).")
parser.add_argument("-g","--ground_truth_file",help="Name of output ground truth file.")
parser.add_argument("-gff","--gff_file",help="GFF file.")
parser.add_argument("-q","--out_file_type",help="Output file type.", default='fastq',type=str)
parser.add_argument("-dist","--expression_distribution",help="Distribution type for expression values.", default='poisson',type=str)
parser.add_argument("-a","--adaptor",help="Adaptor Sequence.", default='TGGAATTCTCGGGTGCCAAGG',type=str)
parser.add_argument("-th","--thread",help="Number of Parallel thread.", default=4,type=int)
parser.add_argument("-r","--replacement",help="Sample RNAs with replacement",default=True)
parser.add_argument("-rna","--rna_type",help="RNA type (miRNA/piRNA/...).", default='miRNA',type=str)

args = parser.parse_args()

print('######################################################################################')
print('miRSim: Seed-based Synthetic Sequence Simualtor')
print('Version: 1.0')
print('Developer: Vivek Ruhela, Ph.d.')
print('Department: Center for Computational Biology')
print('Institute: Indraprastha Institute of Information Technnology, New Delhi, India')
print('######################################################################################')

if args.encoding_quality:
    if args.out_file_type == 'fasta':
        if not args.encoding_quality == 0:
            sys.exit('ERROR : To generate fasta file, the parameter -e should be equal to 0.')

if not args.input:
    sys.exit('Please provide referencce sequence file')

if not args.gff_file:
    print('Please provide referencce GFF file')
    sys.exit(1)

print('--------------------------------------------------------------------------------------')
default_param = {}
if not args.total_seq:
    default_param['Total number of sequences to be generated'] = args.total_seq
if not args.min_depth:
    default_param['Minimum depth'] = args.min_depth
if not args.encoding_quality:
    default_param['Encoding Quality'] = args.encoding_quality
if not args.expression_distribution:
    default_param['Expression Distribution'] = args.expression_distribution
if not args.replacement:
    default_param['Select RNAs with Replacement'] = args.replacement
if not args.adaptor:
    default_param['Adaptor'] = args.adaptor
if not args.thread:
    default_param['Number of thread'] = args.thread
if not args.rna_type:
    default_param['RNA type'] = args.rna_type

print('The following parameters are at default values:')
if bool(default_param):
    for k,v in default_param.items():
        print(k, " = ", v)


print('--------------------------------------------------------------------------------------')



print('--------------------------------------------------------------------------------------')

print('Generating the Synthetic sequences using the following Parameters:')

print('Input reference file: ', args.input)

print('GFF reference file: ',args.gff_file)

print('RNA type: ', args.rna_type)

print('Total number of sequences to be generated: ', args.total_seq)

print('% fraction of standard sequences: ', args.std_seq)

if args.non_rna:
    if not args.seed_error_seq and not args.xseed_error_seq and not args.both_seed_xseed_error_seq:
        print('% fraction of non RNA sequences (approximately equally divided into seed, xseed and both seed xseed error seq): ', args.non_rna)
        args.seed_error_seq = int(args.non_rna/3)
        args.xseed_error_seq = int(args.non_rna/3)
        args.both_seed_xseed_error_seq = args.non_rna - args.seed_error_seq - args.xseed_error_seq
        print('% fraction of sequences having error in seed region: ', args.seed_error_seq)
        print('% fraction of sequences having error in xseed region: ', args.xseed_error_seq)
        print('% fraction of sequences having error in both seed and xseed region: ', args.both_seed_xseed_error_seq)
    else:
        sys.exit("ERROR: Please do not provide both non_RNA fraction and seed,xeed error fractions.")
elif args.seed_error_seq and args.xseed_error_seq and args.both_seed_xseed_error_seq:
    print('% fraction of sequences having error in seed region: ', args.seed_error_seq)
    print('% fraction of sequences having error in xseed region: ', args.xseed_error_seq)
    print('% fraction of sequences having error in both seed and xseed region: ', args.both_seed_xseed_error_seq)
else:
    sys.exit('Please provide either non-mir fraction or seed, xseed error fractions.')

print('Minimum depth: ',args.min_depth)

if not args.out_file_name:
    args.out_file_name = str(args.rna_type) + '.fastq.gz'
    print('Ouput file name: ', args.out_file_name)
else:
    print('Ouput file name: ', args.out_file_name)

if not args.ground_truth_file:
    args.ground_truth_file = str(args.rna_type) + '_ground_truth.csv'
    print('Ouput ground truth file name: ', args.ground_truth_file)
else:
    print('Ouput ground truth file name: ', args.ground_truth_file)

if args.out_file_type == 'fastq':
    print('Output file type: ', args.out_file_type)
    print('Encoding Quality: ', args.encoding_quality)
else:
    print('Output file type: ', args.out_file_type)

if not args.out_path:
    args.out_path = os.getcwd()
    print('Output file location: ', args.out_path)
else:
    print('Output file location: ', args.out_path)

print('Adaptor Sequence: ', args.adaptor)

if not args.seed:
    args.seed = 'Not Provided'
    print('Seed: ', args.seed)
else:
    print('Seed: ', args.seed)

print('Expression distribution: ', args.expression_distribution)

print('The numbper of parallel thread: ', args.thread)

print('Select RNAs with Replacement: ', args.replacement)

print('--------------------------------------------------------------------------------------')    

generate_synthetic_data(args.input, args.total_seq, args.std_seq, args.seed_error_seq,
                        args.xseed_error_seq, args.both_seed_xseed_error_seq, args.min_depth,
                        args.encoding_quality, args.out_path, args.out_file_name, 
                        args.ground_truth_file, args.gff_file, args.rna_type, args.out_file_type,
                        args.replacement, args.expression_distribution, args.seed, args.adaptor,
                        args.thread)
