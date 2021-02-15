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
parser.add_argument("-a","--adaptor",help="Adaptor Sequence.", default='TGGAATTCTCGGGTGCCAAGG',type=str)
parser.add_argument("-b","--both_seed_xseed_error_seq",help="Fraction of sequence having impurity in both seed and xseed region outo of total generated sequence.",type=float)
parser.add_argument("-d","--min_depth",help="Minimum depth of sequence to be generated.", default=5,type=int)
parser.add_argument("-dist","--expression_distribution",help="Distribution type for expression values.", default='poisson',type=str)
parser.add_argument("-e","--encoding_quality",help="Quality score encoding for fastq file (33/64 for fastq, 0 for fasta).", default=33,type=int)
parser.add_argument("-g","--ground_truth_file",help="Name of output ground truth file.")
parser.add_argument("-gff","--gff_file",help="GFF file.")
parser.add_argument("-i","--input",help="Reference fasta file input.",type=str)
parser.add_argument("-ms","--mismatch_seed",help="Maximum number of mismatch in seed region", default=2, type=int)
parser.add_argument("-mxs","--mismatch_xseed",help="Maximum number of mismatch in xseed region", default=2, type=int)
parser.add_argument("-n","--out_file_name",help="Name of output sequence file (fastq/fasta).")
parser.add_argument("-nr","--non_rna",help="Fraction of Non RNA sequence out of total generated sequence.",type=float)
parser.add_argument("-o","--out_path",help="Path of saving output (fastq/fasta) file.")
parser.add_argument("-q","--out_file_type",help="Output file type.", default='fastq',type=str)
parser.add_argument("-r","--replacement",help="Sample RNAs with replacement",default=True)
parser.add_argument("-rna","--rna_type",help="RNA type (miRNA/piRNA/...).", default='miRNA',type=str)
parser.add_argument("-s","--seed_error_seq",help="Fraction of sequence having impurity in seed region out of total generated sequence.",type=float)
parser.add_argument("-se","--seed",help="Seed (random/fixed-prided by user).",type=int)
parser.add_argument("-st","--std_seq",help="Fraction of stadnard sequence out of total generated sequence.",type=float)
parser.add_argument("-t","--total_seq",help="Total number of sequence to be generated.",nargs="?", default=50000,type=int)
parser.add_argument("-th","--thread",help="Number of Parallel thread.", default=4,type=int)
parser.add_argument("-x","--xseed_error_seq",help="Fraction of sequence having impurity in xseed region (extra region outside seed region) out of total generated sequence.",type=float)

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

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
    sys.exit(bcolors.FAIL + 'Please provide referencce sequence file' + bcolors.ENDC)

if not args.gff_file:
    print(bcolors.FAIL + 'Please provide referencce GFF file' + bcolors.ENDC)
    sys.exit(1)
    
if args.mismatch_seed >= 6:
    print(bcolors.FAIL + 'The maximum number of mismatch in seed region should be less than 6.' + bcolors.ENDC)
    sys.exit(1)
    
if args.mismatch_xseed >= 12:
    print(bcolors.FAIL + 'The maximum number of mismatch in xseed region should be less than 12.' + bcolors.ENDC)
    sys.exit(1)

print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
default_param = {}
if args.total_seq == 50000:
    default_param['Total number of sequences to be generated'] = args.total_seq
if args.min_depth == 5:
    default_param['Minimum depth'] = args.min_depth
if args.encoding_quality == 33:
    default_param['Encoding Quality'] = args.encoding_quality
if not args.min_depth:
    default_param['Minimum depth'] = args.min_depth
if args.mismatch_seed == 2:
    default_param['Max Mismatch Seed Region'] = args.mismatch_seed
if args.mismatch_xseed == 2:
    default_param['Max Mismatch Xseed Region'] = args.mismatch_xseed
if not args.encoding_quality:
    default_param['Encoding Quality'] = args.encoding_quality
if args.expression_distribution == 'poisson':
    default_param['Expression Distribution'] = args.expression_distribution
if args.replacement:
    default_param['Select RNAs with Replacement'] = args.replacement
if args.adaptor == 'TGGAATTCTCGGGTGCCAAGG':
    default_param['Adaptor'] = args.adaptor
if args.thread == 4:
    default_param['Number of thread'] = args.thread
if args.rna_type == 'miRNA':
    default_param['RNA type'] = args.rna_type

print(bcolors.HEADER + 'The following parameters are at default values:' + bcolors.ENDC)
if bool(default_param):
    for k,v in default_param.items():
        print(bcolors.WARNING + str(k) + " = " + str(v) + bcolors.ENDC)

print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')



print('--------------------------------------------------------------------------------------')

print(bcolors.BOLD + 'Generating the Synthetic sequences using the following Parameters:' + bcolors.ENDC)

print(bcolors.BOLD + 'Input reference file: ' + bcolors.ENDC, args.input)

print(bcolors.BOLD + 'GFF reference file: ' + bcolors.ENDC,args.gff_file)

print(bcolors.BOLD + 'RNA type: ' + bcolors.ENDC, args.rna_type)

print(bcolors.BOLD + 'Total number of sequences to be generated: ' + bcolors.ENDC, args.total_seq)

print(bcolors.BOLD + '% fraction of standard sequences: ' + bcolors.ENDC, args.std_seq)

if args.non_rna:
    if not args.seed_error_seq and not args.xseed_error_seq and not args.both_seed_xseed_error_seq:
        print(bcolors.BOLD + '% fraction of non RNA sequences (approximately equally divided into seed, xseed and both seed xseed error seq): ' + bcolors.ENDC, args.non_rna)
        args.seed_error_seq = int(args.non_rna/3)
        args.xseed_error_seq = int(args.non_rna/3)
        args.both_seed_xseed_error_seq = args.non_rna - args.seed_error_seq - args.xseed_error_seq
        print(bcolors.BOLD + '% fraction of sequences having error in seed region: ' + bcolors.ENDC, args.seed_error_seq)
        print(bcolors.BOLD + '% fraction of sequences having error in xseed region: ' + bcolors.ENDC, args.xseed_error_seq)
        print(bcolors.BOLD + '% fraction of sequences having error in both seed and xseed region: ' + bcolors.ENDC, args.both_seed_xseed_error_seq)
    else:
        sys.exit(bcolors.FAIL + "ERROR: Please do not provide both non_RNA fraction and seed,xeed error fractions." + bcolors.ENDC)
elif args.seed_error_seq and args.xseed_error_seq and args.both_seed_xseed_error_seq:
    print(bcolors.BOLD + '% fraction of sequences having error in seed region: ' + bcolors.ENDC, args.seed_error_seq)
    print(bcolors.BOLD + '% fraction of sequences having error in xseed region: ' + bcolors.ENDC, args.xseed_error_seq)
    print(bcolors.BOLD + '% fraction of sequences having error in both seed and xseed region: ' + bcolors.ENDC, args.both_seed_xseed_error_seq)
else:
    sys.exit(bcolors.FAIL + 'Please provide either non-mir fraction or seed, xseed error fractions.' + bcolors.ENDC)

print(bcolors.BOLD + 'Minimum depth: ' + bcolors.ENDC,args.min_depth)

print(bcolors.BOLD + 'Max Mismatch in seed region: ' + bcolors.ENDC,args.mismatch_seed)

print(bcolors.BOLD + 'Max Mismatch in xseed region: ' + bcolors.ENDC,args.mismatch_xseed)

if not args.out_file_name:
    args.out_file_name = str(args.rna_type) + '.fastq.gz'
    print(bcolors.BOLD + 'Ouput file name: ' + bcolors.ENDC, args.out_file_name)
else:
    print(bcolors.BOLD + 'Ouput file name: ' + bcolors.ENDC, args.out_file_name)

if not args.ground_truth_file:
    args.ground_truth_file = str(args.rna_type) + '_ground_truth.csv'
    print(bcolors.BOLD + 'Ouput ground truth file name: ' + bcolors.ENDC, args.ground_truth_file)
else:
    print(bcolors.BOLD + 'Ouput ground truth file name: ' + bcolors.ENDC, args.ground_truth_file)

if args.out_file_type == 'fastq':
    print(bcolors.BOLD + 'Output file type: ' + bcolors.ENDC, args.out_file_type)
    print(bcolors.BOLD + 'Encoding Quality: ' + bcolors.ENDC, args.encoding_quality)
else:
    print(bcolors.BOLD + 'Output file type: ' + bcolors.ENDC, args.out_file_type)

if not args.out_path:
    args.out_path = os.getcwd()
    print(bcolors.BOLD + 'Output file location: ' + bcolors.ENDC, args.out_path)
else:
    print(bcolors.BOLD + 'Output file location: ' + bcolors.ENDC, args.out_path)

print(bcolors.BOLD + 'Adaptor Sequence: ' + bcolors.ENDC, args.adaptor)

if not args.seed:
    args.seed = 'Not Provided'
    print(bcolors.BOLD + 'Seed: ' + bcolors.ENDC, args.seed)
else:
    print(bcolors.BOLD + 'Seed: ' + bcolors.ENDC, args.seed)

print(bcolors.BOLD + 'Expression distribution: ' + bcolors.ENDC, args.expression_distribution)

print(bcolors.BOLD + 'The numbper of parallel thread: ' + bcolors.ENDC, args.thread)

print(bcolors.BOLD + 'Select RNAs with Replacement: ' + bcolors.ENDC, args.replacement)

print('--------------------------------------------------------------------------------------')    

start = time.time()

generate_synthetic_data(args.input, args.total_seq, args.std_seq, args.seed_error_seq,
                        args.xseed_error_seq, args.both_seed_xseed_error_seq, args.min_depth,
                        args.encoding_quality, args.out_path, args.out_file_name, 
                        args.ground_truth_file, args.gff_file, args.rna_type, args.out_file_type,
                        args.replacement, args.expression_distribution, args.seed, args.adaptor,
                        args.mismatch_seed, args.mismatch_xseed, args.thread)

end = time.time()
time_taken = '%.2f' % ((end-start)/60)
print(f"Runtime of the program is {time_taken} minutes")