#!/usr/bin/python

import argparse
from collections import defaultdict
import time
import warnings
from generate_synthetic_data import * 

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i","--input",help="Reference fasta file input.",type=str)
parser.add_argument("-t","--total_seq",help="Total number of sequence to be generated.",nargs="?", default=10000,type=int)
parser.add_argument("-p","--pure_seq",help="Fraction of pure sequence out of total generated sequence.",default=20,type=int)
parser.add_argument("-s","--seed_error_seq",help="Fraction of sequence having impurity in seed region out of total generated sequence.",default=10,type=int)
parser.add_argument("-x","--xseed_error_seq",help="Fraction of sequence having impurity in xseed region (extra region outside seed region) out of total generated sequence.",default=10,type=int)
parser.add_argument("-b","--both_seed_xseed_error_seq",help="Fraction of sequence having impurity in both seed and xseed region outo of total generated sequence.",default=10,type=int)
parser.add_argument("-d","--min_depth",help="Minimum depth of sequence to be generated.", default=100,type=int)
parser.add_argument("-e","--encoding_quality",help="Quality score encoding for fastq file (requires -q).", default=64,type=int)
parser.add_argument("-se","--seed",help="Seed.",default=108,type=int)
parser.add_argument("-o","--out_path",help="Path of saving output (fastq/fasta) file.")
parser.add_argument("-n","--out_file_name",help="Name of output sequence file (fastq/fasta).")
parser.add_argument("-g","--ground_truth_file",help="Name of output ground truth file.")
parser.add_argument("-gff","--gff_file",help="GFF file.")
parser.add_argument("-q","--out_file_type",help="Output file type.",nargs="?", default='fastq',type=str)
parser.add_argument("-dist","--expression_distribution",help="Distribution type for expression values.",nargs="?", default='poisson',type=str)
parser.add_argument("-a","--adaptor",help="Adaptor Sequence.",nargs="?", default='TGGAATTCTCGGGTGCCAAGG',type=str)
parser.add_argument("-th","--thread",help="Number of Parallel thread.",nargs="?", default='Number of CPUs',type=int)
parser.add_argument("-r","--replacement",help="Sample RNAs with replacement",nargs="?", default=True)
parser.add_argument("-rna","--rna_type",help="RNA type (miRNA/piRNA/...).",nargs="?", default='mirna',type=str)

args = parser.parse_args()

generate_synthetic_data(args.input, args.total_seq, args.pure_seq, args.seed_error_seq,
                        args.xseed_error_seq, args.both_seed_xseed_error_seq, args.min_depth,
                        args.encoding_quality, args.out_path, args.out_file_name, 
                        args.ground_truth_file, args.gff_file, args.rna_type, args.out_file_type,
                        args.replacement, args.expression_distribution, args.seed, args.adaptor,
                        args.thread)
