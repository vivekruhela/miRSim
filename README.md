# miRSim

<table>
<tr>
  <td>Publications</td>
  <td>
    <a href="http://doi.org/10.5281/zenodo.4560585">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.6546356.svg" alt="zenodo reference">
    </a>
  </td>
</tr>
</table>


![miRSim Synthetic Sequence Simulator](synthetic_reads.png)

This aim of this tool is to generate the synthetic RNA-Seq sequences by mutating the seed and xseed region (extra sequence after removing seed) by following poisson or gamma error distribution. In addition to the synthetic sequences, this tool also generate the ground truth. This tool can be used to evalute the ability of the tool/pipeline to correctly identify the RNA sequences.

## Getting Started

### Prerequisites

This tool has been developed in Python (version 3.6.9). The following libraries are required to run miRSim sequence generation module:

```
-----------------------------------------------------
Packages                        Version
-----------------------------------------------------
pandas                        1.0.1
numpy                         1.18.1
threading                     4.1.0
scipy                         1.4.1
-----------------------------------------------------
```
Rest of the requires packages such as `random,os,sys,gzip,time` are already included into python standard library. This tool has been tested on Linux (Ubuntu) 16.04 and 18.04 version.


### Installing

You can download this tool by cloning the github repository and use directly by switching to the tool directory.

```
git clone https://github.com/vivekruhela/miRSim.git
cd miRSim
pip install -r requirements.txt
```
## Databases
#### Databases Required
[miRBase](http://www.mirbase.org/) and [piRNAdb](https://www.pirnadb.org/)
#### Database used
In refs directory, we have used miRBase (version 22) and piRNAdb (version 1.7.5) for their sequences in fasta files and their genomic location in gff file.

## Project Structure

```
|---- refs
         |---- hsa_high_conf.gff3
         |---- mature_high_conf_hsa.fa
         |---- novel_gff.gff3
         |---- final_novel_seq_filtered.fa
         |---- pirnadb.hg38.gff3
         |---- piRNAdb.hsa.v1_7_5.fa
|---- miRSim.py
|        |
|        |---- generate_synthetic_data.py
|                        |
|                        |---- generate_synthetic_data.py
|                                            |
|                                            |---- gff.py
|                                            |---- update_gff.py
|                                            |---- generate_sequence.py
|                                                          |---- expression_split.py
|                                                          |---- cigar_generation.py
|                                                          |---- mir_location.py
|                                                          |---- sequence_alteration.py
|                                                                        |---- alter_nt.py
|                                            |---- sequence_calculation.py
|                                                          |---- expression_split.py
|                                                          |---- sort_by_chromosome.py
|                                            |---- write_fastq.py
|                                                       |---- execute_parallel_thread_for_file_write.py
|                                                                         |---- write_small_fastq_chunks.py
|---- miRSim.ipynb
|---- synthetic_reads.png
|---- README.md
|---- LICENSE
```

## Arguments

```
usage: miRSim.py [-h] [-a ADAPTOR] [-b BOTH_SEED_XSEED_ERROR_SEQ]
                 [-d MIN_DEPTH] [-dist EXPRESSION_DISTRIBUTION]
                 [-e ENCODING_QUALITY] [-g GROUND_TRUTH_FILE] [-gff GFF_FILE]
                 [-i INPUT] [-ms MISMATCH_SEED] [-mxs MISMATCH_XSEED]
                 [-n OUT_FILE_NAME] [-nr NON_RNA] [-o OUT_PATH]
                 [-q OUT_FILE_TYPE] [-r REPLACEMENT] [-rna RNA_TYPE]
                 [-s SEED_ERROR_SEQ] [-se SEED] [-st STD_SEQ] [-t [TOTAL_SEQ]]
                 [-th THREAD] [-x XSEED_ERROR_SEQ]

optional arguments:
  -h, --help            show this help message and exit
  -a ADAPTOR, --adaptor ADAPTOR
                        Adaptor Sequence. (default: TGGAATTCTCGGGTGCCAAGG)
  -b BOTH_SEED_XSEED_ERROR_SEQ, --both_seed_xseed_error_seq BOTH_SEED_XSEED_ERROR_SEQ
                        Fraction of sequence having impurity in both seed and
                        xseed region outo of total generated sequence.
                        (default: None)
  -d MIN_DEPTH, --min_depth MIN_DEPTH
                        Minimum depth of sequence to be generated. (default:
                        5)
  -dist EXPRESSION_DISTRIBUTION, --expression_distribution EXPRESSION_DISTRIBUTION
                        Distribution type for expression values. (default:
                        poisson)
  -e ENCODING_QUALITY, --encoding_quality ENCODING_QUALITY
                        Quality score encoding for fastq file (33/64 for
                        fastq, 0 for fasta). (default: 33)
  -g GROUND_TRUTH_FILE, --ground_truth_file GROUND_TRUTH_FILE
                        Name of output ground truth file. (default: None)
  -gff GFF_FILE, --gff_file GFF_FILE
                        GFF file. (default: None)
  -i INPUT, --input INPUT fasta file (default: None)
                        Reference fasta file input. (default: None)
  -ms MISMATCH_SEED, --mismatch_seed MISMATCH_SEED
                        Maximum number of mismatch in seed region (default: 2)
  -mxs MISMATCH_XSEED, --mismatch_xseed MISMATCH_XSEED
                        Maximum number of mismatch in xseed region (default:
                        2)
  -n OUT_FILE_NAME, --out_file_name OUT_FILE_NAME
                        Name of output sequence file (fastq/fasta). (default:
                        None)
  -nr NON_RNA, --non_rna NON_RNA
                        Fraction of Non RNA sequence out of total generated
                        sequence. (default: None)
  -o OUT_PATH, --out_path OUT_PATH
                        Path of saving output (fastq/fasta) file. (default:
                        None)
  -q OUT_FILE_TYPE, --out_file_type OUT_FILE_TYPE
                        Output file type. (default: fastq)
  -r REPLACEMENT, --replacement REPLACEMENT
                        Sample RNAs with replacement (default: True)
  -rna RNA_TYPE, --rna_type RNA_TYPE
                        RNA type (miRNA/piRNA/...). (default: miRNA)
  -s SEED_ERROR_SEQ, --seed_error_seq SEED_ERROR_SEQ
                        Fraction of sequence having impurity in seed region
                        out of total generated sequence. (default: None)
  -se SEED, --seed SEED
                        Seed (random/fixed-prided by user). (default: None)
  -st STD_SEQ, --std_seq STD_SEQ
                        Fraction of stadnard sequence out of total generated
                        sequence. (default: None)
  -t [TOTAL_SEQ], --total_seq [TOTAL_SEQ]
                        Total number of sequence to be generated. (default:
                        50000)
  -th THREAD, --thread THREAD
                        Number of Parallel thread. (default: 4)
  -x XSEED_ERROR_SEQ, --xseed_error_seq XSEED_ERROR_SEQ
                        Fraction of sequence having impurity in xseed region
                        (extra region outside seed region) out of total
                        generated sequence. (default: None)
```

## Steps to generate Synthetic data

### Using pre-exieting databases such as miRBase or piRNAdb
You can download the sequence fasta file (.fasta) and their genomic location (.gff) from the miRBase/piRNAdb database to generate the synthetic data. (Example shown below)

### Example
(A) Basic Example: If you want to prepare the synthetic data having total number of sequences = 50000 and contains standard miRNA (let's say 50% i.e. 25000 sequences) and non-miRNAs (let's say 50% i.e. 25000 sequences)

```
python miRSim.py -i refs/mature_high_conf_hsa.fa -gff refs/hsa_high_conf.gff3 -st 50 -nr 50
```

(B) Advance Example: If you want to prepare the synthetic data that contains multiple type of RNAs (let's say miRNA + piRNA + novel miRNA) with following proportion:

```
Total Number of sequence = 500000
% of pure miRNA sequence                                           = 20% (i.e. 20000 sequences)
% of miRNA sequence with error in seed region                      = 10% (i.e. 10000 sequences)
% of miRNA sequence with error in xseed region                     = 10% (i.e. 10000 sequences)
% of miRNA sequence with error in both seed and xseed region       = 5% (i.e. 5000 sequences)
% of pure piRNA sequences                                          = 10% (i.e. 10000 sequences)
% of piRNA sequence with error in seed region                      = 10% (i.e. 10000 sequences)
% of piRNA sequence with error in xseed region                     = 5% (i.e. 5000 sequences)
% of piRNA sequence with error in both seed and xseed region       = 5% (i.e. 5000 sequences)
% of pure novel miRNA                                              = 10% (i.e. 10000 sequences)
% of novel miRNA sequence with error in seed region                = 10% (i.e. 10000 sequences)
% of novel miRNA sequence with error in xseed region               = 3% (i.e. 3000 sequences)
% of novel miRNA sequence with error in both seed and xseed region = 2% (i.e. 2000 sequences)

                                                             Total = 100%
```

In order to generate such data you need to call miRSim module separately for each RNA type using following commands `(assuming minimum depth=default, file_type=fastq, encoding_quality=33, adaptor=default)`:
**For miRNA synthetic data**

```
python miRSim.py -i refs/mature_high_conf_hsa.fa -n mirna_raw_data.fastq.gz -g mirna_ground_truth.csv -gff refs/hsa_high_conf.gff3 -t 500000 -st 20 -s 10 -x 10 -b 5 -se 1001 -th 6 -rna miRNA
```
**For piRNA synthetic data**

```
python miRSim.py -i refs/piRNAdb.hsa.v1_7_5.fa -n pirna_raw_data.fastq.gz -g pirna_ground_truth.csv -gff refs/pirnadb.hg38.gff3 -t 500000 -st 10 -s 10 -x 5 -b 5 -se 1001 -th 6 -rna piRNA
```
**For novel miRNA synthetic data**

```
python miRSim.py -i refs/final_novel_seq_filtered.fa -n novel_mirna_raw_data.fastq.gz -g novel_mirna_ground_truth.csv -gff refs/novel_gff.gff3 -t 500000 -st 10 -s 10 -x 3 -b 2 -se 1001 -th 6 -rna novelRNA
```
After generating synthetic data for each individual RNA, merge these fastq files and their ground truth csv. e.g.
```
zcat mirna_raw_data.fastq.gz pirna_raw_data.fastq.gz novel_mirna_raw_data.fastq.gz | gzip > synthetic_raw_data.fastq.gz
```

### Using manually curated database
1. You can prepare the sequence fasta file that contains all the reference sequences in the standard fasta format.
2. Also, prepare gff file that contains the genomic locations of the sequence mentioned in the fasta file in the standard GFF3 format.
3. After preparing both fasta and gff reference files, you can use them with `-i` and `-gff` argument for synthetic data generation.

## Reference

 Ruhela, V., Gupta, R., Krishnamachari, S., Ahuja, G., Gupta, A.:vivekruhela/miRSim v1.0.0 (Version v1.0.0). Zenodo

## Authors

* **Vivek Ruhela** - *Initial work* - [github](https://github.com/vivekruhela)


## License

See the [LICENSE](LICENSE) file for license rights and limitations (Apache2.0).

## Acknowledgement

1. Authors would like to gratefully acknowledge the support by grant from Department of Biotechnology, Govt. of India [Grant: BT/MED/30/SP11006/2015] and Department of Science and Technology, Govt. of India [Grant: DST/ICPS/CPS-Individual/2018/279(G)].
2. Authors would like to gratefully acknowledge the support of Computational Biology Dept., Indraprastha Institute of Information Technology-Delhi (IIIT-D), India for providing resources for tool development.
3. Authors would like to gratefully acknowledge the support of SBILab, Deptt. of ECE & Centre of Excellence in Healthcare, Indraprastha Institute of Information Technology-Delhi (IIIT-D), India for providing guidance in tool metholody and development.
