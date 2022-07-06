#!/usr/bin/env python3
###########################################################################################################################################################
#
#                               Script to carry out workflow control of microbiome analysis for the Antibiotics Under our Feet project
# 
# Requires: Via conda(fastqc, trim-galore, bmtagger, kraken2, bracken, krona, krakentools)
# Considerations: Remove need for kraken and bracken path arguments. Force installation via conda
# Author: Damilola Oresegun	                                                                                                                                 		              #
###########################################################################################################################################################
import enum
import os
import shutil
import sys
import argparse
import subprocess
from pathlib import Path
####################################################################################################
# set the needed arguments
def get_args():
    parser = argparse.ArgumentParser(description = "Workflow to carry out microbiome analysis " +
                                    "on soil-derived whole-genome metagenomic Illumina sequence " +
                                    "data.")
    ################################################################################################
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-i", "--input", dest = "Input_DIR",
                                action = "store",
                                type = str,
                                help = "The path to the directory holding the demultiplexed " +
                                "FASTQ files. FASTQ files can be gzipped or left uncompressed " +
                                "Note: The FASTQ files have to be named with _1 and _2", 
                                required = True)
    required_args.add_argument("-o", "--output", dest = "Output_DIR",
                                action = "store",
                                type = str,
                                help = "Path to the directory to save analysis outputs",
                                required = True)
    required_args.add_argument("-r", "--reference", dest = "HumanReference",
                                action = "store",
                                type = str,
                                default = "/home/doresegu/scratch/private/CMC_Project/" + 
                                "BMTAGGER_INDEX/hg38.fa",
                                help = "Path to the human reference genome to use to remove " +
                                "human contamination")
    required_args.add_argument("-kr", "--kraken", dest = "Kraken_PATH",
                                action = "store",
                                type = str,
                                help = "Full path to your kraken installation if it is not " +
                                "in the $PATH. If in your $PATH, simply write 'kraken2'",
                                required = True)
    required_args.add_argument("-br", "--bracken", dest = "Bracken_PATH",
                                action = "store",
                                type = str,
                                help = "Full path to your bracken installation if it is not " +
                                "in the $PATH. If in your $PATH, simply write 'bracken'",
                                required = True)
    required_args.add_argument("-kb", "--kraken_DB", dest = "Kraken_DB",
                                action = "store",
                                type = str,
                                help = "Full path to your Kraken database",
                                required = True)
    ################################################################################################
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-t", "--threads", dest = "Threads",
                                action = "store",
                                type = int,
                                default = 24,
                                help = "Number of threads. Default is 24")
    optional_args.add_argument("-m", "--max-memory", dest = "Memory",
                                action = "store",
                                type = int,
                                default = 24000,
                                help = "Maximum available memory to use in megabytes e.g. for " +
                                "100G of memory, enter 10000. Default is '24000' i.e. 24G")
    optional_args.add_argument("-kt", "--kraken_hit_threshold", dest = "Kraken_Hit_Threshold",
                                action = "store",
                                type = int,
                                default = 5,
                                help = "A minimum number of groups that must be matched to place " +
                                "a contig into a taxonomic group. Default is [5]")
    optional_args.add_argument("-bt", "--bracken_hit_threshold", dest = "Bracken_Hit_Threshold",
                                action = "store",
                                type = int,
                                default = 20,
                                help = "A minimum number of kmers that must be matched to place " +
                                "a contig into a taxonomic group by bracken re-estimation. " +
                                "Default is [20]")
    optional_args.add_argument("-bl", "--bracken_kmer_length", dest = "Bracken_Kmer_Length",
                                action = "store",
                                type = int,
                                default = 100,
                                help = "The kmer length used to build your bracken database. " +
                                "Note: Your bracken database is also the kraken database. " +
                                "Default length is [100]")
    ################################################################################################
    args = parser.parse_args()
    return args
####################################################################################################