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
import logging
import logging.handlers
import time
import os
import shutil
import sys
import argparse
import subprocess
from pathlib import Path

from Tools import makeDirectory, prechecks
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
# set global variables
args = get_args()
INPDIR = os.path.abspath(args.Input_DIR)
OUTDIR = os.path.abspath(args.Output_DIR)
THREADS = args.Threads
MXMEMORY = args.Memory
REFRENCE = args.HumanReference
KRAK = args.Kraken_PATH
KRAKDB = args.Kraken_DBPATH
KRAK_THRESH = args.Kraken_Hit_Threshold
BRAK = args.Bracken_PATH
BRAKTHRESH = args.Bracken_Hit_Threshold
BRAKLENGTH = args.Bracken_Kmer_Length
####################################################################################################################################################
''' Run the script and functions '''
####################################################################################################################################################
# Setting up logging
if __name__ == '__main__':
    try:
        logger = logging.getLogger("AuOF_MasterControl.py: %s" % time.asctime())
        logger.setLevel(logging.DEBUG)
        # general information and debugging handler
        d_handler = [logging.StreamHandler(sys.stdout),
                    logging.FileHandler("info.log")]
        #d_handler = logging.StreamHandler(sys.stdout)
        d_handler.setLevel(logging.DEBUG)
        d_formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        d_handler.setFormatter(d_formatter)
        logger.addHandler(d_handler)
        # error handler specifically
        e_handler = [logging.StreamHandler(sys.stderr),
                    logging.FileHandler("error.log")]
        e_handler.setLevel(logging.ERROR)
        e_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        e_handler.setFormatter(e_formatter)
        logger.addHandler(e_handler)
    except BaseException:
        outstr = "Debug log file could not be open for logging"
        logger.error(outstr)
        sys.exit(1)
    # Print the arguments to file
    logger.info("Command line: %s", ' '.join(sys.argv))
    logger.info("Starting: %s", time.asctime())
####################################################################################################################################################
    # check inputs
    checkInp, checkOut = prechecks(INPDIR, OUTDIR)
    if checkInp == "Good":
        logger.info("Input files exists. Files will be checked downstream")
    elif checkInp == "Bad":
        logger.info("Script failed. Check error logs for more information")
        logger.error("Script failed to recognise your input directory. Does it exist? " +
                        "Please check and try again.")
    if checkOut == "Make":
        logger.info("Output directory does not exist. Making this")
        makeDirectory(OUTDIR)
    elif checkOut == "Good":
        logger.info("Output directory already exists. Will be using this.")

