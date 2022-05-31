#!/usr/bin/env python3
# quick layout of the workflow
# pre-qc with fastqc
# trim with trimgalore
# remove host
# post trim fastqc
# kraken
# braken
# krona
# ask clarissa and rachel if they wasnt to assembly the classified and unclassified reads
# functional profiling?
###########################################################################################################################################################
#
# Script to carry out workflow control of microbiome analysis for the Antibiotics Under our Feet project
# 
# RequiresL fastqc, trim-galore, bmtagger, kraken2, bracken,
# Author: Damilola Oresegun	                                                                                                                                 		              #
###########################################################################################################################################################
from curses import flash
import enum
from imaplib import Int2AP
import os
import shutil
import sys
import argparse
import subprocess
from pathlib import Path
##########################################################################
# set the needed arguments
def get_args():
    parser = argparse.ArgumentParser(description="Workflow to carry out" +
                                    "microbiome analysis on soil-derived " +
                                    "whole-genome metagenomic " + 
                                    "Illumina sequence data.")
    ##########################################################################
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-i", "--input", 
                               dest="Input_DIR",
                               action="store",
                               type=str,
                               help="The path to the directory " +
                               "holding the demultiplexed FASTQ " +
                               "files. FASTQ files can be gzipped " +
                               "or left uncompressed. "+ 
                               "Note: The fastq files have to be " +
                               "named with _1 and _2",
                               required=True)
    required_args.add_argument("-o", "--output",
                               dest="Output_DIR",
                               action="store",
                               type=str,
                               help="Path to the directory to save " +
                               "analysis outputs",
                               required=True)
    required_args.add_argument("-r", "--reference",
                                dest="HumanReference",
                                action="store",
                                type=str,
                                default="/home/doresegu/scratch/private/" +
                                "CMC_Project/BMTAGGER_INDEX/hg38.fa",
                                help="Path to the human reference genome" +
                                "to use to remove human contamination.")
    required_args.add_argument("-kr", "--kraken",
                                dest="Kraken_PATH",
                                action="store",
                                type=str,
                                help="Full path to your kraken installation " +
                                "if it is not in your $PATH. If in your $PATH "+
                                "simply write kraken2",
                                required=True)
    required_args.add_argument("-br", "--bracken",
                                dest="Bracken_PATH",
                                action="store",
                                type=str,
                                help="Full path to your bracken installation " +
                                "if it is not in your $PATH. If in your $PATH "+
                                "simply write bracken",
                                required=True)
    required_args.add_argument("-kb", "--kraken_DB",
                                dest="Kraken_DBPATH",
                                action="store",
                                type=str,
                                help="Full path to your kraken database",
                                required=True)
    ##########################################################################
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-p", "--processes",
                               dest='Threads',
                               action="store",
                               type=int,
                               default="24",
                               help="Number of threads. Default is 24")
    optional_args.add_argument("-m","--max-memory",
                                dest="Memory",
                                action="store",
                                type=str,
                                default="24000",
                                help="Maximum available memory to use in " +
                                "megabytes e.g. for 100G of memory, enter" +
                                "100000. Default is [24000] i.e 24G")
    optional_args.add_argument("-kt", "--kraken_hit_threshold",
                                dest='Kraken_Hit_Threshold',
                                action="store",
                                type=int,
                                default=2,
                                help="A minimum number of groups that must" +
                                "be matched to place a contig into a " +
                                "taxonomic group")
    optional_args.add_argument("-bt", "--bracken_hit_threshold",
                                dest='Bracken_Hit_Threshold',
                                action="store",
                                type=int,
                                default=20,
                                help="A minimum number of kmers that must" +
                                "be matched to place a contig into a " +
                                "taxonomic group by bracken reestimation" +
                                "Default is [20]")
    optional_args.add_argument("-bl", "--bracken_kmer_length",
                                dest='Bracken_Kmer_Length',
                                action="store",
                                type=int,
                                default=100,
                                help="The kmer length used to build your " +
                                "bracken index database (from the krakenDB " +
                                "Default is [100]")
    
    """
    Options to add:    
    """
    ##########################################################################
    args = parser.parse_args()
    return args
##############################################################################
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
##############################################################################
# check if the input exists
if os.path.exists(INPDIR):
    pass
else:
    print('Are you sure your input directory exists?')
    print('Input directory not found. Please check again')
    sys.exit(1)
    # check if output exists
if os.path.exists(OUTDIR):
    print('Output folder already exists, outputs will be save here')
    pass
else:
    print('Output folder does not exist, creating it now')
    os.makedirs(OUTDIR)
##############################################################################
'''function for formatting outputs to the user
announcements are for large step changes
tells are for small work steps
'''
def comms(command,text):
    if command == "announce":
        delim = "#"
    elif command == "tell":
        delim = "-"
    print('\n'+delim*120)
    max_len=90
    cut=text.split(" ")
    line=""
    for word in cut:
        if (len(line) + 1 + len(word))>max_len:
            edge1=(120-len(line))/2 - 5
            edge2=120-edge1-len(line) - 10
            print(delim*5 + " "*int(edge1) + line + " "*int(edge2) + delim*5)
            line=word
        else:
            line = line+" "+word
    edge1=(120-len(line))/2 - 5
    edge2=120-edge1-len(line) - 10
    print(delim*5 + " "*int(edge1) + line + " "*int(edge2) + delim*5)
    print(delim*120+'\n')
##############################################################################
''' function to run fastqc for pre and post trimming
usage: fastqc(input read folder, output/preQCm, 36)
can also be output/postQC
'''
def fastqc(phase,inpt,outpt,threads):
    if phase == "pre":
        folName = "PreQC_FastQC"
        for i in Path(os.path.abspath(inpt)).glob('*'):
            fname = os.path.basename(i).split("_")[0]
            fastOut = os.path.join(outpt, fname,folName)
            if os.path.exists(fastOut):
                pass
            else:
                os.makedirs(fastOut)
            # get the reads
            fasfi1 = os.path.join(i, "*_R1*")
            fasfi2 = os.path.join(i, "*_R2*")
            # build the command
            fatq1 = ("fastqc", "-t", str(threads), "-o", fastOut, "-f fastq", fasfi1, fasfi2)
            runFatq = ' '.join(fatq1)
            tellU = "Command for fastqc is: " + runFatq
            comms("tell", tellU)
            print(runFatq)
            # run the fastqc command
            subprocess.call(runFatq, shell=True)
    elif phase == "post":
        #folName = "PostQC_FastQC"
        fname = os.path.dirname(inpt).split("/")[-1]
        #fastOut = os.path.join(outpt, fname,folName)
        print("checc " +fname)
        fasfi1 = os.path.join(inpt, fname+"_R1_clean_reads.fastq")
        fasfi2 = os.path.join(inpt, fname+"_R2_clean_reads.fastq")
        print(fasfi1)
        print(fasfi2)
        fatq1 = ("fastqc", "-t", str(threads), "-o", outpt, "-f fastq", fasfi1, fasfi2)
        runFatq = ' '.join(fatq1)
        tellU = "Command for fastqc is: " + runFatq
        comms("tell", tellU)
        print(runFatq)
        # run the fastqc command
        subprocess.call(runFatq, shell=True)
##############################################################################
'''function to run trim-galore for trimming forward and reverse reads'''
def trimmy(inpt,outpt,threads):
    for i in Path(inpt).glob('*'):
        fname = os.path.basename(i).split("_")[0]
        trmOut = os.path.join(outpt, fname, "TrimmedReads")
        if os.path.exists(trmOut):
            pass
        else:
            os.makedirs(trmOut)
        # get the reads
        fasfi1 = os.path.join(i, "*_R1*")
        fasfi2 = os.path.join(i, "*_R2*")
        trm = ("trim_galore --no_report_file --paired --cores", str(threads), "-o", trmOut, fasfi1, fasfi2)
        runTrm = ' '.join(trm)
        tellU = "Command for trim-galore is: " + runTrm
        comms("tell", tellU)
        # run the command
        subprocess.call(runTrm, shell=True)
        my_message = "Renaming the trimmed reads"
        comms("tell", my_message)
        # get the filenames
        for i in Path(trmOut).glob('*'):
            print(str(i))
            if (str(i).__contains__("R1")):
                outFil = os.path.join(trmOut, fname+"_R1_trimmed.fastq")
                print(outFil)
                shutil.copy2(str(os.path.abspath(i)), os.path.abspath(outFil))
                #os.remove(str(os.path.abspath(i)))
            if (str(i).__contains__("R2")):
                outFil = os.path.join(trmOut, fname+"_R2_trimmed.fastq")
                shutil.copy2(str(os.path.abspath(i)), os.path.abspath(outFil))
                #os.remove(str(os.path.abspath(i)))
##############################################################################
'''function to remove human DNA contamination. Requires the reference to be
indexed already. However the function will check if it has been indexed and
if not, indexing will take place
'''
def bmtagAln(samfol,ref,inp1,inp2,outpt,isolate):
    print(ref)
    chekRef = ref.split(".fa")[0]
    if os.path.exists(os.path.abspath(chekRef)+".bitmask"):
        tellU = "Human reference has already been indexed. Moving on"
        comms("tell", tellU) 
    else:
        tellU = "Human reference has not been indexed. Indexing will begin now"
        comms("tell", tellU)
        # run index
        bmInd = ("bmtool -d", ref, "-o", os.path.abspath(chekRef)+".bitmask")
        bmPsm = ("srprism mkindex -i", ref, "-o", os.path.abspath(chekRef)+".srprism", "-M", MXMEMORY)
        runBmInd = ' '.join(bmInd)
        runBmPsm = ' '.join(bmPsm)
        tellU = "Making index"
        comms("tell",tellU)
        subprocess.call(runBmInd,shell=True)
        subprocess.call(runBmPsm, shell=True)
        tellU = "Index complete. Moving on"
        comms("tell", tellU)
    tellU = "Starting removal of human DNA contamination"
    comms("tell", tellU)
    # make temporary output folder
    tmpOut = os.path.join(outpt,"tempOut")
    if os.path.exists(tmpOut):
        pass
    else:
        os.makedirs(tmpOut)
    refr = os.path.abspath(chekRef)+".bitmask"
    refx = os.path.abspath(chekRef)+".srprism"
    outFile = outpt + "/bmtagger.list"
    bmTag = ("bmtagger.sh -b", refr, "-x", refx, "-T", tmpOut, "-q 1 -1", 
            inp1, "-2", inp2, "-o", outFile)
    runBmTag = ' '.join(bmTag)
    tellU = "Looking for human reads now"
    comms("tell", tellU)
    print(runBmTag)
    subprocess.call(runBmTag, shell=True)
    tellU = "BMTAG complete. Cleaning reads now"
    comms("tell", tellU)
    # get the human reads, skip them and write the non-human reads to file
    cleanOut = os.path.join(samfol, "CleanReads")
    if os.path.exists(cleanOut):
        pass
    else:
        os.makedirs(cleanOut)
    humanRead={}
    for line in open(outFile):
            humanRead[line.strip()]=None
    cleanRead1 = os.path.join(cleanOut, isolate+"_R1_clean_reads.fastq")
    # Print out the fastq file line by line unless the read is human:
    skip=False
    clean1out = open(cleanRead1, "w", newline='')
    for i, line in enumerate(open(inp1)):
        if i%4==0:
            if line[1:].split("/")[0].split()[0] in humanRead:
                skip=True
            else:
                skip=False
        if skip==False:
            #print(line.rstrip())
            clean1out.write(line.rstrip()+"\n")
    cleanRead2 = os.path.join(cleanOut, isolate+"_R2_clean_reads.fastq")
    skip=False
    clean2out = open(cleanRead2, "w", newline='')
    for i, line in enumerate(open(inp2)):
        if i%4==0:
            if line[1:].split("/")[0].split()[0] in humanRead:
                skip=True
            else:
                skip=False
        if skip==False:
            #print(line.rstrip())
            clean2out.write(line.rstrip()+"\n")
    """ humanRead = {}
    for line in open(outFile):
        humanRead[line.strip()]=None
    cleanRead1 = os.path.join(cleanOut, isolate+"_R1_clean_reads.fastq")
    skip=False
    clean1out = open(cleanRead1, "w", newline='')
    for i, line in enumerate(open(inp1)):
        if i%4==0:
            if line[1:].split("/")[0].split()[0] in humanRead:
                skip=True
            else:
                skip=False
        else:
            continue
        if skip==False:
            clean1out.write(line.rstrip()+"\n")
    # do it for read2
    cleanRead2 = os.path.join(cleanOut, isolate+"_R2_clean_reads.fastq")
    skip=False
    clean2out = open(cleanRead2, "w",newline='')
    for i, line in enumerate(open(inp2)):
        if i%4==0:
            if line[1:].split("/")[0].split()[0] in humanRead:
                skip=True
            else:
                skip=False
        if skip==False:
            clean2out.write(line.rstrip()+"\n") """
    os.rmdir(tmpOut)
##############################################################################
'''function to remove human reads and keep non-human reads'''
def remove_human(bmtag,reads):
    humanReads = {}
    for line in open(bmtag):
        humanReads[line.strip()]=None
    skip = False
    for i , line in enumerate(open(reads)):
        if i%4==0:
            if line[1:].split("/")[0].split[0] in humanReads:
                skip = True
            else:
                skip = False
        if skip == False:
            print(line.rstrip())
# check if the reference is indexed
##############################################################################
'''kraken and bracken classification.'''
def krabracken(isolate,outpt, read1, read2):
    #krakOut = os.path.join(outpt,isolate,"Kraken")
    samOut = os.path.join(outpt, isolate)
    Krak = (KRAK, "--db", KRAKDB, "--paired", read1, read2, "--threads", str(THREADS),
            "--output", samOut+"_All_classifications.tsv",
            "--report", samOut+"_fullreport.txt", "--use-names",
            "--unclassified-out", samOut+"_unclassified#.fastq",
            "--classified-out", samOut+"_classified#.fastq", 
            "--minimum-hit-groups", str(KRAK_THRESH), "--report-minimizer-data")
    runKrak = ' '.join(Krak)
    tellU = "Running Kraken using the command: " + runKrak
    comms("tell",tellU)
    print(runKrak)
    #
    subprocess.call(runKrak, shell=True)
    #
    coCol = ("cut -f1-3,6-8", samOut+"_fullreport.txt", ">", samOut+"_classicReport.txt")
    runCoCol = ' '.join(coCol)
    subprocess.call(runCoCol, shell=True)
    #
    tellU = "Kraken completed and saved in: " + samOut
    comms("tell",tellU)
    print(runCoCol)
    #
    tellU = "Starting bracken"
    comms("tell",tellU)
    Brak = (BRAK, "-d", KRAKDB, "-i", samOut+"_fullreport.txt", "-o", 
            samOut+"bracken_fullreport.txt", "-t", str(BRAKTHRESH), "-w",
            samOut+"bracken_classicreport.txt", "-r", str(BRAKLENGTH))
    runBrak = ' '.join(Brak)
    tellU = "Bracken running with the command: " + runBrak
    comms("tell",tellU)
    print(runBrak)
    #
    subprocess.call(runBrak,shell=True)
    #
##############################################################################
if __name__ == '__main__':
    # check if inputs are gzipped
    for i in Path(os.path.abspath(INPDIR)).glob('*'):
        for t in os.listdir(str(i)):
            #print(os.path.basename(t))
            zfile = os.path.join(i, t)
            print(zfile)
            if (zfile.__contains__(".gz")):
                # decompress with pigz
                my_message = "Decompressing the input FASTQ files"
                comms("announce", my_message)
                pggz = ("pigz -d", zfile, "-p", str(THREADS))
                runPggz = ' '.join(pggz)
                tellU = "Decompressing: " + zfile
                comms("tell",tellU)
                print(runPggz)
                subprocess.call(runPggz,shell=True)
            else:
                my_message = "File is not zipped. Proceeding with workflow"
                comms("announce",my_message)
    my_message = "Starting PreQC analysis of raw reads using FastQC"
    comms("announce",my_message)
    #
    fastqc("pre",INPDIR,OUTDIR,THREADS)
    #
    my_message = "PreQC FastQC complete"
    comms("announce",my_message)
    my_message = "Starting trimming using Trim-Galore"
    comms("announce",my_message)
    #
    trimmy(INPDIR,OUTDIR,THREADS)
    #
    my_message = "Trimming completed. Trimmed reads are in the 'TrimmedReads' folder"
    comms("announce", my_message)
    my_message = "Removing Human DNA contamination"
    comms("announce", my_message)
    # removing human reads
    for i in Path(os.path.abspath(OUTDIR)).glob('*'):
        outfol = os.path.join(i, "BMTAG")
        if os.path.exists(outfol):
            pass
        else:
            os.makedirs(outfol)
        sample = os.path.basename(i)
        print(sample)
        read1 = os.path.join(i, "TrimmedReads", sample+"_R1_trimmed.fastq")
        read2 = os.path.join(i, "TrimmedReads", sample+"_R2_trimmed.fastq")
        #
        bmtagAln(i,REFRENCE,read1,read2,outfol,sample)
        #
        my_message = "Human reads removed and clean reads are in: " + os.path.abspath(outfol)
        comms("tell", my_message)
        my_message = "Starting PostQC analysis of raw reads using FastQC"
        comms("announce", my_message)
        # not working yet for some reason!!!!!
        postQC = os.path.join(i, "PostQC_FastQC")
        fastqc("post",outfol,postQC,THREADS)
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tellU = "FastQC PostQC analysis complete. Report are in: " + os.path.abspath(postQC)
        comms("tell", tellU)
        my_message = "Starting taxonomic classification with Kraken"
        comms("announce", my_message)
        cleanedReads = os.path.join(i, "CleanReads")
        krakenOut = os.path.join(i, "Kraken_Bracken")
        if os.path.exists(krakenOut):
            pass
        else:
            os.makedirs(krakenOut)
        karead1 = os.path.join(cleanedReads, sample+"_R1_clean_reads.fastq")
        karead2 = os.path.join(cleanedReads, sample+"_R2_clean_reads.fastq")
        krabracken(sample, krakenOut, karead1, karead2)
    