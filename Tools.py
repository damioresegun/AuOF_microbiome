#!/usr/bin/env python3
'''Script to hold various small tools and recurring commands
Author: Damilola Oresegun
'''

import subprocess
from concurrent.futures import thread
from distutils.archive_util import make_archive
import os
from pathlib import Path
import shutil
from tabnanny import check


def makeDirectory(directory):
    '''Function checks if a directory exists. If so, nothing is done,
    if the folder does not exist, the folder and any parent folder,
    is also made
    Input: Path to directory to be checked/made
    Output: None
    Usage: makeDirectory(path/to/directory)
    '''
    if os.path.exists(directory):
        pass
    else:
        os.makedirs(directory)


def prechecks(input, output):
    if os.path.exists(input):
        inp = "Good"
        pass
    else:
        inp = "Failed"
    if os.path.exists(output):
        out = "Good"
    else:
        out = "Make"
    return inp, out


def zipFile(func, file, threads):
    '''Compress or decompress a file using pigz. Takes in
    a function type (compress or decomp), the file to be zipped
    and the number of threads
    Options:
    func: compress or decomp
    file: unzipped file
    threads: number of threads for pigz to use
    Usage: zipCheck = zipFile("compress", file.fastq, 12)
    '''
    if func == "compress":
        try: 
            pggz = ' '.join("pigz --best", file, "-p", str(threads))
            print(pggz)
            subprocess.call(pggz, shell = True)
            zipped = "Done"
        except (FileNotFoundError, FileExistsError):
            zipped = "Error. File not found or file does not exist. Please check your inputs again"
    if func == "decomp":
        try:
            pggz = ' '.join("pigz -d", file, "-p", str(threads))
            print(pggz)
            subprocess.call(pggz, shell = True)
            zipped = "Done"
        except (FileNotFoundError, FileExistsError):
            zipped = "Error. File not found or file does not exist. Please check your inputs again"
    return zipped


def fastqc(phase, inpt, outpt, threads):
    ''' function to run fastqc for pre and post trimming
    usage: fastqc(input read folder, output/preQCm, 36)
    can also be output/postQC
    '''
    if phase == "pre":
        folName = "PreQC_FastQC"
        for file in Path(inpt).glob('*'):
            fname = os.path.basename(file).split("_")[0]
            fastOut = os.path.join(outpt, fname, folName)
            makeDirectory(fastOut)
            # get the reads
            fasfi1 = os.path.join(file, "*_R1*")
            fasfi2 = os.path.join(file, "*_R2*")
            # build the command
            fatq1 = ' '.join(["fastqc -t", str(threads), "-o", fastOut, "-f fastq", fasfi1, fasfi2])
            print(fatq1)
            subprocess.call(fatq1, shell=True)
    elif phase == "post":
        folName = "PostQC_FastQC"
        for file in Path(inpt).glob('*'):
            fname = os.path.dirname(inpt).split("/")[-1]
            fasfi1 = os.path.join(inpt, fname+"_R1_clean_reads.fastq")
            fasfi2 = os.path.join(inpt, fname+"_R2_clean_reads.fastq")
            fastOut = os.path.join(outpt, fname, folName)
            makeDirectory(fastOut)
            # build the command
            fatq1 = ' '.join(["fastqc -t", str(threads), "-o", fastOut, "-f fastq", fasfi1, fasfi2])
            print(fatq1)
            subprocess.call(fatq1, shell=True)

def trimmy(inpt, outpt, threads):
    '''Function to run trim-galore for trimming forward and reverse reads'''
    for folder in Path(inpt).glob('*'):
        fname = os.path.basename(folder).split("_")[0]
        trmOut = os.path.join(outpt, fname, "TrimmedReads")
        makeDirectory(trmOut)
        # get the reads
        fasfi1 = os.path.join(folder, "*_R1*")
        fasfi2 = os.path.join(folder, "*_R2*")
        # cutadapt only needs 8 cores maximum
        if (threads > 8):
            threads = 8
        # make the command
        runTrm = ' '.join(["trim_galore --no_report_file --paired --cores", str(threads),
                        "-o", trmOut, fasfi1, fasfi2])
        print("Command for trim-galore is: " + runTrm)
        #subprocess.call(runTrm, shell = True)
        print("Trim-galore complete")
        # now rename the files
        for file in Path(trmOut).glob('*'):
            print(str(file))
            if (str(file).__contains__("R1")):
                outFile = os.path.join(trmOut, fname + "_R1_trimmed.fastq")
                shutil.copy2(str(os.path.abspath(file)), os.path.abspath(outFile))
            if (str(file).__contains__("R2")):
                outFile = os.path.join(trmOut, fname + "_R2_trimmed.fastq")
                shutil.copy2(str(os.path.abspath(file)), os.path.abspath(outFile))
    return trmOut


def bmtagAligner(inpDir, reference, memory):
    '''Function to remove human reads. Takes in trimmed reads from TrimmedReads folder'''
    for folder in Path(os.path.abspath(inpDir)).glob('*'):
        if (str(folder).__contains__(".")):
            continue
        else:
            outfol = os.path.join(folder, "BMTAG")
            makeDirectory(outfol)
            sample = os.path.basename(folder)
            read1 = os.path.join(folder, "TrimmedReads", sample + "_R1_trimmed.fastq")
            read2 = os.path.join(folder, "TrimmedReads", sample + "_R2_trimmed.fastq")
            checkRef = reference.split(".fa")[0]
            # check if the reference was indexed
            if os.path.exists(os.path.abspath(checkRef) + ".bitmask"):
                refPres = "Reference has been indexed already"
            else: 
                refPres = "Reference has not been indexed. Indexing will begin now"
                runBmInd = ' '.join(["bmtool -d", reference, "-o", os.path.abspath(checkRef) + ".bitmask"])
                runBmPsm = ' '.join(["sprism mkindex -i", reference, "-o", os.path.abspath(checkRef) + ".srprism", 
                                    "-M", memory])
                print(runBmInd)
                subprocess.call(runBmInd, shell = True)
                print(runBmPsm)
                subprocess.call(runBmPsm)
            # make temporary output folder
            tmpOut = os.path.join(outfol, "tempOut")
            makeDirectory(tmpOut)
            # point to the indexed reference files
            refr = os.path.abspath(checkRef) + ".bitmask"
            refx = os.path.abspath(checkRef) + ".srprism"
            # set output file
            alOutFile = outfol + "/bmtagger.list"
            # start alignment
            bmTag = ' '.join(["bmtagger.sh -b", refr, "-x", refx, "-T", tmpOut, "-q 1 -1", read1, "-2", read2, 
                            "-o", alOutFile])
            print(bmTag)
            subprocess.call(bmTag, shell = True)
            # get the human reads, skip them and write the non-human reads to file
            cleanOut = os.path.join(folder, "CleanReads")
            makeDirectory(cleanOut)
            # make a dictionary to hold human reads
            humanRead = {}
            for line in open(alOutFile):
                humanRead[line.strip()] = None
            cleanReads1 = os.path.join(cleanOut, sample + "_R1_clean_reads.fastq")
            # Print out the fastq file line by line unless the read is human
            skip = False
            cleanOut1 = open(cleanReads1, "w", newline='')
            for i, line in enumerate(open(read1)):
                if i%4 == 0:
                    if line[1:].split("/")[0].split()[0] in humanRead:
                        skip = True
                    else:
                        skip = False
                if skip == False:
                    cleanOut1.write(line.rstrip() + "\n")
            cleanReads2 = os.path.join(cleanOut, sample + "_R2_clean_reads.fastq")
            skip = False
            cleanOut2 = open(cleanReads2, "w", newline='')
            for i, line in enumerate(open(read2)):
                if i%4 == 0:
                    if line[1:].split("/")[0].split()[0] in humanRead:
                        skip = True
                    else:
                        skip = False
                if skip == False:
                    cleanOut2.write(line.rstrip() + "\n")
            os.rmdir(tmpOut)
    return refPres, cleanOut, outfol



        

    
