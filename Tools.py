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
        if not os.listdir(output):
            out = "Good"
        else:
            out = "Failed"
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
        for folder in Path(inpt).glob('*'):
            fname = os.path.basename(folder)
            fasfi1 = os.path.join(inpt, fname, "CleanReads", fname+"_R1_clean_reads.fastq")
            fasfi2 = os.path.join(inpt, fname, "CleanReads", fname+"_R2_clean_reads.fastq")
            fastOut = os.path.join(inpt, fname, folName)
            makeDirectory(fastOut)
            # build the command
            fatq1 = ' '.join(["fastqc -t", str(threads), "-o", fastOut, "-f fastq", fasfi1, fasfi2])
            print(fatq1)
            subprocess.call(fatq1, shell=True)
    return fastOut

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


def krakbracken(isolate, inpt, outpt, krkn, brkn, krkthrs, krkdb, brkthrs, brklen, threads):
    """
    Function to carry out kraken classification then bracken abundance re-estimation
    This is then followed by a conversion of the bracken report to krona visulisation
    files in the form of HTML files. Additionally, the bracken report is also converted
    to the MPA format for use downstream in HUMAnN.

    Parameters:
    - isolate: isolate name
    - inpt: folder containing clean reads
    - outpt: root output folder. Kraken and Bracken folders will be made here
    - krkn: path to kraken2 or the particular kraken2 command in $PATH
    - brkn: path to bracken or the particular bracken command in $PATH
    - krkthrs: the kraken threshold
    - krkdb: path to the kraken database
    - brkthrs: bracken threshold
    - brklen: bracken minimum read/contig length
    - threads: the number of threads to use
    """
    # check the output folder
    krakOut = os.path.join(outpt, "Kraken")
    brakOut = os.path.join(outpt, "Bracken")
    makeDirectory(krakOut)
    makeDirectory(brakOut)
    krakOut = os.path.join(krakOut, isolate)
    brakOut = os.path.join(brakOut, isolate)
    read1 = os.path.join(inpt, isolate+"_R1_clean_reads.fastq")
    read2 = os.path.join(inpt, isolate+"_R1_clean_reads.fastq")
    runKrak = ([krkn, "--db", krkdb, "--paired", read1, read2, "--threads", str(threads),
            "--output", krakOut+"_All_classifications.tsv",
            "--report", krakOut+"_fullreport.txt", "--use-names",
            "--unclassified-out", krakOut+"_unclassified#.fastq",
            "--classified-out", krakOut+"_classified#.fastq", 
            "--minimum-hit-groups", str(krkthrs), "--report-minimizer-data"])
    print(runKrak)
    subprocess.call(runKrak, shell = True)
    print("Kraken completed. Continuing to bracken")
    runBrak = ([brkn, "-d", krkdb, "-i", krakOut+"_fullreport.txt", "-o", 
            brakOut+"_bracken_fullreport.txt", "-t", str(brkthrs), "-w",
            brakOut+"_bracken_classicreport.txt", "-r", str(brklen)])
    print(runBrak)
    subprocess.call(runBrak, shell = True)
    print("Bracken complete. Generating krona plots")
    # generate krona plots
    Runkrnny = (["kreport2krona.py -r", brakOut+"_bracken_classicreport.txt", 
            "-o", brakOut+"_bracken_classicreport.krona"])
    RunKtny = (["ktImportText", brakOut+"_bracken_classicreport.krona", 
            "-o", brakOut+"_bracken_classicreport.html"])
    print(Runkrnny)
    print(RunKtny)
    subprocess.call(Runkrnny, shell=True)
    subprocess.call(RunKtny, shell=True)
    # generate reports in mpa format
    Runmpp = (["kreport2mpa3.py -r", brakOut+"_bracken_classicreport.txt", 
                "-o", brakOut+"_bracken_classicreport.tsv", "-hm",
                "--percentages --display-header"])
    subprocess.call(Runmpp, shell = True)
    # output files
    brakRep = brakOut+"_bracken_classicreport.txt"
    kroFle = brakOut+"_bracken_classicreport.html"
    mpaFle = brakOut+"_bracken_classicreport.tsv"
    return brakRep, kroFle, mpaFle




    
