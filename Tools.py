#!/usr/bin/env python3
'''File to hold various tools and recurring functions for the AuOF pipeline
Author: Damilola Oresegun
'''
# import packages
import subprocess
import os
from pathlib import Path
import shutil
from os import system


# state the functions
def makeDirectory(directory):
    '''Function checks if a directory exists. If so, nothing is done,
    if the folder does not exist, the folder and any parent folder,
    is also made
    Input: Path to directory to be checked/made
    Output: None
    Usage: makeDirectory(path/to/directory)
    Returns: Nothing
    '''
    if os.path.exists(directory):
        print(directory + " already exists")
        pass
    else:
        print(directory + " does not exist. Generating...")
        os.makedirs(directory)


def prechecks(input, output):
    '''
    Functions checks input files and output locations. 
    If input file exists, return 'Good', if not, return 
    'Failed'.
    If output path exists, check if the folder is empty.
    If the folder is empty, return 'Good', or else 'Failed'.
    If the output path does not exist, make the folder.
    Usage: prechecks(path/to/input/file, path/to/output)
    Returns: Two strings: Good, Failed or Make
    '''
    if os.path.exists(input):
        print(input + " exists. Proceeding...")
        inp = "Good"
        pass
    else:
        print(input + " does not exist")
        inp = "Failed"
    if os.path.exists(output):
        print(output + " exists. Proceeding...")
        # if the folder is empty
        if not os.listdir(output):
            print(output + " is empty. Proceeding...")
            out = "Good"
        # if not empty
        else:
            print(output + " exists but is not empty. Cannot proceed.")
            out = "Failed"
    # if the folder doesnt exist
    else:
        print(output + " does not exist. Generating...")
        out = "Make"
    return inp, out


def zipFile(func, file, threads):
    '''Compress or decompress a file using pigz. Takes in
    a function type (compress or decomp), the file to be zipped
    and the number of threads.
    Options =
    func: compress or decomp
    file: unzipped file
    threads: number of threads for pigz to use
    Usage: zipCheck = zipFile("compress", file.fastq, 12)
    Returns: One string: 'Done' or 'File not found' error
    '''
    # to compress files
    if func == "compress":
        try:
            print("Compressing " + file) 
            pggz = ' '.join("pigz --best", file, "-p", str(threads))
            # print the command for log file
            print("Command for compressing is: " + pggz)
            subprocess.call(pggz, shell = True)
            zipped = "Done"
            print("Compression completed.")
        except (FileNotFoundError, FileExistsError):
            print("Error. " + file + " not found or does not exist")
            zipped = "Error. File not found or file does not exist. Please check your inputs again"
    # to decompress files
    if func == "decomp":
        try:
            print("Decompressing " + file)
            pggz = ' '.join("pigz -d", file, "-p", str(threads))
            # print the command for log file
            print("Command for decompressing is: " + pggz)
            subprocess.call(pggz, shell = True)
            zipped = "Done"
        except (FileNotFoundError, FileExistsError):
            print("Error. " + file + " not found or does not exist")
            zipped = "Error. File not found or file does not exist. Please check your inputs again"
    return zipped


def fastqc(phase, inpt, outpt, threads):
    ''' Function to run fastqc for pre and post trimming.
    Requires a 'pre' or 'post' phase description. Pre is to run
    FastQC on raw reads and Post is to run FastQC on quality 
    assessed/trimmed reads. Requires input, output and number of 
    threads.
    Usage: fastqc(pre/post, input/read/folder, output/preQCm, 36)
            can also be output/postQC
    Returns: String: Path/to/isolate/fastqc/output
    '''
    # fastqc raw reads
    if phase == "pre":
        print("Starting FastQC on raw reads...")
        folName = "PreQC_FastQC"
        # for each file in the input folder
        for file in Path(inpt).glob('*'):
            # get the isolate name from the path
            fname = os.path.basename(file).split("_")[0]
            fastOut = os.path.join(outpt, fname, folName)
            # make the isolate output folder
            makeDirectory(fastOut)
            # get the reads
            fasfi1 = os.path.join(file, "*_R1*")
            fasfi2 = os.path.join(file, "*_R2*")
            # build the command
            fatq1 = ' '.join(["fastqc -t", str(threads), "-o", fastOut, "-f fastq", fasfi1, fasfi2])
            # print the command for log file
            print("Command for preQC FastQC is: " + fatq1)
            subprocess.call(fatq1, shell=True)
        print("FastQC completed successfully.")
    # fastqc quality assessed 
    elif phase == "post":
        print("Starting FastQC on cleaned reads...")
        folName = "PostQC_FastQC"
        # for each folder in the input folder path
        for folder in Path(inpt).glob('*'):
            # ignore files in the folder. Only want folders
            if (str(folder).__contains__(".")):
                continue
            else:
                # get the isolate name
                fname = os.path.basename(folder)
                # get the cleaned reads from the CleanReads folder
                fasfi1 = os.path.join(inpt, fname, "CleanReads", fname+"_R1_clean_reads.fastq")
                fasfi2 = os.path.join(inpt, fname, "CleanReads", fname+"_R2_clean_reads.fastq")
                fastOut = os.path.join(inpt, fname, folName)
                # make the output folder
                makeDirectory(fastOut)
                # build the command
                fatq1 = ' '.join(["fastqc -t", str(threads), "-o", fastOut, "-f fastq", fasfi1, fasfi2])
                # print the command for log file
                print("Command for postQC FastQC is: " + fatq1)
                subprocess.call(fatq1, shell=True)
        print("FastQC completed successfully.")
    return fastOut

def trimmy(inpt, outpt, threads):
    '''Function to run trim-galore for trimming forward and reverse reads
    Requires the input folder holding folders from earlier in the pipeline.
    Usage: trimmy(input/folder, where/to/make/output/folder, threads)
    Returns: String: Path/to/isolate/TrimmedReads
    '''
    print("Starting trim-galore for isolates in " + inpt)
    for folder in Path(inpt).glob('*'):
        fname = os.path.basename(folder).split("_")[0]
        trmOut = os.path.join(outpt, fname, "TrimmedReads")
        makeDirectory(trmOut)
        # get the reads
        fasfi1 = os.path.join(folder, "*_R1*")
        fasfi2 = os.path.join(folder, "*_R2*")
        # cutadapt only needs 8 cores maximum
        if (threads > 8):
            print("Your threads are greater than necessary for this step. Limiting to 8")
            threadr = 8
        else:
            threadr = threads
        # make the command
        runTrm = ' '.join(["trim_galore --no_report_file --paired --cores", str(threadr),
                        "-o", trmOut, fasfi1, fasfi2])
        # print the command for log file
        print("Command for trim-galore is: " + runTrm)
        subprocess.call(runTrm, shell = True)
        # print for log file
        print(fname + "trimming complete")
        # now rename the files
        for file in Path(trmOut).glob('*'):
            # rename the trimmed R1
            if (str(file).__contains__("R1")):
                outFile = os.path.join(trmOut, fname + "_R1_trimmed.fastq")
                # copy the file and rename
                shutil.copy2(str(os.path.abspath(file)), os.path.abspath(outFile))
                print("Trimmed forward reads have been renamed and saved as " + outFile)
            if (str(file).__contains__("R2")):
                outFile = os.path.join(trmOut, fname + "_R2_trimmed.fastq")
                # copy the file and rename
                shutil.copy2(str(os.path.abspath(file)), os.path.abspath(outFile))
                print("Trimmed reverse reads have been renamed and saved as " + outFile)
    print("Trim-galore completed successfully.")
    return trmOut


def bmtagAligner(inpDir, reference, memory):
    '''Function to remove human reads from input sequence data.
    Takes in trimmed reads from TrimmedReads folder. Makes the
    output folder in the parent folder of the input folder.
    Requires the reference genome and a maximum memory usage limit

    Usage: bmtagAligner(path/to/trimmed/reads, path/to/reference, max_Memory)
    Returns: Three strings:
        1 - Reference indexed or not indexed
        2 - Path to CleanReads output
        3 - Path to the alignment output
    '''
    print("Carrying out alignment using BMTAGGER...")
    for folder in Path(os.path.abspath(inpDir)).glob('*'):
        # if a file is in the folder, ignore it
        if (str(folder).__contains__(".")):
            continue
        else:
            outfol = os.path.join(folder, "BMTAG")
            # make a folder to hold the alignment outuput
            makeDirectory(outfol)
            sample = os.path.basename(folder)
            # get the trimmed reads
            read1 = os.path.join(folder, "TrimmedReads", sample + "_R1_trimmed.fastq")
            read2 = os.path.join(folder, "TrimmedReads", sample + "_R2_trimmed.fastq")
            checkRef = reference.split(".fa")[0]
            # check if the reference was indexed
            if os.path.exists(os.path.abspath(checkRef) + ".bitmask"):
                refPres = "Reference has been indexed already"
                print(refPres)
            else: 
                refPres = "Reference has not been indexed. Indexing will begin now"
                print(refPres)
                # index the reference using bmtool
                runBmInd = ' '.join(["bmtool -d", reference, "-o", os.path.abspath(checkRef) + ".bitmask"])
                runBmPsm = ' '.join(["srprism mkindex -i", reference, "-o", os.path.abspath(checkRef) + ".srprism", 
                                    "-M", memory])
                # print for log file
                print(runBmInd)
                subprocess.call(runBmInd, shell = True)
                # print for log file
                print(runBmPsm)
                subprocess.call(runBmPsm)
                print("Reference indexed.")
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
            # print for log file
            print(bmTag)
            subprocess.call(bmTag, shell = True)
            print("BMTAGGER complete. Proceeding to remove contaminant reads in read files")
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
        print("Reads are decontaminated in and saved in " + cleanOut)
    return refPres, cleanOut, outfol


def krakbracken(isolate, inpt, outpt, krkn, brkn, krkthrs, krkdb, brkthrs, brklen, threads, ktool):
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
    print("Starting Kraken classification...")
    # check the output folder
    krakOut = os.path.join(outpt, "Kraken")
    brakOut = os.path.join(outpt, "Bracken")
    makeDirectory(krakOut)
    makeDirectory(brakOut)
    krakOut = os.path.join(krakOut, isolate)
    brakOut = os.path.join(brakOut, isolate)
    read1 = os.path.join(inpt, isolate+"_R1_clean_reads.fastq")
    read2 = os.path.join(inpt, isolate+"_R2_clean_reads.fastq")
    # run the kraken command
    runKrak = ' '.join([krkn, "--db", krkdb, "--paired", read1, read2, "--threads", str(threads),
            "--output", krakOut+"_All_classifications.tsv",
            "--report", krakOut+"_fullreport.txt", "--use-names",
            "--unclassified-out", krakOut+"_unclassified#.fastq",
            "--classified-out", krakOut+"_classified#.fastq", 
            "--minimum-hit-groups", str(krkthrs), "--report-minimizer-data"])
    # print for log file
    print(runKrak)
    subprocess.call(runKrak, shell = True)
    print("Kraken completed. Continuing to bracken...")
    # run bracken
    runBrak = ' '.join([brkn, "-d", krkdb, "-i", krakOut+"_fullreport.txt", "-o", 
            brakOut+"_bracken_fullreport.txt", "-t", str(brkthrs), "-w",
            brakOut+"_bracken_classicreport.txt", "-r", str(brklen)])
    # print for log file
    print(runBrak)
    subprocess.call(runBrak, shell = True)
    print("Bracken complete. Generating krona plots...")
    # generate krona plots
    Runkrnny = ' '.join([ktool+"/kreport2krona.py -r", brakOut+"_bracken_classicreport.txt", 
            "-o", brakOut+"_bracken_classicreport.krona"])
    RunKtny = ' '.join(["ktImportText", brakOut+"_bracken_classicreport.krona", 
            "-o", brakOut+"_bracken_classicreport.html"])
    # print for log file
    print(Runkrnny)
    print(RunKtny)
    subprocess.call(Runkrnny, shell=True)
    subprocess.call(RunKtny, shell=True)
    print("Bracken completed. Converting bracken results")
    print("Convert bracken results to mpa style")
    # generate reports in mpa format
    Runmpp = ' '.join([ktool+"/kreport2mpa3.py -r", brakOut+"_bracken_classicreport.txt", 
                "-o", brakOut+"_bracken_classicreport.tsv", "-hm",
                "--percentages --display-header"])
    # print for log file
    print(Runmpp)
    subprocess.call(Runmpp, shell = True)
    print("MPA conversion successful")
    # output files
    brakRep = brakOut+"_bracken_classicreport.txt"
    kroFle = brakOut+"_bracken_classicreport.html"
    mpaFle = brakOut+"_bracken_classicreport.tsv"
    return brakRep, kroFle, mpaFle


def humann3(isolate, inpt, humann, krakmpa, outpt, threads):
    '''
    Function to carry out functional profiling using HUMAnN version 3.1.1
    Takes in the isolate name, the path to the folder containing
    paired-end reads and the path to the human executable. Additionally
    a taxonomic profile file in the metaphlan format is needed, the path 
    to an output folder and the number of threads is also needed.
    Returns the path to the output folder
    '''
    print("Starting HUMAnN functional profiling...")
    reads = [os.path.join(inpt, isolate+"_R1_clean_reads.fastq"), 
            os.path.join(inpt, isolate+"_R2_clean_reads.fastq")]
    # concatenate into one fastq file
    outReads = os.path.join(inpt, isolate+"_combinedReads.fastq")
    run = ' '.join(["cat", reads[0], reads[1], ">", outReads])
    print("Concatenating read files into one...")
    print(run)
    system(run)
    """ with open(outReads, 'wb') as allReads:
        for file in reads:
            with open(file, 'rb') as f:
                shutil.copyfileobj(f, allReads)
                #for line in f:
                #    allReads.write(line)
            allReads.write(b"\n") """
    humanOut = os.path.join(outpt, isolate, "HUMAnN")
    makeDirectory(humanOut)
    runHum = ' '.join([humann, "--input", outReads, "--output", humanOut, 
                "--taxonomic-profile", krakmpa, "--threads", str(threads),
                "--verbose"])
    print(runHum)
    subprocess.call(runHum, shell = True)
    print("HUMAnN completed")
    return humanOut    


def makeBiom(inpt,outpt):
    '''
    Function to take in a folder, search through the folder to find all
    instances of a folder named 'Bracken'. If found, copy the report
    into a new folder and rename using the isolate name. Then convert 
    each report into a biom file AND make a combined biom file of all 
    isolates in the input folder
    Input: tester
    Output: tester/bracken
    Returns: String: The name of the combined Biom file
    '''
    print("Converting files to biom")
    # make the output directory if it doesnt exist
    makeDirectory(outpt)
    for folder in Path(os.path.abspath(inpt)).glob('*'):
        fname = os.path.basename(folder)
        brakFol = os.path.join(folder, "Bracken", 
                        fname + "_bracken_classicreport.txt")
        if os.path.exists(os.path.abspath(brakFol)):
            runErr = "Good"
            print("Copying " + brakFol + " to " + outpt)
            cpfile = os.path.join(outpt, fname + ".txt")
            shutil.copy2(brakFol, cpfile)
            outfile = os.path.join(outpt, fname + ".biom")
            # convert the txt files to biom files
            runBiom = ' '.join(["kraken-biom", cpfile, "-o", outfile])
            # print for log file
            print(runBiom)
            subprocess.call(runBiom, shell = True)
            print("Conversion completed")
        else:
            runErr = "Failed"
            print(str(folder) + " did not contain Bracken files")
            continue
    # now make a combined biom file
    if runErr == "Good":
        print("Now making combined biom file")
        allTfiles = os.path.join(outpt, "*.txt")
        allBiom = os.path.join(outpt, "CombinedIsolate.biom")
        runComBiom = ' '.join(["kraken-biom", allTfiles, "-o", allBiom])
        print(runComBiom)
        subprocess.call(runComBiom, shell = True)
    else:
        print(str(folder) + " did not contain Bracken files")
    return allBiom, runErr