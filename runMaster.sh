#!/bin/bash
# Bash script to run the AuOF master control script during development
# enter the path to the AuOF_MasterControl python script
AuOF=/home/doresegu/scratch/private/CMC_Project/AuOF_MasterControl.py
# enter the path to the folder holding the directories of the reads
input=/home/doresegu/scratch/private/CMC_Project/Demultiplexed
# full path to the output folder
output=/home/doresegu/scratch/private/CMC_Project/AuOF_Revamp2
# path to the reference genome for decontamination
reference=/home/doresegu/scratch/private/CMC_Project/BMTAGGER_INDEX/hg38.fa
# path to the kraken2 package. If in the $PATH then just enter kraken2
kraken="kraken2"
# path to the bracken package. If in the $PATH then just enter bracken
bracken="bracken"
# path to the kraken database
krakenDB=/home/doresegu/scratch/private/Kraken_DB
# number of threads
threads=24
# enter the maximum memory to use in megabytes
memory=150000
# path to the krakentools package directory
krakenTool=/home/doresegu/scratch/private/KrakenTools
# state the kraken hit threshold
krakThres=5
# state the bracken hit threshold
brakThres=20
# state the bracken read length
brakRead=100
# full path to the humann package or if in path enter humann
humann3="humann"
###############################################################
# No need to change below this
###############################################################
# this will automatically make a temporary folder
tempFolder=($(dirname $output)/tempFolder)
# makes the output directoru
mkdir $output
mkdir $tempFolder
# automatically makes a total logfile
touch ${tempFolder}/logFile.txt
# make a run command
run="python3 $AuOF -i $input -o $output -r $reference -kr $kraken -br $bracken -kb $krakenDB -t $threads 
        -m $memory -dk $krakenTool -kt $krakThres -bt $brakThres -bl $brakRead -f $humann3 2>&1 | tee ${tempFolder}/logFile.txt"
# print it to screen
echo $run
# run the command
eval $run
# transfer the logfile to the output directory
mv ${tempFolder}/logFile.txt $output
# delete the temporary folder
rm -r $tempFolder
