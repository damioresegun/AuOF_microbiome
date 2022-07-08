#!/bin/bash
# Bash script to run the AuOF master control script during development
# enter variable
AuOF=/home/doresegu/scratch/private/CMC_Project/AuOF_MasterControl.py
input=/home/doresegu/scratch/private/CMC_Project/Demultiplexed/Pin25307
output=/home/doresegu/scratch/private/CMC_Project/Testing
reference=/home/doresegu/scratch/private/CMC_Project/BMTAGGER_INDEX/hg38.fa
kraken="kraken2"
bracken="bracken"
krakenDB=/home/doresegu/scratch/private/Kraken_DB
threads=24
memory=150000

mkdir $output
touch ${output}/logFile.txt
run="python3 $AuOF -i $input -o $output -r $reference -kr $kraken -br $bracken -kb $krakenDB -t $threads -m $memory 2>&1 | tee ${output}/logFile.txt"
echo $run
eval $run
