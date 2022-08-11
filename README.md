# README
<!-- toc -->
- [AuOF](#AuOF)
	- [Overview](##Overview)
	- [Installation](##Installation)
		- [Beware](###Beware)
		- [Package](###Package)
		- 
<!-- tocstop -->
# AuOF
## Overview
A pipeline to carry out processing and analysis of metagenomic whole genome sequences as part of the Antibiotics under our feet project. While the pipeline is able to automate a majority of the process, downstream analyses still requires further manual intervention. Additionally user input is required to configure the running the pipeline. The pipeline is primarily written in python3 with the running script written in Bash/Shell. The pipeline takes in paired-end Illumina short read sequences from a metagenomic source. It has been extensively tested and developed using Illumina short read sequences from metagenomic soil samples.
## Installation
### Beware
- The pipeline was developed using the full NCBI RefSeq database for taxonomic profiling. This means ~130GB of storage is required for a similar set up.
- Due to the large size of the database used for the taxonomic profiling, the pipeline was tested with 150GB of RAM because the database is loaded into memory for the profiling step
- It is not necessary to use the full RefSeq database for taxonomic profiling however, it IS necessary that the amount of memory used for the pipeline is able to hold your database in memory with overhead to allow other processes to occur at the same time. 
### Package
- The forked KrakenTools package from https://github.com/damioresegun/KrakenTools
- To install:
	- **wget**
		```bash
		# download
		wget https://github.com/damioresegun/KrakenTools/archive/refs/heads/master.zip
		# unzip
		unzip master.zip
		```
	- **git**
		```
		git clone git@github.com:damioresegun/KrakenTools.git
		```
### CONDA/MAMBA

The rest of the tools are easily installed using conda. If you have mamba, it is recommended to use mamba as it is quicker to use.

#### Install MAMBA
`conda install mamba -n base -c conda-forge`

#### Install packages
##### Tools and Versions
| Tools          | Version |
| -------------- | ------- |
| BMTAGGER       | 3.101   |
| Assembly-stats | 0.1.1   | 
| Bracken        | 2.6.2   |
| Cutadapt       | 1.18    |
| FastQC         | 0.11.9  |
| Humann         | 3.1.1   |
| Kraken2        | 2.1.2   |
| Krona          | 2.8.1   |
| Pigz           | 2.6     |
| Pip            | 22.1.1  |
| Python         | 3.7.12  |
| R              | 3.2.2   |
| Samtools       | 0.1.19  |
| Trim-galore    | 0.6.7   |
### Install databases
#### Kraken
- Below is a guide to download and install the full RefSeq nucleotide database for Kraken taxonomic classification
```bash
# make the folder
mkdir Kraken_DB
# download the taxonomy files
kraken-build --download-taxonomy --db Kraken_DB --threads 12
# download the libraries you want for your database
# bacteria
kraken2-build --download-library bacteria --db Kraken_DB --threads 12 --no-masking
# viruses
kraken2-build --download-library viral --db Kraken_DB --threads 12 --no-masking 
# protozoa
kraken2-build --download-library protozoa --db Kraken_DB --threads 12 --no-masking 
# fungi
kraken2-build --download-library fungi --db Kraken_DB --threads 12 --no-masking 
# archaea
kraken2-build --download-library archaea --db Kraken_DB --threads 12 --no-masking 
# plant
kraken2-build --download-library plant --db Kraken_DB --threads 12 --no-masking 
# plasmids
kraken2-build --download-library plasmids --db Kraken_DB --threads 12 --no-masking
```
- In some instances, an error might occur during the installation of the database after running the `kraken2-build` command:
```bash
ERROR: "rsync_from_ncbi.pl: unexpected FTP path (new server?) for https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/128/725/GCF_900128725.1_BCifornacula_v1.0"
```
- This can be solved by opening the `rsync_from_ncbi.pl` script that is installed along with kraken2 and change line below:
```perl
# original line
if (! ($full_path =~ s#^ftp://${qm_server}${qm_server_path}/##)) { 
# changed to:
if (! ($full_path =~ s#^https://${qm_server}${qm_server_path}/##)) {
```
- Other errors that can come up during installation is: `rsync_from_ncbi.pl: unexpected FTP path (new server?) for na`
	- To fix this, the `assembly_summary.txt` that is created from the `kraken2-build` command has to be fixed. This means:
	``` bash
		# change into the library affected. Here it is in the bacteria
		# there should be a file named: assembly_summary.txt
		awk -v FS='\t' '$20 != "na" {print $0}' assembly_summary.txt > new_assembly_summary.txt 
		cp new_assembly_summary.txt assembly_summary.txt
		```
#### Bracken
Once bracken is installed with conda and the kraken database has been downloaded as described above, the bracken index can be built with:
```bash
bracken-build -d Kraken_DB/ -t 24 -k 35 -l 100
# -t is the number of threads used
# -k is the kmer length to use in the re-estimation process of bracken
# -l minimum sequence length in your sequences
```
**Note:** 
## Pipeline Options
### Required Parameters
- `-i/--input`: The path to the directory holding the demultiplexed FASTQ files. FASTQ files can be gzipped or left uncompressed. **Note:** The FASTQ files have to be named with _1 and _2
- `-o\--output`: Full path to the directory to save analysis outputs. This folder should either be empty or not created yet. The pipeline will check if the folder is already created -- if so, the pipeline will continue only IF the folder is empty
- `-r/--reference`: Full path to the human reference genome to use for decontamination. In theory, this can be any other reference genome and the script can be tweaked to accept more than one reference genome for decontamination
- `-kr/--kraken`: Full path to your kraken installation if it is not in the $PATH. If in your $PATH, simply write `-kr kraken2`
- `-br/--bracken`: Full path to your bracken installation if it is not in the $PATH. If in your $PATH, simply write `-br bracken`
- `-kb/--kraken_DB`: Full path to the kraken database you have built
- `-kt/--kraken_tools`: Full path to the KrakenTools package forked by DRO
### Optional Parameters
- `-t/--threads`: Number of threads. Default is `24`
- `-m/--max-memory`:  Maximum available memory to use in megabytes e.g. for `100G` of memory, enter `100000`. Default is `24000` i.e. `24G`
- `-kt/--kraken_hit_threshold`: A minimum number of groups that must be matched to place a contig into a taxonomic group. Default is `5`
- `-bt/--bracken_hit_threshold`: A minimum number of *k-mers* that must be matched to place a contig into a taxonomic group by bracken re-estimation. Default is `20`
- `-bl/--bracken_read_length`: The read length used to build your bracken database. **Note:** Your bracken database is also the kraken database. Default length is `100`
- `-f/--functional`: Full path to the humann v3+ package if not in $PATH. If in $PATH, this parameter is not necessary. Default is `humann`

## Pipeline Running Script
