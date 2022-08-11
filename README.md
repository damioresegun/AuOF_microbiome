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
#### Bracken


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
- `-m/--max-memory`:  Maximum available memory to use in megabytes e.g. for `100G` of memory, enter `100000`. Default is `240000` i.e. `24G`
- `-kt/--kraken_hit_threshold`: A minimum number of groups that must be matched to place a contig into a taxonomic group. Default is `5`