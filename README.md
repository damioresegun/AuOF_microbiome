# README
<!--toc-->
- [AuOF](#AuOF)
	- [Overview](##Overview)
	- [Installation](##Installation)
		- [Beware](###Beware)
		- [Package](###Package)
<!--toc-->

# AuOF
## Overview
A pipeline to carry out processing and analysis of metagenomic whole genome sequences as part of the Antibiotics under our feet project. While the pipeline is able to automate a majority of the process, downstream analyses still requires further manual intervention. Additionally user input is required to configure the running the pipeline. The pipeline is primarily written in python3 with the running script written in Bash/Shell. 
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
	- 
