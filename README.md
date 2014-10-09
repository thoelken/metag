metaG
=====

preprocessing and classification pipline for metagenomic data analysis

*This repo is still a stub for further developement of a broader metagenome project without further working title. In any case, this is probably not the code you are looking for. Anyhow, feel free to 
look at the code or use it otherwise.*

## fastaseq.cpp
The program currently takes a FASTA file as input and prints four random sequences of length 99 somewhere from the entries. For repetitative access, an index with the file extension ".idx" is placed next to the FASTA file.

### Installation
Compile the program with:

	$ g++ -std=c++11 -o fastseq fastaseq.cpp
