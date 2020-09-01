#! /bin/bash

# Loop for removing adaptors using Cutadapt in all fastq file in a directory
# Writen by Carlos Valiente Mullor, 2019

# Usage:
#		# Execute in the directory containing R1 and R2 reads in fastq format
#		./remove_adaptors.sh
#
# "adaptors.fasta": file with the adaptor sequence to remove

ls *.fastq* \
 | awk '{ if ($1 ~ /R2/) {print $0"\n"} else {print $0" "}}' ORS="" \
 | while read fastqR1 fastqR2 
	do
		cutadapt -a file:adaptors.fasta -A file:adaptors.fasta -o $fastqR1.1.cut.fq -p $fastqR2.2.cut.fq $fastqR1 $fastqR2
	done
