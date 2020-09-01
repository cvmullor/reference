#! /bin/bash

# Loop for indexing the reference sequences (fasta) from a species
# 	with BWA and samtools faidx
# Written by Carlos Valiente Mullor, 2019
#
# Usage:
#		./m1_index.sh <path/to/references>
#
# Reference sequences in fna (fasta) format
# Empty arg = actual directory

refsp=$1

mkdir all_indexes
cp $refsp*.fna* all_indexes

ls all_indexes/*.fna | while read file;
do
	bwa index -a bwtsw $file
	samtools faidx $file
done
