#! /bin/bash

# Loop for mapping paired reads from samples in a directory
#	to different reference sequences with BWA MEM and sorting 
#	the resulting files
#	Written by Carlos Valiente Mullor, 2019
#
# Usage:
#		./m2_mapsort <path/to/reads>
#
# R1 and R2 reads in fastq format (not compressed)

readsp=$1 # path/to/reads

# Mapping and sorting
ls $readsp*.fastq \
 | sed -e "s|$readsp||g" \
 | sed 's/.fastq//g' \
 | awk '{ if ($1 ~ /R2/) {print $0"\n"} else {print $0" "}}' ORS="" \
 | while read fastqR1 fastqR2;
do
	ls all_indexes/*.fna \
	 | sed -e 's/all_indexes\///g' \
	 | sed 's/.fna//g' \
	 | while read reference;
	do
		~/tfm_map/softmap/fbwa/bwa mem -t 8 all_indexes/$reference.fna $readsp$fastqR1.fastq $readsp$fastqR2.fastq 2>nohup_bwa.out \
		 | samtools sort -@8 -o $fastqR1.$reference.sorted.bam - 2>nohup_sort.out
	done
done;

# Indexing sorted BAM files
ls *.sorted.bam \
 | while read srtfile;
do 
	samtools index $srtfile
done
