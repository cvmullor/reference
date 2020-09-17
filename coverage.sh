#! /bin/bash

# Loop for calculating mean read coverage on each position of the 
# reference for each mapping
# Written by Carlos Valiente Mullor, 2019
#
# Usage:
#		./s2_coverage.sh
#
# Input: sorted BAM files

ls -1 *.sorted.bam | sed 's/.sorted.bam//g' | while read bamfile;
do
	## Positions with >0 reads mapped
	# n reads mapping against each position
	samtools depth -aa $bamfile.sorted.bam | awk '{ if ($3>0) {print $3} }' > $bamfile.nozero.temp
	# average coverage
	awk -v filename=$bamfile.nozero.temp '{ sum += $0; n++ } END { if (n > 0) print (sum / n)" "filename; }' $bamfile.nozero.temp >> mcov_nozero.csv
	
	## All positions
	samtools depth -aa $bamfile.sorted.bam | awk '{ if ($3>=0) {print $3} }' > $bamfile.all.temp
	awk -v filename=$bamfile.all.temp '{ sum += $0; n++ } END { if (n > 0) print (sum / n)" "filename; }' $bamfile.all.temp >> mcov_all.csv
done

