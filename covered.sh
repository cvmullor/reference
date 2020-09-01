#! /bin/bash

# Loop for calculating percentage of reference covered by reads for each mapping
# Written by Carlos Valiente Mullor, 2019
#
# Usage:
#		./s2_covered.sh
#
# Input: sorted BAM files

ls -1 *.sorted.bam | while read bamfile;
do
	# Positions with >0 reads mapped
	nonzero=$(samtools depth -aa $bamfile | awk '{if($3>0) total+=1}END{print total}')
	# All positions
	all=$(samtools depth -aa $bamfile | awk '{if($3>=0) total+=1}END{print total}')
	percent=$(bc <<< "scale=6; ($nonzero / $all)*100")
	echo "$bamfile;$percent;$all" | sed 's/.sorted.bam//g' >> core_size.csv
done
