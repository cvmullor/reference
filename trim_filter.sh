#! /bin/bash
#
# Loop for quality filtering and trimming on paired reads in a directory with prinseq
# Written by Carlos Valiente Mullor, 2019
#
# Usage:
#		# Execute in the directory containing R1 and R2 reads in fastq format
#		./trim_filter.sh

ls *.fastq* \
 | awk '{ if ($1 ~ /R2/) {print $0"\n"} else {print $0" "}}' ORS="" \
 | while read fastqR1 fastqR2 
	do
		perl prinseq-lite.pl -fastq $fastqR1 -fastq2 $fastqR2 -min_len 50 -min_qual_mean 20 -ns_max_p 10 -trim_qual_right 20 2>>stats.txt
	done
