#! /bin/bash

# Loop for making congruence tests for each MSA in the current directory against 
#	a set of phylogenetic trees
# Written by Carlos Valiente Mullor, 2019
#
# Usage:
#       ./congruence_test <trees_file.nwk>
#
# MSAs files should be fasta

trees=$1

ls *.fas | while read aln;
do
	nohup iqtree -nt AUTO -s $aln -z $trees -m GTR -zb 10000 -zw &
done
