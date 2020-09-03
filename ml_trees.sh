#! /bin/bash

# Loop for building ML trees with IQ-TREE from each of the MSAs in a directory
# Written by Carlos Valiente Mullor, 2019
#
# Usage:
#       ./ml_trees.sh <path/to/alignments.fas>
#
# MSAs files should be fasta

msadir=$1

ls $msadir*.fas \
 | sed 's|$msadir||g' \
 | while read msa
 do
        nohup iqtree -s $msa -m GTR -bb 1000 -nt AUTO &
done
