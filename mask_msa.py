# Remove (mask) MSA columns with gaps in >=1 sequence(s) defined by the user
#
# Written by Carlos Valiente-Mullor, 2019
#
# Usage:
# Remove columns with gaps in ANY sequence:
#       python3 mask_msa.py <FASTA file> all
#
# Remove columns with gaps in specific sequences (indexes start at 0):
#       python3 mask_msa.py <FASTA file> [index1] [index2] ...
#
# Example:
# Remove columns if gaps are present in 1st (index 0), 4th and/or 8th sequence:
#       python3 mask_msa.py example.fas 0 3 7


from Bio import AlignIO
import re
import sys


####  Functions  ####

# Return list of index pairs (start, end) of regions w/o gaps
# Input: DNA sequence (string)
def nogaps_range(s):
	for i in re.finditer('(?!-+)\w+', s):
		yield (i.start(), i.end())

# Concatenate columns w/o gaps (in the reference seq.) from a MSA
def remove_gaps(ref, msa_original):
	msa_edited = msa_original[:, 0:0] # Empty MSA (no sequences, only headers)
	for region in nogaps_range(ref):
		start, end = region[0], region[1]
		msa_edited = msa_edited + msa_original[:, start:end]
	return msa_edited

# Mask MSA
def mask(order_list, alignment):
	for order in order_list:
		reference = str(alignment[order].seq)
		alignment = remove_gaps(reference, alignment)
	return alignment

def main():
	filename = sys.argv[1]
	orders = sys.argv[2:]
	
	msa = AlignIO.read("complete_msa/" + filename, "fasta")

	if orders[0] == 'all':
		outfile = "complete_msa/masked_core." + filename
		orders = range(len(msa))
	else:
		outfile = "complete_msa/masked." + filename
		orders = list(map(int, orders))
	
	masked = mask(orders, msa)
	AlignIO.write(masked, outfile, "fasta")

if __name__ == "__main__":
    main()
