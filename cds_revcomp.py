# Reverse complement of CDSs with strand '-'
#
# Written by Carlos Valiente-Mullor, 2019
#
# Usage:
# 		python3 cds_revcomp.py <CDS multi-fasta>
#
# Input: multi-fasta file with CDSs; last character of header string is
#	'+' or '-' (strand), according to annotation files used for extracting CDSs.
# Output: multi-fasta file with reverse complement of minus strand sequences,
#	the remaining CDS are not modified.


import sys
import Bio
from Bio.Seq import Seq


cds = sys.argv[1] # FASTA file
outfile = "revcomp_" + cds

with open(cds, "r") as myf:
	lines = myf.readlines()
	lall = list(map(str.strip, lines))

# Reverse complement if needed
final_mfa = []
rev = False
for line in lall:
	if line[0] == ">":
		if line[-1] == "+":
			rev = False
		elif line[-1] == "-":
			rev = True
		final_mfa.append(line)
	else:
		if rev == False:
			final_mfa.append(line)
		elif rev == True:
			seq = Seq(line)
			rc = str(seq.reverse_complement())
			final_mfa.append(rc)

with open (outfile, "a") as outf:
	for line in final_mfa:
		outf.write(line + "\n")
