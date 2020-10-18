# Add gaps to regions where homologous sequences were absent in any genome in
#  a XMFA-formatted alignment
#
# Written by Carlos Valiente-Mullor, 2019
#
# Usage:
# 		python3 gaps_xmfa.py <XMFA file>


import sys
import re
import glob
import os


# read XMFA
def read_xmfa(xmfa):
	with open("complete_msa/" + xmfa, "r") as myf:
		lines = myf.readlines()
		lsequences = list(map(str.strip, lines))
	return lsequences

# entry:filename
def sequence_index(lsequences):
	seqname_index = {}
	for entry in lsequences:
		heads = re.search("#Sequence..*File\t[a-zA-Z0-9.]*", entry)
		stop = re.search(">\s\d", entry)
		if heads:
			temp = heads.group(0).split("\t")
			name = temp[1]
			entry = re.search("\d\d*", temp[0])
			seqname_index[entry.group(0)] = name
		elif stop:
			break
	nseq = len(seqname_index) # number of sequences in aln
	return nseq, seqname_index

# read XMFA, fulfill files with separated blocks
def divide_blocks(lsequences, outfile, seqname_index):
	block = 0
	for line in lsequences:
		if "#" in line:
			with open(outfile, "a") as xmfa:
				xmfa.write(line + "\n")
		elif line == "=":
			block += 1
		elif ">" in line:
			blockfile = str(block) + ".fa"
			key = line[2]
			seqname = seqname_index[key]
			with open(blockfile, "a") as bkf:
				bkf.write("\n" + line + "\n")
		else:
			with open(blockfile, "a") as bkf:
				bkf.write(line)

# read .fa files, add gap entries, concat in a new xmfa file
def add_gaps(nseq, seqname_index, outfile):
	fafiles = glob.glob("*.fa") # read files
	fafiles.sort(key=lambda f: int(''.join(filter(str.isdigit, f)))) # sort .fa files (ascending order)
	for fa in fafiles:
		with open(fa, "r") as myf:
			lines = myf.readlines()
			lall2 = list(map(str.strip, lines[1:]))
		blength = len(lall2[1])
		conserved = []
		for line in lall2:
			if ">" in line:
				if line[3] == ":":
					conserved.append(line[2])
				else:
					conserved.append(line[2:4])
		if len(conserved) < nseq:
			for key in seqname_index.keys():
				if key not in conserved:
					header = "> " + key + ":0-0 + " + seqname_index[key]
					gaps = "-" * blength
					with open(fa, "a") as new:
						new.write("\n" + header + "\n" + gaps)
		with open(fa, "a") as new:
			new.write("\n=\n")
		with open(fa, "r") as myf:
			lines = myf.readlines()
			lall2 = list(map(str.strip, lines[1:]))
		with open(outfile, "a") as final:
			for line in lall2:
				final.write(line + "\n")
	# Remove block files
	for fa in glob.glob("*.fa"):
		os.remove(fa)

def main():
	fxmfa = sys.argv[1]
	final = "wgaps." + fxmfa

	lall = read_xmfa(fxmfa)
	seqnumber, seqindex = sequence_index(lall)
	divide_blocks(lall, final, seqindex)
	add_gaps(seqnumber, seqindex, final)

if __name__ == "__main__":
	main()