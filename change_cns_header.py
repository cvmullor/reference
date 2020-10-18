# Change header of fasta sequence
#
# Written by Carlos Valiente-Mullor, 2019
#
# Usage:
#       python3 change_cns_header.py <FASTA file>


import sys


def new_header(fasfile):
	with open(fasfile, "r") as myf:
		lall = myf.readlines()
		seq = "".join(map(str.strip, lall[1:]))

	head = ">" + ".".join(fasfile.split(".")[:-2])
	outf = ".".join(fasfile.split(".")[:-1]) + ".fasta"
	with open(outf, "a") as newf:
		newf.write(head + "\n")
		newf.write(seq + "\n")

def main():
	fas = sys.argv[1]
	new_header(fas)

if __name__ == "__main__":
	main()
