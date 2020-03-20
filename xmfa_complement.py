# Obtain an XMFA alignment complemented by consensus sequences (fasta) associated
#	with different reference sequences included in an XMFA alignment obtained
# 	with progressiveMauve
#
# Written by Iván Ansari Toledano, 2019
#
# Usage:
#	python xmfa_complement.py <XMFA file>
#
# Paths in "def main()" must be modified depending on the reference sequence


import glob, os, re


def xmfa_memory(xmfa):
	fich = open(xmfa, "r")
	cont = 0
	orig_xmfa = []

	block_index = []
	block = []
	seq = ""
	seqid = ""
	ref_num = 0
	reverse = False
	# grab refs in xmfa
	for line in fich:

		if line[0] == "#":
			if "Sequence" in line and "File" in line:
				m = re.search("\\t(.*)\\n", line)
				if m:
					line = m.group(1)
				orig_xmfa.append([line])

		elif line[0] == ">":
			if len(seq) > 0:
				# reverse sequence if necessary
				if reverse:
					seq = seq[::-1]
					reverse = False
				orig_xmfa[ref_num-1].append([seqid, seq])
			m = re.search("> (\d):", line)
			if m:
				ref_num = int(m.group(1))
			n = re.search("(> .*)\s\w*", line)
			if n:
				coords = n.group(1)
				block.append(coords)
			if coords[-1] == "-":
				reverse = True
			seqid = line
			seq = ""

		elif line[0] == "=":
			block_index.append(block)
			block = []

		else:
			line = line.rstrip()
			seq = seq + line
	fich.close()
	orig_xmfa[ref_num-1].append([seqid, seq])

	return orig_xmfa, block_index


def list_fastas():
	"""
	fasta_list = glob.glob("*.sorted.bam")
	for i in range(len(fasta_list)):
		fasta_list[i] = fasta_list[i].replace(".sorted.bam", ".fasta")
	return fasta_list
	"""

	fasta_list = glob.glob("*.fasta")
	return fasta_list


def fasta_dict(fasta, orig_xmfa, index):
	# se crea un diccionario coord vs secuencia
	nucdic = {"A":"T", "T":"A", "a":"t", "t":"a", "G":"C", "C":"G", "g":"c", "c":"g"}
	dic = {}
	seq = ""
	fich = open(fasta, "r")

	for line in fich:
		if line[0] != ">":
			seq = seq + line.rstrip()

	for j in range(len(orig_xmfa[index])):
		if j == 0:
			continue

		m = re.search(":(\d+)", orig_xmfa[index][j][0])
		if m:
			coord_start = int(m.group(1))
		n = re.search("-(\d+)", orig_xmfa[index][j][0])
		if n:
			coord_end = int(n.group(1))

		dicseq = seq[coord_start-1:coord_end]

		# check if it is needed to reverse the sequence
		reverse = re.search(".+\s(-)\s.+", orig_xmfa[index][j][0])

		if reverse:
			dicseq = dicseq[::-1]
			newseq = ""
			for nucl in dicseq:
				if nucl in nucdic:
					newseq = newseq + nucdic[nucl]
				else:
					newseq = newseq + nucl
			dicseq = newseq

		dic[orig_xmfa[index][j][0]] = dicseq

	fich.close()
	return dic


def fasta_mem(ref_list, orig_xmfa):
	fasmem = []
	for index in range(len(orig_xmfa)):
		reffastas = []
		os.chdir(ref_list[index])
		fasta_list = list_fastas()
		for fasta in fasta_list:
			fasdic = fasta_dict(fasta, orig_xmfa, index)
			reffastas.append([fasta.replace(".fasta", ""), fasdic])

		fasmem.append(reffastas)

	return fasmem

def id_numbers_header(orig_xmfa, fasmem, ref_list):
	id_nums = {}
	header = ["#FormatVersion Mauve1\n"]
	cont = 1
	for i in range(len(orig_xmfa)):
		id_nums[orig_xmfa[i][0]] = str(cont)
		idvar = re.sub(r"(\/.+\/)", "", orig_xmfa[i][0])
		header.append("#Sequence"+str(cont)+"File\t"+idvar+"\n")
		header.append("#Sequence"+str(cont)+"Format\tFastA\n")
		cont = cont + 1
		for j in range(len(fasmem[i])):
			id_nums[fasmem[i][j][0]] = str(cont)
			idvar = ref_list[i]+fasmem[i][j][0]+".fasta\n"
			idvar = re.sub(r"(\/.+\/)", "", idvar)
			header.append("#Sequence"+str(cont)+"File\t"+idvar)
			header.append("#Sequence"+str(cont)+"Format\tFastA\n")
			cont = cont + 1
	header.append("#BackboneFile\t101ref.xmfa.bbcols\n")
	return id_nums, header

def folded_seq(seq):
	n = 80
	seq_list = [seq[i:i+n] for i in range(0, len(seq), n)]
	seq = ""
	for i in seq_list:
		seq = seq + i + "\n"
	return seq


def create_gap_index(refseq, fastaseq):
	gapindex = []
	index = 0
	start = 0
	length = 0
	started = False
	for position in range(len(refseq)):
		if not started and refseq[position] == "-":
			start = index
			length = 1
			started = True
		elif started and refseq[position] == "-":
			length = length + 1
		elif started and refseq[position] != "-":
			gapindex.append([index, length])
			index = index + 1
			length = 0
			started = False
		elif not started and refseq[position] != "-":
			index = index + 1
	if started:
		gapindex.append([index, length])
	gappedseq = insert_gaps(gapindex, fastaseq)
	if len(gappedseq) < len(refseq):
		complete_len = len(refseq) - len(gappedseq)
		gappedseq = gappedseq + "-" * complete_len
	return gappedseq


def insert_gaps(gapindex, seq):
	shift = 0
	for gaps in gapindex:
		seq = seq[:(gaps[0] + shift)] + "-" * gaps[1] + seq[(gaps[0] + shift):]
		shift = shift + gaps[1]
	return seq


def write_xmfa(orig_xmfa, fasmem, id_nums, header, block_index, ref_list):
	fich = open("xmfa_final.xmfa", "w")

	for i in header:
		fich.write(i)

	for block in block_index:
		for ref_coords in block:
			m = re.search(">\s(\d):", ref_coords)
			if m:
				ref = int(m.group(1)) - 1
			for i in range(len(orig_xmfa[ref])):
				if i == 0:
					continue
				if ref_coords in orig_xmfa[ref][i][0]:
					idvar = re.sub(r">\s(\d):", "> "+id_nums[orig_xmfa[ref][0]]+":",orig_xmfa[ref][i][0])
					idvar = re.sub(r"(\/.+\/)", "", idvar)
					fich.write(idvar)
					seq = folded_seq(orig_xmfa[ref][i][1])
					fich.write(seq)
					key_dic = orig_xmfa[ref][i][0]
					refseq = orig_xmfa[ref][i][1]
			for i in range(len(fasmem[ref])):
				ref_coordsvar = re.sub(r">\s(\d):", "> "+id_nums[fasmem[ref][i][0]]+":", ref_coords)
				idvar = ref_coordsvar+" "+ref_list[ref]+fasmem[ref][i][0]+".fasta\n"
				idvar = re.sub(r"(\/.+\/)", "", idvar)
				fich.write(idvar)
				fich.write(folded_seq(create_gap_index(refseq, fasmem[ref][i][1][key_dic])))
		fich.write("=\n")

	fich.close()

def main():

	xmfa_directory = "path/to/reference.xmfa"

	# Nº of paths and variables (refX) depend on nº of reference genomes in the XMFA alignment
	# Fasta sequences in "consensus_sequences" directory will be added to the XMFA according to the coordinates on the corresponding reference
	# To include mappings to only one of the references, the remaining paths must lead to empty folders
	ref1 = "path/to/consensus_sequences" # directory of consensus sequences (fasta) obtained from mappings to reference 1
	ref2 = "path/to/empty"
	ref3 = "path/to/empty"

	# List of references may be adapted to the number of reference genomes in the xmfa alignment of references
	# e.g., [ref1, ref2, ref3, ref4, ref5], if there are 5 references sequences
	ref_list = [ref1, ref2, ref3]

	orig_xmfa, block_index = xmfa_memory(xmfa_directory)

	curwd = os.getcwd()

	fasmem = fasta_mem(ref_list, orig_xmfa)
	id_nums, header = id_numbers_header(orig_xmfa, fasmem, ref_list)

	os.chdir(curwd)
	write_xmfa(orig_xmfa, fasmem, id_nums, header, block_index, ref_list)

if __name__ == "__main__":
	main()
