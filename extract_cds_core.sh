#! /bin/bash

# Extract and concatenate core CDS from consensus sequences obtained from
# 	mappings against the same reference genome
# Written by Carlos Valiente-Mullor, 2019
#
# Usage:
#		./cds_extr-concat.sh <strain name>
#
# Required:
#		cds_revcomp.py
#
# CDS coordinates from reference sequences were obtained with Prokka
# Annotation files were named with strain name and gff extension (e.g., FA_1090.gff)
# Consensus sequences obtained from mappings against the same reference were placed
#  in a directory named with the corresponding strain name
# Strict core CDS of the references from the same species were previously obtained
# 	with Proteinortho (file 'core.proteinortho.csv')


strain=$1 # strain name (e.g., FA_1090)

# 'core.proteinortho.csv': proteinortho output file including only strict core of the references
awk -F'\t' -v r=$strain 'NR==1{for (i=1; i<=NF; i++) if ($i==r){c=i; break}; next} {print $c}' core.proteinortho.csv > $strain.core.csv

# CDS coordinates of the references (annotation files in a directory named 'gff')
pos=$(grep -f $strain.core.csv gff/$strain.gff | grep -v "Aragorn" | awk '{print $4"-"$5","}' ORS='')

# Coordinates + strand (+ or -) to generate headers for each CDS
grep -f $strain.core.csv gff/$strain.gff | grep -v "Aragorn" | awk '{print ">"$4"-"$5"_"$7}' > newhead.$strain

# Folders with consensus sequences in a directory named 'cns' (e.g., /cns/FA_1090/<sample_name>.cns.fasta)
ls cns/$strain/ \
 | sed "s/cns\/$strain\///g" \
 | sed 's/.cns.fasta//g' \
 | while read sample;
 do 
	# Extract CDSs from each consensus sequences with EMBOSS (each CDS in a different entry)
	~/Documentos/softw/EMBOSS-6.6.0/emboss/extractseq --sequence cns/$strain/$sample*.cns.fasta -reg $pos -separate -outseq cds.$sample.$strain.fas
	
	# Multi-fasta of CDSs with sequences in one line
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  
		END {printf("\n");}' < cds.$sample.$strain.fas > cds.$sample.$strain.temp
	# Change headers for coordinates + strand
	sed '1{/^$/d}' cds.$sample.$strain.temp \
	 | awk 'NR%2==0' \
	 | paste -d'\n' newhead.$strain - > def.cds.$sample.$strain.fas

	rm cds.$sample.$strain.temp
	
	# Reverse complement of CDS with strand '-' (use custom python script cds_revcomp.py)
	python3 cds_revcomp.py def.cds.$sample.$strain.fas
	
	# Concatenate all CDS of each consensus sequences in a single sequence
	grep -v "^>" revcomp_def.cds.$sample.$strain.fas | awk -v head=$sample"."$strain 'BEGIN { ORS=""; print ">"head".CDS\n" } { print }' > finalcds_core.$sample.$strain.fas
	echo -e "\n" >> finalcds_core.$sample.$strain.fas
	

	rm cds.$sample.$strain.fas def.cds.$sample.$strain.fas revcomp_def.cds.$sample.$strain.fas
done

rm newhead.$strain

# MSA of all the consensus sequences from mappings against the same reference
cat finalcds_core.* > finalcds_core_aln.$strain.fas

mkdir cds.$strain
mv finalcds* cds.$strain
