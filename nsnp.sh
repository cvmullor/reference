#! /bin/bash

# Loop for obtaining nº of SNPs from each mapping
# Written by Carlos Valiente Mullor, 2019
#
# Usage:
#		./s1_nsnp.sh
#
# Input: compressed filtered variants obtained with 'c1_consensus.sh'

ls *.flt.vcf.gz \
 | sed 's/.flt.vcf.gz//g' \
 | while read pass
 do
	bcftools stats $pass.flt.vcf.gz | grep '^SN' | cut -f3- >> allstats.txt
done

ls *.flt.vcf.gz | sed 's/.flt.vcf.gz//g' > fltnames.txt
grep "SNPs:" allstats.txt | cut -f2- > nsnps.txt
paste fltnames.txt nsnps.txt > pass_snps.csv
	
rm allstats.txt fltnames.txt nsnps.txt
