#! /bin/bash

# Pipeline used in the work "One is not enough: on the effects of reference
#	genome for the mapping and subsequent analyses of short-reads". BioRxiv
#	[Preprint]. 2020. Available from: https://doi.org/10.1101/2020.04.14.041004.
#	
#	Written by Carlos Valiente Mullor, 2020
#
# Usage:
#		./refeval_main.sh [-hd] -r <path/to/references> -s <path/to/reads> [-t <threads>] [-q <adapters.fasta>]


## Functions ##

function print_usage {
    echo "Usage: ./refeval_main.sh [-hd] -r <path/to/references> -s <path/to/reads> [-t threads] [-q adapters_file]"
    echo "Options:"
    echo " -r       path to the folder that contains the reference sequences (FASTA) for mapping"
    echo " -s       path to the folder that contains the reads (FASTQ) to be mapped"
    echo " -t       number of threads to be used during mapping"
    echo -e " -q       path to file with adapters to be removed during read quality\n\t  processing. If no file is provided, QC step will not be performed"
    echo " -h       this help message"
    echo " -d       required dependencies and files"
    echo ""
}

function print_required {
    echo "Required software:"
    echo -e "\tbwa"
    echo -e "\tsamtools"
    echo -e "\tbcftools"
    echo -e "\tprogressiveMauve"
    echo -e "\tIQ-tree"
    echo -e "\tTreeCmp"
    echo -e "\tProteinortho"
    echo -e "\textractseq (EMBOSS package)"
    echo -e "\tcodeml (PAML package)"
    echo -e "\tLDJump"
    echo -e "\tLDhat"
    echo -e "\tPhi (PhiPack)"
    echo ""
    echo "Required scripts:"
    echo -e "\tchange_cns_head.py"
    echo -e "\txmfa_complement.py"
    echo -e "\tgaps_xmfa.py"
    echo -e "\tmask_msa.py"
    echo -e "\tcds_revcomp.py"
    echo -e "\trho_LDJump.R"
    echo -e "\tplots.R"
    echo -e "\tstats.R"
    echo -e "\txmfa2fasta.pl"
    echo -e "\n  Python and R scripts are available at https://github.com/cvmullor/reference"
    echo "  'xmfa2fasta.pl' is available at https://github.com/kjolley/seq_scripts"
    echo ""
    echo "Required python libraries:"
    echo -e "\tBiopython"
    echo ""
    echo "Required R libraries:"
    echo -e "\tggplot2"
    echo -e "\tdplyr"
    echo ""
    echo "Required files:"
    echo -e "\tProkka annotation (FFN) of reference sequences"
    echo -e "\tProkka annotation (GFF) of reference sequences"
    echo ""
}

function check_path {
    local path=$1
    if [[ ${path: -1} == "/" ]];
    then
        path=$(echo $path | sed 's/\(.*\)\//\1 /')
    fi
    echo $path
}

function rm_adaptors {
    local readp=$1
    local adapfile=$2

    # Remove adaptors
    ls $readp/*.fastq \
     | sed -e "s|$readp\/||g" \
     | awk '{ if ($1 ~ /R2/) {print $0"\n"} else {print $0" "}}' ORS="" \
     | while read fastqR1 fastqR2;
     do
        cutadapt -a file:$adapfile -A file:$adapfile -o $readp/$fastqR1.1.cut.fq -p $readp/$fastqR2.2.cut.fq $readp/$fastqR1 $readp/$fastqR2
    done
    mkdir $readp/original_fastq
    mv $readp/*.fastq $readp/original_fastq/
}

function quality_control {
    local readp=$1
    local ppath=$2

    # Quality filtering and trimming
    ls $readp/*.fq \
     | sed -e "s|$readp\/||g" \
     | awk '{ if ($1 ~ /R2/) {print $0"\n"} else {print $0" "}}' ORS="" \
     | while read fastqR1 fastqR2;
     do
        perl $ppath/prinseq-lite.pl -fastq $readp/$fastqR1 -fastq2 $readp/$fastqR2 -min_len 50 -min_qual_mean 20 -ns_max_p 10 -trim_qual_right 20 2>>prinseq_stats.out
    done
    mkdir $readp/cut_fastq $readp/bad_QC
    mv $readp/*.cut.fq $readp/cut_fastq/
    mv $readp/*_singletons_* $readp/bad_QC/
    mv $readp/*_bad_* $readp/bad_QC/

    ls -1 $readp/*.fastq > $readp/temp1.txt
    ls -1 $readp/*.fastq | sed 's/_prinseq_good_....//g' > $readp/temp2.txt
    paste $readp/temp* | while read line; do mv $line; done
    rm $readp/temp1.txt $readp/temp2.txt
}

function index_refs {
    local refp=$1
    local ext=$2

    mkdir indexes
    cp $refp/*.$ext indexes

    ls indexes/*.$ext | while read file;
    do
	    bwa index -a bwtsw $file
	    samtools faidx $file
    done
}

function mapsort {
    local readp=$1
    local ext=$2
    local thr=$3

    # Map and sort
    ls $readp/*.fastq \
     | sed -e "s|$readp\/||g" \
     | sed 's/.fastq//g' \
     | awk '{ if ($1 ~ /R2/) {print $0"\n"} else {print $0" "}}' ORS="" \
     | while read fastqR1 fastqR2;
    do
	    ls indexes/*.$ext \
	     | sed -e 's|indexes\/||g' \
	     | sed "s/.$ext//g" \
	     | while read reference;
	    do
		    bwa mem -t $thr indexes/$reference.$ext $readp/$fastqR1.fastq $readp/$fastqR2.fastq 2>bwa.out \
		     | samtools sort -@$thr -o $fastqR1.$reference.sorted.bam - 2>sort.out
	    done
    done;

    # Index sorted BAM files
    ls *.sorted.bam \
     | while read srtfile;
    do 
    	samtools index $srtfile
    done

    mkdir bam_files
    mv *.sorted.bam* bam_files
}

function sample_list {
    local readp=$1

    ls $readp/*.fastq \
     | sed -e "s|$readp\/||g" \
     | sed 's/.fastq//g' \
     | awk '{ if ($1 ~ /R1/) {print} }' > sample_list.temp
}

function varcalling {
    local ext=$1

    ls indexes/*.$ext \
     | sed -e 's/indexes\///g' \
     | sed "s/\.$ext//g" \
     | while read reference;
     do
    	cat sample_list.temp \
    	 | while read sample;
	     do
		    samtools mpileup -u -f indexes/$reference.$ext bam_files/$sample.$reference.sorted.bam \
		     | bcftools call -vmO v --skip-variants indels -o $sample.$reference.raw.vcf
	    done
    done
}

function varfilter {
    ls *.raw.vcf \
     | sed 's/.raw.vcf//g' \
     | while read raw
     do
    	bcftools query -f '%DP\n' $raw.raw.vcf \
	     | awk -v filename=$raw '{ sum += $0; n++ } END { if (n > 0) print (sum / n) * 2" "filename; }' >> avgx2.txt
    done;

    cat avgx2.txt \
     | while read average name
     do
	    bcftools filter -g 10 -e '%QUAL<40 || DP<8 || GT!="1/1" || MQ<30 || DP>'$average $name.raw.vcf -Oz -o $name.flt.vcf.gz
	    bcftools index $name.flt.vcf.gz
    done;

rm avgx2.txt
}

function consensus {
    local ext=$1

    ls indexes/*.$ext \
     | sed -e 's|indexes\/||g' \
     | sed "s/\.$ext//g" \
     | while read reference;
     do
    	cat sample_list.temp \
    	 | while read sample;
	     do
		    cat indexes/$reference.$ext \
		     | bcftools consensus $sample.$reference.flt.vcf.gz > $sample.$reference.cns.fa
	    done
    done

    mkdir vcf_raw vcf_filtered consensus_seqs
    mv *.raw.vcf vcf_raw
    mv *.flt.vcf.gz* vcf_filtered
    mv *.cns.fa consensus_seqs
}

function mapstats {
    # Mean read coverage (depth)
    ls -1 bam_files/*.sorted.bam \
    | sed 's|bam_files\/||g' \
    | sed 's/.sorted.bam//g' \
    | while read bamfile;
    do
        # Only considering positions with al least 1 read mapped
        samtools depth -aa bam_files/$bamfile.sorted.bam | awk '{ if ($3>0) {print $3} }' > $bamfile.nozero.temp
        awk '{ sum += $0; n++ } END { if (n > 0) print (sum / n); }' $bamfile.nozero.temp >> mcov_nozero.csv
        
        # Considering all positions
        samtools depth -aa bam_files/$bamfile.sorted.bam | awk '{ if ($3>=0) {print $3} }' > $bamfile.all.temp
        awk '{ sum += $0; n++ } END { if (n > 0) print (sum / n); }' $bamfile.all.temp >> mcov_all.csv
    done

    rm *.temp

    # % of reference length covered by reads
    ls -1 bam_files/*.sorted.bam \
    | sed 's|bam_files\/||g' \
    | sed 's/.sorted.bam//g' \
    | while read bamfile;
    do
        # number of positions with >0 reads mapped
        nonzero=$(samtools depth -aa bam_files/$bamfile.sorted.bam | awk '{if($3>0) total+=1}END{print total}')
        # All positions (= length of the reference)
        all=$(samtools depth -aa bam_files/$bamfile.sorted.bam | awk '{if($3>=0) total+=1}END{print total}')
        percent=$(bc <<< "scale=6; ($nonzero / $all)*100")
        echo -e "$percent\t$all" >> perc_ref_covered.tsv
    done

    # Number of SNPs
    ls vcf_filtered/*.flt.vcf.gz \
    | sed 's|vcf_filtered\/||g' \
    | sed 's/.flt.vcf.gz//g' \
    | while read pass
    do
        bcftools stats vcf_filtered/$pass.flt.vcf.gz | grep '^SN' | cut -f3- >> vcf_allstats.txt
        refe=$(echo ${pass##*.})
        echo -e $pass"\t"$refe >> fltnames.txt
    done

    grep "SNPs:" vcf_allstats.txt | cut -f2- > nsnps.txt
    paste fltnames.txt nsnps.txt > pass_snps.tsv
        
    rm vcf_allstats.txt fltnames.txt nsnps.txt

    # Merge results of the three parameters
    echo -e "sample\treference\tSNPs\tper.ref.covered\treference_length\tmean.coverage\tmean.coverage.all_pos" > mapstats.tsv
    paste pass_snps.tsv perc_ref_covered.tsv mcov_nozero.csv mcov_all.csv >> mapstats.tsv

    rm pass_snps.tsv perc_ref_covered.tsv mcov_nozero.csv mcov_all.csv
    mkdir mapstats
    mv mapstats.tsv mapstats/
}

function group_cns {
    local ext=$1

    # Change consensus sequence headers (sample and reference names)
    for file in consensus_seqs/*.cns.fa;
    do
        python change_cns_header.py $file
    done
    rm consensus_seqs/*.cns.fa

    # Group consensus from mappings against the same reference in folders
    ls indexes/*.$ext \
     | sed -e 's|indexes\/||g' \
     | sed "s/\.$ext//g" \
     | while read reference;
    do
        mkdir consensus_seqs/$reference
        mv consensus_seqs/*$reference*.fasta consensus_seqs/$reference/
    done
    mkdir consensus_seqs/empty
}

function alignment {
    local ext=$1
    local ppath=$2

    # Alignment of references
    allref=$(ls indexes/*.$ext | awk '{ print $1" " }' ORS='')

    $ppath/progressiveMauve --output=references.aln.xmfa $allref

    mkdir reference_MSA
    mv references.aln.xmfa* reference_MSA

    # Add consensus sequences to MSA of references
    total_refs=$(ls indexes/*.$ext | wc -l)
    ref_idx=0
    full_path=$(pwd)

    ls indexes/*.$ext \
     | sed -e 's|indexes\/||g' \
     | sed "s/\.$ext//g" \
     | while read reference;
    do
        ((ref_idx++))
        python xmfa_complement.py $total_refs $ref_idx $reference $full_path
        mv xmfa_final.xmfa xmfa_final.$reference.xmfa
    done

    mkdir complete_msa
    mv xmfa_final.*.xmfa complete_msa/

    # Add gaps to XMFA alignments
    ls complete_msa \
     | while read file;
    do
        python gaps_xmfa.py $file
    done

    mv wgaps.* complete_msa/

    # Convert XMFA to FASTA
    ls complete_msa/wgaps.* \
     | sed 's|\.xmfa$||g' \
     | while read file;
    do
        perl xmfa2fasta.pl --file $file.xmfa > $file.fasta
    done

    # Write informative MSA headers
    ls indexes/*.$ext \
     | sed -e 's|indexes\/||g' \
     | sed "s/\.$ext//g" \
     | while read reference;
    do
        grep "#Sequence.File" complete_msa/xmfa_final.$reference.xmfa \
         | sed 's|indexes\/||g' \
         | awk '{ print ">"$2 }' > complete_msa/headers.$reference.temp
        
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  
        END {printf("\n");}' < complete_msa/wgaps.xmfa_final.$reference.fasta > complete_msa/aln.$reference.temp

        sed '1{/^$/d}' complete_msa/aln.$reference.temp \
         | awk 'NR%2==0' \
         | paste -d'\n' complete_msa/headers.$reference.temp - > complete_msa/def.wgaps.xmfa_final.$reference.fasta
    done

    rm complete_msa/headers.*.temp complete_msa/aln.*.temp
}

function mask {
    ref_idx=0

    ls complete_msa/def.*.fasta \
     | sed -e 's|complete_msa\/||g' \
     | while read file;
    do
        python mask_msa.py $file all
        python mask_msa.py $file $ref_idx
        ((ref_idx++))
    done
}

function ml_tree {
    # Build ML phylogenies from each MSA
    ls complete_msa/masked*.fasta \
     | sed 's|complete_msa\/||g' \
     | while read msa
     do
        iqtree -s complete_msa/$msa -m GTR -bb 1000 -nt AUTO
    done
    mkdir trees
    mv complete_msa/masked*.fasta.* trees/
}

function topology_test {
    local ext=$1

    # Remove distances from trees (only topology remains)
    cat trees/masked.*.treefile | sed "s/:[0-9\.-]*//g" > masked.trees.nwk
    cat trees/masked_core.*.treefile | sed "s/:[0-9\.-]*//g" > core.trees.nwk

    # Remove reference name
    ls indexes/*.$ext \
     | sed -e 's|indexes\/||g' \
     | sed "s/\.$ext//g" \
     | while read reference;
     do
        sed -i "s|\.$reference\.cns\.fasta||g" masked.trees.nwk
        sed -i "s|\.$reference\.cns\.fasta||g" core.trees.nwk
        sed "s|\.$reference\.cns\.fasta||g" complete_msa/masked.def.wgaps.xmfa_final.$reference.fasta > masked.$reference.aln.fasta
        sed "s|\.$reference\.cns\.fasta||g" complete_msa/masked_core.def.wgaps.xmfa_final.$reference.fasta > core.$reference.aln.fasta
    done

    # Perform topology tests
    ls masked.*.aln.fasta \
     | sed 's/masked\.//g' \
     | while read aln;
     do
        iqtree -nt AUTO -s masked.$aln -z masked.trees.nwk -m GTR -zb 1000 -zw
        iqtree -nt AUTO -s core.$aln -z core.trees.nwk -m GTR -zb 1000 -zw
    done

    mkdir topology_test
    mv *.aln.fasta.* topology_test/
    mv masked.* topology_test/
    mv core.* topology_test/
}

function tree_dist {
    local ppath=$1

    # Compute tree distances
    java -jar $ppath/treeCmp.jar -m -d mc rc -i topology_test/masked.trees.nwk -o masked.trees.summ.tsv -N -I
    java -jar $ppath/treeCmp.jar -m -d mc rc -i topology_test/core.trees.nwk -o core.trees.summ.tsv -N -I
    head -n -5 masked.trees.summ.tsv > masked.trees.dist.tsv # w/o summary
    head -n -5 core.trees.summ.tsv > core.trees.dist.tsv     # w/o summary

    mkdir tree_distance; mv *.tsv tree_distance/
}

function strict_core {
    local ext=$1

    # Run proteinortho
    # Use annotated fasta files (.ffn) of references from Prokka
    proteinortho5.pl -project=refs_orthologs -p=blastn+ -singles -nograph -clean ffn/*.ffn

    ls indexes/*.fna \
     | sed -e 's|indexes\/||g' \
     | sed "s/\.fna//g" \
     | awk '{print}' ORS='\t' \
     | awk '{ print "species\tgenes\talgconn\t"$0}' > core.proteinortho.csv

    # Parse proteinortho output to retain core genes without duplicates
    total_refs=$(ls indexes/*.$ext | wc -l)
    awk -v tr=$total_refs '{if ($1==tr && $2==tr) {print} }' refs_orthologs.proteinortho >> core.proteinortho.csv
}

function extract_cds {
    local strain=$1
    local corext=$2
    local ppath=$3

    # Extract core CDS and concatenate from consensus sequences 
    #	+ reverse complement if necessary
    ls consensus_seqs/$strain/ \
     | sed "s/consensus_seqs\/$strain\///g" \
     | sed 's/.cns.fasta//g' \
     | while read sample;
     do 
        # Extract CDSs from each consensus sequences with EMBOSS (each CDS in a different entry)
        $ppath/extractseq --sequence consensus_seqs/$strain/$sample*.cns.fasta -reg $pos -separate -outseq cds.$sample.fas
        
        # Multi-fasta of CDSs with sequences in one line
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  
            END {printf("\n");}' < cds.$sample.fas > cds.$sample.temp
        # Change headers for coordinates + strand
        sed '1{/^$/d}' cds.$sample.temp \
         | awk 'NR%2==0' \
         | paste -d'\n' newhead.$strain - > def.cds.$sample.fas

        rm cds.$sample.temp
        
        # Reverse complement of CDS with strand '-' (use custom python script cds_revcomp.py)
        python3 cds_revcomp.py def.cds.$sample.fas
        
        # Concatenate all CDS of each consensus sequences in a single sequence
        grep -v "^>" revcomp_def.cds.$sample.fas | awk -v head=$sample 'BEGIN { ORS=""; print ">"head".CDS\n" } { print }' > finalcds$corext.$sample.fas
        echo -e "\n" >> finalcds$corext.$sample.fas

        rm cds.$sample.fas def.cds.$sample.fas revcomp_def.cds.$sample.fas
    done

    rm newhead.$strain

    # MSA of all the consensus sequences from mappings against the same reference
    cat finalcds$corext.* > finalcds$corext.aln.$strain.fas

    mkdir cds$corext.$strain
    mv finalcds* cds$corext.$strain
}

function cds_alignment {
    local ext=$1
    local ppath=$2

    ls indexes/*.$ext \
    | sed -e 's|indexes\/||g' \
    | sed "s/\.$ext//g" \
    | while read reference;
    do
        ## Reference CDSs
        # CDS coordinates of the references
        # Use annotation files (.gff) of references from Prokka
        pos=$(grep "CDS" gff/$reference.gff | awk '{print $4"-"$5","}' ORS='')
        # Coordinates + strand (+ or -) to generate headers for each CDS
        grep "CDS" gff/$reference.gff | awk '{print ">"$4"-"$5"_"$7}' > newhead.$reference
        # Extract CDS
        extract_cds $reference "" $ppath

        ## Core CDSs
        # Extract column with gene ID of a particular reference
        awk -F'\t' -v r=$reference 'NR==1{for (i=1; i<=NF; i++) if ($i==r){c=i; break}; next} {print $c}' core.proteinortho.csv > $reference.core.csv
        # CDS coordinates of the references
        # Use annotation files (.gff) of references from Prokka
        pos=$(grep -f $reference.core.csv gff/$reference.gff | grep -v "Aragorn" | awk '{print $4"-"$5","}' ORS='')
        # Coordinates + strand (+ or -) to generate headers for each CDS
        grep -f $reference.core.csv gff/$reference.gff | grep -v "Aragorn" | awk '{print ">"$4"-"$5"_"$7}' > newhead.$reference
        # Extract core CDS
        extract_cds $reference _core $ppath
    done

    mkdir proteinortho_results
    mv core.proteinortho.csv proteinortho_results/
    mv refs_orthologs.proteinortho proteinortho_results/
    mv refs_orthologs.removed_blast-graph proteinortho_results/

    mkdir cds_msa
    mkdir cds_msa/ref_cds cds_msa/core_cds
    mv cds.* cds_msa/ref_cds
    mv cds_core.* cds_msa/core_cds

    rm *.core.csv
}

function compute_omega {
    local ext=$1

    # Compute dN/dS
    mkdir codeml_results

    ls indexes/*.$ext \
     | sed -e 's|indexes\/||g' \
     | sed "s/\.$ext//g" \
     | while read reference;
     do
        # Compute dN/dS
sed -i "s|path_to_alnfile|cds_msa\/ref_cds\/cds.$reference\/finalcds.aln.$reference.fas|" codeml.ctl
        sed -i "s|codeml_output_file|codeml_ref.$reference.out|" codeml.ctl
        echo -e "\n" | codeml
        mv codeml_ref.*.out codeml_results/
        rm 2ML.* 2NG.* rst rst1 rub
        sed -i "s|cds_msa\/ref_cds\/cds.$reference\/finalcds.aln.$reference.fas|path_to_alnfile|" codeml.ctl
        sed -i "s|codeml_ref.$reference.out|codeml_output_file|" codeml.ctl

        sed -i "s|path_to_alnfile|cds_msa\/core_cds\/cds_core.$reference\/finalcds_core.aln.$reference.fas|" codeml.ctl
        sed -i "s|codeml_output_file|codeml_core.$reference.out|" codeml.ctl
        echo -e "\n" | codeml
        mv codeml_core.*.out codeml_results/
        rm 2ML.* 2NG.* rst rst1 rub
        sed -i "s|cds_msa\/core_cds\/cds_core.$reference\/finalcds_core.aln.$reference.fas|path_to_alnfile|" codeml.ctl
        sed -i "s|codeml_core.$reference.out|codeml_output_file|" codeml.ctl

        # Parse outputs
        grep "dN/dS" codeml_results/codeml_ref.$reference.out | awk -v ref=$reference '{ if (NR>2) {print $8"\t"ref} }' > codeml_results/dnds_ref.$reference.tsv
        grep "dN/dS" codeml_results/codeml_core.$reference.out | awk -v ref=$reference '{ if (NR>2) {print $8"\t"ref} }' > codeml_results/dnds_core.$reference.tsv
    done

    # Replace "99.0000" with NAs
    echo -e "omega\tmsa" > codeml_results/dnds.ref.tsv
    cat codeml_results/dnds_ref.*.tsv | sed 's|99\.0000|NA|g' >> codeml_results/dnds.ref.tsv
    echo -e "omega\tmsa" > codeml_results/dnds.core.tsv
    cat codeml_results/dnds_core.*.tsv | sed 's|99\.0000|NA|g' >> codeml_results/dnds.core.tsv
}

function compute_rho {
    local ext=$1

    # Compute rho
    full_path=$(pwd)
    echo -e "rate\tmsa\tsegment" > rho.tsv

    ls indexes/*.$ext \
    | sed -e 's|indexes\/||g' \
    | sed "s/\.$ext//g" \
    | while read reference;
    do
        Rscript --vanilla rho_LDJump.R $reference $full_path
        rm Phi.inf.list Phi.inf.sites Phi.log Phi.poly.unambig.sites
        tr -s ' ' '\n' < $reference.cte_estimates.txt \
        | awk -v rs=$reference '{ print $1"\t"rs"\t"NR }' >> rho.tsv
    done

    mkdir ldjump_results
    mv *.cte_estimates.txt ldjump_results/
    mv rho.tsv ldjump_results/
}



## Arguments ##

isref='false'
refp=""
issamp='false'
readp=""
isthr='false'
thr=1
qc='false'
adapters=""

pprinseq="/home/carlos/Documentos/tfm/fq_mapeos/sofware_prep/prinseq-lite-0.20.4"
pmauve="/home/carlos/Documentos/softw/mauve/mauve_snapshot_2015-02-13/linux-x64"
ptreecmp="/home/carlos/Documentos/softw/TreeCmp/bin"
pemboss="/home/carlos/Documentos/softw/EMBOSS-6.6.0/emboss"


while getopts 'hdr:s:t:q:' flag; do
  case "${flag}" in
    h)
        print_usage ;;
    d)
        print_required ;;
    r) 
		isref='true'
		refp="${OPTARG}" ;;
    s)
		issamp='true'
		readp="${OPTARG}" ;;
    t)
		isthr='true'
		thr="${OPTARG}" ;;
	q)
		qc='true'
		adapters="${OPTARG}" ;;
  esac
done


## Checks ##

# Path to references
if $isref; then
	refp=$(check_path $refp)
	aref=$( ls $refp | head -n1 )
	ext=$(echo ${aref##*.}) # File extension of the references
else
	echo "error: Path to references should be provided"
fi

# Path to reads
if $issamp; then
	readp=$(check_path $readp)
else
	echo "error: Path to reads should be provided"
fi

# Threads
if $isref && $issamp ; then
	if $isthr; then
		num='^[0-9]+$'
		if ! [[ $thr =~ $num ]] ; then
   			echo "error: -t argument must be integer"
		else
		echo "Mapping with $thr threads"
		fi
	else
		thr=1
		echo "Mapping with 1 thread"
	fi
fi

# QC
if $isref && $issamp ; then
	if $qc; then
		echo -e "Performing QC\n"
		rm_adaptors $readp $adapters
		quality_control $readp $pprinseq
	else
		echo "Skipping QC"
	fi
fi


## Analyses ##

if $isref && $issamp ; then
    # Map each sample to each reference, assumming paired-end reads (2 fastq files
    #   per sample) identified as R1 (forward) and R2 (reverse)
    # Output: indexes references, sorted BAM files and BAM indexes
    index_refs $refp $ext
    mapsort $readp $ext $thr
    
    # Variant calling and filtering. Consensus calling
    # Output: raw and filtered VCF files, consensus sequences in FASTA format
    sample_list $readp
    varcalling $ext
    varfilter
    consensus $ext
   
    # Compute mapping and variant calling statistics
    # Output: number of SNPs, % of reference covered and mean coverage of each mapping
    #	in a TSV file
    mapstats
    
    # Alignment of consensus and reference sequences
    # output: two MSAs (FASTA) per reference (and XMFA intermediate files):
    #	a) "core" MSA (i.e., all columns with gaps in any sequences are removed)
    #	b) masked MSA (i.e., only positions with gaps in the reference genome masked)   
    group_cns $ext
    alignment $ext $pmauve
    mask
    
    # ML tree building, topology tests and RF, MC tree distance metrics
    # Output: tree files (newick), .iqtree files containing topology test results,
    #	TSV files with tree distance comparisons results
    ml_tree
    topology_test $ext
    tree_dist $ptreecmp
    
    # Extract references' CDS from consensus sequences and compute dN/dS in pairwise
    #	comparisons between sequences obtained from mappings to the same reference
    # Output: codeml output files, TSV with dN/dS values per reference
    strict_core $ext
    cds_alignment $ext $pemboss
    compute_omega $ext
    
    # Compute rho from alignments in 1000bp windows
    # Output: TSV file with rho values per reference and MSA region
    compute_rho $ext
    
    # Boxplots of SNPs, % reference covered, coverage and dN/dS per reference;
    #	depict distribution rho along the aligments depending on the reference
    # Output: PNG graphs
    Rscript --vanilla plots.R
    mkdir plots
    mv *.png plots/
    
    # Summary statistics of all the parameters. Perform Kruskal-Wallis and Wilcoxon 
    #	tests (for mapping statsitics and dN/dS), and pairwise Kolmogorov-Smirnov
    #	tests (for rho distributions)
    Rscript --vanilla stats.R
    mkdir stat_test
    mv *.csv stat_test/
fi