# Evaluation of the effect of reference selection
Pipeline for the evaluation of the effect of reference choice on short-read mapping as performed in the work *One is not enough: on the effects of reference genome for the mapping and subsequent analyses of short-reads. BioRxiv [Preprint]. 2020*. Available from: [https://doi.org/10.1101/2020.04.14.041004](https://www.biorxiv.org/content/10.1101/2020.04.14.041004v1).

`refeval_main.sh` executes all the analyses carried out in this work.

The workflow is summarized in Fig 1:
![Fig 1](Fig1_overview.png)

## Options

### `-h` (help)
```
Usage: ./refeval_main.sh [-hd] -r <path/to/references> -s <path/to/reads> [-t threads] [-q adapters_file]
Options:
 -r       path to the folder that contains the reference sequences (FASTA) for mapping
 -s       path to the folder that contains the reads (FASTQ) to be mapped
 -t       number of threads to be used during mapping
 -q       path to file with adapters to be removed during read quality
          processing. If no file is provided, QC step will not be performed
 -h       this help message
 -d       required dependencies and files
```

### `-d` (required dependencies and files)
```
Required software:
        bwa
        samtools
        bcftools
        progressiveMauve
        IQ-tree
        TreeCmp
        Proteinortho
        extractseq (EMBOSS package)
        codeml (PAML package)
        LDJump
        LDhat
        Phi (PhiPack)

Required scripts:
        change_cns_header.py
        xmfa_complement.py
        gaps_xmfa.py
        mask_msa.py
        cds_revcomp.py
        rho_LDJump.R
        plots.R
        stats.R
        xmfa2fasta.pl

  Python and R scripts are available at https://github.com/cvmullor/reference
  'xmfa2fasta.pl' is available at https://github.com/kjolley/seq_scripts

Required python libraries:
        Biopython

Required R libraries:
        ggplot2
        dplyr

Required files:
        Prokka annotation (FFN) of reference sequences
        Prokka annotation (GFF) of reference sequences
```

All required programms should be in $PATH. The required scripts, as well as folders "ffn" and "gff" containing annotation files, should be placed in the same folder where `refeval_main.sh` is executed. 


### `-r` (references)
Path to the folder containing the reference sequences (in fasta/fna format) that will be used for read mapping.

### `-s` (samples)
Path to the folder containing the reads (in fastq format) that will be mapped against the different references.

### `-t` (threads) [optional]
Number of threads that BWA MEM will use for mapping step. Default: 1

### `-q` (quality control) [optional]
Path to the file containing the adapters to be removed during read quality processing. Quality trimming and filtering will also be performed.
If no file is provided, QC step will be skipped.

## Description of scripts

Usage is detailed in the header section of each script.

`refeval_main.sh`: main script. Performs the main steps decribed in the work:
  1) If specified, reads are filtered and trimmed.
  2) Each sample is mapped to each reference, and variants are called and quality-filtered, then a consensus sequence is built from each mapping.
  3) Mapping staistics are computed.
  4) References and consensus sequences are aligned and the final multiple sequence alignments (MSAs) are masked.
  5) Maximum likelihood phylogenies are inferred from each MSA and differences between trees are evaluated with congruence trests and topological distances.
  6) CDS are extracted and pairwise dN/dS are calculates between consensus sequences from mappings against the same reference.
  7) Recombination rates are calculated along the MSAs.
  8) Plots are generated for each analyses/statistics.
  9) To evaluate the significance of the differences depending on the reference genome, different tests are performed for the results of each analyses.

`change_cns_header.py`: change consensus sequence headers to include sample name and reference strain.

`xmfa_complement.py`: complements a XMFA-formatted MSA of reference sequences with consensus sequences associated with one of these references.

`gaps_xmfa.py`: adds gaps to regions of a XMFA-formatted MSA where homologous sequences are absent in any genome. 

`mask_msa.py`: removes columns with gaps in the sequences specified by the user. It was used for (a) removing positions that were absent in the sequence that reads were mapped to and (b) to remove all positions with gaps, thus obtaining a "core" MSA.

`cds_revcomp.py`: makes reverse complement sequence from a CDS when necessary. Executed by the previous one.

`rho_LDJump.R`: computes population recombination rates (rho) from an alignment in 1000pb windows.

`plots.R`: make boxplots for number of SNPs, percentage of reference covered by reads, mean coverage and pairwise dN/dS values, and the distribution of rho along the alignments.

`stats.R`: Get summary statistics of all the parameters. Perform Kruskal-Wallis and Wilcoxon tests (for mapping statsitics and dN/dS), and pairwise Kolmogorov-Smirnov tests (for rho distributions)

