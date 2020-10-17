Pipeline for the evaluation of the effect of reference choice on short-read mapping as performed in the work *One is not enough: on the effects of reference genome for the mapping and subsequent analyses of short-reads. BioRxiv [Preprint]. 2020*. Available from: [https://doi.org/10.1101/2020.04.14.041004](https://www.biorxiv.org/content/10.1101/2020.04.14.041004v1).

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

## `-q` (quality control) [optional]
Path to the file containing the adapters to be removed during read quality processing. Quality trimming and filtering will also be performed.
If no file is provided, QC step will be skipped.
