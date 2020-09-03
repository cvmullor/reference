Code used in the work *One is not enough: on the effects of reference genome for the mapping and subsequent analyses of short-reads. BioRxiv [Preprint]. 2020*. Available from: [https://doi.org/10.1101/2020.04.14.041004](https://www.biorxiv.org/content/10.1101/2020.04.14.041004v1). Usage is detailed in the header section of each script.

The workflow is summarized in Fig 1:
![Fig 1](Fig1_overview.png)


`remove_adaptors.sh`: removes adaptors from raw reads.

`trim_filter.sh`: quality filtering and trimming of raw reads.

`index.sh`: indexes reference sequences with BWA and samtools.

`mapsort.sh`: read mapping against each reference of the same species.

`consensus.sh`: calls variants and computes the consensus sequence for each mapping.

`nsnp.sh`, `covered.sh` and `coverage.sh`: obtain different mapping parameters (number of SNPs, percentage of reference covered by reads and average coverage of the reference, respectively).

`xmfa_complement.py`: complements a XMFA-formatted MSA of reference sequences with consensus sequences associated with one of these references.

`gaps_xmfa.py`: adds gaps to regions of a XMFA-formatted MSA where homologous sequences are absent in any genome. 

`mask_msa.py`: removes columns with gaps in the sequences specified by the user. It was used for (a) removing positions that were absent in the sequence that reads were mapped to and (b) to remove all positions with gaps, thus obtaining a "core" MSA.

`ml_trees.sh`: builds maximum likelihood phylogenies (IQ-TREE) from each MSA.

`congruence_test.sh`: performs ELW tests (IQ-TREE) for MSAs against a set of tree topologies. It was used to assess the congruence between phylogenies that differed only in the genome chosen as reference for mapping.

`extract_cds.sh` and `extract_cds_core.sh`: extract and concatenate CDSs of consensus sequences obtained from mappings against the same reference, considering all CDSs of the reference or only core genome, respectively. Outputs were used for dN/dS analyses.

`cds_revcomp.py`: makes reverse complement sequence from a CDS when necessary. Executed by the previous one.
