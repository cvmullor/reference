Code used in Valiente-Mullor C et al. One is not enough: on the effects of reference genome for the mapping and subsequent analyses of short-reads. BioRxiv:2020.04.14.041004 [Preprint]. 2020. Available from: [https://doi.org/10.1101/2020.04.14.041004](https://www.biorxiv.org/content/10.1101/2020.04.14.041004v1). Usage is detailed in the header section of each script.

The workflow is summarized in the Figure below (Fig 1).
![Fig 1](Fig1_overview.png)


`remove_adaptors.sh`: removes adaptors from raw reads.

`trim_filter.sh`: quality filtering and trimming of raw reads.

`index.sh`: indexes reference sequences with BWA and samtools.

`mapsort.sh`: read mapping against each reference of the same species.

`nsnp.sh`, `covered.sh` and `coverage.sh`: obtain different mapping parameters (number of SNPs, percentage of reference covered by reads and average coverage of the reference, respectively).

`consensus.sh`: calls variants and computes the consensus sequence for each mapping.

`xmfa_complement.py`: Complements an XMFA-formatted MSA of reference sequences with consensus sequences associated with different/one of these references.

`mask_msa.py`: Removes selected/certain columns with gaps (in the sequences specified by the user) from the MSAs. It was used for (a) removing positions that were absent in the sequence that reads were mapped to and also (b) to remove all positions with gaps, thus obtaining a "core" MSA.

`ml_trees.sh`: calls IQ-TREE for building of maximum likelihood trees from each final MSA after removing columns with gaps.

`congruence_test`: performs ELW tests (IQ-TREE) for MSAs against a set of tree topologies. It was used to assess the congruence between phylogenies that differed only in the genome chosen as reference for mapping.

`extract_cds.sh` and `extract_cds_core.sh`: Produces/Obtains an alignment of concatenated CDSs of consensus sequences derived/obtained from mappings against the same reference, including all CDSs of the reference or only core genome CDSs, respectively. Outputs were used for dN/dS analyses.

`cds_revcomp.py`: Makes reverse complement sequence from a CDS when necessary. Executed by the previous one.
