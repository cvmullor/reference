#!/usr/bin/env Rscript

# Compute population recombination rates with LDJump
# Written by Carlos Valiente-Mullor, 2019
# 
# Usage:
#     Rscript --vanilla rho_LDJump.R <reference.strain> <full/actual/path>


# Libraries
require(LDJump)
library(LDJump)
library(ggplot2)

# Arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

strain = args[1] # strain name of the mapping reference
fp = args[2]

msafile=paste0(fp, "/complete_msa/masked_core.def.wgaps.xmfa_final.", strain, ".fasta") # MSA

# Required
ldhat = ""
phipack = ""

# Compute rho with LDJump
rho = LDJump(msafile,
             alpha = 0.05, 
             segLength = 1000, 
             pathLDhat = ldhat, 
             pathPhi = phipack, 
             format = "fasta", 
             refName = NULL, 
             start = NULL, 
             constant = F, 
             status = T, 
             cores = 1, 
             accept = T, 
             demography = F, 
             out = "")

# Save
outfile = paste0(strain, ".cte_estimates.txt")
write(rho$`Constant estimates:`, file = outfile)



