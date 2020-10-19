#!/usr/bin/env Rscript

# Boxplots (mapping parameters and dN/dS) and recombination rates distributions
# Written by Carlos Valiente-Mullor, 2019
# 
# Usage:
#     Rscript --vanilla plots.R


# Libraries
library(ggplot2)


# Mapstats (boxplots)
myfile="mapstats/mapstats.tsv"
mapst = read.table(myfile, header = T, sep = "\t")

png(filename = "SNPs.png", width = 7.5, height = 6.5, units = "in", res = 600)

ggplot (mapst, aes(x=reference, y=SNPs)) + stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(fill="#E44446FF") + 
  theme(axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 20),
        legend.position="none",
        plot.title = element_text(face = "italic", hjust = 0.5, size = 20)) +
  xlab("Reference") + 
  ylab("SNPs")
  
dev.off()

png(filename = "ref_covered.png", width = 7.5, height = 6.5, units = "in", res = 600)

ggplot (mapst, aes(x=reference, y=per.ref.covered)) + stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(fill="#E44446FF") + 
  theme(axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 20),
        legend.position="none",
        plot.title = element_text(face = "italic", hjust = 0.5, size = 20)) +
xlab("Reference") + 
  ylab("% covered")

dev.off()

png(filename = "coverage.png", width = 7.5, height = 6.5, units = "in", res = 600)

ggplot (mapst, aes(x=reference, y=`mean.coverage`)) + stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(fill="#E44446FF") + 
  theme(axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 20),
        legend.position="none",
        plot.title = element_text(face = "italic", hjust = 0.5, size = 20)) +
  xlab("Reference") + 
  ylab("Coverage")

dev.off()


# dN/dS (boxplot)
file.ref=paste0("codeml_results/dnds.ref.tsv")
file.core=paste0("codeml_results/dnds.core.tsv")

dnds.ref = read.table(file.ref, header = T, sep = "\t")
dnds.core = read.table(file.core, header = T, sep = "\t")

# CDSs from 1 ref
png(filename = "dnds.ref.png", width = 7.5, height = 6.5, units = "in", res = 600)

ggplot (dnds.ref, aes(x=msa, y=omega)) + stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(fill="#E44446FF") + 
  theme(axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 20),
        legend.position="none",
        plot.title = element_text(face = "italic", hjust = 0.5, size = 20)) +
  xlab("Reference") + 
  ylab("dN/dS")

dev.off()

# core CDSs
png(filename = "dnds.core.png", width = 7.5, height = 6.5, units = "in", res = 600)

ggplot (dnds.core, aes(x=msa, y=omega)) + stat_boxplot(geom="errorbar", width=0.5) + geom_boxplot(fill="#E44446FF") + 
  theme(axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 20),
        legend.position="none",
        plot.title = element_text(face = "italic", hjust = 0.5, size = 20)) +
  xlab("Reference") + 
  ylab("dN/dS")

dev.off()


# Recombination rates
myfile = "ldjump_results/rho.tsv"
rdata = read.table(myfile, header = T, sep = "\t")

g = ggplot(data = rdata, aes(x=segment, y=rate)) + 
  geom_line(aes(colour=msa), size = 0.3) +
  xlab("genome region") +
  ylab("rho") +
  facet_grid(rows = vars(msa)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 9),
        legend.title = element_blank(),
        strip.text = element_text(size=12),
        legend.key.size = unit(1.5,"line"))

oext = ".png"
outname = paste0("recombination_rates", oext)

png(filename = outname, width = 7.5, height = 6.5, units = "in", res = 600)
g
dev.off()
