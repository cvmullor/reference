#!/usr/bin/env Rscript

# Summary statistics. Kruskal-Wallis + Wilcoxon tests. Pairwise Kolmogorov-Smirnov test
# Written by Carlos Valiente-Mullor, 2019
# 
# Usage:
#     Rscript --vanilla stats.R


# Libraries
library(dplyr)


# Read files (mapstats, dN/dS, recombination rates)
myf.stats = "mapstats/mapstats.tsv"
mstat = read.table(myf.stats, header = T, sep = "\t")

myf.dnds.ref = "codeml_results/dnds.ref.tsv"
myf.dnds.core = "codeml_results/dnds.core.tsv"
mref = read.table(myf.dnds.ref, header = T, sep = "\t")
mcore = read.table(myf.dnds.core, header = T, sep = "\t")

myf.recomb = "ldjump_results/rho.tsv"
mrecomb = read.table(myf.recomb, header = T, sep = "\t")


# Summary statistics
sum.stats = group_by(mstat, reference) %>%
  summarise(
    count = n(),
    mean = mean(SNPs, na.rm = TRUE),
    sd = sd(SNPs, na.rm = TRUE),
    median = median(SNPs, na.rm = TRUE),
    min = min(SNPs, na.rm = TRUE),
    max = max(SNPs, na.rm = TRUE)
  )
write.csv(sum.stats, file = "summary.SNPs.csv", quote=F, row.names=F)

sum.stats = group_by(mstat, reference) %>%
  summarise(
    count = n(),
    mean = mean(per.ref.covered, na.rm = TRUE),
    sd = sd(per.ref.covered, na.rm = TRUE),
    median = median(per.ref.covered, na.rm = TRUE),
    min = min(per.ref.covered, na.rm = TRUE),
    max = max(per.ref.covered, na.rm = TRUE)
  )
write.csv(sum.stats, file = "summary.per_ref_covered.csv", quote=F, row.names=F)

sum.stats = group_by(mstat, reference) %>%
  summarise(
    count = n(),
    mean = mean(mean.coverage, na.rm = TRUE),
    sd = sd(mean.coverage, na.rm = TRUE),
    median = median(mean.coverage, na.rm = TRUE),
    min = min(mean.coverage, na.rm = TRUE),
    max = max(mean.coverage, na.rm = TRUE)
  )
write.csv(sum.stats, file = "summary.mean_coverage.csv", quote=F, row.names=F)

sum.stats = group_by(mref, msa) %>%
  summarise(
    count = n(),
    mean = mean(omega, na.rm = TRUE),
    sd = sd(omega, na.rm = TRUE),
    median = median(omega, na.rm = TRUE),
    min = min(omega, na.rm = TRUE),
    max = max(omega, na.rm = TRUE)
  )
write.csv(sum.stats, file = "summary.dnds_ref.csv", quote=F, row.names=F)

sum.stats = group_by(mcore, msa) %>%
  summarise(
    count = n(),
    mean = mean(omega, na.rm = TRUE),
    sd = sd(omega, na.rm = TRUE),
    median = median(omega, na.rm = TRUE),
    min = min(omega, na.rm = TRUE),
    max = max(omega, na.rm = TRUE)
  )
write.csv(sum.stats, file = "summary.dnds_core.csv", quote=F, row.names=F)

sum.stats = group_by(mrecomb, msa) %>%
  summarise(
    count = n(),
    mean = mean(rate, na.rm = TRUE),
    sd = sd(rate, na.rm = TRUE),
    median = median(rate, na.rm = TRUE),
    min = min(rate, na.rm = TRUE),
    max = max(rate, na.rm = TRUE)
  )
write.csv(sum.stats, file = "summary.recomb_rate.csv", quote=F, row.names=F)


# Kruskal-Wallis + Wilcoxon tests
# dataframe
row.params = c("SNPs", "per.ref.covered", "mean.coverage", "dNdS.ref", "dNdS.core" )
col.test = c("Kruskal_Wallis.chi_squared", "df", "p.value")

kw = data.frame(matrix(ncol = 3, nrow = 5))

colnames(kw) = col.test
rownames(kw) = row.params

# tests
res.kruskal = kruskal.test(SNPs ~ reference, data = mstat)
kw["SNPs",] = c(res.kruskal$statistic[[1]], res.kruskal$parameter[[1]], res.kruskal$p.value)

if(res.kruskal$p.value[1] < 0.05) {
  res.wilcox = pairwise.wilcox.test(mstat$SNPs, mstat$reference, p.adj = "bonf")
  write.csv(res.wilcox$p.value, file="wilcoxon.p_val.SNPs.csv", quote = F)
}

res.kruskal = kruskal.test(per.ref.covered ~ reference, data = mstat)
kw["per.ref.covered",] = c(res.kruskal$statistic[[1]], res.kruskal$parameter[[1]], res.kruskal$p.value)

if(res.kruskal$p.value[1] < 0.05) {
  res.wilcox = pairwise.wilcox.test(mstat$per.ref.covered, mstat$reference, p.adj = "bonf")
  write.csv(res.wilcox$p.value, file="wilcoxon.p_val.per_ref_covered.csv", quote = F)
}

res.kruskal = kruskal.test(mean.coverage ~ reference, data = mstat)
kw["mean.coverage",] = c(res.kruskal$statistic[[1]], res.kruskal$parameter[[1]], res.kruskal$p.value)

if(res.kruskal$p.value[1] < 0.05) {
  res.wilcox = pairwise.wilcox.test(mstat$mean.coverage, mstat$reference, p.adj = "bonf")
  write.csv(res.wilcox$p.value, file="wilcoxon.p_val.mean_coverage.csv", quote = F)
}

res.kruskal = kruskal.test(omega ~ msa, data = mref)
kw["dNdS.ref",] = c(res.kruskal$statistic[[1]], res.kruskal$parameter[[1]], res.kruskal$p.value)

if(res.kruskal$p.value[1] < 0.05) {
  res.wilcox = pairwise.wilcox.test(mref$omega, mref$msa, p.adj = "bonf")
  write.csv(res.wilcox$p.value, file="wilcoxon.p_val.dNdS_ref.csv", quote = F)
}

res.kruskal = kruskal.test(omega ~ msa, data = mcore)
kw["dNdS.core",] = c(res.kruskal$statistic[[1]], res.kruskal$parameter[[1]], res.kruskal$p.value)

if(res.kruskal$p.value[1] < 0.05) {
  res.wilcox = pairwise.wilcox.test(mcore$omega, mcore$msa, p.adj = "bonf")
  write.csv(res.wilcox$p.value, file="wilcoxon.p_val.dNdS_core.csv", quote = F)
}

# save
write.csv(kw, file = "test.kruskal-wallis.csv", quote=F, row.names=T)


# Pairwise Kolmogorov-Smirnov test
# function is available at: https://rdrr.io/github/netlify/NetlifyDS/src/R/pairwise_ks_test.R
pairwise_ks_test <- function(value, group, n_min = 50, warning = 0, alternative = "two.sided" ){
  
  lev <- unique(group)
  
  lst <- lapply( seq_along(lev), function(i) value[group == lev[i]] )
  names(lst)<-lev
  
  if (sum(lengths(lst)< n_min)) {
    lst <- lst [-which(lengths(lst)< n_min)]}
  
  f <- function(x, y){ 
    w <- getOption("warn") 
    options(warn = warning)  # ignore warnings 
    p <- ks.test(x, y, alternative = alternative, exact = 
                   F)$p.value 
    options(warn = w) 
    return(p) 
  } 
  
  res <- lapply(lst, function(x) lapply(lst, function(y) f(x, y))) 
  
  res<-unlist(res)
  res <- matrix(res, nrow = length(lst), ncol = length(lst), byrow = T)
  row.names(res) <- colnames(res) <- names(lst)
  cat("Pairwise Kolmogorov-Smirnov Test p-value Matrix","\n","\n")
  return(res)
}

# test
value = mrecomb$rate
group = mrecomb$msa
ks = pairwise_ks_test(value, group, warning = -1)
write.csv(ks, file = "kolmogorov.p_val.csv")
