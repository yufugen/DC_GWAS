#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

trait <- as.character(args[1])

sex <- as.character(args[2]) # sex_combined, female, or male

chr.type <- as.character(args[3]) #autosome or X


# load GENESIS
library(GENESIS)

# load data
load(sprintf("%s_%s_%s.RData",
             sex, trait, chr.type))

###### fit genesis model

if(chr.type == "autosome"){
  M.snp = 1070777
  gwas.data <- sum.stats.auto
}
if(chr.type == "X"){
  M.snp = 38231
  gwas.data <- sum.stats.X
}


fit2 <- genesis(gwas.data, filter=TRUE, modelcomponents=2, cores=2, 
                LDcutoff=0.1, LDwindow=1, c0=10, M = M.snp,
                qqplot = FALSE, summaryGWASdata.save = TRUE)


###### save the results

save(fit2, file = sprintf("./results/%s_%s_%s_genesis_2_model_est.RData", 
                          sex, trait, chr.type))
