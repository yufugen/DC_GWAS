#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

trait <- as.character(args[1])
sex <- as.character(args[2])
file.name <- as.character(args[3])

library(tidyverse)
library(readr)
library(data.table)

sum.stats <- read_delim(file.name)

###############################
###  effective sample size  ###
###############################

sample.files <- c("ukb22627_imp_chr1-22_v3_s487395.sample",
                  "ukb22627_imp_chrX_v3_s486743.sample",
                  "ukb22627_imp_chrXY_v3_s486429.sample")

incl.files <- sprintf("%s_BOLT_%s.incl", sex, trait)
incl.dat <- fread(incl.files)

##### autosomes #####
sample.dat <- fread(sample.files[1])[-1,]
autosome.sam <- length(intersect(sample.dat$ID_1, incl.dat$V1))

##### X chromosome, non-PAR #####
sample.dat <- fread(sample.files[2])[-1,]
X.sam <- length(intersect(sample.dat$ID_1, incl.dat$V1))

##### X chromosome, PAR #####
sample.dat <- fread(sample.files[3])[-1,]
XY.sam <- length(intersect(sample.dat$ID_1, incl.dat$V1))

###########################
#### Assemble the data ####
###########################

# autosome and X chromoosme separate

sum.stats.auto <- sum.stats %>%
  filter(CHR %in% 1:22) %>%
  mutate(effective_n = autosome.sam,
         z = BETA/SE) %>%
  select(SNP, z, effective_n) %>%
  rename(snp = SNP,
         n = effective_n)


sum.stats.X <- sum.stats %>%
  filter(CHR == 23) %>%
  mutate(effective_n = ifelse(IS_PAR == "yes", XY.sam, X.sam),
         z = BETA/SE) %>%
  select(SNP, z, effective_n) %>%
  rename(snp = SNP,
         n = effective_n)

save(sum.stats.auto, file = sprintf("%s_%s_autosome.RData", trait, sex, trait),
     ascii = FALSE, version = 2,
     compress = "xz", compression_level = 6)

save(sum.stats.X, file = sprintf("%s_%s_X.RData", trait, sex, trait),
     ascii = FALSE, version = 2,
     compress = "xz", compression_level = 6)
