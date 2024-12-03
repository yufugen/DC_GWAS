# extracting GWAS data for STAN model

library(data.table)
library(tidyverse)
library(dplyr)
library(readr)

args = commandArgs(trailingOnly=TRUE)

input.dir <- as.character(args[1])
trait <- as.character(args[2])
output.dir <- as.character(args[3])

# read in BOLT GWAS result, using bgen, as it contain more SNPs
# fread() is fast in importing large data
res.f <- fread(Sys.glob(file.path(input.dir, 
                                  sprintf("BOLT_female*%s.bgen.stats.gz", trait))))
res.m <- fread(Sys.glob(file.path(input.dir, 
                                  sprintf("BOLT_male*%s.bgen.stats.gz", trait))))
res.c <- fread(Sys.glob(file.path(input.dir, 
                                  sprintf("BOLT_sex_combined*%s.bgen.stats.gz", trait))))

select_vars <- read.table("./src_stan/full_chr_qc_ld_table_noMHC.txt", header=FALSE)

female_vars <- subset(res.f, SNP %in% select_vars[,1]) %>%
  select(SNP, CHR, BP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF)%>%
  rename(A1FREQ_FEMALE = A1FREQ,
         BETA_FEMALE = BETA,
         SE_FEMALE = SE,
         P_BOLT_LMM_INF_FEMALE = P_BOLT_LMM_INF)

male_vars <- subset(res.m, SNP %in% select_vars[,1]) %>%
  select(SNP, CHR, BP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF)%>%
  rename(A1FREQ_MALE = A1FREQ,
         BETA_MALE = BETA,
         SE_MALE = SE,
         P_BOLT_LMM_INF_MALE = P_BOLT_LMM_INF)

combine_vars <- subset(res.c, SNP %in% select_vars[,1]) %>%
  select(SNP, CHR, BP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF)%>%
  mutate(MAF = ifelse(A1FREQ < 0.5, A1FREQ, 1 - A1FREQ)) %>%
  rename(A1FREQ_COMB = A1FREQ,
         MAF_COMB = MAF,
         BETA_COMB = BETA,
         SE_COMB = SE,
         P_BOLT_LMM_INF_COMB = P_BOLT_LMM_INF)

if (all(all(female_vars$SNP == male_vars$SNP), all(combine_vars$SNP == male_vars$SNP))){
  
  all.tab <- left_join(combine_vars, left_join(female_vars, male_vars))
  all.tab <- filter(all.tab, MAF_COMB > 0.01)
  
  # Check for multiallelic variants exist
  # If yes, remove them
  if (any(duplicated(all.tab$SNP))){
    dup.snp.list <- all.tab$SNP[duplicated(all.tab$SNP)]
    all.tab <- all.tab[!all.tab$SNP %in% dup.snp.list, ]
    out <- all.tab
  }
  
  out <- all.tab
}

dat.stan <- out %>%
  select(SNP, CHR, BP, ALLELE1, ALLELE0, A1FREQ_FEMALE,A1FREQ_MALE, 
         BETA_FEMALE, BETA_MALE,BETA_COMB,
         SE_FEMALE, SE_MALE,SE_COMB,
         P_BOLT_LMM_INF_FEMALE, P_BOLT_LMM_INF_MALE,P_BOLT_LMM_INF_COMB) %>%
  rename(B.F = BETA_FEMALE,
         B.M = BETA_MALE,
         B.COMB = BETA_COMB,
         SE.F = SE_FEMALE,
         SE.M = SE_MALE,
         SE.COMB = SE_COMB,
         P.F = P_BOLT_LMM_INF_FEMALE,
         P.M = P_BOLT_LMM_INF_MALE,
         P.COMB = P_BOLT_LMM_INF_COMB)

write_delim(dat.stan,
            file = sprintf("%s/%s_gwas_res_selected_vars_STAN.txt", 
                           output.dir,trait))

