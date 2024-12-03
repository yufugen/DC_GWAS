#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

###############################################################################
### SET PARAMETERS

MAIN_PATH = "./gwas_res_STAN"
SRC_PATH =  "./src_stan" 
OUTPUT_PATH = "./res"
USE_THESE_PACKAGES <- c("logger","tidyverse", "rstan", "readr", "ggplot2", "data.table")

TRAIT = args[1]
BIOBANK = args[2] #UKB or FinnGen

NO_COMPO = "M4"
ALPHA = sqrt(2) # expected ratio of sex specific effect size
n.iter <- 4000
n.warmup <- 2000

master_seed = 365
options(mc.cores = parallel::detectCores())

###############################################################################
### IMPORT PACKAGES

new_packages <- USE_THESE_PACKAGES[!(USE_THESE_PACKAGES %in% installed.packages()[,"Package"])]
if(length(new_packages)){
  install.packages(new_packages)
  sapply(USE_THESE_PACKAGES, require, character.only = TRUE)
} else {
  sapply(USE_THESE_PACKAGES, require, character.only = TRUE)
}

###############################################################################
### LOAD THE MODEL
XCI_model <- stan_model(file = sprintf("%s/sex_bias_model.stan", SRC_PATH))

###############################################################################
### LOAD FUNCTIONS
source(sprintf("%s/extract_pi_values.R", SRC_PATH))

###############################################################################
### IMPORT DATA

message(sprintf("Inputting trait: %s", TRAIT))

log_info(sprintf("execute_analysis: Importing GWAS results of %s ...", TRAIT))
gwas_res <- read_table(file = sprintf("%s/%s_gwas_res_selected_vars_STAN.txt",MAIN_PATH,TRAIT))

# subset data, autosomal data
dat <- filter(gwas_res, CHR %in% 1:22)

###############################################################################
### PREPARE DATA 

## Map to LD group data
# load LD group, UKB and FinnGen use different coordination
if(BIOBANK == "UKB"){
  dat$VAR <- paste(dat$CHR, dat$BP, dat$ALLELE0, dat$ALLELE1, dat$SNP,
                   sep = ":")
  LD_group <- read.table(sprintf("%s/chr1-22_noMHC_variant_rsid_pos_info_LDgroup.txt", SRC_PATH),
                         header = TRUE)
  LD_group$VAR1 <- paste(LD_group$chr, LD_group$posBP, LD_group$ref, LD_group$alt, LD_group$rsid,
                         sep = ":")
  LD_group$VAR2 <- paste(LD_group$chr, LD_group$posBP, LD_group$alt, LD_group$ref, LD_group$rsid,
                         sep = ":")
  LD_group_short <- LD_group[,c("VAR1","VAR2","LDgroup")]
  
  dat_1 <- left_join(dat, LD_group_short, by = c("VAR" = "VAR1"))
  dat_1 <- dat_1[complete.cases(dat_1),]
  dat_2 <- left_join(dat, LD_group_short, by = c("VAR" = "VAR2"))
  dat_2 <- dat_2[complete.cases(dat_2),]
  
  dat <- bind_rows(dat_1, dat_2)
  dat <- dat %>% select(-VAR, -VAR2, -VAR1)
  dat <- as.data.frame(dat)
  row.names(dat) <- NULL
  
} else if(BIOBANK == "FinnGen"){
  dat$VAR <- paste(dat$CHR, dat$BP, dat$ALLELE0, dat$ALLELE1,
                   sep = ":")
  LD_group <- read.table(sprintf("%s/chr1-22_noMHC_variant_rsid_pos_info_LDgroup_GRCh38.txt", SRC_PATH),
                         header = TRUE)
  LD_group$VAR1 <- paste(LD_group$chr, LD_group$posBP, LD_group$ref, LD_group$alt,
                         sep = ":")
  LD_group$VAR2 <- paste(LD_group$chr, LD_group$posBP, LD_group$alt, LD_group$ref, 
                         sep = ":")
  LD_group_short <- LD_group[,c("VAR1","VAR2","LDgroup")]
  
  dat_1 <- left_join(dat, LD_group_short, by = c("VAR" = "VAR1"))
  dat_1 <- dat_1[complete.cases(dat_1),]
  dat_2 <- left_join(dat, LD_group_short, by = c("VAR" = "VAR2"))
  dat_2 <- dat_2[complete.cases(dat_2),]
  
  dat <- bind_rows(dat_1, dat_2)
  dat <- dat %>% select(-VAR, -VAR2, -VAR1)
  dat <- as.data.frame(dat)
  row.names(dat) <- NULL
}

### random subset
## randomly select 6 variants per LD block
dat$index <- seq.int(nrow(dat))
random_subset <- function(seed=1, n.variant = 6) {
  set.seed(seed); print(paste0("Seed #: ", seed))
  # METHOD3: LD blocks
  random <- numeric(0)
  for (i in unique(dat$LDgroup)) {
    sample_subset <- dat[dat$LDgroup == i, 'index']
    n.sub <- length(sample_subset)
    if(n.sub < n.variant){
      random[(length(random) + 1): (length(random) + n.sub)] <- sample_subset
    } else {
      random[(length(random) + 1):(length(random) + n.variant)] <- sample(sample_subset,n.variant)
    }
  }
  print(paste0("random subset length: ", length(random)))
  return(random)
}
fit_STAN <- function(random) {
  data.random <- filter(dat, index %in% random)
  
  scale.factor <- cbind(sqrt(2*data.random$A1FREQ_FEMALE*(1-data.random$A1FREQ_FEMALE)),
                        sqrt(2*data.random$A1FREQ_MALE*(1-data.random$A1FREQ_MALE)))
  B <- cbind(data.random$B.F, data.random$B.M)
  # scale the effect size by allele frequency
  B <- B*scale.factor
  SE <- cbind(data.random$SE.F, data.random$SE.M)
  SE <- SE*scale.factor
  SESQ <- apply(SE, c(1,2), function(x) x^2)
  ids <- data.random$SNP
  chr <- data.random$CHR
  bp <- data.random$BP
  PVALS <- cbind(data.random$P.F, data.random$P.M)
  A1FREQ <- cbind(data.random$A1FREQ_FEMALE, data.random$A1FREQ_MALE)
  alpha <- ALPHA
  
  standata <- list("B" = B,
                   "SESQ" = SESQ,
                   "alpha" = alpha,
                   "K" = 4,
                   "N" = nrow(B),
                   "M" = ncol(B),
                   "ID" = ids,
                   "CHR" = chr,
                   "BP" = bp,
                   "PVALS" = PVALS,
                   "A1FREQ" = A1FREQ)
  
  # fit STAN model
  
  fit <- sampling(XCI_model,
                  data = standata,
                  chains = 4,
                  warmup = n.warmup,
                  iter = n.iter,
                  save_warmup = TRUE,
                  seed = master_seed)
  return(fit)
}

###############################################################################
### FIT STAN MODEL

res_list <- list()

random <- random_subset()
# sampling with STAN
results <- fit_STAN(random)
mix <- extract_fit_pi(results)[[1]]


saveRDS(results, file = sprintf("%s/%s_sample_autosomes_unknown_sigma.rds", OUTPUT_PATH, TRAIT))
write.table(mix, file = sprintf("%s/%s_mixprop_sample_autosomes_unknown_sigma.txt", OUTPUT_PATH, TRAIT),
            sep="\t", row.names=TRUE, quote = FALSE)
