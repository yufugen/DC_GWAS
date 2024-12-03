#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

###############################################################################
### SET PARAMETERS

MAIN_PATH = "./gwas_res_STAN"
SRC_PATH =  "./src_stan"
OUTPUT_PATH = "./res"
USE_THESE_PACKAGES <- c("logger","tidyverse", "rstan", "readr", "ggplot2")

SAVE_OUTPUT = TRUE
GENERATE_PLOT = TRUE
CHECK_CONVERGENCE = TRUE

TRAIT = args[1]

CHR_TYPE = "X"

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
# load the stan model
XCI_model <- stan_model(file = sprintf("%s/sex_bias_model.stan", SRC_PATH))

###############################################################################
### LOAD FUNCTIONS
source(sprintf("%s/extract_pi_values.R", SRC_PATH))

###############################################################################
### IMPORT DATA

message(sprintf("Inputting trait: %s", TRAIT))

log_info(sprintf("execute_analysis: Importing GWAS results of %s ...", TRAIT))
gwas_res <- read_table(file = sprintf("%s/%s_gwas_res_selected_vars_STAN.txt",MAIN_PATH,TRAIT))

# subset data, X chromosome data 
dat <- filter(gwas_res, CHR == 23)

log_info("execute_analysis: The dataset includes {nrow(dat)} SNPs for CHR {CHR_TYPE}.")


###############################################################################
### PREPARE DATA

### NOTE! THIS SHOULD BE SCALED DATA!
scale.factor <- cbind(sqrt(2*dat$A1FREQ_FEMALE*(1-dat$A1FREQ_FEMALE)),
                      sqrt(2*dat$A1FREQ_MALE*(1-dat$A1FREQ_MALE)))

B <- cbind(dat$B.F, dat$B.M)
# scale the effect size by allele frequency
B <- B*scale.factor

SE <- cbind(dat$SE.F, dat$SE.M)
# scale the SE by allele frequency
SE <- SE*scale.factor

# the model needs squared SE
SESQ <- apply(SE, c(1,2), function(x) x^2)

# other information to keep
ids <- dat$SNP
chr <- dat$CHR
bp <- dat$BP
PVALS <- cbind(dat$P.F, dat$P.M)
A1FREQ <- cbind(dat$A1FREQ_FEMALE, dat$A1FREQ_MALE)
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

log_info("execute_analysis: standata prepared for modelling.")

###############################################################################
### STAN MODEL

options(mc.cores = parallel::detectCores())

fit <- sampling(XCI_model,
                data = standata,
                chains = 4,
                warmup = n.warmup,
                iter = n.iter,
                save_warmup = TRUE,
                seed = master_seed)

model_time <- get_elapsed_time(fit)/60
log_info("STAN model complete. Run time(min):")
print(model_time)

saveRDS(fit, file = sprintf("%s/%s_chr%s_%s_unknown_sigma.rds",
                            OUTPUT_PATH, TRAIT, CHR_TYPE, NO_COMPO))
###############################################################################
### CHECK CONVERGENCE

if(CHECK_CONVERGENCE){
  # TRACEPLOT <- traceplot(fit,pars = "pi", inc_warmup = TRUE)
  # 
  # ggsave(sprintf("%s/traceplot_%s_%s_%s_%s.png", OUTPUT_PATH, TRAIT, CHR_TYPE, NO_COMPO, MODEL),
  #        TRACEPLOT, device = "png", width = 14, height = 8, dpi = "print")
  # 
  R_hat_sum <- summary(fit, "pi")$summary %>%
    as.data.frame() %>%
    rownames_to_column("parameter") %>%
    select(parameter, Rhat)
  print("Check convergence:")
  print(R_hat_sum)
}

################################################################################
### EXTRACT ESTIMATION OUTPUT & SAVE

fit_pi <- extract_fit_pi(fit)
print(fit_pi)
fit_sigmasq <- extract_fit_vars(fit)
print(fit_sigmasq)
# save the results into two tables
# one contains the estimated proportions
# one contains the posterior probabilities and assignment for each variant
if(SAVE_OUTPUT){
  log_info("Export results into delimited tables")
  write.table(fit_pi[[1]],
              file = sprintf("%s/summary_proportion_%s_chr%s_%s_unknown_sigma.txt", 
                             OUTPUT_PATH, TRAIT, CHR_TYPE, NO_COMPO),
              quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(fit_sigmasq[[1]],
              file = sprintf("%s/summary_sigmasq_%s_chr%s_%s_unknown_sigma.txt", 
                             OUTPUT_PATH, TRAIT, CHR_TYPE, NO_COMPO),
              quote = FALSE, row.names = TRUE, col.names = TRUE)
}

log_info("Script completed...")
