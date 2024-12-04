# script to calculate the posterior probability
args = commandArgs(trailingOnly=TRUE)

trait <- args[1]
set.seed(19)
library(tidyverse)
library(dplyr)
library(data.table)
source("./src_stan/calc_posteior_variant.R ")

####################
#### parameters ####
####################

alpha <- sqrt(2)
threshold <- 0.8
groups_name <- c("null","female-biased","equal","male-biased")

#######################
#### lead variants ####
#######################
lead.snp <- read_csv("./data/list_lead_snp.csv")
sum.stats <- filter(lead.snp, TRAIT == trait)

#####  X chromosome #####

## proportion
proportion <- read.table(sprintf("./res/summary_proportion_%s_chrX_M4_unknown_sigma.txt", trait), check.names = FALSE)

## sigmasq
sigmasq <- read.table(sprintf("./res/summary_sigmasq_%s_chrX_M4_unknown_sigma.txt", trait), check.names = FALSE)

dat <- filter(sum.stats, CHR == 23)
scale.factor <- cbind(sqrt(2*dat$A1FREQ_FEMALE*(1-dat$A1FREQ_FEMALE)),
                      sqrt(2*dat$A1FREQ_MALE*(1-dat$A1FREQ_MALE)))

B <- cbind(dat$BETA_FEMALE, dat$BETA_MALE)
# scale the effect size by allele frequency
B <- B*scale.factor
SE <- cbind(dat$SE_FEMALE, dat$SE_MALE)
# scale the SE by allele frequency
SE <- SE*scale.factor
SESQ <- apply(SE, c(1,2), function(x) x^2)
ids <- dat$SNP
chr <- dat$CHR
bp <- dat$BP
A1FREQ <- cbind(dat$A1FREQ_FEMALE, dat$A1FREQ_MALE)

standata <- list("B" = B,
                 "SESQ" = SESQ,
                 "N" = nrow(B),
                 "M" = ncol(B),
                 "ID" = ids,
                 "CHR" = chr,
                 "BP" = bp,
                 "A1FREQ" = A1FREQ)

####### calculate posterior probability #######
posterior.prob <- calc_posteriors(sum.prop = proportion$mean,
                                  sigmasq = sigmasq$mean,
                                  dat = standata,
                                  alpha = alpha)
posterior.prob <- posterior.prob[!is.na(rowSums(posterior.prob[,1:4])),]

###### assign component ######
components <- apply(posterior.prob, 1, function(post){
  max_index <- which.max(post[1:4])
  max_group <- groups_name[max_index]
  if (as.numeric(post[max_index]) <= threshold){
    max_group <- "uncategorized"
  }
  return(max_group)
})
posterior.prob <- cbind(posterior.prob, components)

dat <- left_join(dat, posterior.prob, by = c("SNP" = "ID"))
dat <- dat %>% relocate(TRAIT, .before = ALLELE1)

res.x <- dat
rm(dat)

#####  autosomes #####

## proportion
proportion <- read.table(sprintf("./res/%s_mixprop_sample_autosomes_unknown_sigma.txt", trait), check.names = FALSE, header = TRUE)

## sigmasq
sigmasq <- read.table(sprintf("./res/summary_sigmasq_%s_sample_autosomes_M4_unknown_sigma.txt", trait), check.names = FALSE)

dat <- filter(sum.stats, CHR %in% 1:22)

scale.factor <- cbind(sqrt(2*dat$A1FREQ_FEMALE*(1-dat$A1FREQ_FEMALE)),
                      sqrt(2*dat$A1FREQ_MALE*(1-dat$A1FREQ_MALE)))

B <- cbind(dat$BETA_FEMALE, dat$BETA_MALE)
# scale the effect size by allele frequency
B <- B*scale.factor
SE <- cbind(dat$SE_FEMALE, dat$SE_MALE)
# scale the SE by allele frequency
SE <- SE*scale.factor
SESQ <- apply(SE, c(1,2), function(x) x^2)
ids <- dat$SNP
chr <- dat$CHR
bp <- dat$BP
A1FREQ <- cbind(dat$A1FREQ_FEMALE, dat$A1FREQ_MALE)
standata <- list("B" = B,
                 "SESQ" = SESQ,
                 "N" = nrow(B),
                 "M" = ncol(B),
                 "ID" = ids,
                 "CHR" = chr,
                 "BP" = bp,
                 "A1FREQ" = A1FREQ)

####### calculate posterior probability #######
posterior.prob <- calc_posteriors(sum.prop = proportion$mean,
                                  sigmasq = sigmasq$summary_sigmasq.mean,
                                  dat = standata,
                                  alpha = alpha)
posterior.prob <- posterior.prob[!is.na(rowSums(posterior.prob[,1:4])),]

###### assign component ######

components <- apply(posterior.prob, 1, function(post){
  max_index <- which.max(post[1:4])
  max_group <- groups_name[max_index]
  if (as.numeric(post[max_index]) <= threshold){
    max_group <- "uncategorized"
  }
  return(max_group)
})
posterior.prob <- cbind(posterior.prob, components)

dat <- left_join(dat, posterior.prob, by = c("SNP" = "ID"))
dat <- dat %>% relocate(TRAIT, .before = ALLELE1)

res.a <- dat

res <- bind_rows(res.x,res.a)

write.table(res,
            sprintf("./res/lead_variant_%s_posterior_prob.txt", trait),
            quote = FALSE, row.names = FALSE, col.names = TRUE)
