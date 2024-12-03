# Role of chromosome X and dosage compensation mechanisms in complex trait genetics

Provided below are resources and scripts used in "Interpreting dosage compensation through complex trait genome-wide associations".

## Resources

+ GWAS summary statistics with UK Biobank and FinnGen R10 data from Zenodo.

## Analyses

### Estimation of SNP-heritability and effect size distribution in chrX and autosomes

We used [GENESIS](https://github.com/yandorazhang/GENESIS) to estimate the SNP-heritability and effect size distribution in autosomes and chrX. For the GENESIS to work on chrX, replace the `LDwindow1MB_cutoff0.1.RData` and `w_hm3.noMHC.snplist.RData` in the original data of GENESIS with the one found in `GENESIS_analysis/data/` in this repository. We only included the default LD reference panel for the X chromosome (i.e. for X chromosome analysis, the `LDcutoff = 0.1` and `LDwindow = 1` when use `genesis()` to perform the model fitting).


### LD blocks in chrX and autosomes

### Four-component sex bias mixture model of genome-wide variants
