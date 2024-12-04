# Role of chromosome X and dosage compensation mechanisms in complex trait genetics

Provided below are resources and scripts used in "Interpreting dosage compensation through complex trait genome-wide associations".

## Resources

+ GWAS summary statistics with UK Biobank and FinnGen R10 data from Zenodo.

## Analyses

### Estimation of SNP-heritability and effect size distribution in the X chromosome and autosomes

We used GENetic Effect-Size distribution Inference from Summary-level data ([GENESIS](https://github.com/yandorazhang/GENESIS)) to estimate the SNP-heritability and effect size distribution in autosomes and chrX. For the GENESIS to work on chrX, replace the `LDwindow1MB_cutoff0.1.RData` and `w_hm3.noMHC.snplist.RData` in the original data of GENESIS with the one found in `GENESIS_analysis/data/` in this repository. We only included the default LD reference panel for the X chromosome (i.e. for X chromosome analysis, the `LDcutoff = 0.1` and `LDwindow = 1` when use `genesis()` to perform the model fitting).


### LD blocks in the X chromosome and autosomes

We partitioned LD blocks for the X chromosome, non-PAR, PAR1, and PAR2, separately. The partition was performed with [LAVA partitioning algorithm](https://github.com/cadeleeuw/lava-partitioning/) using 1000 Genomes (EUR) data. Female data were used when performing LD partition in non-PAR, and both sexes were used for LD partition in PAR. Data in plink format avaiable in `LD_partition/data`.

The autosomal LD blocks were the original blocks created for [LAVA](https://github.com/josefin-werme/lava) and downloaded from [here](https://github.com/cadeleeuw/lava-partitioning/). The full LD blocks (autosomes + X chromosome) can be found in `LD_partition/res/locdef.locfile.eur.complete.par1.par2`.

### Four-component sex bias mixture model of genome-wide variants

Bayesian mixutre model to estimate proportion of sex-biased variants and identify sex-biased variants using GWAS summary statistics of quantitative traits.

The model relies on [`rstan`](https://mc-stan.org/users/interfaces/rstan).
