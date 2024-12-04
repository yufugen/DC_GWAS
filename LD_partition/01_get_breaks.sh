#!/bin/bash

# LD partition for non-PAR, X chromosome
./ldblock ./data/1000G_EUR_full_female_chrX -out ./res/1000G_EUR_full_female_chrX_ldblocks -min-size 2500

# LD partition for PAR1, X chromosome
./ldblock ./data/1000G.sex.combine.eur.PAR1 -out ./res/1000G.sex.combine.eur.full.PAR1_ldblocks -min-size 2500

# LD partition for PAR2, X chromosome
./ldblock ./data/1000G.sex.combine.eur.PAR2 -out ./res/1000G.sex.combine.eur.full.PAR2_ldblocks -min-size 2500
