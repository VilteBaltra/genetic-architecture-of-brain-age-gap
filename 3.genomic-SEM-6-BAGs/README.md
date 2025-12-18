This script runs **genomic structural equation modeling (genomic SEM)** pipeline for a **common factor GWAS of brain age gap** phenotypes derived from six GWASs.

The workflow includes:
1. Munging individual GWAS summary statistics  
2. LD Score Regression 
3. Exploratory and confirmatory factor analysis  
4. Multivariate (common factor) GWAS  
5. Estimation of effective sample size  

---

## Overview of input GWASs

The following Brain Age Gap GWAS were included:

| Trait name (Genomic SEM) | Description |
|-------------------------|-------------|
| Han | Baltramonaityte et al. (current study) BAG GWAS based on Han et al. (2021) ENIGMA model |
| Jawinski | Jawinski et al. BAG GWAS (excluding overlapping UKBB participants) |
| Pyment | Leonardsen et al. BAG |
| Wen | Wen et al. BAG |
| Kaufman | Kaufman et al. BAG |
| Smith_62modes | Smith et al. (62 imaging modes; A1/A2 swapped) |

All traits are treated as **continuous** and analysed using **OLS**.

---

## Software requirements

- R (≥ 4.1 recommended)
- [`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM)
- R packages:
  - `devtools`
  - `dplyr`
  - `Matrix`
  - `stats`
