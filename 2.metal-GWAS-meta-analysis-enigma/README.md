### Scripts overview

`1.METAL-meta-analysis-info0.6-ENIGMAcohorts-noUKBB.sh`
- Runs METAL meta-analysis of ENIGMA cohorts (all smaller sample size cohorts, excluding UKBB).

`2.meta-analysis-two-chunks-ENIGMAcohorts-and-UKBB.sh`
- Runs meta-analysis combining ENIGMA GWAS (from script 1) and UKBB GWAS.

`3.METAL-add-CHR-BP-and-PLOT-samplesize-weighted-20250521.R`
- Checks GWAS summary statistics from script 2 
- Adds rsID from the 1000 Genomes reference file  
- Renames chromosome X to chr 23 for consistency
- Generates Manhattan and QQ plots
