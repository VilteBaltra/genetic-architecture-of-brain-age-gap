  
\# runs metal meta-analysis of ENIGMA cohorts (all smaller sample size
cohorts)  
1.METAL-meta-analysis-info0.6-ENIGMAcohorts-noUKBB.sh  
  
\# run meta-analysis of ENIGMA GWAS (obtained with the first script) and
UKBB GWAS  
2.meta-analysis-two-chunks-ENIGMAcohorts-and-UKBB.sh  
  
\# checks the GWAS sumstats from script 2 and adds rsID from 1000G
reference file  
\# renames chr X to be consistently called chr 23  
\# plots Manhattan + qq plots  
3.METAL-add-CHR-BP-and-PLOT-samplesize-weighted-20250521.R
