# Stage-2 meta-analysis combining:
# 1) Non-UKBB meta (INFO 0.6, GC ON, MAF>0.01)
# 2) UKBB discovery meta (BOLT-LMM)

# META-ANALYSIS OF no UKBB cohorts with UKBB cohorts
./metal
SCHEME SAMPLESIZE  
AVERAGEFREQ ON
MINMAXFREQ ON
ADDFILTER MAF > 0.01
CUSTOMVARIABLE N

MARKER   MarkerName
ALLELE   Allele1 Allele2
FREQ     Freq1
EFFECT   Zscore
WEIGHT TotalSampleSize
PVAL     P-value
LABEL N as TotalSampleSize

# process ENIGMA-only GWAS
PROCESS METAANALYSIS_all_noUKBB_GCon_samplesize_info0.6_20250721_1.tbl

# UKBB combined (discov1 and discov2)
MARKER   SNP
WEIGHT   N
STDERRLABEL SE
ALLELE   ALLELE1 ALLELE0
FREQ     A1FREQ
FREQLABEL MAF
EFFECT   BETA
STDERR   SE
PVAL     P_BOLT_LMM_INF
LABEL N as N

# process UKBB-only GWAS derived with BOLT-LMM
PROCESS UKBB_discov1to2_sumstats_main_rsid_info0.6.txt.gz

# no MAF filter as no MAF column found (its the freq1 one) but no need to add as already filtered for MAF 0.01 in previous steps
OUTFILE METAANALYSIS_ALL_two_stage_GCon_samplesize_info0.6_20250721 .tbl
MINWEIGHT 5000
ANALYZE HETEROGENEITY
QUIT

