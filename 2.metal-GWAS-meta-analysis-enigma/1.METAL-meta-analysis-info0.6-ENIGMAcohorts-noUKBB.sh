# BrainAge GWAS meta-analysis of ENIGMA cohorts (no UKBB)
# Excludes: UKBB, GenR
# MAF: prefiltered to 0.01
# INFO threshold: 0.6

#!/bin/bash
# this script runs METAL meta-analysis for BrainAge GWAS
# download precompiled binary (https://csg.sph.umich.edu/abecasis/Metal/download/)
cd ${PROJECT_ROOT}/generic-metal
wget https://csg.sph.umich.edu/abecasis/metal/download/Linux-metal.tar.gz 
tar -xvzf Linux-metal.tar.gz


### STAGE 1: META-ANALYSIS OF ALL COHORTS (but not UKBB)

./metal
#GENOMICCONTROL ON # automatically correct test statistics to account for small amounts of population stratification or unaccounted for relatedness
#AVERAGEFREQ ON # alspac lacks MAF, need to add, then can run this
# MINMAXFREQ ON
# VERBOSE OFF
#Selecting an Analysis Scheme (below list available options, in "MY RUN" I use 'SAMPLESIZE')
# SCHEME SAMPLESIZE        - default approach, uses p-value and direction of effect, weighted according to sample size
# SCHEME STDERR            - classical approach, uses effect size estimates and standard errors
# STDERR SE                - specify the label for the standard error column.

# MY RUN (useful guide for Metal code: https://github.com/huw-morris-lab/meta-analysis)
#SCHEME STDERR  
GENOMICCONTROL ON 
SCHEME SAMPLESIZE  
AVERAGEFREQ ON
MINMAXFREQ ON
ADDFILTER MAF > 0.01
#VERBOSE ON
CUSTOMVARIABLE TotalSampleSize

## ## Meta-analysis will be based on sample sizes, p-values and direction of effect ...

MARKER   SNP
WEIGHT   N
STDERRLABEL SE
ALLELE   A1 A2
FREQ     EAF
FREQLABEL MAF
EFFECT   BETA
STDERR   SE
PVAL     P
LABEL TotalSampleSize as N

#CUSTOMVARIABLE MAF

# define project dir
PROJECT_ROOT=/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024

# RMW
PROCESS ${PROJECT_ROOT}/gwas-results/ACP/ACP_sumstats_main_chr1to22_rsid_keycolumns_noinfo.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/BHRCS/BHRCS_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/BILGIN/BILGIN_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/brazilunicamp/brazilunicamp_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/COBRE/COBRE_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/COMPIMP/COMPIMP_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/DBIS/DBIS_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/DCHS/DCHS_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/DNS/DNS_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/FBIRN/FBIRN_sumstats_main_rsid_keycolumns_info0.6.txt.gz 
PROCESS ${PROJECT_ROOT}/gwas-results/FIDMAG/FIDMAG_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/GOBS/GOBS_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/IGP/IGP_sumstats_main_chr1to22_rsid_keycolumns.txt.gz 
PROCESS ${PROJECT_ROOT}/gwas-results/imagen/imagen_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/IMH_study/125samples/IMH_study_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/LBC1936/LBC1936_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/mprc/mprc_sumstats_main_chr1to22_rsid_keycolumns.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/PAFIP_3/PAFIP_3_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/VOX/VOX_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/SHIP-2/SHIP-2_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/SHIP-T/SHIP-T_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/SHIP-T_B2/SHIP-T_B2_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/GIG/GIG_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/OSLO/OSLO_sumstats_main_chr1to22_rsid_keycolumns_noinfo.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/SYS/SYS_sumstats_main_rsid_keycolumns_info0.6.txt.gz
PROCESS ${PROJECT_ROOT}/gwas-results/NESDA/NESDA_sumstats_main_rsid_keycolumns_info0.6.txt.gz

# REGENIE
# ID	CHR	POS	A0	A1	FreqA1	Rsq	Beta	SE	p	N	MAF	SNP
MARKER   SNP
WEIGHT   N
STDERRLABEL SE
ALLELE   A1 A0
FREQ     FreqA1
FREQLABEL MAF
EFFECT   Beta
STDERR   SE
PVAL     p
PROCESS ${PROJECT_ROOT}/gwas-results/ISHARE/ISHARE_sumstats_main_rsid_keycolumns_info0.6.txt.gz

# ALSPAC (plink)
# ID #CHROM POS REF ALT A1 TEST OBS_CT BETA SE T_STAT P ERRCODE ALT_FREQS MAF
MARKER   ID
WEIGHT   OBS_CT
STDERRLABEL SE
ALLELE   ALT REF
FREQ     ALT_FREQS
FREQLABEL MAF
EFFECT   BETA
STDERR   SE
PVAL     P
LABEL TotalSampleSize as OBS_CT
PROCESS ${PROJECT_ROOT}/gwas-results/ALSPAC/alspac_sumstats_main_chr1to22_rsid_maf05.txt

# Execute meta-analysis
OUTFILE METAANALYSIS_all_noUKBB_GCon_samplesize_info0.6_20250721_ .tbl
MINWEIGHT 5000
ANALYZE HETEROGENEITY
QUIT

