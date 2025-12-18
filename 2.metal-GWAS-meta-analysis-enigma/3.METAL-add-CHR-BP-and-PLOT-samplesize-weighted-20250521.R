# MANHATTAN PLOT USING SAMPLESIZE WEIGHING FOR BOTH STEPS OF TWO-STAGE META-ANALYSIS + THE OVERALL META-ANALYSIS

########## METAL ##########
library(data.table)
library(tidyverse)
library(stringr)
library(qqman)

PROJECT_ROOT <- "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024"
METAL_DIR <- file.path(PROJECT_ROOT, "generic-metal")
setwd(METAL_DIR)

metal <- fread('METAANALYSIS_ALL_two_stage_GCon_samplesize_info0.6_202507211.tbl')
metal <- as.data.frame(metal)
# explore
head(metal)
dim(metal)
#11380448

metal$`P-value` <- as.numeric(metal$`P-value`)
#metal_10k <- metal %>% filter(TotalSampleSize > 10000)
metal_10k <- metal %>% filter(N > 10000)
nrow(metal_10k)

# check output
gwas.hits <- metal %>% filter(`P-value` < 5e-8) 
nrow(gwas.hits)
gwas.hits2 <- metal_10k %>% filter(`P-value` < 5e-8) 
nrow(gwas.hits2)

# filter out heterogeneous snps
metal_het <- metal_10k[!(metal_10k$HetISq > 0.75 & metal_10k$HetPVal < 0.001), ] # 8477 snps removed
gwas.hits3 <- metal_het %>% filter(`P-value` < 5e-8) 
nrow(gwas.hits3)

# read in ref to get chr and pos
ref <- fread("../reference.1000G.maf.0.005.txt", sep = " ") # which uses build GRCh37 as doubled checked via dbSNP here https://www.ncbi.nlm.nih.gov/snp/

ref$MAF_1000G <- ref$MAF
ref$A1_1000G <- ref$A1
ref$A2_1000G <- ref$A2
# add CHR and BP column from 1000G reference 
metal2 <- merge(metal_het, ref[, c('SNP', 'CHR', 'BP', 'MAF_1000G', 'A1_1000G', 'A2_1000G')], by.x = 'MarkerName', by.y = 'SNP', all.x = TRUE) 

# add CHR and BP for missing entries (extracting it from Markername format that do not have snp ids)
metal3 <- metal2 %>% mutate(CHR = ifelse(is.na(CHR), str_extract(MarkerName, "^[^:]+"), CHR), # extracts everything before the colon to fill the CHR column
                            BP = ifelse(is.na(BP), str_extract(MarkerName, "(?<=:)(\\d+)"), BP) ) # extracts the digits after the colon for the BP column

# ***!!!! differs for stderr vs samplesize schemes ***!!!
# select key columns 
sumstats <- metal3[, c("MarkerName", "CHR", "BP", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Weight", "Zscore", "P-value", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "N")] 
dim(sumstats) # 9844732      18

# see rows where CHR is not in the range 1 to 23
rows_with_others <- sumstats %>% filter(!(CHR %in% as.character(c(1:23, 'X')))) # snps here seem to be mostly only available in the ukbb sumstats
dim(rows_with_others)

# remove chr outside fo 1:23 (or X)
sumstats_clean <- sumstats %>% filter(CHR %in% c(1:23, 'X'))
summary(as.factor((sumstats_clean$CHR))) 

# use same name for chr 23 and X ---> PROBLEMATIC!!! --> converts string to numeric and changed chr nr
# ensure CHR is character to handle mixed types
sumstats_clean$CHR <- as.character(sumstats_clean$CHR)

# replace "X" with "23"
sumstats_clean$CHR[sumstats_clean$CHR == "X"] <- "23"

# convert CHR to numeric
sumstats_clean$CHR <- as.integer(sumstats_clean$CHR)
sumstats_clean$BP <- as.integer(sumstats_clean$BP)

# order by chr and bp
sumstats_ordered <- sumstats_clean[order(sumstats_clean$CHR, sumstats_clean$BP), ]
# 8878894 snps

gwas.hits.clean <- sumstats_ordered %>% filter(`P-value` < 5e-8) 
nrow(gwas.hits.clean)

# rename "P-value" to "P"
names(sumstats_ordered)[names(sumstats_ordered) == "P-value"] <- "P"

### ----- SAVE CLEAN SUMSTATS FOR FUMA ----- ###
fwrite(sumstats_ordered, file = paste0(PROJECT_ROOT, '/METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_', Sys.Date(), '.txt'),
       sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

### ----- MANHATTAN PLOT ----- ###

# # CHR and BP should be numeric
sumstats_ordered$BP <- as.numeric(sumstats_ordered$BP) 

# two stage 
options(bitmapType='cairo') # to disable no X11 warning on HPC
jpeg(paste0(PROJECT_ROOT, "/metal_manhattan_ENIGMA_combinedUKBB_GCon_only-samplesize_", Sys.Date(), ".jpeg"), width = 10, height = 6, units = 'in', res = 300)
manhattan(sumstats_ordered, chr="CHR",bp="BP",p="P", snp='MarkerName', main = "Manhattan plot: METAL two stage, all samplesize weighted ", col = c("slategray2", "blue4"), suggestiveline = -log10(1e-05))
dev.off()

### ----- QQ PLOT ----- ###
# QQ plot of beta p-values
jpeg(paste0(PROJECT_ROOT, "/metal_qq_ENIGMA_combinedUKBB_GCon_only-samplesize_", Sys.Date(), ".jpeg"), width = 6, height = 4, units = 'in', res = 300)
qq(sumstats_ordered$P, main = "Q-Q plot of GWAS p-values: METAL two stage, all samplesize weighted")
dev.off()

