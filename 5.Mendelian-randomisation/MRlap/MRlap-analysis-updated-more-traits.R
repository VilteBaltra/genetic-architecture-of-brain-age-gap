

# Directly install the package from github
# install.packages("remotes")
#remotes::install_github("n-mounier/MRlap")
library(MRlap)
library(data.table)
library(writexl)

#setwd

# read in gwas files
brainage <- fread("brainage_factor_nogenr_excl2k_jawinski_noGC_renamed.txt.gz")
head(brainage)

SmokingInitiation <- fread("SmokingInitiation_renamed.txt.gz")
head(SmokingInitiation)

# SmokingInitiation
B = MRlap(exposure = SmokingInitiation,
          exposure_name = "SmokingInitiation",
          outcome = brainage,
          outcome_name = "brainageFactor",
          ld = "eur_w_ld_chr",
          hm3 = "w_hm3.snplist",
          MR_plink = genetics.binaRies::get_plink_binary(),
          MR_bfile = "/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/EUR/EUR",
          MR_threshold = 5e-8,
          MR_pruning_LD = 0.001,
          MR_pruning_dist = 1000)

# save corrected results but exclude rsid names (rids saved in harmonised data)
MRcorrection_df <- as.data.frame(B$MRcorrection[setdiff(names(B$MRcorrection), "IVs")])

# combine into list
df_list <- list(
  MRcorrection = MRcorrection_df,
  LDSC = as.data.frame(B$LDSC),
  GeneticArchitecture = as.data.frame(B$GeneticArchitecture),
  HarmonisedMRData = B$harmonised_mr_data
)

# save to excel
write_xlsx(df_list, path = "MRlap_results_SmokingInitiation.xlsx")

# reverse MRlap for smoking initiation
# SmokingInitiation
B = MRlap(exposure = brainage,
          exposure_name = "brainageFactor",
          outcome = SmokingInitiation,
          outcome_name = "SmokingInitiation",
          ld = "eur_w_ld_chr",
          hm3 = "w_hm3.snplist",
          MR_plink = genetics.binaRies::get_plink_binary(),
          MR_bfile = "/Users/vb506/Documents/EUR/EUR",
          MR_threshold = 5e-8,
          MR_pruning_LD = 0.001,
          MR_pruning_dist = 1000)

# save corrected results but exclude rsid names (rids saved in harmonised data)
MRcorrection_df <- as.data.frame(B$MRcorrection[setdiff(names(B$MRcorrection), "IVs")])

# combine into list
df_list <- list(
  MRcorrection = MRcorrection_df,
  LDSC = as.data.frame(B$LDSC),
  GeneticArchitecture = as.data.frame(B$GeneticArchitecture),
  HarmonisedMRData = B$harmonised_mr_data
)

# save to excel
write_xlsx(df_list, path = "MRlap_results_SmokingInitiation_reverse.xlsx")

## SBP
sbp <- fread("SBP_GCST90310294_renamed.h.tsv.gz")
head(sbp)

# SBP
sbp_results = MRlap(exposure = sbp,
          exposure_name = "SBP",
          outcome = brainage,
          outcome_name = "brainageFactor",
          ld = "eur_w_ld_chr",
          hm3 = "w_hm3.snplist",
          MR_plink = genetics.binaRies::get_plink_binary(),
          MR_bfile = "/Users/vb506/Documents/EUR/EUR",
          MR_threshold = 5e-8,
          MR_pruning_LD = 0.001,
          MR_pruning_dist = 1000)

# save corrected results but exclude rsid names (rids saved in harmonised data)
MRcorrection_df <- as.data.frame(sbp_results$MRcorrection[setdiff(names(sbp_results$MRcorrection), "IVs")])


# DBP
dbp <- fread("DBP_GCST90310295_renamed.h.tsv.gz")
head(dbp)

# DBP
dbp_results = MRlap(exposure = dbp,
                    exposure_name = "DBP",
                    outcome = brainage,
                    outcome_name = "brainageFactor",
                    ld = "eur_w_ld_chr",
                    hm3 = "w_hm3.snplist",
                    MR_plink = genetics.binaRies::get_plink_binary(),
                    MR_bfile = "/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/EUR/EUR",
                    MR_threshold = 5e-8,
                    MR_pruning_LD = 0.001,
                    MR_pruning_dist = 1000)

# save corrected results but exclude rsid names (rids saved in harmonised data)
MRcorrection_df <- as.data.frame(dbp_results$MRcorrection[setdiff(names(dbp_results$MRcorrection), "IVs")])

# combine into list
df_list <- list(
  MRcorrection = MRcorrection_df,
  LDSC = as.data.frame(dbp_results$LDSC),
  GeneticArchitecture = as.data.frame(dbp_results$GeneticArchitecture),
  HarmonisedMRData = dbp_results$harmonised_mr_data
)

# save to excel
write_xlsx(df_list, path = "MRlap_results_DBP.xlsx")



# Longevity 90th
longevity <- fread("longevity_90th_percentile_rsid_renamed2.txt.gz")
head(longevity)

# capitalise alleles
longevity$alt <- toupper(longevity$alt)
longevity$ref <- toupper(longevity$ref)

longevity_results = MRlap(exposure = longevity,
                    exposure_name = "Longevity_90th",
                    outcome = brainage,
                    outcome_name = "brainageFactor",
                    ld = "eur_w_ld_chr",
                    hm3 = "w_hm3.snplist",
                    MR_plink = genetics.binaRies::get_plink_binary(),
                    MR_bfile = "/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/EUR/EUR",
                    MR_threshold = 5e-6,
                    MR_pruning_LD = 0.001, # with 0.001 no IV left after thresholding
                    MR_pruning_dist = 1000,
                    verbose = TRUE)

# save corrected results but exclude rsid names (rids saved in harmonised data)
MRcorrection_df <- as.data.frame(longevity_results$MRcorrection[setdiff(names(longevity_results$MRcorrection), "IVs")])

# combine into list
df_list <- list(
  MRcorrection = MRcorrection_df,
  LDSC = as.data.frame(longevity_results$LDSC),
  GeneticArchitecture = as.data.frame(longevity_results$GeneticArchitecture),
  HarmonisedMRData = longevity_results$harmonised_mr_data
)

# save to excel
write_xlsx(df_list, path = "MRlap_results_longevity.xlsx")




### REVERSE MRLap MR (BAG on smoking initiation)

# read in gwas files
brainage <- fread("brainage_factor_nogenr_excl2k_jawinski_noGC_renamed.txt.gz")
head(brainage)

SmokingInitiation <- fread("SmokingInitiation_renamed.txt.gz")
head(SmokingInitiation)

# SmokingInitiation
B = MRlap(exposure = brainage,
          exposure_name = "brainageFactor",
          outcome = SmokingInitiation,
          outcome_name = "SmokingInitiation",
          ld = "eur_w_ld_chr",
          hm3 = "w_hm3.snplist",
          MR_plink = genetics.binaRies::get_plink_binary(),
          MR_bfile = "/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/EUR/EUR",
          MR_threshold = 5e-8,
          MR_pruning_LD = 0.001,
          MR_pruning_dist = 1000)

# save corrected results but exclude rsid names (rids saved in harmonised data)
MRcorrection_df <- as.data.frame(B$MRcorrection[setdiff(names(B$MRcorrection), "IVs")])

# combine into list
df_list <- list(
  MRcorrection = MRcorrection_df,
  LDSC = as.data.frame(B$LDSC),
  GeneticArchitecture = as.data.frame(B$GeneticArchitecture),
  HarmonisedMRData = B$harmonised_mr_data
)

# save to excel
write_xlsx(df_list, path = "MRlap_results_BAG_on_SmokingInitiation.xlsx")



# ASD
asd <- fread("sumstats-check/iPSYCH-PGC_ASD_Nov2017_logbeta.gz")
asd$N <- 46351 # 18,382 cases, 27,969 controls as detailed here https://ipsych.dk/fileadmin/iPSYCH/Dokumenter/iPSYCH-PGC_ASD_Nov2017_readme.pdf 
setnames(asd, c("OR", "logbeta"), c("ignore", "beta"))
head(asd)
asd <- as.data.frame(asd)

# # convert OR to beta before func as there is an error with OR, i.e., in MRlap:::tidy_inputGWAS function below
# GWASData_clean <- dplyr::mutate(Z = log(.data$beta)/.data$se, 
#                                 beta = NULL, se = NULL)
# # "GWASData_clean %>%" is missing, should be fixed to below:  
# GWASData_clean <- GWASData_clean %>% dplyr::mutate(Z = log(.data$beta)/.data$se, 
#                                 beta = NULL, se = NULL)

# ASD reverse
asd_results = MRlap(exposure = brainage,
                    exposure_name = "brainageFactor",
                    outcome = asd,
                    outcome_name = "ASD",
                    ld = "eur_w_ld_chr",
                    hm3 = "w_hm3.snplist",
                    MR_plink = genetics.binaRies::get_plink_binary(),
                    MR_bfile = "/Users/vb506/Documents/EUR/EUR",
                    MR_threshold = 5e-8,
                    MR_pruning_LD = 0.001,
                    MR_pruning_dist = 1000)

# save corrected results but exclude rsid names (rids saved in harmonised data)
MRcorrection_df <- as.data.frame(asd_results$MRcorrection[setdiff(names(asd_results$MRcorrection), "IVs")])


# combine into list
df_list <- list(
  MRcorrection = MRcorrection_df,
  LDSC = as.data.frame(asd_results$LDSC),
  GeneticArchitecture = as.data.frame(asd_results$GeneticArchitecture),
  HarmonisedMRData = asd_results$harmonised_mr_data
)

# save to excel
write_xlsx(df_list, path = "MRlap_results_ASD_reverse.xlsx")

# Sleep Duration
sleep <- fread("sumstats-check/sleepdurationsumstats.txt")
sleep$N <- 446118 # from https://pmc.ncbi.nlm.nih.gov/articles/PMC6405943/
sleep <- as.data.frame(sleep)
write.table(sleep, file = gzfile("sleepdurationsumstats.txt.gz"), sep = "\t", row.names = FALSE, quote = FALSE)
setnames(sleep, c('ALLELE1', 'ALLELE0', 'BETA_SLEEPDURATION', 'SE_SLEEPDURATION', 'P_SLEEPDURATION'), c("A1", "A2", "beta", "se","p"))
head(sleep)

# sleep forward
sleep_results = MRlap(exposure = sleep,
                    exposure_name = "sleep_duration",
                    outcome = brainage,
                    outcome_name = "brainageFactor",
                    ld = "eur_w_ld_chr",
                    hm3 = "w_hm3.snplist",
                    MR_plink = genetics.binaRies::get_plink_binary(),
                    MR_bfile = "/Users/vb506/Documents/EUR/EUR",
                    MR_threshold = 5e-8,
                    MR_pruning_LD = 0.001,
                    MR_pruning_dist = 1000)

# save corrected results but exclude rsid names (rids saved in harmonised data)
MRcorrection_df <- as.data.frame(sleep_results$MRcorrection[setdiff(names(sleep_results$MRcorrection), "IVs")])


# combine into list
df_list <- list(
  MRcorrection = MRcorrection_df,
  LDSC = as.data.frame(sleep_results$LDSC),
  GeneticArchitecture = as.data.frame(sleep_results$GeneticArchitecture),
  HarmonisedMRData = sleep_results$harmonised_mr_data
)

# save to excel
write_xlsx(df_list, path = "MRlap_results_sleep_forward.xlsx")

# sleep reverse
sleep <- fread("sumstats-check/sleepdurationsumstats.txt")
sleep$N <- 446118 # from https://pmc.ncbi.nlm.nih.gov/articles/PMC6405943/
sleep <- as.data.frame(sleep)
setnames(sleep, c('ALLELE1', 'ALLELE0', 'BETA_SLEEPDURATION', 'SE_SLEEPDURATION', 'P_SLEEPDURATION'), c("A1", "A2", "beta", "se","p"))
head(sleep)

# sleep reverse
sleep_results2 = MRlap(exposure = brainage,
                      exposure_name = "brainageFactor",
                      outcome = sleep,
                      outcome_name = "sleep_duration",
                      ld = "eur_w_ld_chr",
                      hm3 = "w_hm3.snplist",
                      MR_plink = genetics.binaRies::get_plink_binary(),
                      MR_bfile = "/Users/vb506/Documents/EUR/EUR",
                      MR_threshold = 5e-8,
                      MR_pruning_LD = 0.001,
                      MR_pruning_dist = 1000)

# save corrected results but exclude rsid names (rids saved in harmonised data)
MRcorrection_df <- as.data.frame(sleep_results2$MRcorrection[setdiff(names(sleep_results2$MRcorrection), "IVs")])


# combine into list
df_list <- list(
  MRcorrection = MRcorrection_df,
  LDSC = as.data.frame(sleep_results2$LDSC),
  GeneticArchitecture = as.data.frame(sleep_results2$GeneticArchitecture),
  HarmonisedMRData = sleep_results2$harmonised_mr_data
)

# save to excel
write_xlsx(df_list, path = "MRlap_results_sleep_reverse.xlsx")


# PD
# PD =list(outcome_file="GP2_PD_with_Kaviar_IDs_unique_N.gz",
#          outcome_name="PD", snp_col="ID", beta_col="beta", se_col="standard_error",
#          ea_col="effect_allele", nea_col="other_allele", pval_col="p_value",
#          eaf_col="effect_allele_frequency",
#          output_prefix="PD")

pd <- fread("GP2_PD_with_Kaviar_IDs_unique_N.gz")
pd$N <- sum(pd$ncases, pd$ncontrols) 
setnames(pd, c("ID", "standard_error", "effect_allele", "other_allele", "effect_allele_frequency"), c("SNP", "se", "A1", "A2", "EAF"))
head(pd)
# pd <- as.data.frame(pd)

# list 15 SNPs that are kept in Steiger analysis for PD (i.e., 'rs62065444' is dropped)
IVs_15_steiger <- c("rs10494988", "rs10753232","rs12146713", "rs12263364", "rs185726277", "rs28364628",  
  "rs28520337", "rs337637","rs429358", "rs5743091", "rs6442411",
  "rs674243", "rs765724","rs7704770", "rs8067545")
# select these 15

# # # for PD remove 1 SNP that is dropped in Steiger analysis to compare results
# pd=pd[!(pd$SNP %in% 'rs62065444'),]

# PD reverse
pd_results = MRlap(exposure = brainage,
                    exposure_name = "brainageFactor",
                    outcome = pd,
                    outcome_name = "PD",
                    ld = "eur_w_ld_chr",
                    hm3 = "w_hm3.snplist",
                    MR_plink = genetics.binaRies::get_plink_binary(),
                    MR_bfile = "/Users/vb506/Documents/EUR/EUR",
                    MR_threshold = 5e-8,
                    MR_pruning_LD = 0.001,
                    MR_pruning_dist = 1000,
                    user_SNPsToKeep = IVs_15_steiger) # adding IVs_15_steiger here did not work, so had to modify the source function to run MR on the correct 15 SNPs

# save corrected results but exclude rsid names (rids saved in harmonised data)
MRcorrection_df <- as.data.frame(pd_results$MRcorrection[setdiff(names(pd_results$MRcorrection), "IVs")])

# combine into list
df_list <- list(
  MRcorrection = MRcorrection_df,
  LDSC = as.data.frame(pd_results$LDSC),
  GeneticArchitecture = as.data.frame(pd_results$GeneticArchitecture),
  HarmonisedMRData = pd_results$harmonised_mr_data
)

# save to excel
write_xlsx(df_list, path = "MRlap_results_PD_reverse.xlsx")


# import pandas as pd
# import os
# from pathlib import Path
# 
# # Set the path to your directory
# data_dir = Path("/campaign/VB-FM5HPC-001/Vilte/Documents/MRlap-results")
#
# # List all Excel files (.xls or .xlsx)
# file_list = list(data_dir.glob("*.xls*"))
# 
# df_combined = pd.concat(
#   [pd.read_excel(file).assign(source_file=file.name) for file in file_list],
#   ignore_index=True
# )
# 
# # Save the combined DataFrame to a new Excel file
# output_path = data_dir / "combined_data.xlsx"
# df_combined.to_excel(output_path, index=False)
# 




