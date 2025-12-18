# setwd

############## add EAF ###################
data_name = "Immune_age_gap_PLINK_EUR.glm.linear"
data_name = "Musculoskeletal_age_gap_PLINK_EUR.glm.linear"
data_name = "Pulmonary_age_gap_PLINK_EUR.glm.linear"
data_name = "Renal_age_gap_PLINK_EUR.glm.linear"
data_name = "Cardiovascular_age_gap_PLINK_EUR.glm.linear"
  
# read in ref to get chr and pos
ref <- fread("../reference.1000G.maf.0.005.txt", sep = " ") # which uses build GRCh37 as doubled checked via dbSNP here https://www.ncbi.nlm.nih.gov/snp/

# read in sumstats that need eaf added
sumstats <- fread(paste0("sumstats-check/9-PhenoBAG-NA-NM/", data_name))
sumstats <- sumstats %>% rename(A1other = A1)
sumstats <- sumstats %>% rename(A1 = ALT, A2 = REF)

# First merge as above
sumstats_eaf <- merge(sumstats, ref[, c("SNP", "A1", "A2", "MAF")], by.x="ID", by.y = "SNP", all.x = TRUE)
#sumstats_merged <- merge(sumstats, ref[, c("SNP", "A1", "A2", "MAF")], by.x="ID", by.y = "SNP")

# Rename columns for clarity
names(sumstats_eaf)[names(sumstats_eaf) == "A1.x"] <- "A1"
names(sumstats_eaf)[names(sumstats_eaf) == "A2.x"] <- "A2"
names(sumstats_eaf)[names(sumstats_eaf) == "A1.y"] <- "A1_ref"
names(sumstats_eaf)[names(sumstats_eaf) == "A2.y"] <- "A2_ref"

# derive EAF based on alignment
sumstats_eaf$EAF <- with(sumstats_eaf, ifelse(A1 == A1_ref, MAF,
                                              ifelse(A1 == A2_ref, 1 - MAF, NA)))

# save
write.table(sumstats_eaf, file = gzfile(paste0(data_name, "_with_eaf.gz")), sep = "\t", row.names = FALSE, quote = FALSE)
rm(sumstats,sumstats_eaf)


############## add EAF to ASD ###################
sumstats <- fread("sumstats-check/iPSYCH-PGC_ASD_Nov2017_logbeta.gz")
data_name="iPSYCH-PGC_ASD_Nov2017_logbeta"
# First merge as above
sumstats_eaf <- merge(sumstats, ref[, c("SNP", "A1", "A2", "MAF")], by = "SNP", all.x = TRUE)

# Rename columns for clarity
names(sumstats_eaf)[names(sumstats_eaf) == "A1.x"] <- "A1"
names(sumstats_eaf)[names(sumstats_eaf) == "A2.x"] <- "A2"
names(sumstats_eaf)[names(sumstats_eaf) == "A1.y"] <- "A1_ref"
names(sumstats_eaf)[names(sumstats_eaf) == "A2.y"] <- "A2_ref"

# derive EAF based on alignment
sumstats_eaf$EAF <- with(sumstats_eaf, ifelse(A1 == A1_ref, MAF,
                                              ifelse(A1 == A2_ref, 1 - MAF, NA)))
# add N 
# 18,382 cases, 27,969 controls as detailed here https://ipsych.dk/fileadmin/iPSYCH/Dokumenter/iPSYCH-PGC_ASD_Nov2017_readme.pdf 
sumstats_eaf$ncase <- 18382
sumstats_eaf$ncontrol <- 27969
head(sumstats_eaf)

# save
write.table(sumstats_eaf, file = gzfile(paste0(data_name, "_with_eaf.gz")), sep = "\t", row.names = FALSE, quote = FALSE)
rm(sumstats,sumstats_eaf)

############## add N to PD ###################
df <- fread("GP2_PD_with_Kaviar_IDs_unique.txt.gz")
df$ncases <- 81255  # from https://www.medrxiv.org/content/10.1101/2025.03.14.24319455v1
df$ncontrols <- 1746386 
write.table(df, file = gzfile("GP2_PD_with_Kaviar_IDs_unique_N.gz"), sep = "\t", row.names = FALSE, quote = FALSE)

############## STANDARDISATION OF GWASs ###################
library(data.table)

# -----------------------------
# Standardise all datasets
# -----------------------------

standardize_beta_se <- function(df, beta_col, se_col, eaf_col, n_col) {
  required_cols <- c(beta_col, se_col, eaf_col, n_col)
  if (!all(required_cols %in% names(df))) {
    stop("Missing required columns: ", paste(setdiff(required_cols, names(df)), collapse = ", "))
  }
  
  A1_FREQ <- df[[eaf_col]]
  A1_FREQ[A1_FREQ >= 1] <- max(A1_FREQ[A1_FREQ < 1], na.rm = TRUE)
  A1_FREQ[A1_FREQ <= 0] <- min(A1_FREQ[A1_FREQ > 0], na.rm = TRUE)
  
  Z <- df[[beta_col]] / df[[se_col]]
  denom <- sqrt(2 * A1_FREQ * (1 - A1_FREQ) * (df[[n_col]] + Z^2))
  
  df$beta.std <- Z / denom
  df$se.std   <- 1 / denom
  
  return(df)
}


### traits to standardise ###
outcome_info <- list(
  IEAA = list(
    outcome_file = "sumstats-check/IEAA_EUR_summary_statistics.txt",  # compress if needed
    sep = "\t",
    beta_col = c("Effect"),
    se_col = c("SE"),
    eaf_col = c("Freq1"),
    samplesize_col = c("N"),
    output_prefix = "IEAA_EUR"
  ),
  GrimAge = list(
    outcome_file = "sumstats-check/GrimAge_EUR_summary_statistics.txt",
    sep = "\t",
    beta_col = c("Effect"),
    se_col = c("SE"),
    eaf_col = c("Freq1"),
    samplesize_col = c("N"),
    output_prefix = "GrimAge_EUR"
  )
  ,
  PAI1 = list(
    outcome_file = "sumstats-check/PAI1_EUR_summary_statistics.txt",
    sep = "\t",
    beta_col = c("Effect"),
    se_col = c("SE"),
    eaf_col = c("Freq1"),
    samplesize_col = c("N"),
    output_prefix = "PAI1_EUR"
  )
  ,
  PhenoAge = list(
    outcome_file = "sumstats-check/PhenoAge_EUR_summary_statistics.txt",
    sep = "\t",
    beta_col = c("Effect"),
    se_col = c("SE"),
    eaf_col = c("Freq1"),
    samplesize_col = c("N"),
    output_prefix = "PhenoAge_EUR"
  ),
  BMI = list(
    outcome_file = "sumstats-check/bmi-25673413-GCST002783-EFO_0004340.h.tsv.gz",
    sep = "\t",
    beta_col = c("beta"),
    se_col = c("standard_error"),
    eaf_col = c("eaf_ref"),
    samplesize_col = c("n"),
    output_prefix = "BMI"
  ),
  Cardiovascular_BAG = list(
    outcome_file = "Cardiovascular_age_gap_PLINK_EUR.glm.linear_with_eaf.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "BETA", se_col = "SE",
    eaf_col = "EAF",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    samplesize_col = "OBS_CT",
    output_prefix = "Cardiovascular_BAG"
  ),
  Immune_BAG = list(
    outcome_file = "Immune_age_gap_PLINK_EUR.glm.linear_with_eaf.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "BETA", se_col = "SE",
    eaf_col = "EAF",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    samplesize_col = "OBS_CT",
    output_prefix = "Immune_BAG"
  ),
  Renal_BAG = list(
    outcome_file = "Renal_age_gap_PLINK_EUR.glm.linear_with_eaf.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "BETA", se_col = "SE",
    eaf_col = "EAF",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    samplesize_col = "OBS_CT",
    output_prefix = "Renal_BAG"
  ),
  Pulmonary_BAG = list(
    outcome_file = "Pulmonary_age_gap_PLINK_EUR.glm.linear_with_eaf.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "BETA", se_col = "SE",
    eaf_col = "EAF",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    samplesize_col = "OBS_CT",
    output_prefix = "Pulmonary_BAG"
  ),
  Musculoskeletal_BAG = list(
    outcome_file = "Musculoskeletal_age_gap_PLINK_EUR.glm.linear_with_eaf.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "BETA", se_col = "SE",
    eaf_col = "EAF",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    samplesize_col = "OBS_CT",
    output_prefix = "Musculoskeletal_BAG"
  ),
  PP = list(
    outcome_file = "sumstats-check/PP_GCST90310296.h.tsv.gz",
    sep = "\t",
    snp_col = "rsid", beta_col = "beta", se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    ea_col = "effect_allele", nea_col = "other_allele", pval_col = "p_value",
    chr_col = "chromosome", pos_col = "base_pair_location",
    samplesize_col = "n",
    output_prefix = "PP"),
  RHR = list(
    outcome_file = "sumstats-check/RHR_UKBplusICRHR_without23andme.gz",
    sep = "\t",
    snp_col = "SNP_UBK", beta_col = "Effect", se_col = "StdErr",
    eaf_col = "A1Freq",
    ea_col = "A1", nea_col = "A2", pval_col = "P-value",
    chr_col = "CHR", pos_col = "BP",
    samplesize_col = "TotalSampleSize",
    output_prefix = "RHR"),
  HbA1c = list(
    outcome_file = "sumstats-check/MAGIC1000G_HbA1c_EUR.tsv.gz",
    sep = "\t",
    snp_col = "variant", beta_col = "beta", se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    ea_col = "effect_allele", nea_col = "other_allele", pval_col = "p_value",
    chr_col = "chromosome", pos_col = "base_pair_location",
    samplesize_col = "sample_size",
    output_prefix = "HbA1c")
)

for (name in names(outcome_info)) {
  info <- outcome_info[[name]]
  
  message("Processing: ", name)
  
  file_path <- info$outcome_file
  sep <- info$sep
  
  df <- fread(file_path, sep = sep, data.table = FALSE)
  
  # Helper function: select first column name that exists in df
  select_col <- function(candidates, df_cols) {
    if (is.null(candidates)) return(NULL)
    match <- candidates[candidates %in% df_cols]
    if (length(match) > 0) return(match[1]) else return(NULL)
  }
  
  # Allow multiple possible names
  beta_col <- select_col(info$beta_col, names(df))
  se_col   <- select_col(info$se_col, names(df))
  eaf_col  <- select_col(info$eaf_col, names(df))
  
  # Check if sample size is constant or a column
  n_col <- NULL
  if (!is.null(info$samplesize_col)) {
    if (is.numeric(info$samplesize_col) && length(info$samplesize_col) == 1) {
      df$N <- info$samplesize_col  # add constant column
      n_col <- "N"
    } else {
      n_col <- select_col(info$samplesize_col, names(df))
    }
  }
  
  # Standardize if all needed cols exist
  if (!is.null(beta_col) && !is.null(se_col) &&
      !is.null(eaf_col) && !is.null(n_col)) {
    
    df <- standardize_beta_se(df,
                              beta_col = beta_col,
                              se_col = se_col,
                              eaf_col = eaf_col,
                              n_col = n_col)
    
    message("Standardized: ", name)
    
  } else {
    message("Skipped standardization for: ", name, " (missing beta/se/eaf/n)")
  }
  
  # view
  print(head(df))
  
  # save standardised data
  outfile <- paste0(info$output_prefix, "_standardised.txt.gz")
  message("Saving: ", outfile)
  write.table(df, file = gzfile(outfile), sep = "\t", row.names = FALSE, quote = FALSE)
  
}



############### MR PART #########################
#setwd("/Users/vb506/Documents/Documents-Nav-PC/Vilte-2025-NaveenPC/2025/3.MR/scripts/input-GWASs/scripts/source/")
setwd("/Users/vb506/Documents/MR-analysis-local/")

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(genetics.binaRies)
library(ggplot2)
library(openxlsx)
library(ieugwasr)
library(readxl)

source("my_MR_tests.R")
plink_bin <- genetics.binaRies::get_plink_binary()
bfile <- "/Users/vb506/Documents/EUR/EUR"


# create an output folder if it does not exist
if (!dir.exists("MR-output-forward-std-check")) dir.create("MR-output-forward-std-check")


# -----------------------------
# Helper to format exposure data
# -----------------------------

prepare_exposure_data <- function(exposure_info, exposure_name, clump_p=5e-8, filter_pval=1) {
  cat("\n\n🟢 Preparing exposure:", exposure_name, "\n"); flush.console()
  
  gwas <- fread(exposure_info$file, sep = exposure_info$sep) |> as.data.frame()
  gwas[[exposure_info$pval_col]] <- as.numeric(gwas[[exposure_info$pval_col]])
  if (!is.null(filter_pval)) gwas <- gwas[gwas[[exposure_info$pval_col]] < filter_pval, ]
  
  format_args <- list(
    gwas, type="exposure",
    snp_col = exposure_info$snp_col,
    beta_col = exposure_info$beta_col,
    se_col = exposure_info$se_col,
    effect_allele_col = exposure_info$effect_allele_col,
    other_allele_col = exposure_info$other_allele_col,
    pval_col = exposure_info$pval_col,
    min_pval = 1e-200
  )
  if (!is.null(exposure_info$eaf_col)) format_args$eaf_col <- exposure_info$eaf_col
  if (!is.null(exposure_info$chr_col)) format_args$chr_col <- exposure_info$chr_col
  if (!is.null(exposure_info$pos_col)) format_args$pos_col <- exposure_info$pos_col
  if (!is.null(exposure_info$samplesize_col)) format_args$samplesize_col <- exposure_info$samplesize_col
  
  exp_dat <- do.call(format_data, format_args)
  exp_dat$exposure <- exposure_name
  
  # Clump
  clumped_snps <- ld_clump(
    tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id=exp_dat$id.exposure),
    plink_bin = plink_bin,
    bfile = bfile,
    clump_kb = 10000, clump_r2 = 0.001,
    clump_p = clump_p, pop = "EUR"
  )
  exp_dat_clumped <- subset(exp_dat, SNP %in% clumped_snps$rsid)
  cat("  ✔ Clumped SNPs:", nrow(exp_dat_clumped), "\n"); flush.console()
  
  return(exp_dat_clumped)
}


# -----------------------------
# Function to run MR for an outcome
# -----------------------------

run_MR_for_outcome <- function(exp_dat_clumped, exposure_name, outcome_info) {
  
  cat("\n\n🔵 Running MR:", exposure_name, " -> ", outcome_info$outcome_name, "\n"); flush.console()
  
  # Load outcome data
  outcome_data <- read_outcome_data(
    snps = exp_dat_clumped$SNP, 
    filename = outcome_info$outcome_file,
    sep = outcome_info$sep, #"\t",
    snp_col = outcome_info$snp_col,
    beta_col = outcome_info$beta_col, 
    se_col = outcome_info$se_col,
    effect_allele_col = outcome_info$ea_col, 
    other_allele_col = outcome_info$nea_col,
    pval_col = outcome_info$pval_col,
    eaf_col = if (!is.null(outcome_info$eaf_col)) outcome_info$eaf_col else NULL
  )
  outcome_data$outcome <- outcome_info$outcome_name
  
  # Harmonise
  dat <- harmonise_data(exp_dat_clumped, outcome_data)
  
  # Run MR
  res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
  print(res)
  final_results <- my_MR_tests(res, dat)
  res_single <- mr_singlesnp(dat)
  res_loo <- mr_leaveoneout(dat)
  
  # Plots
  p1 <- mr_scatter_plot(res, dat)
  p2 <- mr_forest_plot(res_single)
  p3 <- mr_leaveoneout_plot(res_loo)
  p4 <- mr_funnel_plot(res_single)
  
  # Save plots
  dir.create("MR-output-forward-std-check", showWarnings = FALSE)
  pdf(paste0("MR-output-forward-std-check/", exposure_name, "_to_", outcome_info$outcome_name, "_MR_plots.pdf"), width=8, height=8)
  print(p1[[1]]); print(p2[[1]]); print(p3[[1]]); print(p4[[1]])
  dev.off()
  
  ggsave(paste0("MR-output-forward-std-check/scatter_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p1[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-forward-std-check/forest_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p2[[1]], width=6, height=8, dpi=300)
  ggsave(paste0("MR-output-forward-std-check/leaveoneout_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p3[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-forward-std-check/funnel_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p4[[1]], width=6, height=6, dpi=300)
  
  # Save results
  output_file <- paste0("MR-output-forward-std-check/MR_results_", exposure_name, "_to_", outcome_info$outcome_name, "_", Sys.Date(), ".xlsx")
  write.xlsx(
    list(main_MR_results = final_results,
         single_SNP_results = res_single,
         leave_one_out_results = res_loo),
    file = output_file,
    rowNames = FALSE, overwrite = TRUE)
  
  cat("✅ Done:", exposure_name, " -> ", outcome_info$outcome_name, "\n")
  cat("📁 Results saved to:", output_file, "\n"); flush.console()
  
}


# -----------------------------
# Define outcome datasets 
# -----------------------------

# # first convert zscore to beta in ENIGMA GWAS
# gwas_enigma <- fread("METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22.txt.gz")
# gwas_enigma <- as.data.frame(gwas_enigma)
# head(gwas_enigma)
# # get beta and se (as detailed here https://www.biostars.org/p/319584/)
# gwas_enigma <- gwas_enigma %>% mutate(beta = Zscore / sqrt(2 * Freq1 * (1 - Freq1) * (N + Zscore^2)),
#                                       SE   = 1 / sqrt(2 * Freq1 * (1 - Freq1) * (N + Zscore^2)))
# # save
# write.table(gwas_enigma, file = gzfile("METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22_withBeta.txt.gz"), sep = "\t", row.names = FALSE, quote = FALSE)

# OUTCOMES
outcome_info <- list(
  Kaufman = list(
    outcome_file = "Kaufman_standardised.txt.gz",
    outcome_name = "Kaufman",
    sep = "\t",
    snp_col = "SNP", beta_col = "beta.std", se_col = "se.std",
    ea_col = "A1_kaufman", nea_col = "A2_kaufman", pval_col = "PVAL",
    output_prefix = "Kaufman"
  ),
  Smith = list(
    outcome_file = "Smith_standardised.txt.gz",
    outcome_name = "Smith",
    sep = "\t",
    snp_col = "rsid", beta_col = "beta.std", se_col = "se.std",
    ea_col = "a1", nea_col = "a2", pval_col = "P",
    output_prefix = "Smith"
  ),
  brainageFactor = list(
    outcome_file = "brainage_factor_nogenr_excl2k_jawinski_noGC_withEAF.txt.gz", # brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz", #brainage_factor_withSmith_withN.txt.gz",
    outcome_name = "brainageFactor",
    sep = "\t",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    output_prefix = "brainageFactor"
  )
  ,
  Han = list(
    outcome_file = "METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22_withBeta.txt.gz", # was "METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-05-27_withBeta.txt.gz"
    outcome_name = "Han",
    sep = "\t",
    snp_col = "MarkerName", beta_col = "beta", se_col = "SE",
    ea_col = "Allele1", nea_col = "Allele2", pval_col = "P",
    output_prefix = "Han"
  ),
  Jawinski = list(
    outcome_file = "Jawinski_standardised.txt.gz", # brainage2025.full.eur.gwm.gz",
    outcome_name = "Jawinski",
    sep =  "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    output_prefix = "Jawinski"
  ),
  Leonardsen = list(
    outcome_file = "Leonardsen_standardised.txt.gz",
    outcome_name = "Leonardsen",
    sep = "\t",
    snp_col = "SNP", beta_col = "beta.std", se_col = "se.std",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    output_prefix = "Leonardsen"
  ),
  Wen = list(
    outcome_file = "Wen_standardised.txt.gz",
    outcome_name = "Wen",
    sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    ea_col = "A1_wen", nea_col = "A2_wen", pval_col = "P",
    output_prefix = "Wen"
  )
)


# EPI CLOCK EXPOSURES 
exposure_list2 <- list(
  PhenoAge = list(
    file = "PhenoAge_EUR_standardised.txt.gz",
    sep = "\t",
    snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P"
  ),
  PAI1 = list(
    file = "PAI1_EUR_standardised.txt.gz",
    sep = "\t",
    snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P"
  )
  ,
  GrimAge = list(
    file = "GrimAge_EUR_standardised.txt.gz",
    sep = "\t",
    snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P"
  )
  ,
  IEAA = list(
    file = "IEAA_EUR_standardised.txt.gz",
    sep = "\t",
    snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P"
  ),
  ASD = list( # adding here as only 2 genome-wide significant snps
    file = "sumstats-check/iPSYCH-PGC_ASD_Nov2017_logbeta.gz", 
    sep = "\t",
    snp_col = "SNP", beta_col = "logbeta", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P"
    #eaf_col = "MAF" # none 
  )
)

# ADDITIONAL EXPOSURES REQUESTED BY CO-AUTHORS
## binary traits not standardised
exposure_list <- list(
  PD = list(
    file = "GP2_PD_with_Kaviar_IDs_unique.txt.gz", # "GP2_ALL_EUR_ALL_DATASET_HG38_12162024_rsid.txt.gz", # "TRACTOR_EUR_Meta-analysis_random_LARGE-PD_P1_P2_hg38.txt.gz",  # "GP2_EUR_ONLY_HG38_12162024.txt.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "beta", se_col = "standard_error",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    eaf_col = "effect_allele_frequency"
  ),
  Schizophrenia = list(
    file = "sumstats-check/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "PVAL",
    eaf_col = "MAF"
  ),
  Alzheimer = list(
    file = "sumstats-check/alzheimer-35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
    sep = "\t",
    snp_col = "variant_id", beta_col = "beta", se_col = "standard_error",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    eaf_col = "effect_allele_frequency"
  ),
  BMI = list(
    file = "BMI_standardised.txt.gz",
    sep = "\t",
    snp_col = "variant_id", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    eaf_col = "eaf_ref"
  ),
  Cardiovascular_BAG = list(
    file = "Cardiovascular_BAG_standardised.txt.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    eaf_col = "EAF"
  ),
  Immune_BAG = list(
    file = "Immune_BAG_standardised.txt.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    eaf_col = "EAF"
  ),
  Renal_BAG = list(
    file = "Renal_BAG_standardised.txt.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    eaf_col = "EAF"
  ),
  Pulmonary_BAG = list(
    file = "Pulmonary_BAG_standardised.txt.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    eaf_col = "EAF"
  ),
  Musculoskeletal_BAG = list(
    file = "Musculoskeletal_BAG_standardised.txt.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    eaf_col = "EAF"
  ),
  PP = list(
    file = "PP_standardised.txt.gz",
    sep = "\t",
    snp_col = "rsid", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    eaf_col = "effect_allele_frequency", chr_col = "chromosome", pos_col = "base_pair_location"
  ),
  RHR = list(
    file = "RHR_standardised.txt.gz",
    sep = "\t",
    snp_col = "SNP_UBK", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P-value",
    eaf_col = "A1Freq", chr_col = "CHR", pos_col = "BP"
  ),
  HbA1c = list(
    file = "HbA1c_standardised.txt.gz",
    sep = "\t",
    snp_col = "variant", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    eaf_col = "effect_allele_frequency", chr_col = "chromosome", pos_col = "base_pair_location"
  ),
  CAD = list(
    file = "sumstats-check/cad.add.160614.website.txt",
    sep = "\t",
    snp_col = "markername", beta_col = "beta", se_col = "se_dgc",
    effect_allele_col = "effect_allele", other_allele_col = "noneffect_allele", pval_col = "p_dgc",
    eaf_col = "effect_allele_freq", chr_col = "chr", pos_col = "bp_hg19"
  ),
  # MDD = list(
  #   file = "sumstats-check/MDD2018_ex23andMe_logbeta.gz",
  #   sep = "\t",
  #   snp_col = "SNP", beta_col = "logbeta", se_col = "SE",
  #   effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
  #   eaf_col = "FRQ_A_59851", chr_col = "CHR", pos_col = "BP"
  # ),
  MDD = list(
    file = "sumstats-check/pgc-mdd2025_no23andMe_eur_v3-49-24-11_formatted.tsv.gz",
    sep = "\t",
    snp_col = "rsid", beta_col = "beta", se_col = "standard_error",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    eaf_col = "effect_allele_frequency", chr_col = "chromosome", pos_col = "base_pair_location"
  ),
  Sleep_duration = list(
    file = "sumstats-check/sleepdurationsumstats.txt", # seems to be standardised already
    sep = "\t",
    snp_col = "SNP", beta_col = "BETA_SLEEPDURATION", se_col = "SE_SLEEPDURATION",
    effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", pval_col = "P_SLEEPDURATION",
    eaf_col = "A1FREQ", chr_col = "CHR", pos_col = "BP"
  ),
  Loneliness = list(
    file = "sumstats-check/Day_2018_NatComs_Social/MTAG_results.txt.gz", # seems to be standardised already 
    sep = "\t",
    snp_col = "snpid", beta_col = "mtag_beta", se_col = "mtag_se",
    effect_allele_col = "a1", other_allele_col = "a2", pval_col = "mtag_pval",
    eaf_col = "freq", chr_col = "chr", pos_col = "bpos"
  )
)

# -----------------------------
# Run MR analyses (additional, PD, SCZ, ALZ, BMI, Cardiovascular_BAG, Renal_BAG, Pulmonary_BAG, Musculoskeletal_BAG, PP, RHR, HbA1c, CAD
# -----------------------------
for (exposure_name in names(exposure_list)) { 
  exposure_info <- exposure_list[[exposure_name]]
  
  exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-8, filter_pval=5e-6)
  
  for (outcome_name in names(outcome_info)) { 
    run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]])
  }
}


# -----------------------------
# Run MR analyses (epi clocks + ASD)
# -----------------------------
# using clump_p=5e-6 as Gran only has 2 ind. snps and Hannum 2-6 depending on exposure-outcome combination
for (exposure_name in names(exposure_list2)) {
  exposure_info <- exposure_list2[[exposure_name]]
  
  exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-6, filter_pval=5e-6)
  
  for (outcome_name in names(outcome_info)) { 
    run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]])
  }
}


# -----------------------------
# Combine output to spreadsheet
# -----------------------------

# get list of all your Excel files in the directory
setwd("MR-output-forward-std-check")
file_list <- list.files(pattern = "\\.xlsx$")

# read first sheet of each file and store in a list
data_list <- lapply(file_list, function(file) { 
  read_excel(file, sheet = 1) %>% mutate(SourceFile = file) # to keep track of source
})

# combine all into one data frame
combined_data <- bind_rows(data_list)

# view result
print(combined_data)
print(combined_data %>% filter(method == "Inverse variance weighted") %>% filter(pval < 0.05)) # none

# save results
write.xlsx(combined_data, file = paste0("MR_results_remaining_traits_to_BAG.std_", Sys.Date(), ".xlsx"), rowNames = FALSE, overwrite = TRUE)

# # tidy output in terminal
# cd MR-output-forward-std-check
# mkdir individual-plots combined-plots forward-MR-results
# mv *.png individual-plots
# mv *.pdf combined-plots
# mkdir forward-MR-results/combined
# mv MR_results_combined_*.xlsx forward-MR-results/combined
# mkdir forward-MR-results/per-trait
# mv MR_results*.xlsx forward-MR-results/per-trait
