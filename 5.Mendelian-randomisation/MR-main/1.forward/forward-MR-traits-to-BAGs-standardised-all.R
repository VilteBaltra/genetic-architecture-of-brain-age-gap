# ================================
# Setup
# ================================
.libPaths("/campaign/VB-FM5HPC-001/Vilte/R/x86_64-pc-linux-gnu-library/4.3")
setwd("/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/MR-analysis")

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
bfile <- "/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/EUR/EUR"

# create an output folder if it does not exist
if (!dir.exists("MR-output-forward-std")) dir.create("MR-output-forward-std")


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
  dir.create("MR-output-forward-std", showWarnings = FALSE)
  pdf(paste0("MR-output-forward-std/", exposure_name, "_to_", outcome_info$outcome_name, "_MR_plots.pdf"), width=8, height=8)
  print(p1[[1]]); print(p2[[1]]); print(p3[[1]]); print(p4[[1]])
  dev.off()
  
  ggsave(paste0("MR-output-forward-std/scatter_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p1[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-forward-std/forest_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p2[[1]], width=6, height=8, dpi=300)
  ggsave(paste0("MR-output-forward-std/leaveoneout_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p3[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-forward-std/funnel_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p4[[1]], width=6, height=6, dpi=300)
  
  # Save results
  output_file <- paste0("MR-output-forward-std/MR_results_", exposure_name, "_to_", outcome_info$outcome_name, "_", Sys.Date(), ".xlsx")
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
  # brainageFactor = list(
  #   outcome_file = "brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz", #brainage_factor_withSmith_withN.txt.gz",
  #   outcome_name = "brainageFactor",
  #   sep = "\t",  
  #   snp_col = "SNP", beta_col = "est", se_col = "se_c",
  #   ea_col = "A1", nea_col = "A2", pval_col = "Pval_Estimate",
  #   output_prefix = "brainageFactor"
  # ),
  brainageFactor = list(
    outcome_file = "brainage_factor_nogenr_excl2k_jawinski_noGC_withEAF.txt.gz", # brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz", #brainage_factor_withSmith_withN.txt.gz",
    outcome_name = "brainageFactor",
    sep = "\t",  
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    output_prefix = "brainageFactor"
  ),
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


# EXPOSURES 
exposure_list <- list(
  Metabolic_BAG = list(
    file = "Metabolic_BAG_standardised.txt.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P"
  ),
  Longevity_90th = list(
    file = "Longevity_90th_standardised.txt.gz",
    sep = "\t",
    snp_col = "SNP", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "EA", other_allele_col = "NEA", pval_col = "P-value",
    eaf_col = "EAF", chr_col = "Chr", pos_col = "Position"
  ),
  T2D = list(
    file = "Mahajan.NatGenet2018b.T2D.European_mapped_CHR_ALL.txt.gz",
    sep = "\t",
    snp_col = "RSID", beta_col = "Beta", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P-value",
    eaf_col = "EAF", chr_col = "Chr", pos_col = "Pos"
  ),
  DBP = list(
    file = "DBP_standardised.txt.gz",
    sep = "\t",
    snp_col = "rsid", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    eaf_col = "effect_allele_frequency", chr_col = "chromosome", pos_col = "base_pair_location"
  ),
  SBP = list(
    file = "SBP_standardised.txt.gz",
    sep = "\t",
    snp_col = "rsid", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    eaf_col = "effect_allele_frequency", chr_col = "chromosome", pos_col = "base_pair_location"
  ),
  SmokingInitiation = list(
    file = "sumstats2/SmokingInitiation.txt.gz",
    sep = "\t",
    snp_col = "RSID", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "PVALUE",
    eaf_col = "AF", chr_col = "CHROM", pos_col = "POS"
  ),
  CigarettesPerDay = list(
    file = "CigarettesPerDay_standardised.txt.gz",
    sep = "\t",
    snp_col = "RSID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "PVALUE",
    eaf_col = "AF", chr_col = "CHROM", pos_col = "POS"
  ),
  ADHD = list(
    file = "ADHD2022_iPSYCH_deCODE_PGC.meta.logBeta.gz",
    sep = "\t",
    snp_col = "SNP", beta_col = "logBeta", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    eaf_col = "FRQ_A_38691", chr_col = "CHR", pos_col = "BP"
  ),
  Bipolar = list(
    file = "daner_bip_pgc3_nm_noukbiobank_beta.gz",
    sep = "\t",
    snp_col = "SNP", beta_col = "logBeta", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    eaf_col = "FRQ_A_40463", chr_col = "CHR", pos_col = "BP"
  ),
  CRP = list(
    file = "CRP_standardised.txt.gz",
    sep = "\t",
    snp_col = "variant_id", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "p_value",
    chr_col = "chromosome", pos_col = "base_pair_location"
  )
)

# EPI CLOCK EXPOSURES 
exposure_list2 <- list(
  Gran = list(
    file = "Gran_EUR_standardised.txt.gz",
    sep = "\t",
    snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P"
  ),
  Hannum = list(
    file = "Hannum_EUR_standardised.txt.gz",
    sep = "\t",
    snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P"
  )
)

# ADDITIONAL EXPOSURES REQUESTED BY CO-AUTHORS
exposure_list <- list(
  PD_BAG = list(
    file = "GP2_PD_with_Kaviar_IDs_unique.txt.gz", # "GP2_ALL_EUR_ALL_DATASET_HG38_12162024_rsid.txt.gz", # "TRACTOR_EUR_Meta-analysis_random_LARGE-PD_P1_P2_hg38.txt.gz",  # "GP2_EUR_ONLY_HG38_12162024.txt.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "beta", se_col = "standard_error",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value"
  ))
# exposure_list <- list(
#   PD_BAG = list(
#     file =  "GP2_EUR_ONLY_HG38_12162024.txt.gz",
#     sep = "\t",
#     snp_col = "SNP_ID", beta_col = "beta", se_col = "standard_error",
#     effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value"
#   ))
# # chromosome	base_pair_position	SNP_ID	effect_allele	other_allele	effect_allele_frequency	N_datasets	p_value	beta	standard_error	p_value(random)	beta(random)

# -----------------------------
# Run MR analyses (clump_p=5e-8)
# -----------------------------
for (exposure_name in names(exposure_list)[names(exposure_list) != "Longevity_90th"][9]) {
  exposure_info <- exposure_list[[exposure_name]]
  
  exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-8, filter_pval=5e-6)
  
  for (outcome_name in names(outcome_info)) {
    run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]])
  }
}

# -----------------------------
# Run MR analyses (clump_p=5e-6)
# -----------------------------
# longevity 90th run separately using clump_p=5e-6 as only ~1 SNP after clumping
for (exposure_name in "Longevity_90th") {
  exposure_info <- exposure_list[[exposure_name]]
  
  exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-6, filter_pval=5e-6)
  
  for (outcome_name in names(outcome_info)) {
    run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]])
  }
}

# -----------------------------
# Run MR analyses (epi clocks)
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
# Run MR analyses (additional, PD)
# -----------------------------
for (exposure_name in names(exposure_list)) {
  exposure_info <- exposure_list[[exposure_name]]
  
  exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-8, filter_pval=5e-6)
  
  for (outcome_name in "brainageFactor") {
    run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]])
  }
}

# -----------------------------
# Combine output to spreadsheet
# -----------------------------

# get list of all your Excel files in the directory
setwd("MR-output-forward-std")
file_list <- list.files(pattern = "\\.xlsx$")

# read first sheet of each file and store in a list
data_list <- lapply(file_list, function(file) { 
  read_excel(file, sheet = 1) %>% mutate(SourceFile = file) # to keep track of source
})

# combine all into one data frame
combined_data <- bind_rows(data_list)

# view result
print(combined_data)

# save results
write.xlsx(combined_data, file = paste0("MR_results_trait_to_BAG.std_", Sys.Date(), ".xlsx"), rowNames = FALSE, overwrite = TRUE)

# # tidy output in terminal
# cd MR-output-forward-std
# mkdir individual-plots combined-plots forward-MR-results
# mv *.png individual-plots
# mv *.pdf combined-plots
# mkdir forward-MR-results/combined
# mv MR_results_combined_*.xlsx forward-MR-results/combined
# mkdir forward-MR-results/per-trait
# mv MR_results*.xlsx forward-MR-results/per-trait
