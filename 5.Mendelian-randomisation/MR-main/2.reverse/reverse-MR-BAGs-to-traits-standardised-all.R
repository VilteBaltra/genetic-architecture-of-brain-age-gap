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
if (!dir.exists("MR-output-reverse-std")) dir.create("MR-output-reverse-std")


# -----------------------------
# Helper to format exposure data
# -----------------------------

prepare_exposure_data <- function(exposure_info, exposure_name, clump_p=5e-8, filter_pval=1) {
  cat("\n\n🟢 Preparing exposure:", exposure_name, "\n"); flush.console()
  
  gwas <- fread(exposure_info$file, sep = exposure_info$sep) |> as.data.frame()
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
    sep = "\t",
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
  final_results <- my_MR_tests(res, dat)
  res_single <- mr_singlesnp(dat)
  res_loo <- mr_leaveoneout(dat)
  
  # Plots
  p1 <- mr_scatter_plot(res, dat)
  p2 <- mr_forest_plot(res_single)
  p3 <- mr_leaveoneout_plot(res_loo)
  p4 <- mr_funnel_plot(res_single)
  
  # Save plots
  dir.create("MR-output-reverse-std", showWarnings = FALSE)
  pdf(paste0("MR-output-reverse-std/", exposure_name, "_to_", outcome_info$outcome_name, "_MR_plots.pdf"), width=8, height=8)
  print(p1[[1]]); print(p2[[1]]); print(p3[[1]]); print(p4[[1]])
  dev.off()
  
  ggsave(paste0("MR-output-reverse-std/scatter_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p1[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-reverse-std/forest_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p2[[1]], width=6, height=8, dpi=300)
  ggsave(paste0("MR-output-reverse-std/leaveoneout_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p3[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-reverse-std/funnel_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p4[[1]], width=6, height=6, dpi=300)
  
  # Save results
  output_file <- paste0("MR-output-reverse-std/MR_results_", exposure_name, "_to_", outcome_info$outcome_name, "_", Sys.Date(), ".xlsx")
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
# Define exposure datasets
# -----------------------------

# define exposure list
exposure_list <- list(
  # brainageFactor = list(
  #   file = "brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz", 
  #   sep = "\t",
  #   snp_col = "SNP", beta_col = "est", se_col = "se_c",
  #   effect_allele_col = "A1", other_allele_col = "A2",
  #   pval_col = "Pval_Estimate", samplesize_col = "N"
  # ),
  brainageFactor = list(
    file = "brainage_factor_nogenr_excl2k_jawinski_noGC_withEAF.txt.gz", 
    sep = "\t",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "A1", other_allele_col = "A2",
    pval_col = "P", samplesize_col = "N"
  ),
  Han = list(
    file = "METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22_withBeta.txt.gz",
    sep = "\t",
    snp_col = "MarkerName", beta_col = "beta", se_col = "SE",
    effect_allele_col = "Allele1", other_allele_col = "Allele2",
    pval_col = "P", samplesize_col = "N"
  ),
  Jawinski = list(
    file =  "Jawinski_standardised.txt.gz", 
    sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2",
    pval_col = "P", samplesize_col = "N"
  ),
  Leonardsen = list(
    file = "Leonardsen_standardised.txt.gz",
    sep = "\t",
    snp_col = "SNP", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2",
    pval_col = "P", samplesize_col = "N"
  ),
  Wen = list(
    file = "Wen_standardised.txt.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1_wen", other_allele_col = "A2_wen",
    pval_col = "P", samplesize_col = "OBS_CT"
  )
)

# Kaufman and Smith only have 1 clumped snp, so using clump_p = 5e-6 instead: 
exposure_list2 <- list(
  Kaufman = list(
    file = "Kaufman_standardised.txt.gz",
    sep = "\t",
    snp_col = "SNP",  beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1_kaufman", other_allele_col = "A2_kaufman", pval_col = "PVAL"
  ),
  Smith = list(
    file = "Smith_standardised.txt.gz",
    sep = "\t",
    snp_col = "rsid", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "a1", other_allele_col = "a2",
    pval_col = "P"
  )
)


# -----------------------------
# Define outcome datasets
# -----------------------------
outcome_info <- list(
  # PD
  PD =list(outcome_file="GP2_PD_with_Kaviar_IDs_unique.txt.gz",
                      outcome_name="PD", snp_col="ID", beta_col="beta", se_col="standard_error",
                      ea_col="effect_allele", nea_col="other_allele", pval_col="p_value",
                      output_prefix="PD"),
  # Metabolic-BAG
  Metabolic_BAG =list(outcome_file="Metabolic_BAG_standardised.txt.gz",
                      outcome_name="Metabolic-BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
                      ea_col="A1", nea_col="A2", pval_col="P",
                      output_prefix="Metabolic_BAG"),
  # Longevity-90th
  Longevity_90th = list(outcome_file="Longevity_90th_standardised.txt.gz",
                        outcome_name="Longevity-90th", snp_col="SNP",beta_col = "beta.std", se_col = "se.std",
                        ea_col="EA", nea_col="NEA", pval_col="P-value", eaf_col="EAF",
                        chr_col="Chr", pos_col="Position", 
                        output_prefix="Longevity_90th"),
  # T2D (Mahajan)
  T2D = list(outcome_file="Mahajan.NatGenet2018b.T2D.European_mapped_CHR_ALL.txt.gz",
             outcome_name="T2D", snp_col="RSID", beta_col="Beta", se_col="SE",
             ea_col="A1", nea_col="A2", pval_col="P-value", eaf_col="EAF",
             chr_col="Chr", pos_col="Pos", 
             output_prefix="T2D"),
  # DBP
  DBP = list(outcome_file="DBP_standardised.txt.gz",
             outcome_name="DBP", snp_col="rsid", beta_col = "beta.std", se_col = "se.std",
             ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col="effect_allele_frequency",
             chr_col="chromosome", pos_col="base_pair_location", 
             output_prefix="DBP"),
  # SBP
  SBP = list(outcome_file="SBP_standardised.txt.gz",
             outcome_name="SBP", snp_col="rsid", beta_col = "beta.std", se_col = "se.std",
             ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col="effect_allele_frequency", 
             chr_col="chromosome", pos_col="base_pair_location", 
             output_prefix="SBP"),
  # Smoking Initiation
  SmokingInitiation = list(outcome_file="sumstats2/SmokingInitiation.txt.gz",
                           outcome_name="SmokingInitiation", snp_col="RSID", beta_col="BETA", se_col="SE",
                           ea_col="ALT", nea_col="REF", pval_col="PVALUE", eaf_col="AF", 
                           chr_col="CHROM", pos_col="POS",
                           output_prefix="SmokingInitiation"),
  # Cigarettes per Day
  CigarettesPerDay = list(outcome_file="CigarettesPerDay_standardised.txt.gz",
                          outcome_name="CigarettesPerDay", snp_col="RSID", beta_col = "beta.std", se_col = "se.std",
                          ea_col="ALT", nea_col="REF", pval_col="PVALUE", eaf_col="AF", 
                          chr_col="CHROM", pos_col="POS",
                          output_prefix="CigarettesPerDay"),
  # ADHD
  ADHD = list(outcome_file="ADHD2022_iPSYCH_deCODE_PGC.meta.logBeta.gz",
              outcome_name="ADHD", snp_col="SNP", beta_col="logBeta", se_col="SE",
              ea_col="A1", nea_col="A2", pval_col="P", eaf_col="FRQ_A_38691",
              chr_col="CHR", pos_col="BP",
              output_prefix="ADHD"),
  # Bipolar
  Bipolar = list(outcome_file="daner_bip_pgc3_nm_noukbiobank_beta.gz",
                 outcome_name="Bipolar", snp_col="SNP", beta_col="logBeta", se_col="SE",
                 ea_col="A1", nea_col="A2", pval_col="P", eaf_col="FRQ_A_40463",
                 chr_col="CHR", pos_col="BP",
                 output_prefix="Bipolar"),
  # CRP
  CRP = list(outcome_file="CRP_standardised.txt.gz",
             outcome_name="CRP", snp_col="variant_id", beta_col = "beta.std", se_col = "se.std",
             ea_col="A1", nea_col="A2", pval_col="p_value",
             chr_col="chromosome", pos_col="base_pair_location",
             output_prefix="CRP"),
  # Gran
  Gran = list(outcome_file="Gran_EUR_standardised.txt.gz",
             outcome_name="Gran", snp_col="rsID", beta_col="beta.std", se_col="se.std",
             ea_col="A1", nea_col="A2", pval_col="P",
             chr_col="chr", pos_col="bp",
             output_prefix="Gran"),
  # Hannum
  Hannum = list(outcome_file="Hannum_EUR_standardised.txt.gz",
             outcome_name="Hannum", snp_col="rsID", beta_col="beta.std", se_col="se.std",
             ea_col="A1", nea_col="A2", pval_col="P",
             chr_col="chr", pos_col="bp",
             output_prefix="Hannum")
  
)




# -----------------------------
# Run MR analyses for exposure_list (clump_p=5e-8)
# -----------------------------

for (exposure_name in names(exposure_list)) {
  exposure_info <- exposure_list[[exposure_name]]
  
  exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-8, filter_pval=5e-6)
  
  for (outcome_name in names(outcome_info)) {
    run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]])
  }
}

# -----------------------------
# Run MR analyses for exposure_list2 (clump_p=5e-6)
# -----------------------------

for (exposure_name in names(exposure_list2)[2]) {
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
setwd("MR-output-reverse-std")
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
write.xlsx(combined_data, file = paste0("MR_results_BAG_to_trait.std_", Sys.Date(), ".xlsx"), rowNames = FALSE, overwrite = TRUE)

# # tidy output in terminal
# cd MR-output-reverse-std
# mkdir individual-plots combined-plots reverse-MR-results
# mv *.png individual-plots
# mv *.pdf combined-plots
# mkdir reverse-MR-results/combined
# mv MR_results_combined_*.xlsx reverse-MR-results/combined
# mkdir reverse-MR-results/per-trait
# mv MR_results*.xlsx reverse-MR-results/per-trait
