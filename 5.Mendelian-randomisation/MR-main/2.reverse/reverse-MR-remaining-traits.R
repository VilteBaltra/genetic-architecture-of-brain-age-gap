# ================================
# Setup
# ================================
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
if (!dir.exists("MR-output-reverse-std-check")) dir.create("MR-output-reverse-std-check")
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
  dir.create("MR-output-reverse-std-check", showWarnings = FALSE)
  pdf(paste0("MR-output-reverse-std-check/", exposure_name, "_to_", outcome_info$outcome_name, "_MR_plots.pdf"), width=8, height=8)
  print(p1[[1]]); print(p2[[1]]); print(p3[[1]]); print(p4[[1]])
  dev.off()
  
  ggsave(paste0("MR-output-reverse-std-check/scatter_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p1[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-reverse-std-check/forest_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p2[[1]], width=6, height=8, dpi=300)
  ggsave(paste0("MR-output-reverse-std-check/leaveoneout_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p3[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-reverse-std-check/funnel_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p4[[1]], width=6, height=6, dpi=300)
  
  # Save results
  output_file <- paste0("MR-output-reverse-std-check/MR_results_", exposure_name, "_to_", outcome_info$outcome_name, "_", Sys.Date(), ".xlsx")
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
  PD =list(outcome_file="GP2_PD_with_Kaviar_IDs_unique_N.gz",
           outcome_name="PD", snp_col="ID", beta_col="beta", se_col="standard_error",
           ea_col="effect_allele", nea_col="other_allele", pval_col="p_value",
           eaf_col="effect_allele_frequency",
           output_prefix="PD"),
  # PAI1
  PAI1 = list(outcome_file="PAI1_EUR_standardised.txt.gz",
              outcome_name="PAI1", snp_col="rsID", beta_col="beta.std", se_col="se.std",
              ea_col="A1", nea_col="A2", pval_col="P",
              chr_col="chr", pos_col="bp",
              output_prefix="PAI1"),
  # GrimAge
  GrimAge = list(outcome_file="GrimAge_EUR_standardised.txt.gz",
              outcome_name="GrimAge", snp_col="rsID", beta_col="beta.std", se_col="se.std",
              ea_col="A1", nea_col="A2", pval_col="P",
              chr_col="chr", pos_col="bp",
              output_prefix="GrimAge"),
  # IEAA
  IEAA = list(outcome_file="IEAA_EUR_standardised.txt.gz",
              outcome_name="IEAA", snp_col="rsID", beta_col="beta.std", se_col="se.std",
              ea_col="A1", nea_col="A2", pval_col="P",
              chr_col="chr", pos_col="bp",
              output_prefix="IEAA"),
  # PhenoAge
  PhenoAge = list(outcome_file="PhenoAge_EUR_standardised.txt.gz",
                outcome_name="PhenoAge", snp_col="rsID", beta_col="beta.std", se_col="se.std",
                ea_col="A1", nea_col="A2", pval_col="P",
                chr_col="chr", pos_col="bp",
                output_prefix="PhenoAge"),
  # ASD
  ASD = list(outcome_file="sumstats-check/iPSYCH-PGC_ASD_Nov2017_logbeta.gz",
                  outcome_name="ASD", snp_col="SNP", beta_col="logbeta", se_col="SE",
                  ea_col="A1", nea_col="A2", pval_col="P",
                  chr_col="CHR", pos_col="BP",
                  output_prefix="ASD"),
  # SCZ
  Schizophrenia = list(outcome_file="sumstats-check/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",
             outcome_name="Schizophrenia", snp_col="ID", beta_col="BETA", se_col="SE",
             ea_col="A1", nea_col="A2", pval_col="PVAL",
             chr_col="CHROM", pos_col="POS",
             output_prefix="Schizophrenia"),
  # Alzheimer
  Alzheimer = list(outcome_file="sumstats-check/alzheimer-35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
                       outcome_name="Alzheimer", snp_col="variant_id", beta_col="beta", se_col="standard_error",
                       ea_col="effect_allele", nea_col="other_allele", pval_col="p_value",
                       chr_col="chromosome", pos_col="base_pair_location",
                       output_prefix="Alzheimer"),
  # BMI
  BMI = list(outcome_file="BMI_standardised.txt.gz",
                   outcome_name="BMI", snp_col="variant_id", beta_col="beta.std", se_col="se.std",
                   ea_col="effect_allele", nea_col="other_allele", pval_col="p_value",
                   chr_col="chromosome", pos_col="base_pair_location",
                   output_prefix="BMI"),
  # Cardiovascular_BAG
  Cardiovascular_BAG = list(outcome_file="Cardiovascular_BAG_standardised.txt.gz",
                   outcome_name="Cardiovascular_BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
                   ea_col="A1", nea_col="A2", pval_col="P",
                   chr_col="#CHROM", pos_col="POS",
                   output_prefix="Cardiovascular_BAG"),
  #Immune_BAG
  Immune_BAG = list(outcome_file="Immune_BAG_standardised.txt.gz",
                   outcome_name="Immune_BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
                   ea_col="A1", nea_col="A2", pval_col="P",
                   chr_col="#CHROM", pos_col="POS",
                   output_prefix="Immune_BAG"),
  # Renal_BAG
  Renal_BAG = list(outcome_file="Renal_BAG_standardised.txt.gz",
             outcome_name="Renal_BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
             ea_col="A1", nea_col="A2", pval_col="P",
             chr_col="#CHROM", pos_col="POS",
             output_prefix="Renal_BAG"),
  # Pulmonary_BAG
  Pulmonary_BAG = list(outcome_file="Pulmonary_BAG_standardised.txt.gz",
                   outcome_name="Pulmonary_BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
                   ea_col="A1", nea_col="A2", pval_col="P",
                   chr_col="#CHROM", pos_col="POS",
                   output_prefix="Pulmonary_BAG"),
  # Musculoskeletal_BAG
  Musculoskeletal_BAG = list(outcome_file="Musculoskeletal_BAG_standardised.txt.gz",
                   outcome_name="Musculoskeletal_BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
                   ea_col="A1", nea_col="A2", pval_col="P",
                   chr_col="#CHROM", pos_col="POS",
                   output_prefix="Musculoskeletal_BAG"),
  # PP
  PP = list(outcome_file="PP_standardised.txt.gz",
             outcome_name="PP", snp_col="rsid", beta_col = "beta.std", se_col = "se.std",
             ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col="effect_allele_frequency",
             chr_col="chromosome", pos_col="base_pair_location",
             output_prefix="PP"),
  # RHR
  RHR = list(outcome_file="RHR_standardised.txt.gz",
            outcome_name="RHR", snp_col="SNP_UBK", beta_col = "beta.std", se_col = "se.std",
            ea_col="A1", nea_col="A2", pval_col="P-value", eaf_col="A1Freq",
            chr_col="CHR", pos_col="BP",
            output_prefix="RHR"),
  # HbA1c
  HbA1c = list(outcome_file="HbA1c_standardised.txt.gz",
             outcome_name="HbA1c", snp_col="variant", beta_col = "beta.std", se_col = "se.std",
             ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col="effect_allele_frequency",
             chr_col="chromosome", pos_col="base_pair_location",
             output_prefix="HbA1c"),
  # CAD
  CAD = list(outcome_file="sumstats-check/cad.add.160614.website.txt",
             outcome_name="CAD", snp_col="markername", beta_col = "beta", se_col = "se_dgc",
             ea_col="effect_allele", nea_col="noneffect_allele", pval_col="p_dgc", eaf_col="effect_allele_freq",
             chr_col="chr", pos_col="bp_hg19",
             output_prefix="CAD"),
  # MDD
  MDD = list(outcome_file="sumstats-check/pgc-mdd2025_no23andMe_eur_v3-49-24-11_formatted.tsv.gz",
             outcome_name="MDD", snp_col="rsid", beta_col = "beta", se_col = "standard_error",
             ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col="effect_allele_frequency",
             chr_col="chromosome", pos_col="base_pair_location",
             output_prefix="MDD"),
  # Sleep_duration
  Sleep_duration = list(outcome_file="sumstats-check/sleepdurationsumstats.txt",
             outcome_name="Sleep_duration", snp_col="SNP", beta_col = "BETA_SLEEPDURATION", se_col = "SE_SLEEPDURATION",
             ea_col="ALLELE1", nea_col="ALLELE0", pval_col="P_SLEEPDURATION", eaf_col="A1FREQ",
             chr_col="CHR", pos_col="BP",
             output_prefix="Sleep_duration"),
  # Loneliness
  Loneliness = list(outcome_file="sumstats-check/Day_2018_NatComs_Social/MTAG_results.txt.gz",
                        outcome_name="Loneliness", snp_col="snpid", beta_col = "mtag_beta", se_col = "mtag_se",
                        ea_col="a1", nea_col="a2", pval_col="mtag_pval", eaf_col="freq",
                        chr_col="chr", pos_col="bpos",
                        output_prefix="Loneliness")
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

# # for PD remove 1 SNP that is dropped in Steiger analysis to compare results
# exp_dat_clumped=(exp_dat_clumped[!(exp_dat_clumped$SNP %in% 'rs62065444'),])

# -----------------------------
# Run MR analyses for exposure_list2 (clump_p=5e-6)
# -----------------------------

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
setwd("MR-output-reverse-std-check")
file_list <- list.files(pattern = "\\.xlsx$")

# read first sheet of each file and store in a list
data_list <- lapply(file_list, function(file) { 
  read_excel(file, sheet = 1) %>% mutate(SourceFile = file) # to keep track of source
})

# combine all into one data frame
combined_data <- bind_rows(data_list)

# view result
print(combined_data %>% filter(method == "Inverse variance weighted") %>% filter(pval < 0.05))

# save results
write.xlsx(combined_data, file = paste0("MR_results_BAG_to_remaining_traits.std_", Sys.Date(), ".xlsx"), rowNames = FALSE, overwrite = TRUE)

# # tidy output in terminal
# cd MR-output-reverse-std-check
# mkdir individual-plots combined-plots reverse-MR-results
# mv *.png individual-plots
# mv *.pdf combined-plots
# mkdir reverse-MR-results/combined
# mv MR_results_combined_*.xlsx reverse-MR-results/combined
# mkdir reverse-MR-results/per-trait
# mv MR_results*.xlsx reverse-MR-results/per-trait
