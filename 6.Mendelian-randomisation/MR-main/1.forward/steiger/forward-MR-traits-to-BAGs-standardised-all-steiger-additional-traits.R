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
#bfile <- "/campaign/VB-FM5HPC-001/Vilte/Projects/MR-mediation/metabolites/EUR/EUR"
bfile <- "/Users/vb506/Documents/EUR/EUR"

# create an output folder if it does not exist
if (!dir.exists("MR-output-forward-std-steiger")) dir.create("MR-output-forward-std-steiger")


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
  if (!is.null(exposure_info$ncase_col)) format_args$ncase_col <- exposure_info$ncase_col
  if (!is.null(exposure_info$ncontrol_col)) format_args$ncontrol_col <- exposure_info$ncontrol_col
  if (!is.null(exposure_info$units_col)) format_args$units_col <- exposure_info$units_col
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

run_MR_for_outcome <- function(exp_dat_clumped, exposure_name, outcome_info, binary="NULL") {
  
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
    samplesize_col = outcome_info$samplesize_col,
    eaf_col = if (!is.null(outcome_info$eaf_col)) outcome_info$eaf_col else NULL
  )
  outcome_data$outcome <- outcome_info$outcome_name
  
  # Harmonise
  dat <- harmonise_data(exp_dat_clumped, outcome_data)
  
  # Set prevalence if binary
  if (!is.null(exposure_info$prevalence.exposure) && isTRUE(binary)) {
    dat$exposure_info <- exposure_info$prevalence.exposure
    dat$units.exposure <- "log odds"
    cat("Setting prevalence.exposure =", exposure_info$prevalence.exposure, "for binary exposure\n")
  }
  
  # Steiger filtering to ensure correct causal direction
  if (isTRUE(binary)) {
    cat("Treating", exposure_name, "as binary exposure (log odds)\n"); flush.console()
    dat$units.exposure <- "log odds"
  }
  

  dat <- steiger_filtering(dat)
  
  # Count SNPs before filtering
  n_before <- nrow(dat)
  n_pass <- sum(dat$steiger_dir, na.rm = TRUE)
  
  # Keep only SNPs supporting exposure → outcome direction
  dat <- subset(dat, steiger_dir == TRUE)
  
  cat("SNPs retained after Steiger filtering:", n_pass, "of", n_before, "\n")
  
  # Skip if too few remain
  if (nrow(dat) < 3) {
    cat("Too few SNPs after Steiger filtering for MR. Skipping this pair.\n")
    return(NULL)
  }
  
  print(table(dat$steiger_dir))
  print(summary(dat$steiger_pval))
  
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
  dir.create("MR-output-forward-std-steiger", showWarnings = FALSE)
  pdf(paste0("MR-output-forward-std-steiger/", exposure_name, "_to_", outcome_info$outcome_name, "_MR_plots.pdf"), width=8, height=8)
  print(p1[[1]]); print(p2[[1]]); print(p3[[1]]); print(p4[[1]])
  dev.off()
  
  ggsave(paste0("MR-output-forward-std-steiger/scatter_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p1[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-forward-std-steiger/forest_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p2[[1]], width=6, height=8, dpi=300)
  ggsave(paste0("MR-output-forward-std-steiger/leaveoneout_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p3[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-forward-std-steiger/funnel_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p4[[1]], width=6, height=6, dpi=300)
  
  # Save results
  output_file <- paste0("MR-output-forward-std-steiger/MR_results_", exposure_name, "_to_", outcome_info$outcome_name, "_", Sys.Date(), ".xlsx")
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
  # Kaufman = list(
  #   outcome_file = "Kaufman_standardised.txt.gz",
  #   outcome_name = "Kaufman",
  #   sep = "\t",  
  #   snp_col = "SNP", beta_col = "beta.std", se_col = "se.std",
  #   ea_col = "A1_kaufman", nea_col = "A2_kaufman", pval_col = "PVAL",
  #   output_prefix = "Kaufman"
  # ),
  # Smith = list(
  #   outcome_file = "Smith_standardised.txt.gz",
  #   outcome_name = "Smith",
  #   sep = "\t",  
  #   snp_col = "rsid", beta_col = "beta.std", se_col = "se.std",
  #   ea_col = "a1", nea_col = "a2", pval_col = "P",
  #   output_prefix = "Smith"
  # ),
  brainageFactor = list(
    outcome_file = "brainage_factor_nogenr_excl2k_jawinski_noGC_withEAF.txt.gz", # brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz", #brainage_factor_withSmith_withN.txt.gz",
    outcome_name = "brainageFactor",
    sep = "\t",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    eaf_col = "eaf",
    samplesize_col = "N",
    output_prefix = "brainageFactor"
  )
  # Han = list(
  #   outcome_file = "METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22_withBeta.txt.gz", # was "METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-05-27_withBeta.txt.gz"
  #   outcome_name = "Han",
  #   sep = "\t",  
  #   snp_col = "MarkerName", beta_col = "beta", se_col = "SE",
  #   ea_col = "Allele1", nea_col = "Allele2", pval_col = "P",
  #   output_prefix = "Han"
  #  ),
  # Jawinski = list(
  #   outcome_file = "Jawinski_standardised.txt.gz", # brainage2025.full.eur.gwm.gz",
  #   outcome_name = "Jawinski",
  #   sep =  "\t",  
  #   snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
  #   ea_col = "A1", nea_col = "A2", pval_col = "P",
  #   output_prefix = "Jawinski"
  # ),
  # Leonardsen = list(
  #   outcome_file = "Leonardsen_standardised.txt.gz",
  #   outcome_name = "Leonardsen",
  #   sep = "\t",  
  #   snp_col = "SNP", beta_col = "beta.std", se_col = "se.std",
  #   ea_col = "A1", nea_col = "A2", pval_col = "P",
  #   output_prefix = "Leonardsen"
  # ),
  # Wen = list(
  #   outcome_file = "Wen_standardised.txt.gz",
  #   outcome_name = "Wen",
  #   sep = "\t",  
  #   snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
  #   ea_col = "A1_wen", nea_col = "A2_wen", pval_col = "P",
  #   output_prefix = "Wen"
  # )
)


# EXPOSURES 
exposure_list <- list(
  # Metabolic_BAG = list(
  #   file = "Metabolic_BAG_standardised.txt.gz",
  #   sep = "\t",
  #   snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
  #   eaf_col = "EAF", 
  #   samplesize_col = "OBS_CT"
  # ),
  # Longevity_90th = list(
  #   file = "Longevity_90th_standardised.txt.gz",
  #   sep = "\t",
  #   snp_col = "SNP", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "EA", other_allele_col = "NEA", pval_col = "P-value",
  #   eaf_col = "EAF", chr_col = "Chr", pos_col = "Position",
  #   samplesize_col = "Effective_N"
  # ),
  # T2D = list(
  #   file = "Mahajan.NatGenet2018b.T2D.European_mapped_CHR_ALL.txt.gz",
  #   sep = "\t",
  #   snp_col = "RSID", beta_col = "Beta", se_col = "SE",
  #   effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P-value",
  #   eaf_col = "EAF", chr_col = "Chr", pos_col = "Pos",
  #   samplesize_col = "Neff"
  # ),
  # DBP = list(
  #   file = "DBP_standardised.txt.gz",
  #   sep = "\t",
  #   snp_col = "rsid", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
  #   eaf_col = "effect_allele_frequency", chr_col = "chromosome", pos_col = "base_pair_location",
  #   samplesize_col = "n"
  # ),
  # SBP = list(
  #   file = "SBP_standardised.txt.gz",
  #   sep = "\t",
  #   snp_col = "rsid", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
  #   eaf_col = "effect_allele_frequency", chr_col = "chromosome", pos_col = "base_pair_location",
  #   samplesize_col = "n"
  # ),
  # SmokingInitiation = list(
  #   file = "sumstats2/SmokingInitiation.txt.gz",
  #   sep = "\t",
  #   snp_col = "RSID", beta_col = "BETA", se_col = "SE",
  #   effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "PVALUE",
  #   eaf_col = "AF", chr_col = "CHROM", pos_col = "POS",
  #   samplesize_col = "EFFECTIVE_N"
  # ),
  # CigarettesPerDay = list(
  #   file = "CigarettesPerDay_standardised.txt.gz",
  #   sep = "\t",
  #   snp_col = "RSID", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "PVALUE",
  #   eaf_col = "AF", chr_col = "CHROM", pos_col = "POS",
  #   samplesize_col = "EFFECTIVE_N"
  # ),
  # ADHD = list(
  #   file = "ADHD2022_iPSYCH_deCODE_PGC.meta.logBeta.gz",
  #   sep = "\t",
  #   snp_col = "SNP", beta_col = "logBeta", se_col = "SE",
  #   effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
  #   eaf_col = "FRQ_A_38691", chr_col = "CHR", pos_col = "BP",
  #   ncase_col = "Nca",
  #   ncontrol_col = "Nco",
  #   units_col = "log odds",
  #   binary = TRUE
  # ),
  # Bipolar = list(
  #   file = "daner_bip_pgc3_nm_noukbiobank_beta.gz",
  #   sep = "\t",
  #   snp_col = "SNP", beta_col = "logBeta", se_col = "SE",
  #   effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
  #   eaf_col = "FRQ_A_40463", chr_col = "CHR", pos_col = "BP",
  #   ncase_col = "Nca",
  #   ncontrol_col = "Nco",
  #   units_col = "log odds",
  #   binary = TRUE
  # ),
  # CRP = list(
  #   file = "CRP_standardised.txt.gz",
  #   sep = "\t",
  #   snp_col = "variant_id", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "A1", other_allele_col = "A2", pval_col = "p_value",
  #   chr_col = "chromosome", pos_col = "base_pair_location",
  #   samplesize_col = "N"
  # ),
  ASD = list(
    file = "iPSYCH-PGC_ASD_Nov2017_logbeta_with_eaf.gz",
    sep = "\t", snp_col = "SNP", beta_col = "logbeta",
    se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2",
    pval_col = "P", chr_col = "CHR", pos_col = "BP", eaf_col = "EAF",
    ncase_col = "ncase", ncontrol_col = "ncontrol", units_col = "log odds",
    prevalence.exposure = 0.012, binary = TRUE
  ),
  Sleep_duration = list(
    file = "sleepdurationsumstats.txt.gz", sep = "\t", snp_col = "SNP",
    beta_col = "BETA_SLEEPDURATION", se_col = "SE_SLEEPDURATION",
    effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0",
    pval_col = "P_SLEEPDURATION",chr_col = "CHR", pos_col = "BP",
    samplesize_col = "N", eaf_col = "A1FREQ",binary = FALSE
  ),
  PD = list(
    file = "GP2_PD_with_Kaviar_IDs_unique_N.gz",
    sep = "\t",
    snp_col = "ID",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "p_value",
    chr_col = "chromosome",
    pos_col = "base_pair_position",
    eaf_col = "effect_allele_frequency",
    ncase_col = "ncases",
    ncontrol_col = "ncontrols",
    units_col = "log odds",
    prevalence.exposure = 0.014,
    binary = TRUE
  ),
  # Schizophrenia
  Schizophrenia = list(file = "sumstats-check/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz", sep = "\t",
                       snp_col = "ID", beta_col = "BETA", se_col = "SE",
                       effect_allele_col = "A1", other_allele_col = "A2",
                       pval_col = "PVAL", chr_col = "CHROM", pos_col = "POS",
                       eaf_col = "FCAS",
                       ncase_col = "NCAS", ncontrol_col = "NCON",
                       units_col = "log odds", prevalence.exposure = 0.01, binary = TRUE),

  # Alzheimer
  Alzheimer = list(file = "sumstats-check/alzheimer-35379992-GCST90027158-MONDO_0004975.h.tsv.gz", sep = "\t",
                   snp_col = "variant_id", beta_col = "beta", se_col = "standard_error",
                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                   pval_col = "p_value", chr_col = "chromosome", pos_col = "base_pair_location",
                   eaf_col = "effect_allele_frequency",
                   ncase_col = "n_cas", ncontrol_col = "n_con",
                   units_col = "log odds", prevalence.exposure = 0.03, binary = TRUE),
  # BMI
  BMI = list(file = "BMI_standardised.txt.gz", sep = "\t",
             snp_col = "variant_id", beta_col = "beta.std", se_col = "se.std",
             effect_allele_col = "effect_allele", other_allele_col = "other_allele",
             pval_col = "p_value", chr_col = "chromosome", pos_col = "base_pair_location",
             samplesize_col = "n", eaf_col = "effect_allele_frequency",
             binary = FALSE),
  # Cardiovascular_BAG
  Cardiovascular_BAG = list(
    file = "Cardiovascular_BAG_standardised.txt.gz", sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    chr_col = "#CHROM", pos_col = "POS", samplesize_col = "OBS_CT",
    eaf_col = "EAF", binary = FALSE
  ),
  # Immune_BAG
  Immune_BAG = list(
    file = "Immune_BAG_standardised.txt.gz", sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    chr_col = "#CHROM", pos_col = "POS", samplesize_col = "OBS_CT",
    eaf_col = "EAF", binary = FALSE
  ),
  # Renal_BAG
  Renal_BAG = list(
    file = "Renal_BAG_standardised.txt.gz", sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    chr_col = "#CHROM", pos_col = "POS", samplesize_col = "OBS_CT",
    eaf_col = "EAF", binary = FALSE
  ),
  # Pulmonary_BAG
  Pulmonary_BAG = list(
    file = "Pulmonary_BAG_standardised.txt.gz", sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    chr_col = "#CHROM", pos_col = "POS", samplesize_col = "OBS_CT",
    eaf_col = "EAF", binary = FALSE
  ),
  # Musculoskeletal_BAG
  Musculoskeletal_BAG = list(
    file = "Musculoskeletal_BAG_standardised.txt.gz", sep = "\t",
    snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
    chr_col = "#CHROM", pos_col = "POS", samplesize_col = "OBS_CT",
    eaf_col = "EAF", binary = FALSE
  ),
  # PP
  PP = list(
    file = "PP_standardised.txt.gz", sep = "\t",
    snp_col = "rsid", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    chr_col = "chromosome", pos_col = "base_pair_location", samplesize_col = "n",
    eaf_col = "effect_allele_frequency", binary = FALSE
  ),
  # RHR
  RHR = list(
    file = "RHR_standardised.txt.gz", sep = "\t",
    snp_col = "SNP_UBK", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P-value",
    chr_col = "CHR", pos_col = "BP", samplesize_col = "TotalSampleSize",
    eaf_col = "A1Freq", binary = FALSE
  ),
  # HbA1c
  HbA1c = list(
    file = "HbA1c_standardised.txt.gz", sep = "\t",
    snp_col = "variant", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    chr_col = "chromosome", pos_col = "base_pair_location", samplesize_col = "sample_size",
    eaf_col = "effect_allele_frequency", ncase_col = "", ncontrol_col = "", units_col = "",
    prevalence.exposure = "", binary = FALSE
  ),
  #CAD
  CAD = list(
    file = "cad.add.160614.website_N.txt.gz", sep = "\t",
    snp_col = "markername", beta_col = "beta", se_col = "se_dgc",
    effect_allele_col = "effect_allele", other_allele_col = "noneffect_allele", pval_col = "p_dgc",
    chr_col = "chr", pos_col = "bp_hg19",
    eaf_col = "effect_allele_freq", ncase_col = "ncases", ncontrol_col = "ncontrols",
    units_col = "log odds", prevalence.exposure = 0.06, binary = TRUE
  ),
  # MDD
  MDD = list(
    file = "sumstats-check/pgc-mdd2025_no23andMe_eur_v3-49-24-11_formatted.tsv.gz", sep = "\t",
    snp_col = "rsid", beta_col = "beta", se_col = "standard_error",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "p_value",
    chr_col = "chromosome", pos_col = "base_pair_location", samplesize_col = "",
    eaf_col = "effect_allele_frequency", ncase_col = "ncases", ncontrol_col = "ncontrols",
    units_col = "log odds", prevalence.exposure = 0.057, binary = TRUE
  ),
  # Loneliness
  Loneliness = list(
    file = "sumstats-check/Day_2018_NatComs_Social/MTAG_results.txt.gz", sep = "\t",
    snp_col = "snpid", beta_col = "mtag_beta", se_col = "mtag_se",
    effect_allele_col = "a1", other_allele_col = "a2", pval_col = "mtag_pval",
    chr_col = "chr", pos_col = "bpos", samplesize_col = "n",
    eaf_col = "freq", binary = FALSE
  )
)

# EPI CLOCK EXPOSURES, ASD, Longevity 
exposure_list2 <- list(
  # Gran = list(
  #   file = "Gran_EUR_standardised.txt.gz",
  #   sep = "\t",
  #   snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
  #   samplesize_col = "N"
  # ),
  # Hannum = list(
  #   file = "Hannum_EUR_standardised.txt.gz",
  #   sep = "\t",
  #   snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
  #   samplesize_col = "N"
  # ),
  # PAI1
  PAI1 = list(file = "PAI1_EUR_standardised.txt.gz", sep = "\t",
              snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
              effect_allele_col = "A1", other_allele_col = "A2",
              pval_col = "P", chr_col = "chr", pos_col = "bp",
              samplesize_col = "N", eaf_col = "Freq1", binary = FALSE),
  # GrimAge
  GrimAge = list(file = "GrimAge_EUR_standardised.txt.gz", sep = "\t",
                 snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
                 effect_allele_col = "A1", other_allele_col = "A2",
                 pval_col = "P", chr_col = "chr", pos_col = "bp",
                 samplesize_col = "N", eaf_col = "Freq1", binary = FALSE),

  # IEAA
  IEAA = list(file = "IEAA_EUR_standardised.txt.gz", sep = "\t",
              snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
              effect_allele_col = "A1", other_allele_col = "A2",
              pval_col = "P", chr_col = "chr", pos_col = "bp",
              samplesize_col = "N", eaf_col = "Freq1", binary = FALSE),

  # PhenoAge
  PhenoAge = list(file = "PhenoAge_EUR_standardised.txt.gz", sep = "\t",
                  snp_col = "rsID", beta_col = "beta.std", se_col = "se.std",
                  effect_allele_col = "A1", other_allele_col = "A2",
                  pval_col = "P", chr_col = "chr", pos_col = "bp",
                  samplesize_col = "N", eaf_col = "Freq1", binary = FALSE),
  # ASD
  ASD = list( # adding here as only 2 genome-wide significant snps
    file = "iPSYCH-PGC_ASD_Nov2017_logbeta_with_eaf.gz",
    sep = "\t", snp_col = "SNP", beta_col = "logbeta",
    se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2",
    pval_col = "P", chr_col = "CHR", pos_col = "BP", eaf_col = "EAF",
    ncase_col = "ncase", ncontrol_col = "ncontrol", units_col = "log odds",
    prevalence.exposure = 0.012, binary = TRUE
  ),
  Longevity_90th = list(
    file = "Longevity_90th_standardised.txt.gz",
    sep = "\t",
    snp_col = "SNP", beta_col = "beta.std", se_col = "se.std",
    effect_allele_col = "EA", other_allele_col = "NEA", pval_col = "P-value",
    eaf_col = "EAF", chr_col = "Chr", pos_col = "Position",
    samplesize_col = "Effective_N"
  )
)


# -----------------------------
# Run MR analyses (clump_p=5e-8)
# -----------------------------
for (exposure_name in names(exposure_list)) { 
  exposure_info <- exposure_list[[exposure_name]]
  
  exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-8, filter_pval=5e-6)
  
  binary_flag <- ifelse(!is.null(exposure_info$binary) && exposure_info$binary, TRUE, FALSE)
  message("Binary flag set to: ", binary_flag)
  
  for (outcome_name in names(outcome_info)) {
    run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]], binary=binary_flag)
  }
}


# -----------------------------
# Run MR analyses (clump_p=5e-6)
# -----------------------------
# longevity 90th run separately using clump_p=5e-6 as only ~1 SNP after clumping
# ASD as well, as only 2 SNPs after clumping
for (exposure_name in names(exposure_list2)) {
  exposure_info <- exposure_list2[[exposure_name]]
  
  exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-6, filter_pval=5e-6)
  
  binary_flag <- ifelse(!is.null(exposure_info$binary) && exposure_info$binary, TRUE, FALSE)
  message("Binary flag set to: ", binary_flag)
  
  for (outcome_name in names(outcome_info)) {
    run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]], binary=binary_flag)
  }
}


# -----------------------------
# Run MR analyses (epi clocks)
# -----------------------------
# using clump_p=5e-6 as Gran only has 2 ind. snps and Hannum 2-6 depending on exposure-outcome combination
for (exposure_name in names(exposure_list2)) {
  exposure_info <- exposure_list2[[exposure_name]]
  
  exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-6, filter_pval=5e-6)
  
  binary_flag <- ifelse(!is.null(exposure_info$binary) && exposure_info$binary, TRUE, FALSE)
  message("Binary flag set to: ", binary_flag)
  
  for (outcome_name in names(outcome_info)[3]) {
    run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]])
  }
}


# -----------------------------
# Combine output to spreadsheet
# -----------------------------

# get list of all your Excel files in the directory
setwd("MR-setwdoutput-forward-std-steiger")
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
write.xlsx(combined_data, file = paste0("MR_results_trait_to_BAG.std_steiger_", Sys.Date(), ".xlsx"), rowNames = FALSE, overwrite = TRUE)

# # tidy output in terminal
# cd MR-output-forward-std-steiger
# mkdir individual-plots combined-plots forward-MR-results
# mv *.png individual-plots
# mv *.pdf combined-plots
# mkdir forward-MR-results/combined
# mv MR_results_combined_*.xlsx forward-MR-results/combined
# mkdir forward-MR-results/per-trait
# mv MR_results*.xlsx forward-MR-results/per-trait
