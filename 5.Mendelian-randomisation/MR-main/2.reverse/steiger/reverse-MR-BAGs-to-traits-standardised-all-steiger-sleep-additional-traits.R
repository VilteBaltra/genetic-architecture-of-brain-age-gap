# ================================
# Setup
# ================================
#setwd("/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/MR-analysis")
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
if (!dir.exists("MR-output-reverse-std-steiger")) dir.create("MR-output-reverse-std-steiger")


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

run_MR_for_outcome <- function(exp_dat_clumped, exposure_name, outcome_info, binary="NULL") {
  
  cat("\n\n🔵 Running MR:", exposure_name, " -> ", outcome_info$outcome_name, "\n"); flush.console()
  
  args <- list(
    snps = exp_dat_clumped$SNP,
    filename = outcome_info$outcome_file,
    sep = "\t",
    snp_col = outcome_info$snp_col,
    beta_col = outcome_info$beta_col,
    se_col = outcome_info$se_col,
    effect_allele_col = outcome_info$ea_col,
    other_allele_col = outcome_info$nea_col,
    pval_col = outcome_info$pval_col
  )
  
  # Add optional columns only if not NULL
  if (!is.null(outcome_info$eaf_col)) args$eaf_col <- outcome_info$eaf_col
  if (!is.null(outcome_info$samplesize_col)) args$samplesize_col <- outcome_info$samplesize_col
  if (!is.null(outcome_info$ncase_col)) args$ncase_col <- outcome_info$ncase_col
  if (!is.null(outcome_info$ncontrol_col)) args$ncontrol_col <- outcome_info$ncontrol_col
  if (!is.null(outcome_info$units_col)) args$units_col <- outcome_info$units_col
  
  outcome_data <- do.call(read_outcome_data, args)
  
  # Harmonise
  dat <- harmonise_data(exp_dat_clumped, outcome_data)
  
  # Set prevalence if binary
  if (!is.null(outcome_info$prevalence.outcome) && isTRUE(binary)) {
    dat$prevalence.outcome <- outcome_info$prevalence.outcome
    dat$units.outcome <- "log odds"
    cat("Setting prevalence.outcome =", outcome_info$prevalence.outcome, "for binary outcome\n")
  }
  
  # # Steiger filtering to ensure correct causal direction
  # if (isTRUE(binary)) {
  #   cat("Treating", outcome_name, "as binary exposure (log odds)\n"); flush.console()
  #   dat$units.outcome <- "log odds"
  # }
  
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
  final_results <- my_MR_tests(res, dat)
  res_single <- mr_singlesnp(dat)
  res_loo <- mr_leaveoneout(dat)
  
  # Plots
  p1 <- mr_scatter_plot(res, dat)
  p2 <- mr_forest_plot(res_single)
  p3 <- mr_leaveoneout_plot(res_loo)
  p4 <- mr_funnel_plot(res_single)
  
  # Save plots
  dir.create("MR-output-reverse-std-steiger", showWarnings = FALSE)
  pdf(paste0("MR-output-reverse-std-steiger/", exposure_name, "_to_", outcome_info$outcome_name, "_MR_plots.pdf"), width=8, height=8)
  print(p1[[1]]); print(p2[[1]]); print(p3[[1]]); print(p4[[1]])
  dev.off()
  
  ggsave(paste0("MR-output-reverse-std-steiger/scatter_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p1[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-reverse-std-steiger/forest_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p2[[1]], width=6, height=8, dpi=300)
  ggsave(paste0("MR-output-reverse-std-steiger/leaveoneout_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p3[[1]], width=6, height=6, dpi=300)
  ggsave(paste0("MR-output-reverse-std-steiger/funnel_", exposure_name, "_to_", outcome_info$outcome_name, ".png"), p4[[1]], width=6, height=6, dpi=300)
  
  # Save results
  output_file <- paste0("MR-output-reverse-std-steiger/MR_results_", exposure_name, "_to_", outcome_info$outcome_name, "_", Sys.Date(), ".xlsx")
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
    file = "brainage_factor_nogenr_excl2k_jawinski_noGC_withEAF.txt.gz", # "brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz", 
    sep = "\t",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "A1", other_allele_col = "A2",
    pval_col = "P", samplesize_col = "N",
    eaf_col = "eaf"
  )
  #,
  # Han = list(
  #   file = "METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22_withBeta.txt.gz",
  #   sep = "\t",
  #   snp_col = "MarkerName", beta_col = "beta", se_col = "SE",
  #   effect_allele_col = "Allele1", other_allele_col = "Allele2",
  #   pval_col = "P", samplesize_col = "N"
  # ),
  # Jawinski = list(
  #   file =  "Jawinski_standardised.txt.gz", 
  #   sep = "\t",
  #   snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "A1", other_allele_col = "A2",
  #   pval_col = "P", samplesize_col = "N"
  # ),
  # Leonardsen = list(
  #   file = "Leonardsen_standardised.txt.gz",
  #   sep = "\t",
  #   snp_col = "SNP", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "A1", other_allele_col = "A2",
  #   pval_col = "P", samplesize_col = "N"
  # ),
  # Wen = list(
  #   file = "Wen_standardised.txt.gz",
  #   sep = "\t",
  #   snp_col = "ID", beta_col = "beta.std", se_col = "se.std",
  #   effect_allele_col = "A1_wen", other_allele_col = "A2_wen",
  #   pval_col = "P", samplesize_col = "OBS_CT"
  # )
)

# # Kaufman and Smith only have 1 clumped snp, so using clump_p = 5e-6 instead: 
# exposure_list2 <- list(
#   Kaufman = list(
#     file = "Kaufman_standardised.txt.gz",
#     sep = "\t",
#     snp_col = "SNP",  beta_col = "beta.std", se_col = "se.std",
#     effect_allele_col = "A1_kaufman", other_allele_col = "A2_kaufman", pval_col = "PVAL"
#   ),
#   Smith = list(
#     file = "Smith_standardised.txt.gz",
#     sep = "\t",
#     snp_col = "rsid", beta_col = "beta.std", se_col = "se.std",
#     effect_allele_col = "a1", other_allele_col = "a2",
#     pval_col = "P"
#   )
# )


# -----------------------------
# Define outcome datasets
# -----------------------------
outcome_info <- list(
  # # Metabolic-BAG
  # Metabolic_BAG =list(outcome_file="Metabolic_BAG_standardised.txt.gz",
  #                     outcome_name="Metabolic-BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
  #                     ea_col="A1", nea_col="A2", pval_col="P",
  #                     output_prefix="Metabolic_BAG",
  #                     eaf_col = "EAF", 
  #                     samplesize_col = "OBS_CT"),
  # 
  # # Longevity-90th
  # Longevity_90th = list(outcome_file="Longevity_90th_standardised.txt.gz",
  #                       outcome_name="Longevity-90th", snp_col="SNP",beta_col = "beta.std", se_col = "se.std",
  #                       ea_col="EA", nea_col="NEA", pval_col="P-value", eaf_col="EAF",
  #                       chr_col="Chr", pos_col="Position", 
  #                       output_prefix="Longevity_90th", 
  #                       samplesize_col = "Effective_N"),
  # # T2D (Mahajan)
  # T2D = list(outcome_file="Mahajan.NatGenet2018b.T2D.European_mapped_CHR_ALL.txt.gz",
  #            outcome_name="T2D", snp_col="RSID", beta_col="Beta", se_col="SE",
  #            ea_col="A1", nea_col="A2", pval_col="P-value", eaf_col="EAF",
  #            chr_col="Chr", pos_col="Pos", 
  #            output_prefix="T2D",
  #            samplesize_col = "Neff"),
  # # DBP
  # DBP = list(outcome_file="DBP_standardised.txt.gz",
  #            outcome_name="DBP", snp_col="rsid", beta_col = "beta.std", se_col = "se.std",
  #            ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col="effect_allele_frequency",
  #            chr_col="chromosome", pos_col="base_pair_location", 
  #            output_prefix="DBP",
  #            samplesize_col = "n"),
  # # SBP
  # SBP = list(outcome_file="SBP_standardised.txt.gz",
  #            outcome_name="SBP", snp_col="rsid", beta_col = "beta.std", se_col = "se.std",
  #            ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col="effect_allele_frequency", 
  #            chr_col="chromosome", pos_col="base_pair_location", 
  #            output_prefix="SBP",
  #            samplesize_col = "n"),
  # # Smoking Initiation
  # SmokingInitiation = list(outcome_file="sumstats2/SmokingInitiation.txt.gz",
  #                          outcome_name="SmokingInitiation", snp_col="RSID", beta_col="BETA", se_col="SE",
  #                          ea_col="ALT", nea_col="REF", pval_col="PVALUE", eaf_col="AF", 
  #                          chr_col="CHROM", pos_col="POS",
  #                          output_prefix="SmokingInitiation",
  #                          samplesize_col = "EFFECTIVE_N"),
  # # Cigarettes per Day
  # CigarettesPerDay = list(outcome_file="CigarettesPerDay_standardised.txt.gz",
  #                         outcome_name="CigarettesPerDay", snp_col="RSID", beta_col = "beta.std", se_col = "se.std",
  #                         ea_col="ALT", nea_col="REF", pval_col="PVALUE", eaf_col="AF", 
  #                         chr_col="CHROM", pos_col="POS",
  #                         output_prefix="CigarettesPerDay",
  #                         samplesize_col = "EFFECTIVE_N"),
  # # ADHD
  # ADHD = list(outcome_file="ADHD2022_iPSYCH_deCODE_PGC.meta.logBeta.gz",
  #             outcome_name="ADHD", snp_col="SNP", beta_col="logBeta", se_col="SE",
  #             ea_col="A1", nea_col="A2", pval_col="P", eaf_col="FRQ_A_38691",
  #             chr_col="CHR", pos_col="BP",
  #             output_prefix="ADHD",
  #             ncase_col = "Nca",
  #             ncontrol_col = "Nco",
  #             units_col = "log odds",
  #             binary = TRUE),
  # # Bipolar
  # Bipolar = list(outcome_file="daner_bip_pgc3_nm_noukbiobank_beta.gz",
  #                outcome_name="Bipolar", snp_col="SNP", beta_col="logBeta", se_col="SE",
  #                ea_col="A1", nea_col="A2", pval_col="P", eaf_col="FRQ_A_40463",
  #                chr_col="CHR", pos_col="BP",
  #                output_prefix="Bipolar",
  #                ncase_col = "Nca",
  #                ncontrol_col = "Nco",
  #                units_col = "log odds",
  #                binary = TRUE),
  # # CRP
  # CRP = list(outcome_file="CRP_standardised.txt.gz",
  #            outcome_name="CRP", snp_col="variant_id", beta_col = "beta.std", se_col = "se.std",
  #            ea_col="A1", nea_col="A2", pval_col="p_value",
  #            chr_col="chromosome", pos_col="base_pair_location",
  #            output_prefix="CRP",
  #            samplesize_col = "N"),
  # # Gran
  # Gran = list(outcome_file="Gran_EUR_standardised.txt.gz",
  #            outcome_name="Gran", snp_col="rsID", beta_col="beta.std", se_col="se.std",
  #            ea_col="A1", nea_col="A2", pval_col="P",
  #            chr_col="chr", pos_col="bp",
  #            output_prefix="Gran",
  #            samplesize_col = "N"),
  # # Hannum
  # Hannum = list(outcome_file="Hannum_EUR_standardised.txt.gz",
  #            outcome_name="Hannum", snp_col="rsID", beta_col="beta.std", se_col="se.std",
  #            ea_col="A1", nea_col="A2", pval_col="P",
  #            chr_col="chr", pos_col="bp",
  #            output_prefix="Hannum",
  #            samplesize_col = "N"),
  # # ASD
  ASD = list(outcome_file="iPSYCH-PGC_ASD_Nov2017_logbeta_with_eaf.gz",
             outcome_name="ASD", snp_col="SNP", beta_col="logbeta", se_col="SE",
             ea_col="A1", nea_col="A2", pval_col="P", eaf_col="EAF",
             chr_col="CHR", pos_col="BP",
             output_prefix="ASD",
             ncase_col = "ncase",
             ncontrol_col = "ncontrol",
             units_col = "log odds",
             prevalence.outcome = 0.012, # says 1.2% in https://www.intechopen.com/chapters/84388
             binary = TRUE),
  #  # Sleep duration
  #  Sleep_duration = list(outcome_file="sleepdurationsumstats.txt.gz",
  #            outcome_name="Sleep_duration", snp_col="SNP", beta_col="BETA_SLEEPDURATION", se_col="SE_SLEEPDURATION",
  #            ea_col="ALLELE1", nea_col="ALLELE0", pval_col="P_SLEEPDURATION", eaf_col = "A1FREQ",
  #            chr_col="CHR", pos_col="BP",
  #            output_prefix="Sleep_duration",
  #            samplesize_col = "N",
  #            binary = FALSE),
  # # PD
  # PD =list(outcome_file="GP2_PD_with_Kaviar_IDs_unique_N.gz",
  #          outcome_name="PD", snp_col="ID", beta_col="beta", se_col="standard_error",
  #          ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col = "effect_allele_frequency",
  #          output_prefix="PD",
  #          ncase_col = "ncases",
  #          ncontrol_col = "ncontrols",
  #          units_col = "log odds",
  #          prevalence.outcome = 0.014,
  #          binary = TRUE),
  # # PAI1
  # PAI1 = list(outcome_file="PAI1_EUR_standardised.txt.gz",
  #             outcome_name="PAI1", snp_col="rsID", beta_col="beta.std", se_col="se.std",
  #             ea_col="A1", nea_col="A2", pval_col="P",
  #             chr_col="chr", pos_col="bp",
  #             output_prefix="PAI1",
  #             eaf_col = "Freq1", samplesize_col = "N", binary = FALSE),
  # #GrimAge
  # GrimAge = list(outcome_file="GrimAge_EUR_standardised.txt.gz",
  #                outcome_name="GrimAge", snp_col="rsID", beta_col="beta.std", se_col="se.std",
  #                ea_col="A1", nea_col="A2", pval_col="P",
  #                chr_col="chr", pos_col="bp",
  #                output_prefix="GrimAge",
  #                eaf_col = "Freq1", samplesize_col = "N", binary = FALSE),
  # # IEAA
  # IEAA = list(outcome_file="IEAA_EUR_standardised.txt.gz",
  #             outcome_name="IEAA", snp_col="rsID", beta_col="beta.std", se_col="se.std",
  #             ea_col="A1", nea_col="A2", pval_col="P",
  #             chr_col="chr", pos_col="bp",
  #             output_prefix="IEAA",
  #             eaf_col = "Freq1", samplesize_col = "N", binary = FALSE),
  # # PhenoAge
  # PhenoAge = list(outcome_file="PhenoAge_EUR_standardised.txt.gz",
  #                 outcome_name="PhenoAge", snp_col="rsID", beta_col="beta.std", se_col="se.std",
  #                 ea_col="A1", nea_col="A2", pval_col="P",
  #                 chr_col="chr", pos_col="bp",
  #                 output_prefix="PhenoAge",
  #                 eaf_col = "Freq1", samplesize_col = "N", binary = FALSE),
  # # SCZ
  # Schizophrenia = list(outcome_file="sumstats-check/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",
  #                      outcome_name="Schizophrenia", snp_col="ID", beta_col="BETA", se_col="SE",
  #                      ea_col="A1", nea_col="A2", pval_col="PVAL",
  #                      chr_col="CHROM", pos_col="POS",
  #                      output_prefix="Schizophrenia",
  #                      eaf_col = "FCAS", binary = TRUE,
  #                      ncase_col = "NCAS",
  #                      ncontrol_col = "NCON",
  #                      units_col = "log odds",
  #                      prevalence.outcome = 0.01),
  # Alzheimer
  Alzheimer = list(outcome_file="sumstats-check/alzheimer-35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
                   outcome_name="Alzheimer", snp_col="variant_id", beta_col="beta", se_col="standard_error",
                   ea_col="effect_allele", nea_col="other_allele", pval_col="p_value",
                   chr_col="chromosome", pos_col="base_pair_location",
                   output_prefix="Alzheimer",
                   eaf_col = "effect_allele_frequency", binary = TRUE,
                   ncase_col = "n_cas",
                   ncontrol_col = "n_con",
                   units_col = "log odds",
                   prevalence.outcome = 0.03), # 65–74 years: About 3.0% have probable Alzheimer's. 75–84 years: Approximately 18.7% have probable Alzheimer's. 85+ years: Around 47.2% have probable Alzheimer's.
  # # BMI
  # BMI = list(outcome_file="BMI_standardised.txt.gz",
  #            outcome_name="BMI", snp_col="variant_id", beta_col="beta.std", se_col="se.std",
  #            ea_col="effect_allele", nea_col="other_allele", pval_col="p_value",
  #            chr_col="chromosome", pos_col="base_pair_location",
  #            output_prefix="BMI",
  #            ea_col="A1", nea_col="A2", pval_col="P",
  #            chr_col="chr", pos_col="bp",
  #            eaf_col = "effect_allele_frequency", samplesize_col = "n", binary = FALSE),
  # # Cardiovascular_BAG
  # Cardiovascular_BAG = list(outcome_file="Cardiovascular_BAG_standardised.txt.gz",
  #                           outcome_name="Cardiovascular_BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
  #                           ea_col="A1", nea_col="A2", pval_col="P",
  #                           chr_col="#CHROM", pos_col="POS",
  #                           output_prefix="Cardiovascular_BAG",
  #                           ea_col="A1", nea_col="A2", pval_col="P",
  #                           chr_col="chr", pos_col="bp",
  #                           eaf_col = "EAF", samplesize_col = "OBS_CT", binary = FALSE),
  # #Immune_BAG
  # Immune_BAG = list(outcome_file="Immune_BAG_standardised.txt.gz",
  #                   outcome_name="Immune_BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
  #                   ea_col="A1", nea_col="A2", pval_col="P",
  #                   chr_col="#CHROM", pos_col="POS",
  #                   output_prefix="Immune_BAG",
  #                   eaf_col = "EAF", samplesize_col = "OBS_CT", binary = FALSE),
  # # Renal_BAG
  # Renal_BAG = list(outcome_file="Renal_BAG_standardised.txt.gz",
  #                  outcome_name="Renal_BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
  #                  ea_col="A1", nea_col="A2", pval_col="P",
  #                  chr_col="#CHROM", pos_col="POS",
  #                  output_prefix="Renal_BAG",
  #                  eaf_col = "EAF", samplesize_col = "OBS_CT", binary = FALSE),
  # # Pulmonary_BAG
  # Pulmonary_BAG = list(outcome_file="Pulmonary_BAG_standardised.txt.gz",
  #                      outcome_name="Pulmonary_BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
  #                      ea_col="A1", nea_col="A2", pval_col="P",
  #                      chr_col="#CHROM", pos_col="POS",
  #                      output_prefix="Pulmonary_BAG",
  #                      eaf_col = "EAF", samplesize_col = "OBS_CT", binary = FALSE),
  # # Musculoskeletal_BAG
  # Musculoskeletal_BAG = list(outcome_file="Musculoskeletal_BAG_standardised.txt.gz",
  #                            outcome_name="Musculoskeletal_BAG", snp_col="ID", beta_col="beta.std", se_col="se.std",
  #                            ea_col="A1", nea_col="A2", pval_col="P",
  #                            chr_col="#CHROM", pos_col="POS",
  #                            output_prefix="Musculoskeletal_BAG",
  #                            eaf_col = "EAF", samplesize_col = "OBS_CT", binary = FALSE),
  # # PP
  # PP = list(outcome_file="PP_standardised.txt.gz",
  #           outcome_name="PP", snp_col="rsid", beta_col = "beta.std", se_col = "se.std",
  #           ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col="effect_allele_frequency",
  #           chr_col="chromosome", pos_col="base_pair_location",
  #           output_prefix="PP",
  #           eaf_col = "effect_allele_frequency", samplesize_col = "n", binary = FALSE),
  # # RHR
  # RHR = list(outcome_file="RHR_standardised.txt.gz",
  #            outcome_name="RHR", snp_col="SNP_UBK", beta_col = "beta.std", se_col = "se.std",
  #            ea_col="A1", nea_col="A2", pval_col="P-value", eaf_col="A1Freq",
  #            chr_col="CHR", pos_col="BP",
  #            output_prefix="RHR",
  #            eaf_col = "A1Freq", samplesize_col = "TotalSampleSize", binary = FALSE),
  # # HbA1c
  # HbA1c = list(outcome_file="HbA1c_standardised.txt.gz",
  #              outcome_name="HbA1c", snp_col="variant", beta_col = "beta.std", se_col = "se.std",
  #              ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col="effect_allele_frequency",
  #              chr_col="chromosome", pos_col="base_pair_location",
  #              output_prefix="HbA1c",
  #              eaf_col = "effect_allele_frequency", samplesize_col = "sample_size", binary = FALSE),
  # CAD
  CAD = list(outcome_file="cad.add.160614.website_N.txt.gz",
             outcome_name="CAD", snp_col="markername", beta_col = "beta", se_col = "se_dgc",
             ea_col="effect_allele", nea_col="noneffect_allele", pval_col="p_dgc", eaf_col="effect_allele_freq",
             chr_col="chr", pos_col="bp_hg19",
             output_prefix="CAD",
             eaf_col = "effect_allele_freq",
             ncase_col = "ncases",
             ncontrol_col = "ncontrols", # 60,801 CAD cases and 123,504 controls from https://cardiogramplusc4d.org/data-downloads/
             units_col = "log odds",
             binary = TRUE,
             prevalence.outcome = 0.06), # overall CAD prevalence of 5.8%-6.2% in the US population from https://www.researchprotocols.org/2025/1/e72451
  # MDD
  MDD = list(outcome_file="sumstats-check/pgc-mdd2025_no23andMe_eur_v3-49-24-11_formatted.tsv.gz",
             outcome_name="MDD", snp_col="rsid", beta_col = "beta", se_col = "standard_error",
             ea_col="effect_allele", nea_col="other_allele", pval_col="p_value", eaf_col="effect_allele_frequency",
             chr_col="chromosome", pos_col="base_pair_location",
             output_prefix="MDD",
             ncase_col = "ncases",
             ncontrol_col = "ncontrols",
             units_col = "log odds",
             binary = TRUE,
             prevalence.outcome = 0.057) # 5.7% in https://www.who.int/news-room/fact-sheets/detail/depression
  # # Loneliness
  # Loneliness = list(outcome_file="sumstats-check/Day_2018_NatComs_Social/MTAG_results.txt.gz",
  #                   outcome_name="Loneliness", snp_col="snpid", beta_col = "mtag_beta", se_col = "mtag_se",
  #                   ea_col="a1", nea_col="a2", pval_col="mtag_pval", eaf_col="freq",
  #                   chr_col="chr", pos_col="bpos",
  #                   output_prefix="Loneliness",
  #                   eaf_col = "freq", samplesize_col = "n", binary = FALSE)
)

# -----------------------------
# Run MR analyses for exposure_list (clump_p=5e-8)
# -----------------------------

for (exposure_name in names(exposure_list)) {
  exposure_info <- exposure_list[[exposure_name]]
  
  exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-8, filter_pval=5e-6)
  
  for (outcome_name in names(outcome_info)) {
    binary_flag <- ifelse(!is.null(outcome_info[[outcome_name]]$binary) && outcome_info[[outcome_name]]$binary, TRUE, FALSE)
    message("Binary flag set to: ", binary_flag)
    run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]], binary = binary_flag)
  }
}

# # -----------------------------
# # Run MR analyses for exposure_list2 (clump_p=5e-6)
# # -----------------------------
# 
# for (exposure_name in names(exposure_list2)[2]) {
#   exposure_info <- exposure_list2[[exposure_name]]
#   
#   exp_dat_clumped <- prepare_exposure_data(exposure_info, exposure_name, clump_p=5e-6, filter_pval=5e-6)
#   
#   for (outcome_name in names(outcome_info)) {
#     run_MR_for_outcome(exp_dat_clumped, exposure_name, outcome_info[[outcome_name]])
#   }
# }

# -----------------------------
# Combine output to spreadsheet
# -----------------------------

# get list of all your Excel files in the directory
setwd("MR-output-reverse-std-steiger")
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
write.xlsx(combined_data, file = paste0("MR_results_BAG_to_trait.std_steiger_additional_", Sys.Date(), ".xlsx"), rowNames = FALSE, overwrite = TRUE)

# # tidy output in terminal
# cd MR-output-reverse-std-steiger
# mkdir individual-plots combined-plots reverse-MR-results
# mv *.png individual-plots
# mv *.pdf combined-plots
# mkdir reverse-MR-results/combined
# mv MR_results_combined_*.xlsx reverse-MR-results/combined
# mkdir reverse-MR-results/per-trait
# mv MR_results*.xlsx reverse-MR-results/per-trait
