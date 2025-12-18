.libPaths("/campaign/VB-FM5HPC-001/Vilte/R/x86_64-pc-linux-gnu-library/4.3")
setwd("/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/MR-analysis")

library(data.table)

# -----------------------------
# add EAF to Kaufmann and Wen
# -----------------------------

### KAUFMANN ###
# add EAF from Jawinski (or 1000 genomes EUR) for oag_pheno_normalized_residualized.Brain_age_gap.glm.linear.gz (Wen) and Kaufman_sumstats.gz (Kaufmann)

kaufman <- fread("Kaufman_sumstats.gz",  data.table = FALSE)

# read in ref to get chr and pos
ref <- fread("../reference.1000G.maf.0.005.txt", sep = " ") # which uses build GRCh37 as doubled checked via dbSNP here https://www.ncbi.nlm.nih.gov/snp/

# First merge as above
kaufman_eaf <- merge(kaufman, ref[, c("SNP", "A1", "A2", "MAF")], by = "SNP", all.x = TRUE)
#kaufman_merged <- merge(kaufman, ref[, c("SNP", "A1", "A2", "MAF")], by = "SNP")

# Rename columns for clarity
names(kaufman_eaf)[names(kaufman_eaf) == "A1.x"] <- "A1_kaufman"
names(kaufman_eaf)[names(kaufman_eaf) == "A2.x"] <- "A2_kaufman"
names(kaufman_eaf)[names(kaufman_eaf) == "A1.y"] <- "A1_ref"
names(kaufman_eaf)[names(kaufman_eaf) == "A2.y"] <- "A2_ref"

# derive EAF based on alignment
kaufman_eaf$EAF <- with(kaufman_eaf, ifelse(A1_kaufman == A1_ref, MAF,
                                            ifelse(A1_kaufman == A2_ref, 1 - MAF, NA)))

# save
write.table(kaufman_eaf, file = gzfile("Kaufman_sumstats_with_eaf.gz"), sep = "\t", row.names = FALSE, quote = FALSE)

### WEN ###

wen <- fread("oag_pheno_normalized_residualized.Brain_age_gap.glm.linear.gz",  data.table = FALSE)

# first merge as above
wen_eaf <- merge(wen, ref[, c("SNP", "A1", "A2", "MAF")], by.x="ID", by.y = "SNP", all.x = TRUE)
#wen_merged <- merge(wen, ref[, c("SNP", "A1", "A2", "MAF")], by = "SNP")

# Rename columns for clarity
names(wen_eaf)[names(wen_eaf) == "A1.x"] <- "A1_wen"
names(wen_eaf)[names(wen_eaf) == "A2.x"] <- "A2_wen"
names(wen_eaf)[names(wen_eaf) == "A1.y"] <- "A1_ref"
names(wen_eaf)[names(wen_eaf) == "A2.y"] <- "A2_ref"

# derive EAF based on alignment
wen_eaf$EAF <- with(wen_eaf, ifelse(A1_wen == A1_ref, MAF,
                                            ifelse(A1_wen == A2_ref, 1 - MAF, NA)))

# save
write.table(wen_eaf, file = gzfile("Wen_sumstats_with_eaf.gz"), sep = "\t", row.names = FALSE, quote = FALSE)


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


# Load outcome_info here

outcome_info <- list(
  Smith = list(
    outcome_file = "Smith_V0140_sumstats_A2swap.txt.gz",
    outcome_name = "Smith",
    sep = "\t",  
    snp_col = "rsid",
    beta_col = c("BETA", "beta", "est"),       # ← vector of possible names
    se_col = c("SE", "se", "se_c"),
    eaf_col = c("af", "eaf", "EAF", "FRQ", "A1_FREQ"),
    samplesize_col = 15952,
    output_prefix = "Smith"
  ),
  Leonardsen = list(
    outcome_file = "leonardsen.txt.gz",
    outcome_name = "Leonardsen",
    sep = "\t",  
    snp_col = "SNP", 
    beta_col = c("BETA", "beta", "est"),    
    se_col = c("SE", "se", "se_c"),
    eaf_col = c("af", "eaf", "EAF", "FRQ", "A1_FREQ"),
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    samplesize_col = c("N"),
    output_prefix = "Leonardsen"
  ),
  Jawinski = list(
    outcome_file = "brainage2025.full.eur.excl2k.gwm.gz",
    outcome_name = "Jawinski",
    sep = " ",  
    snp_col = "ID", 
    beta_col = c("BETA", "beta", "est"),    
    se_col = c("SE", "se", "se_c"),
    eaf_col = c("af", "eaf", "EAF", "FRQ", "A1_FREQ"),
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    samplesize_col = "N",
    output_prefix = "Jawinski"
  ),
  Kaufman = list(
    outcome_file = "Kaufman_sumstats_with_eaf.gz",
    outcome_name = "Kaufman",
    sep = "\t",  
    snp_col = "SNP", beta_col = "BETA", se_col = "SE",
    eaf_col = c("af", "eaf", "EAF", "FRQ", "A1_FREQ"),
    ea_col = "A1_kaufman", nea_col = "A2_kaufman", pval_col = "PVAL",
    samplesize_col = 20170,
    output_prefix = "Kaufman"
  ),
  Wen = list(
    outcome_file = "Wen_sumstats_with_eaf.gz",
    outcome_name = "Wen",
    sep = "\t",  
    snp_col = "ID", beta_col = "BETA", se_col = "SE",
    eaf_col = c("af", "eaf", "EAF", "FRQ", "A1_FREQ"),
    ea_col = "A1_wen", nea_col = "A2_wen", pval_col = "P",
    samplesize_col = 30108,
    output_prefix = "Wen"
  )
)

for (name in names(outcome_info)) {
  info <- outcome_info[[name]]
  
  message("Processing: ", name)
  
  file_path <- info$outcome_file
  sep <- info$sep
  
  df <- fread(cmd = paste("gunzip -c", file_path), sep = sep, data.table = FALSE)
  
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


### epi traits ###
outcome_info <- list(
  Hannum = list(
    outcome_file = "Hannum_EUR_summary_statistics.txt",  # compress if needed
    sep = "\t",
    beta_col = c("Effect"),
    se_col = c("SE"),
    eaf_col = c("Freq1"),
    samplesize_col = c("N"),
    output_prefix = "Hannum_EUR"
  ),
  Gran = list(
    outcome_file = "Gran_EUR_summary_statistics.txt",
    sep = "\t",
    beta_col = c("Effect"),
    se_col = c("SE"),
    eaf_col = c("Freq1"),
    samplesize_col = c("N"),
    output_prefix = "Gran_EUR"
  )
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



# ### other traits below (not done) ###


# -----------------------------
# add EAF to OTHER TRAITS
# -----------------------------

# read in ref to get chr and pos
ref <- fread("../reference.1000G.maf.0.005.txt", sep = " ") # which uses build GRCh37 as doubled checked via dbSNP here https://www.ncbi.nlm.nih.gov/snp/

# read in sumstats that need eaf added
sumstats <- fread("oag_pheno_normalized_residualized.Metabolic_age_gap.glm.linear.gz")
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
write.table(sumstats_eaf, file = gzfile("oag_pheno_normalized_residualized.Metabolic_age_gap_with_eaf.gz"), sep = "\t", row.names = FALSE, quote = FALSE)

## SAME FOR CRP 
# (GCST90029070_buildGRCh37.tsv.gz)

# read in sumstats that need eaf added
sumstats <- fread("GCST90029070_buildGRCh37.tsv.gz")
sumstats$N <- 575531 # from https://www.ebi.ac.uk/gwas/studies/GCST90029070
sumstats <- sumstats %>% rename(A1 = effect_allele, A2 = other_allele)

# First merge as above
sumstats_eaf <- merge(sumstats, ref[, c("SNP", "A1", "A2", "MAF")], by.x="variant_id", by.y = "SNP", all.x = TRUE)

# Rename columns for clarity
names(sumstats_eaf)[names(sumstats_eaf) == "A1.x"] <- "A1"
names(sumstats_eaf)[names(sumstats_eaf) == "A2.x"] <- "A2"
names(sumstats_eaf)[names(sumstats_eaf) == "A1.y"] <- "A1_ref"
names(sumstats_eaf)[names(sumstats_eaf) == "A2.y"] <- "A2_ref"

# derive EAF based on alignment
sumstats_eaf$EAF <- with(sumstats_eaf, ifelse(A1 == A1_ref, MAF,
                                              ifelse(A1 == A2_ref, 1 - MAF, NA)))

# save
write.table(sumstats_eaf, file = gzfile("GCST90029070_buildGRCh37_with_eaf.tsv.gz"), sep = "\t", row.names = FALSE, quote = FALSE)


# define datasets that will be standardised 
outcome_info <- list(
  Metabolic_BAG = list(
    outcome_file = "oag_pheno_normalized_residualized.Metabolic_age_gap_with_eaf.gz",
    sep = "\t",
    snp_col = "ID", beta_col = "BETA", se_col = "SE",
    eaf_col = "EAF",
    ea_col = "A1", nea_col = "A2", pval_col = "P",
    samplesize_col = "OBS_CT",
    output_prefix = "Metabolic_BAG"
  ),
  
  Longevity_90th = list(
    outcome_file = "longevity_90th_percentile_rsid.txt.gz",
    sep = "\t",
    snp_col = "SNP", beta_col = "Beta", se_col = "SE",
    eaf_col = "EAF",
    ea_col = "EA", nea_col = "NEA", pval_col = "P-value",
    chr_col = "Chr", pos_col = "Position",
    samplesize_col = "Effective_N",
    output_prefix = "Longevity_90th"
  ),
  
  # T2D = list( --> no need to standardise as binary phenotype
  #   outcome_file = "Mahajan.NatGenet2018b.T2D.European_mapped_CHR_ALL.txt.gz",
  #   sep = "\t",
  #   snp_col = "RSID", beta_col = "Beta", se_col = "SE",
  #   eaf_col = "EAF",
  #   ea_col = "A1", nea_col = "A2", pval_col = "P-value",
  #   chr_col = "Chr", pos_col = "Pos",
  #   samplesize_col = "Neff",
  #   output_prefix = "T2D"
  # ),
  
  DBP = list(
    outcome_file = "DBP_GCST90310295.h.tsv.gz",
    sep = "\t",
    snp_col = "rsid", beta_col = "beta", se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    ea_col = "effect_allele", nea_col = "other_allele", pval_col = "p_value",
    chr_col = "chromosome", pos_col = "base_pair_location",
    samplesize_col = "n",
    output_prefix = "DBP"
  ),
  
  SBP = list(
    outcome_file = "SBP_GCST90310294.h.tsv.gz",
    sep = "\t",
    snp_col = "rsid", beta_col = "beta", se_col = "standard_error",
    eaf_col = "effect_allele_frequency",
    ea_col = "effect_allele", nea_col = "other_allele", pval_col = "p_value",
    chr_col = "chromosome", pos_col = "base_pair_location",
    samplesize_col = "n",
    output_prefix = "SBP"
  ),
  
  # SmokingInitiation = list( --> no need to standardise as binary 
  #   outcome_file = "sumstats2/SmokingInitiation.txt.gz",
  #   sep = "\t",
  #   snp_col = "RSID", beta_col = "BETA", se_col = "SE",
  #   eaf_col = "AF",
  #   ea_col = "ALT", nea_col = "REF", pval_col = "PVALUE",
  #   chr_col = "CHROM", pos_col = "POS",
  #   samplesize_col = "EFFECTIVE_N",
  #   output_prefix = "SmokingInitiation"
  # ),
  
  CigarettesPerDay = list(
    outcome_file = "sumstats2/CigarettesPerDay.txt.gz",
    sep = "\t",
    snp_col = "RSID", beta_col = "BETA", se_col = "SE",
    eaf_col = "AF",
    ea_col = "ALT", nea_col = "REF", pval_col = "PVALUE",
    chr_col = "CHROM", pos_col = "POS",
    samplesize_col = "EFFECTIVE_N",
    output_prefix = "CigarettesPerDay"
  ),
  
  # ADHD and bipolar --> no need to standardise as binary (OR)
  
  CRP = list(
    outcome_file = "GCST90029070_buildGRCh37_with_eaf.tsv.gz",
    sep = "\t",
    snp_col = "variant_id", beta_col = "beta", se_col = "standard_error",
    eaf_col = "EAF",
    ea_col = "A1", nea_col = "A2", pval_col = "p_value",
    chr_col = "chromosome", pos_col = "base_pair_location",
    samplesize_col = "N",
    output_prefix = "CRP"
  )
)


# Standardising: "Metabolic_BAG" "Longevity_90th" "DBP" "SBP" "CigarettesPerDay" "CRP"    
for (name in names(outcome_info)) {
  info <- outcome_info[[name]]
  
  message("Processing: ", name)
  
  file_path <- info$outcome_file
  sep <- info$sep
  
  df <- fread(cmd = paste("gunzip -c", file_path), sep = sep, data.table = FALSE)
  
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

