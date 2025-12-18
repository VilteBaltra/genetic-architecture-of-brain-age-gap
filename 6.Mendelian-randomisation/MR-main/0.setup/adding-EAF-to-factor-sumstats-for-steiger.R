# merge in a proxy EAF from one of the original GWASs
setwd("/Users/naveennainwal/Documents/Vilte-2025-NaveenPC/six-brainage-sumstats/")

library(data.table)

# ---------------------------
# Step 1: Read your data
# ---------------------------
# Genomic SEM summary stats
#sem_sumstats <- fread("~/Documents/brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz",header = TRUE)
sem_sumstats <- fread("brainage_factor_nogenr_excl2k_jawinski_noGC_renamed.txt.gz",header = TRUE)
head(sem_sumstats)

# Original GWAS (one of the 6 used in SEM)
gwas <- fread("brainage2025.full.eur.gwm.gz", header = TRUE)
head(gwas)

# ---------------------------
# Step 2: Merge on SNP
# ---------------------------
dat <- merge(sem_sumstats, gwas[, c("ID", "A1", "A2", "A1_FREQ", "N")], by.x = "SNP", by.y = "ID", all.x = TRUE, suffixes = c(".sem", ".gwas"))

# # ---------------------------
# # Step 3: Align alleles
# # ---------------------------
# # If effect alleles do not match, flip the beta
# dat$flip <- dat$A1.sem != dat$A1.gwas
# dat$beta <- ifelse(dat$flip, -dat$beta, dat$beta)
# 
# # Assign harmonized effect and other alleles
# dat$effect_allele <- dat$A1.gwas
# dat$other_allele  <- dat$A2.gwas

# ---------------------------
# Step 3: Align alleles and assign EAF
# ---------------------------
# Determine if effect allele in SEM matches GWAS A1
dat[, flip := A1.sem != A1.gwas]

# Assign effect allele frequency for SEM summary stats
dat[, eaf := ifelse(flip, 1 - A1_FREQ, A1_FREQ)]

# ---------------------------
# Step 4: Optional cleanup
# ---------------------------
# Keep only original SEM columns plus new eaf
# Replace ".sem" with nothing in column names
names(dat) <- gsub("\\.sem$", "", names(dat))
sem_sumstats_eaf <- dat[, c(names(sem_sumstats), "eaf"), with = FALSE]

# ---------------------------
# Step 5: Save updated file
# ---------------------------
## Aligning with MAF 

# define a tolerance for "very close"
tolerance <- 0.08

# Align eaf with MAF --> seems to work well!
sem_sumstats_eaf[, eaf2 := ifelse(abs(eaf - MAF) <= tolerance, 
                                 MAF,          # keep MAF
                                 1 - MAF)]     # flip to match

# keep the eaf2 column and rename to eaf as its based on our BAG factor MAF frequencies rather than Jawinski's
sem_sumstats_eaf <- sem_sumstats_eaf[, !("eaf"), with = FALSE]      # drop eaf
setnames(sem_sumstats_eaf, "eaf2", "eaf")                          # rename eaf2 → eaf


#fwrite(sem_sumstats_eaf, "~/Documents/brainage_factor_nogenr_excl2k_jawinski_noGC_withEAF.txt.gz", sep = "\t")
fwrite(sem_sumstats_eaf, "brainage_factor_nogenr_excl2k_jawinski_noGC_withEAF.txt.gz", sep = "\t")


