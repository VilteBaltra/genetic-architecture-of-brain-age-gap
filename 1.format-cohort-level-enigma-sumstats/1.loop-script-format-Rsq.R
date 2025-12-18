setwd("/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/")

library(qqman)
library(data.table)
library(tidyverse)
library(dplyr) # for format_loop
library(stringr) # for format_loop
library(data.table) # for format_loop

# # read 1000G reference (for MAF info) from: https://utexas.app.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v
ref <- fread("reference.1000G.maf.0.005.txt", sep = " ")
ref <- as.data.frame(ref)
ref$ID <- paste0(ref$CHR, ":", ref$BP)
dim(ref)
# 9667224

# Base directory
base_dir <- "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results"

# List all subdirectories
subdirs <- dir(base_dir, full.names = TRUE)

# source script for gwas sumstats formatting and basic QC
source("functions/format_loop_Rsq_no_info_filter.R")

# # name the log file
# output_file <- "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/format_info_DCHS_DNS.txt"
# output_file <- "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/format_info_DNS_FBIRN_IMAGEN.txt" # FBIRN was repated manually, as initially no output
# output_file <- "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/format_info_VOX.txt"
#output_file <- "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/format_info_DBIS_DCHS_DNS_IMAGEN_LBC1936_VOX_GIG_SHIP.txt"
output_file <- "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/format_info_DNS_IMAGEN_LBC1936_VOX_GIG_SHIP.txt"
output_file <- "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/format_info_VOX_GIG_SHIP.txt"
output_file <- "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/format_info_SHIP.txt"

# open the file to capture output
sink(output_file) 

#final run that redoes cobre, no MRISHARE as single text files!
# can try this fixed log output as above records nothing
for (subdir in subdirs) {
  print(subdir)
  # Run function
  format_loop_Rsq(subdir=subdir)
}

sink()

