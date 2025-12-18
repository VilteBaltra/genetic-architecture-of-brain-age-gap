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
source("functions/format_loop_no_info_filter.R")

# name the log file
output_file <- "format_noinfo_COBRE_FIDMAG_IMH_PAFIP.txt"
output_file <- "format_noinfo_GOBS.txt"

# open the file to capture output
sink(output_file) # 6,11,16,21

# run loop for selected cohorts
for (subdir in subdirs) {
  print(subdir)
  # Run function
  format_loop(subdir=subdir)
}

sink()

