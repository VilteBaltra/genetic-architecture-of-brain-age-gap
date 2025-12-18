# Load required library
library(data.table)
library(qqman)

OUT_PATH="/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/generic-metal/output-2025/"

# list of file paths to formatted GWAS sumstats
file_paths <- c(
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/NESDA/NESDA_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/COMPIMP/COMPIMP_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/GOBS/GOBS_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/ACP/ACP_sumstats_main_chr1to22_rsid_keycolumns.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/BHRCS/BHRCS_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/BILGIN/BILGIN_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/brazilunicamp/brazilunicamp_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/COBRE/COBRE_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/DBIS/DBIS_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/DCHS/DCHS_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/DNS/DNS_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/FBIRN/FBIRN_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/FIDMAG/FIDMAG_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/GenRparents/GenerationRparents_sumstats_main_chr1to22_rsid_keycolumns.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/IGP/IGP_sumstats_main_chr1to22_rsid_keycolumns.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/imagen/imagen_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/IMH_study/IMH_study_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/LBC1936/LBC1936_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/mprc/mprc_sumstats_main_chr1to22_rsid_keycolumns.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/PAFIP_3/PAFIP_3_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/VOX/VOX_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/SHIP-2/SHIP-2_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/SHIP-T/SHIP-T_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/SHIP-T_B2/SHIP-T_B2_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/GIG/GIG_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
  "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/OSLO/OSLO_sumstats_main_chr1to22_rsid_keycolumns_noinfo.txt.gz"
)

# loop and plot each file
for(cohort_path in file_paths){
  print(cohort_path)
  gwas <- fread(cohort_path)
  print(head(gwas))
  # gwas$CHR <- ifelse(gwas$CHR == "X", 23, gwas$CHR)
  # gwas$CHR <- as.numeric(gwas$CHR)
  gwas$CHR[gwas$CHR == "X"] <- 23
  gwas$CHR <- as.integer(gwas$CHR)
  
  # Manhattan plot
  cat("Plotting Manhattan plot... \n")
  options(bitmapType='cairo') # to disable no X11 warning on HPC
  jpeg(paste0(OUT_PATH, basename(dirname(cohort_path)), "_manhattan_info0.6_", Sys.Date(), ".jpeg"), width = 10, height = 6, units = 'in', res = 300)
  manhattan(gwas, chr="CHR",bp="BP",p="P", snp='SNP', main = paste0("Manhattan plot: ", basename(dirname(cohort_path))), col = c("slategray2", "blue4"), suggestiveline = -log10(1e-05))
  dev.off()

  # QQ plot
  cat("Plotting QQ plot... \n")
  jpeg(paste0(OUT_PATH, basename(dirname(cohort_path)), "_qq_info0.6_", Sys.Date(), ".jpeg"), width = 6, height = 4, units = 'in', res = 300)
  qq(gwas$P, main = paste0("Q-Q plot of GWAS p-values for ", basename(dirname(cohort_path))))
  dev.off()
}

# paths to files that may have different format
file_paths <- c("/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/ISHARE/ISHARE_sumstats_main_rsid_keycolumns_info0.6.txt.gz",
                "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/ALSPAC/alspac_sumstats_main_chr1to22_rsid_maf05.txt", 
                "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/UKBB/UKBB_discov1_sumstats_main_rsid_info0.6.txt.gz", 
                "/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results/UKBB/UKBB_discov2_sumstats_main_rsid_info0.6.txt.gz")


# ISHARE
# loop and plot each file
for(cohort_path in file_paths[1]){
  print(cohort_path)
  gwas <- fread(cohort_path)
  print(head(gwas))
  # gwas$CHR <- ifelse(gwas$CHR == "X", 23, gwas$CHR)
  # gwas$CHR <- as.numeric(gwas$CHR)
  gwas$CHR[gwas$CHR == "X"] <- 23
  gwas$CHR <- as.integer(gwas$CHR)
  
  # Manhattan plot
  cat("Plotting Manhattan plot... \n")
  options(bitmapType='cairo') # to disable no X11 warning on HPC
  jpeg(paste0(OUT_PATH, basename(dirname(cohort_path)), "_manhattan_info0.6_", Sys.Date(), ".jpeg"), width = 10, height = 6, units = 'in', res = 300)
  manhattan(gwas, chr="CHR",bp="POS",p="p", snp='SNP', main = paste0("Manhattan plot: ", basename(dirname(cohort_path))), col = c("slategray2", "blue4"), suggestiveline = -log10(1e-05))
  dev.off()
  
  # QQ plot
  cat("Plotting QQ plot... \n")
  jpeg(paste0(OUT_PATH, basename(dirname(cohort_path)), "_qq_info0.6_", Sys.Date(), ".jpeg"), width = 6, height = 4, units = 'in', res = 300)
  qq(gwas$p, main = paste0("Q-Q plot of GWAS p-values for ", basename(dirname(cohort_path))))
  dev.off()
}

# ALSPAC
# loop and plot each file
for(cohort_path in file_paths[2]){
  print(cohort_path)
  gwas <- fread(cohort_path)
  print(head(gwas))
  #gwas$CHR[gwas$CHR == "X"] <- 23
  gwas$CHR <- as.integer(gwas$CHR)
  
  # Manhattan plot
  cat("Plotting Manhattan plot... \n")
  options(bitmapType='cairo') # to disable no X11 warning on HPC
  jpeg(paste0(OUT_PATH, basename(dirname(cohort_path)), "_manhattan_info0.6_", Sys.Date(), ".jpeg"), width = 10, height = 6, units = 'in', res = 300)
  manhattan(gwas, chr="#CHROM",bp="POS",p="P", snp='ID', main = paste0("Manhattan plot: ", basename(dirname(cohort_path))), col = c("slategray2", "blue4"), suggestiveline = -log10(1e-05))
  dev.off()
  
  # QQ plot
  cat("Plotting QQ plot... \n")
  jpeg(paste0(OUT_PATH, basename(dirname(cohort_path)), "_qq_info0.6_", Sys.Date(), ".jpeg"), width = 6, height = 4, units = 'in', res = 300)
  qq(gwas$P, main = paste0("Q-Q plot of GWAS p-values for ", basename(dirname(cohort_path))))
  dev.off()
}

# UKBB
# loop and plot each file
for(cohort_path in file_paths[3]){
  print(cohort_path)
  gwas <- fread(cohort_path)
  print(head(gwas))
  
  # Manhattan plot
  cat("Plotting Manhattan plot... \n")
  options(bitmapType='cairo') # to disable no X11 warning on HPC
  jpeg(paste0(OUT_PATH, basename(dirname(cohort_path)), "discov1_manhattan_info0.6_", Sys.Date(), ".jpeg"), width = 10, height = 6, units = 'in', res = 300)
  manhattan(gwas, chr="CHR",bp="BP",p="P_BOLT_LMM_INF", snp='SNP', main = paste0("Manhattan plot: ", basename(dirname(cohort_path))), col = c("slategray2", "blue4"), suggestiveline = -log10(1e-05))
  dev.off()
  
  # QQ plot
  cat("Plotting QQ plot... \n")
  jpeg(paste0(OUT_PATH, basename(dirname(cohort_path)), "_qq_info0.6_", Sys.Date(), ".jpeg"), width = 6, height = 4, units = 'in', res = 300)
  qq(gwas$P_BOLT_LMM_INF, main = paste0("Q-Q plot of GWAS p-values for ", basename(dirname(cohort_path))))
  dev.off()
}

for(cohort_path in file_paths[4]){
  print(cohort_path)
  gwas <- fread(cohort_path)
  print(head(gwas))
  
  # Manhattan plot
  cat("Plotting Manhattan plot... \n")
  options(bitmapType='cairo') # to disable no X11 warning on HPC
  jpeg(paste0(OUT_PATH, basename(dirname(cohort_path)), "_discov2_manhattan_info0.6_", Sys.Date(), ".jpeg"), width = 10, height = 6, units = 'in', res = 300)
  manhattan(gwas, chr="CHR",bp="BP",p="P_BOLT_LMM_INF", snp='SNP', main = paste0("Manhattan plot: ", basename(dirname(cohort_path))), col = c("slategray2", "blue4"), suggestiveline = -log10(1e-05))
  dev.off()
  
  # QQ plot
  cat("Plotting QQ plot... \n")
  jpeg(paste0(OUT_PATH, basename(dirname(cohort_path)), "_discov2_qq_info0.6_", Sys.Date(), ".jpeg"), width = 6, height = 4, units = 'in', res = 300)
  qq(gwas$P_BOLT_LMM_INF, main = paste0("Q-Q plot of GWAS p-values for ", basename(dirname(cohort_path))))
  dev.off()
}

