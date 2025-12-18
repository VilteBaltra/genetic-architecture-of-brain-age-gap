format_loop_Rsq <- function(subdir){
  
  # get cohort name
  cohort_name <- basename(subdir)
  
  # set output dir
  output_dir = subdir
  
  # check for subdirectories that match my pattern
  gwas_paths <- list.dirs(subdir, recursive = TRUE, full.names = TRUE)
  
  # filter for paths ending with "GWAS-output/main"
  target_dirs <- grep("GWAS-output/main$", gwas_paths, value = TRUE)
  
  # change working to gwas output
  cat("Setting working directory to: ", target_dirs, "\n")
  setwd(target_dirs)
  
  # GWAS processing code here:
  ### ----- READ IN AND COMBINE SEPARATE CHR FOR ACP ----- ###
  
  # specify header
  column_header <- c('CHROM','POS','REF','ALT','N_INFORMATIVE','FOUNDER_AF','ALL_AF','INFORMATIVE_ALT_AC','CALL_RATE','HWE_PVALUE','N_REF','N_HET','N_ALT','U_STAT','SQRT_V_STAT','ALT_EFFSIZE','PVALUE')
  
  # specify function to read in data 
  my_fread <- function(chr_number, trait, model){
    chr <- fread(cmd = paste0("zcat ", cohort_name, "_", model, "_", trait, "_chr", chr_number, ".", trait, ".singlevar.score.txt.gz | grep -v '^#'"), header = FALSE)
    cat("Reading in", paste0("zcat ", cohort_name, "_", model, "_", trait, "_chr", chr_number, ".", trait, ".singlevar.score.txt.gz\n"))
    chr <- as.data.frame(chr)
    colnames(chr) <- column_header
    return(chr)
  }
  
  cat("Reading in chromosome files.\n")
  # define empty df to store output 
  chr_all <- data.frame() 
  # check if the "X" directory is present
  # and read in all chromosomes
  print(paste0("Reading in all chr for ", cohort_name, " from ", target_dirs))
  if (any(grepl("X", list.files(target_dirs, recursive = TRUE, full.names = TRUE)))) {
    for (i in c(1:22, "X")) {
      out <- my_fread(chr_number = i, trait = "devage", model = "main")
      chr_all <- rbind(chr_all, out)
    }
  } else {
    message("No X chromosome data found. Reading only chromosomes 1 to 22.\n")
    for (i in c(1:22)) {
      out <- my_fread(chr_number = i, trait = "devage", model = "main")
      chr_all <- rbind(chr_all, out)
    }
  }
  cat("All chromosomes have been read.\n")
  
  # ### ----- ADD INFO ----- ###

  # initiate empty df
  all_chr_info = data.frame()
  
  cat("Reading in info files.\n")
  # read in all chr for info
  if (any(grepl("X", list.files(paste0(subdir, "/INFO"), recursive = TRUE, full.names = TRUE)))) {
    for(i in c(1:22, "X")){ # !think how to handle X!
      print(i)
      chr <- fread(paste0(subdir, "/INFO/chr", i, ".info.gz"))
      all_chr_info <- rbind(all_chr_info, chr)
    }
  } else {
    message("No X chromosome data found. Reading only chromosomes 1 to 22 for INFO files.\n")
    for(i in c(1:22)){ # !think how to handle X!
      print(i)
      chr <- fread(paste0(subdir, "/INFO/chr", i, ".info.gz"))
      all_chr_info <- rbind(all_chr_info, chr)
    }
  }
  cat("All info files have been read.\n")

  # merge (will create more rows than in input dataframes due to duplicates, which we will then remove)
  cat("Merging chromosomes with info.\n")
  all_chr_info$ID <- all_chr_info$SNP
  chr_all$ID <- paste0(chr_all$CHROM, ":", chr_all$POS, ":", chr_all$REF, ":", chr_all$ALT)
 # chr_all$ID <- paste0(chr_all$CHROM, ":", chr_all$POS)
  merged <- merge(chr_all, all_chr_info[, c("ID", "Rsq", "MAF")], by = "ID") 
  head(merged)
  dim(merged)
  
  # remove rows with exact duplicates
  cat("Removing rows with exact duplicates.\n")
  merged_noduplicates <- merged %>%
    distinct()
  
  # # filter on info 0.8
  # merged.filtered = merged %>% filter(Rsq >= 0.80)
  
  ### ----- BASIC QC ----- ###
  # MAF and HWE filtering
  
  # define QC function
  my_QC <- function(data){
    # starting SNP number
    cat("SNP number before QC:", nrow(data), "\n")
    # remove NA
    data_noNA <- data %>% filter(!is.na(PVALUE))
    cat("SNP number after filtering out rows with NA in PVALUE column:", nrow(data_noNA), "\n")
    
    # remove MAF < 0.01 (using this threshold as sample size small)
    data_noNA_maf05 <- data_noNA[data_noNA$MAF >= 0.01,] 
    cat("SNP number after MAF >= 0.01 filtering:", nrow(data_noNA_maf05), "\n")
    
    # exclude snps with HWE_PVALUE < 1e−6
    data_noNA_maf05_hwe <- data_noNA_maf05 %>% filter(HWE_PVALUE >= 1e-6)
    cat("SNP number after HWE < 1e−6 filtering:", nrow(data_noNA_maf05_hwe), "\n")
    
    return(data_noNA_maf05_hwe)
  }
  
  #chr_all_QC <- my_QC(chr_all)
  cat("Running QC.\n")
  chr_all_QC <- my_QC(merged_noduplicates)
  
  ### ----- ADD SNP COLUMN ----- ###
  # ## ADD SNP COLUMN from 1000G ref we read in before the loop
  
  # add SNP column from 1000G reference
  # Convert to data.table
  setDT(chr_all_QC)
  setDT(ref)
  
  # Perform the merge
  cat("Adding SNP column.\n")
  # Remove the last two characters after the second colon in the 'ID' column
  chr_all_QC$ID <- sub("(:[^:]+:[^:]+)$", "", chr_all_QC$ID)
  chr_all_QC_snp <- merge(chr_all_QC, ref[, .(SNP, ID)], by = "ID", all.x = TRUE)
  chr_all_QC_snp$SNP <- ifelse(is.na(chr_all_QC_snp$SNP), chr_all_QC_snp$ID, chr_all_QC_snp$SNP) # if SNP is NA and then CHR:BP_REF_ALT string
  
  # print duplicate nr
  cat("Duplicate SNPs.\n")
  cat(nrow(chr_all_QC_snp), "-", length(unique(chr_all_QC_snp$SNP)), "=", nrow(chr_all_QC_snp) - length(unique(chr_all_QC_snp$SNP)), "duplicate SNPs\n")
  
  # # for duplicates, retain the SNP with the most significant p-value and highest MAF 
  cat("Removing duplicate SNPs based on p-value and MAF.\n")
  chr_all_QC_snp2 <- chr_all_QC_snp %>% group_by(SNP) %>% filter(PVALUE == min(PVALUE)) %>% filter(MAF == max(MAF)) %>% ungroup()
  dim(chr_all_QC_snp2)
  
  # still duplicate snps remain
  # so keeping only one (i.e., only unique SNPs)
  chr_all_QC_snp2 <- chr_all_QC_snp2 %>% distinct(SNP, .keep_all = TRUE)
  dim(chr_all_QC_snp2)
  
  # subset 
  cat("Subsetting to info between 0 and 1.\n")
  chr_all_QC_snp2$Rsq <- as.numeric(chr_all_QC_snp2$Rsq)
  chr_all_QC_snp2_info <- chr_all_QC_snp2[chr_all_QC_snp2$Rsq >= 0 & chr_all_QC_snp2$Rsq <= 1, ]
  cat("Viewing excluded SNPs.\n")
  excluded_snps <- chr_all_QC_snp2[!(chr_all_QC_snp2$SNP %in% chr_all_QC_snp2_info$SNP), ]
  print(excluded_snps)
  
  # select key columns and rename them
  chr_all_QC_snp2_info_sub <- chr_all_QC_snp2_info %>% select(SNP, CHROM, POS, REF, ALT, N_INFORMATIVE, ALL_AF, MAF, ALT_EFFSIZE, SQRT_V_STAT, PVALUE, Rsq)
  names(chr_all_QC_snp2_info_sub) <- c("SNP", "CHR","BP","A2","A1", "N", "EAF", "MAF", "BETA", "SE","P", "Rsq")
  
  # plot R2
  cat("Plotting histogram of Rsq.\n")
  pdf(paste0(output_dir, "/histogram_Rsq_", cohort_name, ".pdf"))
  
  # Create the histogram plot
  hist(as.numeric(chr_all_QC_snp2_info_sub$Rsq), 
       main = "Histogram of Rsq", 
       xlab = "Rsq", 
       ylab = "Frequency", 
       col = "skyblue", 
       border = "white", 
       breaks = 10)  # Adjust 'breaks' for more or fewer bins
  
  # Close the PDF device to save the file
  dev.off()
  
  ### ----- COUNT SIG SNPS ----- ###
  
  # define the function to count SNPs at given p-value thresholds
  cat("Counting SNPs at different thresholds.\n")
  count_sig_snps <- function(data, thresholds) {
    # loop through each threshold
    for (pvalue in thresholds) {
      gwas_hits <- data %>% filter(P < pvalue)
      cat(sprintf("Number of SNPs with P < %g: %d\n", pvalue, nrow(gwas_hits)))
    }
  }
  
  # count snps at defined thresholds
  print(count_sig_snps(chr_all_QC_snp2_info_sub, thresholds = c(5e-8, 5e-7, 5e-6, 5e-5)))
  
  ### ----- SAVE SUMSTATS ----- ###
  cat("Saving summary statistics....\n")
  # save joined and minimally QC'd sumstats
  write.table(chr_all_QC_snp2_info, file=gzfile(paste0(output_dir, "/", cohort_name, '_sumstats_main_chr1to22_rsid_noinfo.txt.gz')),sep="\t",col.names=T,row.names=F,quote=F) 
  # save sumstats with key columns only
  write.table(chr_all_QC_snp2_info_sub, file=gzfile(paste0(output_dir, "/", cohort_name, '_sumstats_main_chr1to22_rsid_keycolumns_noinfo.txt.gz')),sep="\t",col.names=T,row.names=F,quote=F) 
}
