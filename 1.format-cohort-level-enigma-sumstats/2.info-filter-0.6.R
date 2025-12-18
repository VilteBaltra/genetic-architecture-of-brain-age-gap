

library(data.table)

# set wd
setwd("/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024/gwas-results")

# read in gwas file

# Loop through cohorts and filter for info 0.60
cohort_names=c('GOBS', 'BHRCS', 'BILGIN','brazilunicamp', 'COBRE', 'FIDMAG', 'IMH_study', 'PAFIP_3')

for (cohort_name in cohort_names) {
  print(cohort_name)
  cat("Reading in ", paste0(cohort_name, "/", cohort_name, '_sumstats_main_chr1to22_rsid_keycolumns_noinfo.txt.gz', "\n"))
  gwas <- fread(paste0(cohort_name, "/", cohort_name, '_sumstats_main_chr1to22_rsid_keycolumns_noinfo.txt.gz')) 
  
  # filter by info 0.6
  gwas_info <- gwas[gwas$R2 >= 0.6, ]
  
  # print the number of SNPs before and after filtering
  cat("Number of SNPs before filtering:", nrow(gwas), "\n")
  cat("Number of SNPs after filtering:", nrow(gwas_info), "\n")
  
  # sanity check maf
  cat("Sanity checking maf, all should be TRUE \n")
  min(gwas_info$MAF) >= 0.01 # should be TRUE
  max(gwas_info$MAF) <= 0.5 # should be TRUE
  
  cat("Sorting sumstats by CHR and BP \n")
  gwas_info_sorted <- gwas_info[order(gwas_info$CHR, gwas_info$BP), ]
  
  cat("Final SNP nr: \n")
  nrow(gwas_info_sorted)
  
  cat("Saving sumstats filtered for info 0.6 \n")
  write.table(gwas_info_sorted, file=gzfile(paste0(cohort_name, "/", cohort_name, '_sumstats_main_rsid_keycolumns_info0.6.txt.gz')),sep="\t",col.names=T,row.names=F,quote=F) 
  
}

# Rsq
# same loop as before but this time for cohorts that have Rsq column 
#cohort_names2=c('DBIS','DCHS','DNS', 'FBIRN','imagen', 'LBC1936', 'VOX', 'GIG', 'SHIP-2', 'SHIP-T', 'SHIP-T_B2')
cohort_names2=c('GIG', 'SHIP-2', 'SHIP-T', 'SHIP-T_B2', 'NESDA')
cohort_names2=c('COMPIMP')

for (cohort_name in cohort_names2) {
  print(cohort_name)
  cat("Reading in ", paste0(cohort_name, "/", cohort_name, '_sumstats_main_chr1to22_rsid_keycolumns_noinfo.txt.gz', "\n"))
  gwas <- fread(paste0(cohort_name, "/", cohort_name, '_sumstats_main_chr1to22_rsid_keycolumns_noinfo.txt.gz')) 
  
  # filter by info 0.6
  gwas_info <- gwas[gwas$Rsq >= 0.6, ]
  
  # print the number of SNPs before and after filtering
  cat("Number of SNPs before filtering:", nrow(gwas), "\n")
  cat("Number of SNPs after filtering:", nrow(gwas_info), "\n")
  
  # sanity check maf
  cat("Sanity checking maf, all should be TRUE \n")
  print(min(gwas_info$MAF) >= 0.01) # should be TRUE
  print(max(gwas_info$MAF) <= 0.5) # should be TRUE
  
  cat("Sorting sumstats by CHR and BP \n")
  gwas_info_sorted <- gwas_info[order(gwas_info$CHR, gwas_info$BP), ]
  
  cat("Final SNP nr: \n")
  nrow(gwas_info_sorted)
  
  cat("Saving sumstats filtered for info 0.6 as", paste0(cohort_name, "/", cohort_name, '_sumstats_main_rsid_keycolumns_info0.6.txt.gz'), "\n")
  write.table(gwas_info_sorted, file=gzfile(paste0(cohort_name, "/", cohort_name, '_sumstats_main_rsid_keycolumns_info0.6.txt.gz')),sep="\t",col.names=T,row.names=F,quote=F) 
  
}


## separate for ISHARE

cohort_name='ISHARE'
print(cohort_name)
cat("Reading in ", paste0(cohort_name, "/", cohort_name, '_sumstats_main_chr1to22_rsid_noinfo.txt.gz', "\n"))
gwas <- fread(paste0(cohort_name, "/", cohort_name, '_sumstats_main_chr1to22_rsid_noinfo.txt.gz')) 

# filter by info 0.6
gwas_info <- gwas[gwas$Rsq >= 0.6, ]

# print the number of SNPs before and after filtering
cat("Number of SNPs before filtering:", nrow(gwas), "\n")
cat("Number of SNPs after filtering:", nrow(gwas_info), "\n")

# sanity check maf
cat("Sanity checking maf, all should be TRUE \n")
print(min(gwas_info$MAF) >= 0.01) # should be TRUE
print(max(gwas_info$MAF) <= 0.5) # should be TRUE

cat("Sorting sumstats by CHR and BP \n")
gwas_info_sorted <- gwas_info[order(gwas_info$CHR, gwas_info$POS), ]

cat("Final SNP nr: \n")
nrow(gwas_info_sorted)

cat("Saving sumstats filtered for info 0.6 \n")
write.table(gwas_info_sorted, file=gzfile(paste0(cohort_name, "/", cohort_name, '_sumstats_main_rsid_keycolumns_info0.6.txt.gz')),sep="\t",col.names=T,row.names=F,quote=F) 

                

## separate for UKBB and adding N


#python
#import pandas as pd

# UKBB_discov1
# read the gzipped file into a DataFrame
file_path = 'UKBB_discov1_sumstats_main_chr1to23_rsid_noinfo.txt.gz'
df = pd.read_csv(file_path, sep='\t', compression='gzip')

# filter rows where 'INFO' column > 0.6
filtered_df = df[df['INFO'] > 0.6]

# add sample size for metal
filtered_df['N'] = 30809

# save the filtered DataFrame to a gzipped file with the same options as in R
output_file_path = 'UKBB_discov1_sumstats_main_rsid_info0.6.txt.gz'
filtered_df.to_csv(output_file_path, sep='\t', index=False, header=True, quoting=False, compression='gzip')

# UKBB_discov1
# read the gzipped file into a DataFrame
file_path = 'UKBB_discov2_sumstats_main_chr1to23_rsid_noinfo.gz'
df = pd.read_csv(file_path, sep='\t', compression='gzip')

# filter rows where 'INFO' column > 0.6
filtered_df = df[df['INFO'] > 0.6]

# add sample size for metal
filtered_df['N'] = 5170

# save the filtered DataFrame to a gzipped file with the same options as in R
output_file_path = 'UKBB_discov2_sumstats_main_rsid_info0.6.txt.gz'
filtered_df.to_csv(output_file_path, sep='\t', index=False, header=True, quoting=False, compression='gzip')








                
                
