###### Brain age project: Correlations of FreeSurfer ROIs with brain age #######
# Authors: C.Constantinides, E. Walton
# Date: 06/05/2021

# Set working directory 

#create log file
messages=file("3_BrainAGE.log", open="wt")
#rest=file("rest.Rout", open="wt")
sink(messages, type="message")
sink(messages, type="output")


# install/load libraries
cat("Prep: installing/loading libraries\n")
library("Hmisc")
library("dplyr")

# source required functions
cat("Prep: sourcing functions\n")
source("functions/get.means.R")
source("functions/prepare.files.R")

# load df from script 1
df <- readRDS("df_covs_tsv.rds")
males=df$males
females=df$females
rm(df)

cat("Step 5: Run Freesurfer ROIs-BrainAge correlation analysis\n")

# males 
output <- read.csv("enigma-predictions/males_raw_out.csv", header=T, sep="\t")
males$predAge <- output$age_prediction
rm(output)

# females
output <- read.csv("enigma-predictions/females_raw_out.csv", header=T, sep="\t")
females$predAge <- output$age_prediction
rm(output)

ROI_BA=rbind(males,females)

# str(ROI_BA)

# check that no sample duplication
cat("Duplicated IDs:",which(duplicated(ROI_BA$SUBJID)),"\n")

cat('Checking correlations\n')
load("data_covs_tsv.Rdata") # loads data, Covs and TSV saved in script 1

if (all(is.na(data$DX))){ # Check correlations in total sample 
  
  # Check correlations in total sample 
  ROI_BA <- select(ROI_BA, -c('DX'))
  dat_matrix=data.matrix(ROI_BA)
  correlations=rcorr(dat_matrix)
  
  save(correlations,
       file="ROI_BA_correlations.RData")
  
  } else { #  this will only run for case/control cohorts
    
    # Check correlations in total sample (i.e. across HCs and cases)
    dat_matrix=data.matrix(ROI_BA)
    correlations=rcorr(dat_matrix)
    
    # Repeat for HCs and cases only
    # Extract HCs
    ROI_BA_HC <- ROI_BA %>% filter(DX==0)
    # Extract cases
    ROI_BA_CASES <- ROI_BA %>% filter(DX==1)
    
    # Check correlation for HCs and cases separately
    # HCs
    dat_matrix=data.matrix(ROI_BA_HC)
    correlations.HC=rcorr(dat_matrix)
    # Cases
    dat_matrix=data.matrix(ROI_BA_CASES)
    correlations.cases=rcorr(dat_matrix)
    
    save(correlations,correlations.HC, correlations.cases,
         file="ROI_BA_correlations_case_control.RData")
  }


cat('Finished with step 5!\n')

sink()
closeAllConnections()
print(readLines("3_BrainAGE.log"))

