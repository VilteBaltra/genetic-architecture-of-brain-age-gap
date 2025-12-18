#### Brain Age 
# written by LAura HAn (L.Han@ggzingeest.nl) and Esther Walton (waltonesther@gmail.com)
# Edited by Constantinos Constantinides (cc2557@bath.ac.uk) on 19/08/2021 
# Adapted for the brain age GWAS by Vilte Baltramonaityte (vb505@bath.ac.uk) 11/12/2021

# setwd

#create log file
log_file <- file("2_BrainAGE.log", open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

# install/load libraries
cat("Prep: installing/loading libraries\n")

load.lib <- function(x){
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      #require( i , character.only = TRUE )
      library(i)
    }
  }
}

#  Then try/install packages...
load.lib( c("ppcor" , "lsmeans" , "multcomp","data.table","plyr","ModelMetrics",
       "caret","gridExtra","Hmisc","pastecs","psych","ggplot2", "dplyr", "testit") )


# source functions
cat("Prep: sourcing functions\n")

source("functions/get.means.R")
source("functions/prepare.files.R")
source("functions/model.fits.brainAge.R")
source("functions/model.fits.brainAge.dx.R")
source("functions/winsorize.R")

# read in Covs and TSV df derived in script 1
load("data_covs_tsv.Rdata") # loads data, Covs and TSV
#names(Covs) <- toupper(names(Covs))
df <- readRDS("df_covs_tsv.rds")
males=df$males
females=df$females
rm(df)

# Select PCs
PCs <- grep("PC", names(Covs), value = TRUE)

cat("Step 3: QC and basic stats for your brainAge measure \n")

model.fits.brainAge(males,"enigma-predictions/males_raw_out.csv")
model.fits.brainAge(females,"enigma-predictions/females_raw_out.csv")

# only for case/conrtol studies
if (!all(is.na(data$DX))){  
model.fits.brainAge.dx(males,"enigma-predictions/males_raw_out.csv")
model.fits.brainAge.dx(females,"enigma-predictions/females_raw_out.csv")
}

cat("Step 4: Run brainAge analysis\n") 

# read in brainAge for males and females and merge together

# males 
output <- read.csv("enigma-predictions/males_raw_out.csv", header=T, sep="\t")
males$predAge <- output$age_prediction
males$devAge <- (males$predAge-males$AGE)
males$AE <- abs(males$devAge)
rm(output)

# females
output <- read.csv("enigma-predictions/females_raw_out.csv", header=T, sep="\t")
females$predAge <- output$age_prediction
females$devAge <- (females$predAge-females$AGE)
females$AE <- abs(females$devAge)
rm(output)

BA=rbind(males,females)

# check that no sample duplication
cat("Duplicated IDs: ",which(duplicated(BA$SUBJID)),"\n")

# merge covariate file with brainAge measure
data=merge(data[,names(Covs),],BA[,c("SUBJID","predAge","devAge","AE")],by="SUBJID")

cat("Saving data\n") 
# save the phenotypic file with covariates 
save(data, file="data_pheno_covariates.Rdata")

# select continuous variables
continuous_covs <- Covs[!(names(Covs) %in% c("DX","SEX", "ETHNICITY","ANCESTRY", "SES", "HAND","DISEASE_CATEGORY", "SUBTYPE"))]

# merge continuous variable in the covariate file with brainAge measure
data_cont=merge(data[,names(continuous_covs),],BA[,c("SUBJID","predAge","devAge","AE")],by="SUBJID")

# get sum stats and correlations for all continuous variables in covariate file
cat('Basic descriptives before winsorizing \n')
if (has_warning(stat.desc(data_cont[,2:length(data_cont)]))) {
  message("NOTE: If you see a warning message with 'no non-missing arguments to min/max' and 'NaNs produced' below, this is likely due to missing values in your Covariates.csv file. You can insepct the 'sumstats.HC' dataframe, if some NAs are expected, the warning messages are okay to ignore.")
}
sumstats=stat.desc(data_cont[,2:length(data_cont)])
tmp=apply(as.matrix(data_cont[,2:length(data_cont)]),2,quantile, na.rm=T)
sumstats=rbind(sumstats[c(1:5,8:10,12:13),],tmp[c(2,4),]) 
  
# this section is only relevant for case/control studies
if (!all(is.na(data$DX))){  # this will only not run if DX is set to NA for all vars
  if (has_warning(stat.desc(data_cont[data$DX==1,2:length(data_cont)]))) {
    message("NOTE: If you see a warning message with 'no non-missing arguments to min/max' and 'NaNs produced' below, this is likely due to missing values in your Covariates.csv file. You can insepct the 'sumstats.cases' dataframe, if some NAs are expected, the warning messages are okay to ignore.")
  }
  sumstats.cases=stat.desc(data_cont[data$DX==1,2:length(data_cont)])
  tmp=apply(as.matrix(data_cont[data$DX==1,2:length(data_cont)]),2,quantile, na.rm=T)
  sumstats.cases=rbind(sumstats.cases[c(1:5,8:10,12:13),],tmp[c(2,4),]) 
  
  if (has_warning(stat.desc(data_cont[data$DX==0,2:length(data_cont)]))) {
    message("NOTE: If you see a warning message with 'no non-missing arguments to min/max' and 'NaNs produced' below, this is likely due to missing values in your Covariates.csv file. You can insepct the 'sumstats.HC' dataframe, if some NAs are expected, the warning messages are okay to ignore.")
  }
  sumstats.HC=stat.desc(data_cont[data$DX==0,2:length(data_cont)])
  tmp=apply(as.matrix(data_cont[data$DX==0,2:length(data_cont)]),2,quantile, na.rm=T)
  sumstats.HC=rbind(sumstats.HC[c(1:5,8:10,12:13),],tmp[c(2,4),])
  rm(tmp) 
}
# if some variables have NA, the code will produce the following warning message (which is expected):

# Warning messages:
# 1: In min(x) : no non-missing arguments to min; returning Inf
# 2: In max(x) : no non-missing arguments to max; returning -Inf
# 3: In qt((0.5 + p/2), (Nbrval - 1)) : NaNs produced

#  obtain counts for categorical variables
categ=c("DX","SEX","ETHNICITY","ANCESTRY", "SES", "HAND","DISEASE_CATEGORY", "SUBTYPE")
data[,categ] <- as.data.frame(lapply(data[,categ], as.factor))
counts=sapply(data[,categ], summary)

# this section is only relevant for case/control studies
if (!all(is.na(data$DX))){ 
  counts.cases=sapply(data[data$DX==1,categ], summary)
  counts.HC=sapply(data[data$DX==0,categ], summary)
}


cat('Checking correlations\n')

dat_matrix=data.matrix(data)
if (has_warning(rcorr(dat_matrix))) {
  message("NOTE: If a warning message appears below stating that NaNs were produced, it is safe to ignore provided some of your variables in the Covariates.csv only contain NA.")
}
correlations=rcorr(dat_matrix) 
# Warning "In sqrt(npair - 2) : NaNs produced" will appear if some variable are NA (it is expected)

if (!all(is.na(data$DX))){  
  save(sumstats,sumstats.HC,sumstats.cases,
       counts,counts.HC,counts.cases,correlations,
       file="sumstats.bf.wins.Rdata")
} else {
  save(sumstats,counts,correlations,
       file="sumstats.bf.wins.Rdata")
}

data_noNA <- data[,colSums(is.na(data))==0] # this removes any columns that have NA
dat_matrix=data.matrix(data_noNA) # will asign numeric values to variables with strings

# if any variable is homogeneous (e.g., all values for Ancestry = 'European') then exclude it from the correlation matrix
# variable is homogeneous if its SD = 0 
Covs_names2 <- colnames(dat_matrix)
for(i in 1:ncol(dat_matrix)) {
  if( sd(dat_matrix[, Covs_names2[i]]) == 0 ) { 
    message("The variable below has standard deviation = 0 and hence is removed from the correlation matrix:")
    print(Covs_names2[i])
    remove <- as.character(Covs_names2[i]) 
    dat_matrix = subset(dat_matrix, select = -c( i ) ) 
  }
}

# save correlation coefficients 
correlations=cor(dat_matrix, use="pairwise.complete.obs")
saveRDS(correlations, file="correlations_noNA.bf.wins.rds")

# get sum stats and correlations for all variables in covariate file
cat('Winsorizing brain data \n')

# winsorize brain data
data$devAge=winsorize(data$devAge)
data$AE <- abs(data$devAge)

# redefine the 'data_cont' variable using winsorized data
data_cont=data[,c(names(continuous_covs),"predAge","devAge","AE")]

# get sum stats and correlations for all continuous variables in covariate file as before
cat('Basic descriptives after winsorizing \n')
if (has_warning(stat.desc(data_cont[,2:length(data_cont)]))) {
  message("NOTE: If you see a warning message with 'no non-missing arguments to min/max' and 'NaNs produced' below, this is likely due to missing values in your Covariates.csv file. You can insepct the 'sumstats.HC' dataframe, if some NAs are expected, the warning messages are okay to ignore.")
}
sumstats=stat.desc(data_cont[,2:length(data_cont)])
tmp=apply(as.matrix(data_cont[,2:length(data_cont)]),2,quantile, na.rm=T)
sumstats=rbind(sumstats[c(1:5,8:10,12:13),],tmp[c(2,4),])


# this section is only relevant for case/control studies
if (!all(is.na(data$DX))){  # this will only not run if DX is set to NA for all vars
  if (has_warning(stat.desc(data_cont[data$DX==1,2:length(data_cont)]))) {
    message("NOTE: If you see a warning message with 'no non-missing arguments to min/max' and 'NaNs produced' below, this is likely due to missing values in your Covariates.csv file. You can insepct the 'sumstats.cases' dataframe, if some NAs are expected, the warning messages are okay to ignore.")
  }
  sumstats.cases=stat.desc(data_cont[data$DX==1,2:length(data_cont)])
  tmp=apply(as.matrix(data_cont[data$DX==1,2:length(data_cont)]),2,quantile, na.rm=T)
  sumstats.cases=rbind(sumstats.cases[c(1:5,8:10,12:13),],tmp[c(2,4),]) 
  
  if (has_warning(stat.desc(data_cont[data$DX==0,2:length(data_cont)]))) {
    message("NOTE: If you see a warning message with 'no non-missing arguments to min/max' and 'NaNs produced' below, this is likely due to missing values in your Covariates.csv file. You can insepct the 'sumstats.HC' dataframe, if some NAs are expected, the warning messages are okay to ignore.")
  }
  sumstats.HC=stat.desc(data_cont[data$DX==0,2:length(data_cont)])
  tmp=apply(as.matrix(data_cont[data$DX==0,2:length(data_cont)]),2,quantile, na.rm=T)
  sumstats.HC=rbind(sumstats.HC[c(1:5,8:10,12:13),],tmp[c(2,4),])
  rm(tmp) 
}
# if some variables are all NA, the code will produce the following warning message (which is expected):

# Warning messages:
# 1: In min(x) : no non-missing arguments to min; returning Inf
# 2: In max(x) : no non-missing arguments to max; returning -Inf
# 3: In qt((0.5 + p/2), (Nbrval - 1)) : NaNs produced

#  obtain counts for categorical variables
categ=c("DX","SEX","ETHNICITY","ANCESTRY", "SES", "HAND","DISEASE_CATEGORY", "SUBTYPE")
data[,categ] <- as.data.frame(lapply(data[,categ], as.factor))
counts=sapply(data[,categ], summary)

# this section is only relevant for case/control studies
if (!all(is.na(data$DX))){ 
  counts.cases=sapply(data[data$DX==1,categ], summary)
  counts.HC=sapply(data[data$DX==0,categ], summary)
}

cat('Checking correlations after winsorizing\n')

dat_matrix=data.matrix(data)
if (has_warning(rcorr(dat_matrix))) {
  message("NOTE: If a warning message appears below stating that NaNs were produced, it is safe to ignore provided some of your variables in the Covariates.csv only contain NA.")
}
correlations=rcorr(dat_matrix) 
# Warning "In sqrt(npair - 2) : NaNs produced" will appear if some variable are NA (it is expected)

if (!all(is.na(data$DX))){  
  save(sumstats,sumstats.HC,sumstats.cases,
     counts,counts.HC,counts.cases,correlations,
     file="sumstats.after.wins.Rdata")
} else {
  save(sumstats,counts,correlations,
       file="sumstats.after.wins.Rdata")
}

# save a correlation plot for a visual overview for winsorised data
data_noNA <- data[,colSums(is.na(data))==0] # this removes any columns that have NA
dat_matrix=data.matrix(data_noNA) # will asign numeric values to variables with strings

# if any variable is homogeneous (e.g., all values for Ancestry = 'European') then exclude it from the correlation matrix
# variable is homogeneous if its SD = 0 
Covs_names2 <- colnames(dat_matrix)
for(i in 1:ncol(dat_matrix)) {
  if( sd(dat_matrix[, Covs_names2[i]]) == 0 ) { 
    message("The variable below has standard deviation = 0 and hence is removed from correlation matrix:")
    print(Covs_names2[i])
    remove <- as.character(Covs_names2[i]) 
    dat_matrix = subset(dat_matrix, select = -c( i ) ) 
  }
}

# save correlation coefficients 
correlations=cor(dat_matrix, use="pairwise.complete.obs")
saveRDS(correlations, file="correlations_noNA.after.wins.rds")

############ REGRESSION MODELS ##############
# These models are not part of the main analysis but may be useful to keep in case reviewers ask for it (e.g., show that there are significant differences btw cases/controls) 

cat("Running the following regression models:\n")

# Model 1.1: Age association with brain age, controlling for site (if needed) 

if (exists("site_regr")){
  cat('model 1.1: age predicting brain-PAD controlling for site\n')
  form <- as.formula(paste("devAge~AGE+AGE2+",site_regr))
  out=lm(formula = form, data = data)
} else {
  cat('model 1.1: age predicting brain-PAD\n')
  out=lm(devAge ~ AGE+AGE2, data = data)
}

#saveRDS(out,file="devAge_Age_ifSite.rds")
capture.output(summary(out), file ="devAge_Age_ifSite.txt")

# Model 1.2: Age association with brain age, controlling for site (if needed) and hand

if ("HAND" %in% names(data) & length(table(data$HAND))>1){
  if (exists("site_regr")){
    cat('model 1.2: age predicting brain-PAD controlling for site and hand\n')
    form <- as.formula(paste("devAge~as.factor(HAND)+AGE+AGE2+",site_regr))
    out=lm(formula = form, data = data)
  } else {
    cat('model 1.2: age predicting brain-PAD controlling for hand\n')
    out=lm(devAge ~ as.factor(HAND)+AGE+AGE2, data = data)
  }
  #saveRDS(out,file="devAge_Age_HAND_ifSite.rds")
  capture.output(summary(out), file ="devAge_Age_HAND_ifSite.txt")
}

# Model 1.3: Age association with brain age, controlling for sex, site (if needed) and genetic PCs

if ("PC1" %in% names(data)){
  PCs_regr = paste(grep("PC",colnames(Covs),value=T),collapse="+")
  if(!all(is.na(data[, PCs]))){
    if (exists("site_regr")){
      cat('model 1.3: age predicting brain-PAD controlling for sex, site and genetic PCs\n')
      form <- as.formula(paste("devAge~AGE+AGE2+as.factor(SEX)",site_regr, PCs_regr, sep = " + "))
      out=lm(formula = form, data = data)
    } else {
      cat('model 1.3: age predicting brain-PAD controlling for sex and genetic PCs\n')
      form2 <- as.formula(paste("devAge~AGE+AGE2+as.factor(SEX)+", PCs_regr))
      out=lm(formula = form2, data = data)
    }
    #saveRDS(out,file="devAge_Age_Sex_PCs_ifSite.rds")
    capture.output(summary(out), file ="devAge_Age_Sex_PCs_ifSite.txt")
  }
}

############ ACCELERATED / DECELERATED AGING
# Model 2.1: Age assocition with accelerated aging, controlling for sex, site (if needed) and genetic PCs

accelerated_data <- data %>% filter(devAge > 0) # select individuals with accelerated aging

if (exists("site_regr")){
  cat('model 2.1: age predicting accelerated aging controlling for sex, age, site and genetic PCs\n')
  form <- as.formula(paste("devAge~AGE+AGE2+as.factor(SEX)", site_regr, PCs_regr, sep = " + "))
  out=lm(formula = form, data = accelerated_data)
} else {
  cat('model 2.1: age predicting accelerated aging controlling for sex, age, and genetic PCs\n')
  form2 <- as.formula(paste("devAge~AGE+AGE2+as.factor(SEX)+", PCs_regr))
  out=lm(formula = form2, data = accelerated_data)
}
#saveRDS(out,file="acceleratedAge_Age_Sex_PCs_ifSite.rds")
capture.output(summary(out), file ="acceleratedAge_Age_Sex_PCs_ifSite.txt")

# Model 2.2: Age assocition with decelerated aging, controlling for sex, site (if needed) and genetic PCs

decelerated_data <- data %>% filter(devAge < 0) # select individuals with decelerated aging
  
if (exists("site_regr")){
  cat('model 2.2: age predicting decelerated aging controlling for sex, age, site and genetic PCs\n')
  form <- as.formula(paste("devAge~AGE+AGE2+as.factor(SEX)",site_regr, PCs_regr, sep = " + "))
  out=lm(formula = form, data = decelerated_data)
} else {
  cat('model 2.2: age predicting decelerated aging controlling for sex, age, and genetic PCs\n')
  form2 <- as.formula(paste("devAge~AGE+AGE2+as.factor(SEX)+", PCs_regr))
  out=lm(formula = form2, data = decelerated_data)
}
#saveRDS(out,file="deceleratedAge_Age_Sex_PCs_ifSite.rds")
capture.output(summary(out), file ="deceleratedAge_Age_Sex_PCs_ifSite.txt")

############ REPEAT FOR CASES/CONTROLS
# below models are only relevant for case/control studies

# Model 3.1: Brain age diff between cases and HC, controlling for key covariates: sex, age and site (if applicable)

if (!all(is.na(data$DX))){ 
  if (exists("site_regr")){
    cat('model 3.1: case status predicting brain-PAD controlling for sex, age, and site\n')
    form <- as.formula(paste("devAge~as.factor(DX)+as.factor(SEX)+AGE+AGE2+",site_regr))
    out=lm(formula = form, data = data)
    } else {
      cat('model 3.1: case status predicting brain-PAD controlling for sex, and age\n')
      out=lm(devAge ~ as.factor(DX) + as.factor(SEX) + AGE + AGE2, data = data)
      }
  #saveRDS(out,file="devAge_Dx_Sex_ifSite.rds")
  capture.output(summary(out), file ="devAge_Dx_Sex_ifSite.txt")
}


# Model 3.2: Brain age diff between cases and HC, controlling sex, age, site (if applicable) and hand

if (!all(is.na(data$DX))){ #  this will only run for case/control cohorts
  if ("HAND" %in% names(data) & length(table(data$HAND))>1){
    if(!all(is.na(data[,"HAND"]))){
      
      if (exists("site_regr")){
        cat('model 3.2: case status predicting brain-PAD controlling for sex, age, site and hand\n')
        form <- as.formula(paste("devAge~as.factor(DX)+as.factor(SEX)+AGE+AGE2+as.factor(HAND)+",site_regr))
        out=lm(formula = form, data = data)
      } else {
        cat('model 3.2: case status predicting brain-PAD controlling for sex, age, and hand\n')
        out=lm(devAge ~ as.factor(DX) + as.factor(SEX)+AGE+AGE2+as.factor(HAND), data = data)
      }
      
      #saveRDS(out,file="devAge_Dx_Sex_HAND_ifSite.rds")
      capture.output(summary(out), file ="devAge_Dx_Sex_HAND_ifSite.txt")
    }
  }
}


# Model 3.3: Brain age diff between cases and HC, controlling sex, age, site (if applicable) and genetic PCs

if (!all(is.na(data$DX))){ #  this will only run for case/control cohorts
  if ("PC1" %in% names(data)){
    if(!all(is.na(data[, PCs]))){
      
      if (exists("site_regr")){
        cat('model 3.3: case status predicting brain-PAD controlling for sex, age, site, and genetic PCs\n')
        form <- as.formula(paste("devAge~as.factor(DX)+as.factor(SEX)+AGE+AGE2+", site_regr, PCs_regr, sep = " + "))
        out=lm(formula = form, data = data)
      } else {
        cat('model 3.3: case status predicting brain-PAD controlling for sex, age, and genetic PCs\n')
        form2 <- as.formula(paste("devAge~as.factor(DX)+as.factor(SEX)+AGE+AGE2+", PCs_regr))
        out=lm(formula = form2, data = data)
      }
      
      #saveRDS(out,file="devAge_Dx_Age_Sex_PCs_ifSite.rds")
      capture.output(summary(out), file ="devAge_Dx_Age_Sex_PCs_ifSite.txt")
    } 
  }
}

############ ACCELERATED / DECELERATED AGING FOR CASE/CONTROLS STUDIES
# Model 4.1: Age assocition with accelerated aging, controlling for sex, site (if needed) and genetic PCs

if (!all(is.na(data$DX))){ #  this will only run for case/control cohorts
  if (exists("site_regr")){
  cat('model 4.1: age predicting accelerated aging controlling for case status, sex, age, site and genetic PCs\n')
  form <- as.formula(paste("devAge~as.factor(DX)+AGE+AGE2+as.factor(SEX)", site_regr, PCs_regr, sep = " + "))
  out=lm(formula = form, data = accelerated_data)
  } else {
  cat('model 4.1: age predicting accelerated aging controlling for case status, sex, age, and genetic PCs\n')
  form2 <- as.formula(paste("devAge~as.factor(DX)+AGE+AGE2+as.factor(SEX)+", PCs_regr))
  out=lm(formula = form2, data = accelerated_data)
  }
  #saveRDS(out,file="acceleratedAge_DX_Age_Sex_PCs_ifSite.rds")
  capture.output(summary(out), file ="acceleratedAge_DX_Age_Sex_PCs_ifSite.txt")
}

# Model 4.2: Age assocition with decelerated aging, controlling for sex, site (if needed) and genetic PCs


if (!all(is.na(data$DX))){ #  this will only run for case/control cohorts
  if (exists("site_regr")){
    cat('model 4.2: age predicting decelerated aging controlling for case status, sex, age, site and genetic PCs\n')
    form <- as.formula(paste("devAge~as.factor(DX)+AGE+AGE2+as.factor(SEX)",site_regr, PCs_regr, sep = " + "))
    out=lm(formula = form, data = decelerated_data)
  } else {
    cat('model 4.2: age predicting decelerated aging controlling for case status, sex, age, and genetic PCs\n')
    form2 <- as.formula(paste("devAge~as.factor(DX)+AGE+AGE2+as.factor(SEX)+", PCs_regr))
    out=lm(formula = form2, data = decelerated_data)
  }
  #saveRDS(out,file="deceleratedAge_DX_Age_Sex_PCs_ifSite.rds")
  capture.output(summary(out), file ="deceleratedAge_DX_Age_Sex_PCs_ifSite.txt")
}


cat('Finished with step 4!\n')

sink()
closeAllConnections()
print(readLines("2_BrainAGE.log"))



