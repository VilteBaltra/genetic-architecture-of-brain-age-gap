# Clear workspace
rm(list = ls())
cat("\014")

# written by Esther Walton, inspired by https://github.com/AndreAllegrini/IRISK-p/tree/main/scripts
install.packages("gtsummary")
install.packages("flextable")

library(dplyr)
library(gtsummary) 
library(flextable)
library(haven)
library(psych) 

list.files()

# read in GenR data with phenotypic variables, PGSs, PCs, and any other relevant technical covariates
data <- read_sav("prelim_data_2025-12-08.sav")

# rename GenR variables
data <- data %>% rename(
  # Outcomes
  Insulin.child = InsulineChild9_log,
  Glucose.child = GlucoseChild9_clean,
  HDLchol.child = HDLCholesterolChild9_clean,
  Trigl.child = Triglyceriden9_clean,
  fatmass.child = fat_mass_totalCHILD9,
  CRP.child = CRP9_clean,
  CBCLatt.child = Attention_Problems_TScore_9m,
  CBCLint.child = Internalizing_Problems_TScore_9m,
  IQ.child = WISC13_FSIQ,
  
  # Covariates
  age_metab.child = age_serum, 
  age_fatmass.child = age_serum, 
  age_CRP.child = age_serum,  
  age_CBCL.child = age_cbcl,
  age_IQ.child = age_wisc,
  age_MRI.child = age_mri,
  sex.child = sex,
  education.maternal = education, 
  academicachievement_12.child = academicachievement_12, # this is for descriptives only
  age_academicachievement_12.child = age_academicachievement_12,  # this is for descriptives only
  
  # genetic PCs and technical covariates 
  PC1.child3 = C1_genr3,
  PC2.child3 = C2_genr3,
  PC3.child3 = C3_genr3,
  PC4.child3 = C4_genr3,
  PC5.child3 = C5_genr3,
  
  PC1.child4 = C1_genr4,
  PC2.child4 = C2_genr4,
  PC3.child4 = C3_genr4,
  PC4.child4 = C4_genr4,
  PC5.child4 = C5_genr4,
  
  PC1.mother = C1_mother,
  PC2.mother = C2_mother,
  PC3.mother = C3_mother,
  PC4.mother = C4_mother,
  PC5.mother = C5_mother,
  
  PC1.father = C1_partner,
  PC2.father = C2_partner,
  PC3.father = C3_partner,
  PC4.father = C4_partner, 
  PC5.father = C5_partner,
  
  # Polygenic scores
  PRS_BA.child3 = pgs_child_genr3,
  PRS_BA.child4 = pgs_child_genr4,
  PRS_BA.mother = pgs_mother,
  PRS_BA.father = pgs_partner,
  
  # three brain PADs
  BrainPAD.ENIGMA.child = enigma_brainAge,
  BrainPAD.pyment.child = pyment_brainAge,
  BrainPAD.PyBrainAge.child = pybrain_brainAge
)

#Exclude kids without MRI consent, with braces or incidental findings                
library(dplyr) 
data <- data %>%
  filter(
    mri_consent_f09        == 1,  # consent = yes
    has_braces_mri_f09     == 1,  # braces = no
    exclude_incidental_f09 == 1   # incidental = no
  )


outcomes = c(# phenotypes
  'Insulin.child','Glucose.child','HDLchol.child','Trigl.child','fatmass.child',
  'CRP.child','CBCLatt.child','CBCLint.child','IQ.child',
  # brain ages (diff scores)
  'BrainPAD.pyment.child','BrainPAD.PyBrainAge.child','BrainPAD.ENIGMA.child'
)

# View first few rows
head(data)
data <- zap_labels(data)

# descriptives table
summary_tbl <- data %>% 
  tbl_summary(
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} / {N} ({p}%)"
    ),
    digits = all_continuous() ~ 3,  
    missing_text = "(Missing)",
    missing = 'always'
  ) %>%
  modify_caption("**Sample Characteristics**")

# save word doc

summary_tbl %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "GenR_descriptives_BrainAge.GWAS_trio.docx")

# save csv
summary_df <- as.data.frame(summary_tbl)
write.csv(summary_df, "GenR_descriptives_BrainAge.GWAS_trio.csv", row.names = FALSE)

# Define covariates 
covariates.child3 = c(paste0("PC", 1:5,'.child3'))

covariates.child4 = c(paste0("PC", 1:5,'.child4'))

covariates.mother = c(paste0("PC", 1:5,'.mother'))

covariates.father = c(paste0("PC", 1:5,'.father'))

# Residualize outcomes for covariates
data[,'PRS_BA.child_std3'] <- scale(data[,'PRS_BA.child3']) #standardize PGS first
data[,'PRS_BA.child_std_res3'] <- resid(lm(as.formula(paste('PRS_BA.child_std3', paste(c(covariates.child3), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude))

data[,'PRS_BA.child_std4'] <- scale(data[,'PRS_BA.child4']) #standardize PGS first
data[,'PRS_BA.child_std_res4'] <- resid(lm(as.formula(paste('PRS_BA.child_std4', paste(c(covariates.child4), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude))

data$PRS_BA.child_std <- ifelse(is.na(data$PRS_BA.child_std3), data$PRS_BA.child_std4, data$PRS_BA.child_std3)
data$PRS_BA.child_std_res <- ifelse(is.na(data$PRS_BA.child_std_res3), data$PRS_BA.child_std_res4, data$PRS_BA.child_std_res3)

data[,'PRS_BA.mother_std'] <- scale(data[,'PRS_BA.mother']) #standardize PGS first
data[,'PRS_BA.mother_std_res'] <- resid(lm(as.formula(paste('PRS_BA.mother_std', paste(c(covariates.mother), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude))

data[,'PRS_BA.father_std'] <- scale(data[,'PRS_BA.father']) #standardize PGS first
data[,'PRS_BA.father_std_res'] <- resid(lm(as.formula(paste('PRS_BA.father_std', paste(c(covariates.father), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude))


# Check a few correlations (off-diagonal should be ~0)

data %>% select(PRS_BA.child_std_res, PC1.child3) %>% cor(., use = 'p')
data %>% select(PRS_BA.mother_std_res, PC1.mother) %>% cor(., use = 'p')
data %>% select(PRS_BA.father_std_res, PC1.father) %>% cor(., use = 'p')

# covariance table
data.cor <- data.frame(lapply(data , as.numeric))
c = round(cor(data.cor, use='pairwise.complete.obs'),2)
write.csv(c, file="GenR_corr_BrainAge.GWAS_trio.csv")

#########################################  Regressions

# #1: only covariates
outcomes = c(# phenotypes
  'Insulin.child','Glucose.child','HDLchol.child','Trigl.child','fatmass.child',
  'CRP.child','CBCLatt.child','CBCLint.child','IQ.child',
  # brain ages (diff scores)
  'BrainPAD.pyment.child','BrainPAD.PyBrainAge.child','BrainPAD.ENIGMA.child'
)

# define different covariate sets:

covars_metab = c('age_CRP.child','sex.child')
covars_fatmass =  c('age_CRP.child','sex.child')
covars_CRP =  c('age_CRP.child','sex.child')
covars_CBCL =  c('age_CBCL.child','sex.child')
covars_IQ =  c('age_IQ.child','sex.child')
covars_MRI =  c('age_MRI.child','sex.child')

results_trio_covars.only=list()

for (j in outcomes) {
  if (grepl("BrainPAD", j)) {
    covariates <- covars_MRI
  } else if (grepl("CBCL", j)) {
    covariates <- covars_CBCL
  } else if (grepl("IQ", j)) {
    covariates <- covars_IQ
  } else if (grepl("fatmass", j)) {
    covariates <- covars_fatmass
  } else if (grepl("CRP", j)) {
    covariates <- covars_CRP
  } else {
    covariates <- covars_metab
  }

results_trio_covars.only[[j]]=list()
results_trio_covars.only[[j]][['lm']]=lm(as.formula(paste(j, paste(c(covariates), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude)
results_trio_covars.only[[j]][['rho']]=NA
}

# #2: children only

results_trio_childPRS.covars=list()

for (j in outcomes) {
  if (grepl("BrainPAD", j)) {
    covariates <- covars_MRI
  } else if (grepl("CBCL", j)) {
    covariates <- covars_CBCL
  } else if (grepl("IQ", j)) {
    covariates <- covars_IQ
  } else if (grepl("fatmass", j)) {
    covariates <- covars_fatmass
  } else if (grepl("CRP", j)) {
    covariates <- covars_CRP
  } else {
    covariates <- covars_metab
  }
  
results_trio_childPRS.covars[[j]]=list()
results_trio_childPRS.covars[[j]][['lm']]=lm(as.formula(paste(j, paste(c('PRS_BA.child_std_res',covariates), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude)
results_trio_childPRS.covars[[j]][['rho']]['PRS_BA.child_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.child_std_res'), y = covariates, use = 'pairwise', method = 'pearson')[[2]]
}

# #3: parents only

results_trio_parentsPRS.covars=list()

for (j in outcomes) {
  if (grepl("BrainPAD", j)) {
    covariates <- covars_MRI
  } else if (grepl("CBCL", j)) {
    covariates <- covars_CBCL
  } else if (grepl("IQ", j)) {
    covariates <- covars_IQ
  } else if (grepl("fatmass", j)) {
    covariates <- covars_fatmass
  } else if (grepl("CRP", j)) {
    covariates <- covars_CRP
  } else {
    covariates <- covars_metab
  }

results_trio_parentsPRS.covars[[j]]=list() 
results_trio_parentsPRS.covars[[j]][['lm']]=lm(as.formula(paste(j, paste(c('PRS_BA.mother_std_res','PRS_BA.father_std_res',covariates), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude)
results_trio_parentsPRS.covars[[j]][['rho']]['PRS_BA.father_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.father_std_res'), y = c(covariates,'PRS_BA.mother_std_res'), use = 'pairwise', method = 'pearson')[[2]]
results_trio_parentsPRS.covars[[j]][['rho']]['PRS_BA.mother_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.mother_std_res'), y = c(covariates,'PRS_BA.father_std_res'), use = 'pairwise', method = 'pearson')[[2]]
}


# #4: trio

results_trio_allPRS.covars=list()

for (j in outcomes) {
  if (grepl("BrainPAD", j)) {
    covariates <- covars_MRI
  } else if (grepl("CBCL", j)) {
    covariates <- covars_CBCL
  } else if (grepl("IQ", j)) {
    covariates <- covars_IQ
  } else if (grepl("fatmass", j)) {
    covariates <- covars_fatmass
  } else if (grepl("CRP", j)) {
    covariates <- covars_CRP
  } else {
    covariates <- covars_metab
  }

results_trio_allPRS.covars[[j]]=list() 
results_trio_allPRS.covars[[j]][['lm']]=lm(as.formula(paste(j, paste(c('PRS_BA.child_std_res','PRS_BA.mother_std_res','PRS_BA.father_std_res',covariates), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude) 
results_trio_allPRS.covars[[j]][['rho']]['PRS_BA.child_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.child_std_res'), y = c(covariates,'PRS_BA.mother_std_res', 'PRS_BA.father_std_res'), use = 'pairwise', method = 'pearson')[[2]] 
results_trio_allPRS.covars[[j]][['rho']]['PRS_BA.mother_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.mother_std_res'), y = c(covariates, 'PRS_BA.child_std_res','PRS_BA.father_std_res'), use = 'pairwise', method = 'pearson')[[2]] 
results_trio_allPRS.covars[[j]][['rho']]['PRS_BA.father_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.father_std_res'), y = c(covariates, 'PRS_BA.mother_std_res', 'PRS_BA.child_std_res'), use = 'pairwise', method = 'pearson')[[2]] 
}
  
# save tmp
save(results_trio_covars.only,
     results_trio_childPRS.covars,
     results_trio_parentsPRS.covars,
     results_trio_allPRS.covars,
     file= 'tmp_results.Rdata')

# save model output into .csv files 
# first define objects
objects <- list(
  covars_only = results_trio_covars.only,
  childPRS = results_trio_childPRS.covars,
  parentsPRS = results_trio_parentsPRS.covars,
  allPRS = results_trio_allPRS.covars
)


# loop over each object
for (obj_name in names(objects)) {
  
  model_data <- objects[[obj_name]]
  summary_list <- list()
  
  # Loop over models within each outcome object
  for (model_name in names(model_data)) {
    
    # Extract model object
    model_obj <- model_data[[model_name]][["lm"]]
    model_summary <- summary(model_obj)
    
    # Grab coefficients table
    model_df <- as.data.frame(model_summary$coefficients)
    model_df$model <- model_name
    
    # Add R-squared
    model_df$r.squared      <- model_summary$r.squared
    model_df$adj.r.squared  <- model_summary$adj.r.squared
    
    # Add partial correlations (default NA)
    model_df$rho <- NA
    
    # RHO values 
    rho_vals  <- model_data[[model_name]][["rho"]]
    rho_names <- names(rho_vals)
    
    # Fill rho values for matching predictors
    for (i in seq_along(rho_vals)) {
      nm <- rho_names[i]
      
      if (is.null(nm) || length(nm) == 0 || nm == "") next
      
      if (nm %in% rownames(model_df)) {
        model_df[grep(nm, rownames(model_df)), "rho"] <- rho_vals[i]
      }
    }
    
    # Extract F-statistic and p-value for the model
    f_stat_value <- model_summary$fstatistic[1]
    numdf        <- model_summary$fstatistic[2]
    dendf        <- model_summary$fstatistic[3]
    
    model_df$model_pval <- pf(f_stat_value, numdf, dendf, lower.tail = FALSE)
    
    # Add sample size N
    model_df$N <- numdf + dendf + 1
    
    # store this model's summary
    summary_list[[model_name]] <- model_df
  }
  
  # Combine results for all models in this object
  summary_all <- do.call(rbind, summary_list)
  
  # Save output
  file_name <- paste0("BrainPAD_models_summary_", obj_name, ".csv")
  write.csv(summary_all, file_name, row.names = TRUE)
}

############################################### stand betas

# scale all outcomes and covariates(predictors already scaled), other than sex
scale_var=unique(c(outcomes,
                   covars_metab,
                   covars_fatmass,
                   covars_CRP,
                   covars_CBCL,
                   covars_IQ,
                   covars_MRI))

scale_var=setdiff(scale_var,c('sex.child'))
data[scale_var] <- scale(data[scale_var])

results_trio_covars.only_std=list()

for (j in outcomes) {
  if (grepl("BrainPAD", j)) {
    covariates <- covars_MRI
  } else if (grepl("CBCL", j)) {
    covariates <- covars_CBCL
  } else if (grepl("IQ", j)) {
    covariates <- covars_IQ
  } else if (grepl("fatmass", j)) {
    covariates <- covars_fatmass
  } else if (grepl("CRP", j)) {
    covariates <- covars_CRP
  } else {
    covariates <- covars_metab
  }
  
  
  results_trio_covars.only_std[[j]]=list()
  results_trio_covars.only_std[[j]][['lm']]=lm(as.formula(paste(j, paste(c(covariates), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude)
  results_trio_covars.only_std[[j]][['rho']]=NA
  
}

results_trio_childPRS.covars_std=list()

for (j in outcomes) {
  if (grepl("BrainPAD", j)) {
    covariates <- covars_MRI
  } else if (grepl("CBCL", j)) {
    covariates <- covars_CBCL
  } else if (grepl("IQ", j)) {
    covariates <- covars_IQ
  } else if (grepl("fatmass", j)) {
    covariates <- covars_fatmass
  } else if (grepl("CRP", j)) {
    covariates <- covars_CRP
  } else {
    covariates <- covars_metab
  }
 
  results_trio_childPRS.covars_std[[j]]=list()
  results_trio_childPRS.covars_std[[j]][['lm']]=lm(as.formula(paste(j, paste(c('PRS_BA.child_std_res',covariates), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude)
  results_trio_childPRS.covars_std[[j]][['rho']]['PRS_BA.child_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.child_std_res'), y = covariates, use = 'pairwise', method = 'pearson')[[2]]
}

results_trio_parentsPRS.covars_std=list()

for (j in c(outcomes,'PRS_BA.child_std_res')) {
  if (grepl("BrainPAD|PRS", j)) {
    covariates <- covars_MRI
  } else if (grepl("CBCL", j)) {
    covariates <- covars_CBCL
  } else if (grepl("IQ", j)) {
    covariates <- covars_IQ
  } else if (grepl("fatmass", j)) {
    covariates <- covars_fatmass
  } else if (grepl("CRP", j)) {
    covariates <- covars_CRP
  } else {
    covariates <- covars_metab
  }
  # EW: 4 new lines 20251411
  results_trio_parentsPRS.covars_std[[j]]=list() 
  results_trio_parentsPRS.covars_std[[j]][['lm']]=lm(as.formula(paste(j, paste(c('PRS_BA.mother_std_res','PRS_BA.father_std_res',covariates), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude)
  results_trio_parentsPRS.covars_std[[j]][['rho']]['PRS_BA.father_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.father_std_res'), y = c(covariates,'PRS_BA.mother_std_res'), use = 'pairwise', method = 'pearson')[[2]]
  results_trio_parentsPRS.covars_std[[j]][['rho']]['PRS_BA.mother_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.mother_std_res'), y = c(covariates,'PRS_BA.father_std_res'), use = 'pairwise', method = 'pearson')[[2]]
}

results_trio_allPRS.covars_std=list()

for (j in outcomes) {
  if (grepl("BrainPAD", j)) {
    covariates <- covars_MRI
  } else if (grepl("CBCL", j)) {
    covariates <- covars_CBCL
  } else if (grepl("IQ", j)) {
    covariates <- covars_IQ
  } else if (grepl("fatmass", j)) {
    covariates <- covars_fatmass
  } else if (grepl("CRP", j)) {
    covariates <- covars_CRP
  } else {
    covariates <- covars_metab
  }
 
  results_trio_allPRS.covars_std[[j]]=list() 
  results_trio_allPRS.covars_std[[j]][['lm']]=lm(as.formula(paste(j, paste(c('PRS_BA.child_std_res','PRS_BA.mother_std_res','PRS_BA.father_std_res',covariates), collapse=" + "), sep = " ~ ")), data = data, na.action=na.exclude) 
  results_trio_allPRS.covars_std[[j]][['rho']]['PRS_BA.child_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.child_std_res'), y = c(covariates,'PRS_BA.mother_std_res', 'PRS_BA.father_std_res'), use = 'pairwise', method = 'pearson')[[2]] 
  results_trio_allPRS.covars_std[[j]][['rho']]['PRS_BA.mother_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.mother_std_res'), y = c(covariates, 'PRS_BA.child_std_res','PRS_BA.father_std_res'), use = 'pairwise', method = 'pearson')[[2]] 
  results_trio_allPRS.covars_std[[j]][['rho']]['PRS_BA.father_std_res']= partial.r(data = data, x = c(j, 'PRS_BA.father_std_res'), y = c(covariates, 'PRS_BA.mother_std_res', 'PRS_BA.child_std_res'), use = 'pairwise', method = 'pearson')[[2]] 
}

save(results_trio_covars.only_std,
     results_trio_childPRS.covars_std,
     results_trio_parentsPRS.covars_std,
     results_trio_allPRS.covars_std,
     file= 'tmp_results_std.Rdata')

# save model output into .csv files 
# first define objects
objects <- list(
  covars_only_std = results_trio_covars.only_std,
  childPRS_std = results_trio_childPRS.covars_std,
  parentsPRS_std = results_trio_parentsPRS.covars_std,
  allPRS_std = results_trio_allPRS.covars_std
)


# loop over each object
for (obj_name in names(objects)) {
  
  model_data <- objects[[obj_name]]
  summary_list <- list()
  
  # Loop over models within each outcome object
  for (model_name in names(model_data)) {
    
    # Extract model object
    model_obj <- model_data[[model_name]][["lm"]]
    model_summary <- summary(model_obj)
    
    # Grab coefficients table
    model_df <- as.data.frame(model_summary$coefficients)
    model_df$model <- model_name
    
    # Add R-squared
    model_df$r.squared      <- model_summary$r.squared
    model_df$adj.r.squared  <- model_summary$adj.r.squared
    
    # Add partial correlations (default NA)
    model_df$rho <- NA
    
    # RHO values 
    rho_vals  <- model_data[[model_name]][["rho"]]
    rho_names <- names(rho_vals)
    
    # Fill rho values for matching predictors
    for (i in seq_along(rho_vals)) {
      nm <- rho_names[i]
      
      if (is.null(nm) || length(nm) == 0 || nm == "") next
      
      if (nm %in% rownames(model_df)) {
        model_df[grep(nm, rownames(model_df)), "rho"] <- rho_vals[i]
      }
    }
    
    # Extract F-statistic and p-value for the model
    f_stat_value <- model_summary$fstatistic[1]
    numdf        <- model_summary$fstatistic[2]
    dendf        <- model_summary$fstatistic[3]
    
    model_df$model_pval <- pf(f_stat_value, numdf, dendf, lower.tail = FALSE)
    
    # Add sample size N
    model_df$N <- numdf + dendf + 1
    
    # store this model's summary
    summary_list[[model_name]] <- model_df
  }
  
  # Combine results for all models in this object
  summary_all <- do.call(rbind, summary_list)
  
  # Save output
  file_name <- paste0("BrainPAD_models_summary_", obj_name, "_std.csv")
  write.csv(summary_all, file_name, row.names = TRUE)
}

