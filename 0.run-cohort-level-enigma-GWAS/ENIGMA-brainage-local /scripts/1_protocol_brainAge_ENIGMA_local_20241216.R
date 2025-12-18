#### Brain Age Project
# Authors: Laura Han (L.Han@ggzingeest.nl) and Esther Walton (waltonesther@gmail.com)

# --- Setup ---

# Set working directory (optional, uncomment if needed)
# setwd("your/directory/path")

# Create a log file
log_file <- file("1_BrainAGE.log", open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

# --- Install and Load Libraries ---

cat("Prep: Installing/Loading Libraries\n")

load_libraries <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

required_packages <- c(
  "ppcor", "lsmeans", "multcomp", "data.table", "plyr", "ModelMetrics",
  "caret", "gridExtra", "Hmisc", "pastecs", "psych", "ggplot2", "dplyr", "testit"
)

load_libraries(required_packages)

# --- Source Functions ---

cat("Prep: Sourcing Functions\n")
source("functions/get.means.R")
source("functions/prepare.files.R")
source("functions/model.fits.brainAge.R")
source("functions/model.fits.brainAge.dx.R")

# --- Step 2: Deriving Brain Measures ---

cat("Step 2: Obtaining Measures of Brain Age\n")

# Load mean volumes, thickness, and surface area data
Thick <- get.means("CorticalMeasuresENIGMA_ThickAvg.csv")
Thick$ICV <- NULL

Surf <- get.means("CorticalMeasuresENIGMA_SurfAvg.csv")
Surf$ICV <- NULL

Vol <- get.means("SubcorticalMeasuresENIGMA_VolAvg.csv")

# Merge data
TS <- merge(Thick, Surf, by = "row.names")
TSV <- merge(TS, Vol, by.x = "Row.names", by.y = "row.names")

# --- Read and Validate Covariates ---

Covs <- read.csv("Covariates.csv")
names(Covs) <- toupper(names(Covs))

# Select PCs
PCs <- grep("PC", names(Covs), value = TRUE)

# Check for missing data
if (all(is.na(Covs$DX))) {
  cat("No case/control data detected. Skipping case/control analysis.\n")
  if (any(is.na(Covs[, c("SEX", "AGE", PCs)]))) {
    stop("Missing data in SEX, AGE, or PCs. Check Covariates.csv file.\n")
  }
} else {
  cat("Case/control data detected. Applying relevant checks.\n")
  if (any(is.na(Covs[, c("DX", "SEX", "AGE", PCs)]))) {
    stop("Missing data in DX, SEX, AGE, or PCs. Check Covariates.csv file.\n")
  }
}

# Validate required columns
required_columns <- c(
  "SUBJID", "DX", "SEX", "AGE", PCs, "ETHNICITY", "ANCESTRY", "SES", "IQ", 
  "HAND", "LENGTH_OF_ILLNESS", "DISEASE_CATEGORY", "SUBTYPE"
)
missing_columns <- required_columns[!required_columns %in% names(Covs)]
if (length(missing_columns) > 0) {
  stop("Missing required columns in Covariates.csv: ", paste(missing_columns, collapse = ", "), "\n")
}

# Check for case/control data consistency
if (!all(is.na(Covs$DX))) {
  HC <- Covs[Covs$DX == 0, ]
  if (!all(is.na(HC[, c("LENGTH_OF_ILLNESS", "DISEASE_CATEGORY", "SUBTYPE")]))) {
    stop("Control subjects have entries for patient-only variables. Please code as 'NA'.\n")
  }
  rm(HC)
}

# create dummy variable for SITE
if ("SITE" %in% names(Covs)){
  SITE.f = factor(Covs$SITE)
  #Covs$SITE=NULL
  
  if (length(grep("SITE.f",colnames(Covs)))>0){
    site_regr=paste(grep("SITE.f",colnames(Covs),value=T),collapse="+")
  }
  
  Covs=cbind(Covs,model.matrix(~SITE.f)[,-1])
}
# Add AGE^2 for modeling
Covs$AGE2 <- Covs$AGE^2

# Exclude subjects outside age range
age_filter <- Covs$AGE >= 18 & Covs$AGE <= 75
if (!all(age_filter)) {
  cat("Excluding subjects younger than 18 or older than 75 years.\n")
  Covs <- Covs[age_filter, ]
}

# Check for duplicate SUBJIDs
if (anyDuplicated(Covs$SUBJID)) {
  stop("Duplicate SUBJIDs detected in Covariates.csv. Ensure unique IDs.\n")
}

# --- Combine Data ---

data <- merge(Covs, TSV, by.x = "SUBJID", by.y = "Row.names")
save(data, Covs, TSV, file = "data_covs_tsv.Rdata") # save as will be needed in script 3

cat("Creating CSV files for brainAge estimation\n")

# Prepare data for export
df <- prepare.files(data, names(Covs))
saveRDS(df, file="df_covs_tsv.rds") # save as will be needed in script 3
males <- df$males
females <- df$females

# Clean up temporary data
rm(df)

# Prompt user for next steps
invisible(readline(prompt = paste(
  "You should see 2 CSV files in your working directory:\n",
  "- females_raw.csv\n",
  "- males_raw.csv\n",
  "Upload/download these files as instructed, then run the next script.\n"
)))

cat("Finished with Step 1!\n")

# Close log and print contents
sink()
closeAllConnections()
print(readLines("1_BrainAGE.log"))
