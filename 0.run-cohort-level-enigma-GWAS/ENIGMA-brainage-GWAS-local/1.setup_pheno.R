# this script is adapted from Dirk Smit
# ___________________________________________________#
#           SET PATHS
# ___________________________________________________#

PHENOPATH = "/your/path/to/BrainAgeGWAS/folder/" # replace with your path to 'data_pheno_covariates.Rdata' file
setwd(PHENOPATH)

# ___________________________________________________#
#           FORMAT FILES FOR GWAS STAGE
# ___________________________________________________#

# read in the dataset with the brain age phenotype and the covariates 
load("data_pheno_covariates.Rdata") # this file was created with '1.pheno-covars-preparation.R' script

# add a first column with FID_IID (genetic family and individual identifiers) 
# note that:
# (1) FID (family ID) and IID (individuals ID) might be two distinct columns in your cohort, and
# (2) the IID might be identical to SUBJID in the Covariates.csv file. If it differs from SUBJID, 
# you will most likely need to read in a linking file which will permit you to match phenoypic SUBJID
# with genetic IID.

# for example, SUBJID for ALSPAC cohort is identical to the FID and IID columns in the genetic file
# so we can simply copy SUBJID to FID and IID and paste them into one variable separated by an underscore (e.g., S12_S12)
data$FID <- data$SUBJID
data$IID <- data$SUBJID
merged <- data.frame(FID_IID = paste0(data$FID, "_", data$IID), data) # ONLY run this line if your family id (FID) and individual id (IID) are separate columns
# if FID and IID are already joined in your dataset but have a different name simply rename it to FID_IID

# Save pheno.txt file with FID_IID as first column 
write.table(merged[, c("FID_IID", "devAge", "predAge")], file = "pheno.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Save covars.txt file with FID_IID as first column 
# first recode SEX into pedsex so that 1 = male, 2 = female (before males were 0, females were 1)
merged$pedsex <- ifelse(merged$SEX == 1,2,1)

cat("SELECT RELEVANT COVARIATES:
You may need to remove 'SITE' (scanning site) and 'DX' (case-control status), if not applicable to your cohort.")
# select all covariates: sex, age, age2, total intracranial volume, SITE (scanning site, if applicable), DX (case-control status, if applicable)
# if 'SITE' and/or 'DX' are not applicable to your cohort, simply delete them from the line below 
write.table(merged[, c("FID_IID", 'pedsex', 'SEX', 'AGE', 'AGE2', 'ICV', 'SITE', 'DX')], file = "covars.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Save PCs.txt file with FID_IID as first column and 4 PCs 
write.table(merged[, c("FID_IID", "PC1", "PC2","PC3","PC4")], file = "PCs.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)

covar = read.table('covars.txt',header=T,stringsAsFactors=F)

colnames(covar) = tolower(colnames(covar))

# Check if sex is the first covariate, otherwise move it to the first position. This is 
# because the Sex variable takes a special position in the PED file.
idx_id = colnames(covar) %in% c('fid_iid')
idx_sex = colnames(covar) %in% c('sex')
if (sum(idx_id) != 1 || sum(idx_sex) != 1) {
  cat('Wrong number of variables named "fid_iid" or "sex" in the data\n')
  stop()
}
covar = cbind(covar[,which(idx_id)],covar[,which(idx_sex)], covar[,which(!(idx_id | idx_sex))])
colnames(covar)[1] = 'fid_iid'
colnames(covar)[2] = 'sex'

# genoPCs may be larger or smaller (or just different) 
genoPCs = read.table('PCs.txt',header=T,stringsAsFactors=F)
colnames(genoPCs) = tolower(colnames(genoPCs))

# genetic IDs
genoID = genoPCs[,1]

# select phenotype file
cPhenoStr = 'pheno.txt'

pheno = read.table(cPhenoStr,header=T,stringsAsFactors=F);
colnames(pheno) = tolower(colnames(pheno))
# first create a G matrix with genoIDs and start adding stuff
G = as.data.frame(cbind(genoID,genoID))
colnames(G) = c('fid_iid','fid_iid2')
G$parent1 = 0
G$parent2 = 0

# The use of join instead of merge makes sure the original order remains. Standard
# join is a left join (on G) adding NAs for the nonmatching entries. First add the sex 
# variable as column 5. THIS WILL NOT BE INCLUDED AS A COVARITE so sex needs to appear TWICE
# Use install.packages('plyr') if unavailable

M = plyr::join(G,as.data.frame(covar[,c('fid_iid', 'pedsex')]),by='fid_iid')

# subsequently add phenotypes
nvar = dim(pheno)[2]-1
M2 = plyr::join(M,as.data.frame(pheno),by='fid_iid')
ncov = ncol(covar)-1 # minus 'fid_iid'
M3 = plyr::join(M2,as.data.frame(covar[,c(1,2, 4:(ncov+1))]),by='fid_iid')
PED = plyr::join(M3,genoPCs,by='fid_iid')

# Create a list of individuals that have genetic and phenotypic data available
# and save them in 'KeepList.txt'
ndx = PED$fid_iid %in% covar$fid_iid
tmp = t(as.data.frame(strsplit(PED$fid_iid[ndx],"_")))
colnames(tmp) = c('FID','IID')
rownames(tmp) = c()
write.table(tmp,file='KeepList.txt',col.names=F,row.names=F,quote=F)

# write PED file with devage phenotype
ped_nopredage = subset(PED[ndx,], select = -c(predage))
write.table(ped_nopredage,file='devage_main.ped',sep=" ",col.names=F,row.names=F,quote=F)
# write PED file with predage phenotype
ped_nodevage = subset(PED[ndx,], select = -c(devage))
write.table(ped_nodevage,file='predage_main.ped',sep=" ",col.names=F,row.names=F,quote=F)

# create the DAT file data structure
DAT = as.data.frame(c(rep('T',1),rep('C',ncov+ncol(genoPCs)-2)))
colnames(DAT) = 'c1'
DAT$c2 = colnames(PED)[c(6, 8:ncol(PED))] # selects devage + covars (skips predage)
write.table(DAT,file='devage_main.dat',sep=" ",col.names=F,row.names=F,quote=F)
DAT$c2 = colnames(PED)[c(7:dim(PED)[2])] # selects predage + covars
write.table(DAT,file='predage_main.dat',sep=" ",col.names=F,row.names=F,quote=F)

# create another PED and DAT file without ICV
ped_noICV = subset(ped_nopredage, select = -c(icv))
write.table(ped_noICV,file='devage_noICV.ped',sep=" ",col.names=F,row.names=F,quote=F)
DAT[1, 2] = 'devage' # replace predage with devage
write.table(DAT[-5,],file='devage_noICV.dat',sep=" ",col.names=F,row.names=F,quote=F) # -5 removes ICV row from file

## This process will result in seven files: 
# devage_main.ped, devage_main.dat, predage_main.ped, predage_main.dat, devage_noICV.ped, devage_noICV.dat, and Keeplist.txt. 

## devage_main.dat file should list the trait (T) and the covariates (C) as below:
# T devage
# C sex
# C dx
# C icv
# C age
# C age2
# C pc1
# C pc2
# C pc3
# C pc4

## phenotype.ped file contains the values for the traits and covariates:
# 5626A_5626A 5626A_5626A 0 0 1 0.26 1 0 23.58 556.17 0.01 0.02 -0.02 0.02

## KeepList.txt stores IDs:
# S16_S16 S16_S16
# S17_S17 S17_S17

# ___________________________________________________#
#     FOR COHORTS THAT HAVE AT LEAST 50 FEMALES
# ___________________________________________________#

## Create separate .ped and .dat files for females
# Subset to females (females should be coded as 1 for sex and 2 for pedsex)
ped = PED[ndx,]
ped_females = ped[ped$sex == 1, ]
ped_females = subset(ped_females, select = -c(sex)) # remove second instance of sex

# Write the PED file with devage phenotype 
ped_females = subset(ped_females, select = -c(predage))
write.table(ped_females, file='devage_females.ped',sep=" ",col.names=F,row.names=F,quote=F)

# create the DAT file data structure
DAT = as.data.frame(c(rep('T',1),rep('C',ncov-1+ncol(genoPCs)-2)))
colnames(DAT) = 'c1'
DAT$c2 = colnames(ped_females)[6:ncol(ped_females)]
write.table(DAT,file='devage_females.dat',sep=" ",col.names=F,row.names=F,quote=F)

# ___________________________________________________#
#     FOR COHORTS THAT HAVE AT LEAST 50 MALES
# ___________________________________________________#

## Create separate .ped and .dat files for males
# Subset to males (males should be coded as 0)
ped = PED[ndx,]
ped_males = ped[ped$sex == 0, ]
ped_males = subset(ped_males, select = -c(sex))

# Write the PED file with devage and predage phenotypes (controls only)
ped_males = subset(ped_males, select = -c(predage))
write.table(ped_males, file='devage_males.ped',sep=" ",col.names=F,row.names=F,quote=F)

# create the DAT file data structure
DAT = as.data.frame(c(rep('T',1),rep('C',ncov-1+ncol(genoPCs)-2)))
colnames(DAT) = 'c1'
DAT$c2 = colnames(ped_males)[6:ncol(ped_males)]
write.table(DAT,file='devage_males.dat',sep=" ",col.names=F,row.names=F,quote=F)


# ___________________________________________________#
#          FOR CASE-CONRTOL COHORTS ONLY
# ___________________________________________________#

## The section below is for case-control cohorts only 
## it will create the same three files as above, but this time for controls only

# Subset to controls 
ped = PED[ndx,]
ped_controls = ped[ped$dx == 0, ]
ped_controls = subset(ped_controls, select = -c(dx))

# Write the PED file with devage and predage phenotypes (controls only)
ped_controls = subset(ped_controls, select = -c(predage))
write.table(ped_controls, file='devage_controls.ped',sep=" ",col.names=F,row.names=F,quote=F)

# create the DAT file data structure
DAT = as.data.frame(c(rep('T',1),rep('C',ncov-1+ncol(genoPCs)-2)))
colnames(DAT) = 'c1'
DAT$c2 = colnames(ped_controls)[6:ncol(ped_controls)]
write.table(DAT,file='devage_controls.dat',sep=" ",col.names=F,row.names=F,quote=F)

cat("*.ped and *.dat have been derived.\n") 
# end script #
