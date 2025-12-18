# This script runs genomic SEM on 6 BAG GWASs to get the BAG factor

# before opening R in linux run below line (will make parallel analysis work and therefore reduce run time):
# export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# Wen et al. (2024; brain BAG), Smith et al. (2020; all-in-one IDPs), Kaufmann et al. (2019; brainAge),
# Wen: https://labs-laboratory.com/medicine/multimodal_brain_bag
# downloaded GWAS sumstats from https://labs-laboratory.com/medicine/multimodal_brain_bag (from the 9 organs paper, as more recent one than brain BAG only paper!)

# Kaufman
#git clone https://github.com/tobias-kaufmann/brainage.git
#cd brainage
#git checkout 51e650d109347eddb420de6242fabac0485695b4
# concatanate the files
#cat Brainage_GWAS_sumstat_final.gz.* > Brainage_GWAS_sumstat_final.gz
#These .001, .002, .003 files aren’t split text files (where each part might have its own header).
#They’re binary chunks of a single compressed .gz archive — created by a tool like 7-Zip split archive.
#When you concatenate them, you’re essentially reconstructing the original .gz file — and then decompressing it as a whole.
#The decompression happens after rejoining, so the final .txt file will be exactly as it was before splitting — with a single header at the top.
#mv Brainage_GWAS_sumstat_final.gz Kaufman_sumstats.gz

# FORMAT NOTE
# had to remove "#" from header to make this run for oag_pheno_normalized_residualized.Brain_age_gap.glm.linear
# 0 snps were left in Wen.sumstats.gz, ALT was treated as A2, so renamed REF to A2 and ALT to ALT_original (A1 is also present, which is the same as ALT)


#install.packages("devtools")
library(devtools)
#install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
library(dplyr)

# GenomicSEM for Common Factor GWAS
# Step 1: Munge the summary statistics
# setwd

# Specify individual BAG GWAS sumstats
files<-c("METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22.txt.gz", "brainage2025.full.eur.excl2k.gwm.gz", "leonardsen.txt.gz","oag_pheno_normalized_residualized.Brain_age_gap.glm.linear.gz", "Kaufman_sumstats.gz")
ref= "reference.1000G.maf.0.005.txt"
trait.names<-c("ENIGMAsamplesize", "Jawinski_excl2k", "Pyment", "Wen", "Kaufman")
se.logit=c(F,F,F,F,F)
linprob=c(F,F,F,F,F)
OLS=c(T,T,T,T,T)
N=c(60735,52890,28104,30108,20170)
betas=NULL
#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3
hm3<-"eur_w_ld_chr/w_hm3.snplist"
#definte the imputation quality filter
info.filter=0.9
#define the MAF filter
maf.filter=0.01

#run munge
# remember to format the files so munge finds the right columns
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)


# Smith62modes with A2 swapped for A1 (as effect seems to be reversed)

#gunzip -cd Smith/Smith_V0140_sumstats.gz | \
#awk 'NR==1{for(i=1;i<=NF;i++){if($i=="a1")$i="a2"; else if($i=="a2")$i="a1"}; print; next} 1' \
#> Smith/Smith_V0140_sumstats_A2swap.txt

files<-c("Smith/Smith_V0140_sumstats_A2swap.txt") 
ref= "reference.1000G.maf.0.005.txt"
trait.names<-c("Smith_62modes_A2swap") # "Smith" to be added
info.filter=.6
maf.filter=0.01
OLS=T
N=15952 
betas=NULL
hm3<-"eur_w_ld_chr/w_hm3.snplist"
#definte the imputation quality filter
info.filter=0.9
#define the MAF filter
maf.filter=0.01

munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)



# ldsc

#vector of munged summary statisitcs
traits<-c("ENIGMAsamplesize.sumstats.gz","Pyment.sumstats.gz","Jawinski_excl2k.sumstats.gz", "Wen.sumstats.gz", "Kaufman.sumstats.gz", "Smith_62modes_A2swap.sumstats.gz")

#enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
sample.prev<-NA # value should be NA for continuous traits

#vector of population prevalences
population.prev<-NA # same as above

#the folder of LD scores
ld<-"eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld<-"eur_w_ld_chr/"

#name the traits
trait.names<-c("ENIGMAsamplesize","Pyment","Jawinski_excl2k", "Wen", "Kaufman", "Smith_62modes")

#run LDSC
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)

#optional command to save the output as a .RData file for later use
# updated, no GenR:
# save(LDSCoutput,file="LDSCoutput_brainage_samplesize_20250723.RData")

# Step 3: Specify and Estimate a Structural Equation Model

#To run using DWLS estimation#
CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")

#print CommonFactor_DWLs output#
CommonFactor_DWLS

# EFA for 2 factors
#smooth the S matrix for EFA using the nearPD function in the Matrix package.
require(Matrix)
Ssmooth<-as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)

#run EFA with promax rotation and 2 factors using the factanal function in the stats package
require(stats)
EFA<-factanal(covmat = Ssmooth, factors = 2, rotation = "promax")
EFA

#EFA<-factanal(covmat = Ssmooth, factors = 1, rotation = "promax")
#EFA


#Specify the Genomic confirmatory factor model
CFAofEFA <- 'F1 =~ NA*Pyment + Jawinski_excl2k + Wen + Smith_62modes
             F2 =~ NA*ENIGMAsamplesize + Kaufman 
F1~~F2
Kaufman ~~ a*Kaufman
a > .001' 
# Kaufman had negative residual variance, so constraining to be positive

#run the model
out<-usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)

#print the results
out


# common factor GWAS 
# for parallel computing: export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

files<-c("METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22.txt.gz", "leonardsen.txt.gz", "brainage2025.full.eur.excl2k.gwm.gz","oag_pheno_normalized_residualized.Brain_age_gap.glm.linear.gz", "Kaufman_sumstats.gz", "Smith/Smith_V0140_sumstats_A2swap.txt")
ref= "reference.1000G.maf.0.005.txt.gz" 
trait.names<-c("ENIGMAsamplesize","Pyment","Jawinski_excl2k", "Wen", "Kaufman", "Smith_62modes") # "Smith" to be added
se.logit=c(F,F,F,F,F,F)
linprob=c(F,F,F,F,F,F)
info.filter=.6
maf.filter=0.01
OLS=c(T,T,T,T,T,T)
N=c(60735,28104,52890,30108,20170,15952)
betas=NULL

brainage_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=OLS,linprob=linprob,N=N,betas=NULL,info.filter=info.filter,maf.filter=maf.filter,keep.indel=FALSE,parallel=F,cores=6)

dim(brainage_sumstats)

saveRDS(brainage_sumstats, file="brainage_sumstats_20250727.rds")

#run the multivariate GWAS using parallel processing (might need to handle that heywood case)
# export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

brainage_sumstats <- readRDS("brainage_sumstats_20250727.rds")
load("LDSCoutput_brainage_samplesize_20250723.RData")

# run GWAS
start_time <- Sys.time()
brainage_factor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = brainage_sumstats, estimation = "DWLS", cores = 30, toler = FALSE, SNPSE = FALSE, parallel = TRUE,GC="none",MPI=FALSE)
end_time <- Sys.time()

# calculate runtime in minutes
runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
print(paste("Runtime:", round(runtime_minutes, 2), "minutes"))

# ensure its df
brainage_factor_df <- as.data.frame(brainage_factor)

# save summary statistics
write.table(brainage_factor_df, "brainage_factor_withSmith_part1_noGenR_noGC.txt", row.names=FALSE, sep="\t", quote=FALSE)


# ensure its df
brainage_factor_df <- as.data.frame(brainage_factor)

# get effective sample size

#restrict to MAF >= 10%
brainage_factor_df2<-subset(brainage_factor_df, brainage_factor_df$MAF >= .1)
#calculate expected sample size (N_hat)
N_hat<-mean(1/((2*brainage_factor_df2$MAF*(1-brainage_factor_df2$MAF))*brainage_factor_df2$se_c^2))
print(N_hat)
# 160656.1

# add N
brainage_factor_df$N <- 160656
head(brainage_factor_df)

# save dataset
write.table(brainage_factor_df, file = gzfile("brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz"), row.names = FALSE, sep = "\t", quote = FALSE)
