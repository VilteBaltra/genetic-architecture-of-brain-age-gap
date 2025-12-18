#!/bin/bash
#==================================================================
# LD Score Regression pipeline: munge sumstats and run genetic rg
#==================================================================

# Activate LDSC conda environment
source activate /shared/home/vb506/.conda/envs/ldsc/

# create traits.tsv file that details relevant column names for each GWAS sumstats
# traits.tsv — input for munge # defined in first script of 3.LDSC-brainageFactor

#-------------------------------
# Define directories
#-------------------------------
BASE_DIR="/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024"
SUMSTATS_DIR="${BASE_DIR}/genetic-corr/sumstats"
OUTPUT_DIR="${BASE_DIR}/genetic-corr/output2"
REF_LD="eur_w_ld_chr"

#-------------------------------
# Munge enigma
#-------------------------------
./munge_sumstats.py \
  --sumstats "${SUMSTATS_DIR}/METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22.txt.gz" \
  --p P \
  --a1 Allele1 \
  --a2 Allele2 \
  --signed-sumstats Zscore,0 \
  --ignore N \
  --snp MarkerName \
  --out "${SUMSTATS_DIR}/enigma" \
  --merge-alleles w_hm3.snplist

##-------------------------------
## Function to munge traits from TSV --> already done in 3.LDSC-brainageFactor, so commenting it out here
##-------------------------------
#munge_traits() {
#  local tsv_file="$1"
#  tail -n +2 "$tsv_file" | while IFS=$'\t' read -r infile N_col pval a1 a2 snp outname; do
#    echo "[MUNGE] $outname..."
#    args=(--sumstats "$infile" --p "$pval" --a1 "$a1" --a2 "$a2" --snp "$snp" \
#          --out "${SUMSTATS_DIR}/${outname}" --merge-alleles w_hm3.snplist)
#
#    [[ "$N_col" != "-" && -n "$N_col" ]] && args+=(--N "$N_col")
#    ./munge_sumstats.py "${args[@]}"
#  done
#}
#
## munge traits
#munge_traits traits.tsv

#-------------------------------
# Function to run LDSC rg
#-------------------------------
run_ldsc() {
  local tsv_file="$1"
  tail -n +2 "$tsv_file" | while IFS=$'\t' read -r _ _ _ _ _ _ outname; do
    echo "[LDSC] $outname..."
    ./ldsc.py \
      --rg "${SUMSTATS_DIR}/enigma.sumstats.gz,${SUMSTATS_DIR}/${outname}.sumstats.gz" \
      --ref-ld-chr "${REF_LD}/" \
      --w-ld-chr "${REF_LD}/" \
      --out "${OUTPUT_DIR}/enigma-cor-${outname}"
  done
}

# run LDSC on traits
run_ldsc traits.tsv


#-------------------------------
# Additional traits (PD, MDD)
#-------------------------------

# Parkinson's disease (already munged before)
./ldsc.py \
  --rg "${SUMSTATS_DIR}/enigma.sumstats.gz,${SUMSTATS_DIR}/Parkinsons.sumstats.gz" \
  --ref-ld-chr "${REF_LD}/" \
  --w-ld-chr "${REF_LD}/" \
  --out "${OUTPUT_DIR}/enigma-cor-PD"

# MDD (already munged before)
./ldsc.py \
  --rg "${SUMSTATS_DIR}/enigma.sumstats.gz,${SUMSTATS_DIR}/pgc-mdd2025_no23andMe_eur.sumstats.gz" \
  --ref-ld-chr "${REF_LD}/" \
  --w-ld-chr "${REF_LD}/" \
  --out "${OUTPUT_DIR}/enigma-cor-mdd2025"
