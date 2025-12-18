#!/bin/bash
#==================================================================
# LD Score Regression pipeline: munge sumstats and run genetic rg
#==================================================================

# Activate LDSC conda environment
source activate /shared/home/vb506/.conda/envs/ldsc/

# create traits.tsv file that details relevant column names for each GWAS sumstats
# traits.tsv — input for munge
infile_path    N    pval    a1    a2    snp    outname
${SUMSTATS_DIR}/oag_pheno_normalized_residualized.Cardiovascular_age_gap.glm.linear.gz    111374    P    ALT    REF    ID    wen-cardiovascular
${SUMSTATS_DIR}/oag_pheno_normalized_residualized.Immune_age_gap.glm.linear.gz    111276    P    ALT    REF    ID    wen-immune
${SUMSTATS_DIR}/oag_pheno_normalized_residualized.Metabolic_age_gap.glm.linear.gz    111175    P    ALT    REF    ID    wen-metabolic
${SUMSTATS_DIR}/oag_pheno_normalized_residualized.Musculoskeletal_age_gap.glm.linear.gz    111354    P    ALT    REF    ID    wen-musculoskeletal
${SUMSTATS_DIR}/oag_pheno_normalized_residualized.Pulmonary_age_gap.glm.linear.gz    111354    P    ALT    REF    ID    wen-pulmonary
${SUMSTATS_DIR}/oag_pheno_normalized_residualized.Renal_age_gap.glm.linear.gz    111354    P    ALT    REF    ID    wen-renal
${SUMSTATS_DIR}/alzheimer-35379992-GCST90027158-MONDO_0004975.h.tsv.gz    487511    p_value    effect_allele    other_allele    variant_id    Bellenguez-alzheimer
${SUMSTATS_DIR}/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz    N    PVAL    A1    A2    SNP    PGC3_SCZ_wave3.european.autosome
${SUMSTATS_DIR}/daner_bip_pgc3_nm_noukbiobank_beta.txt    353899    p_value    A1    A2    SNP    daner_bip_pgc3_nm_noukbb
${SUMSTATS_DIR}/iPSYCH-PGC_ASD_Nov2017.gz    225534    P    A1    A2    SNP    ASD
${SUMSTATS_DIR}/ADHD2022_iPSYCH_deCODE_PGC.meta.gz    36351    P    A1    A2    SNP    ADHD
${SUMSTATS_DIR}/bmi-25673413-GCST002783-EFO_0004340.h.tsv.gz    238944    p_value    effect_allele    other_allele    hm_rsid    bmi-GCST002783
${SUMSTATS_DIR}/longevity_90th_percentile_rsid.txt.gz    36745    P-value    EA    NEA    rsid    longevity_90th_percentile
${SUMSTATS_DIR}/crp_filtered_sumstats.tsv    500000    P    effect_allele    other_allele    variant_id    crp
${SUMSTATS_DIR}/Mahajan.NatGenet2018b.T2D.European_mapped_CHR_ALL.txt.gz    226311    P    A1    A2    RSID    t2d
${SUMSTATS_DIR}/cad.add.160614.website.txt    185000    p_dgc    effect_allele    noneffect_allele    markername    cad
${SUMSTATS_DIR}/sleepdurationsumstats.txt    446118    P    ALLELE1    ALLELE0    SNP    sleep-duration
${SUMSTATS_DIR}/RHR_UKBplusICRHR_without23andme.gz    835465    p_value    A1    A2    SNP_1000G    resting-heart-rate
${SUMSTATS_DIR}/SmokingInitiation.txt.gz    1232091    PVALUE    ALT    REF    RSID    smoking-initiation
${SUMSTATS_DIR}/MAGIC1000G_HbA1c_EUR.tsv    281416    p_value    effect_allele    other_allele    variant    hba1c
${SUMSTATS_DIR}/CigarettesPerDay.txt.gz    337334    PVALUE    ALT    REF    RSID    CigarettesPerDay
${SUMSTATS_DIR}/SBP_GCST90310294.h.tsv.gz    -    p_value    effect_allele    other_allele    rsid    sbp
${SUMSTATS_DIR}/DBP_GCST90310295.h.tsv.gz    -    p_value    effect_allele    other_allele    rsid    dbp
${SUMSTATS_DIR}/PP_GCST90310296.h.tsv.gz    -    p_value    effect_allele    other_allele    rsid    pp
${SUMSTATS_DIR}/Gran_EUR_summary_statistics.txt    N    P    A1    A2    rsID    Gran_EUR
${SUMSTATS_DIR}/GrimAge_EUR_summary_statistics.txt    N    P    A1    A2    rsID    GrimAge_EUR
${SUMSTATS_DIR}/Hannum_EUR_summary_statistics.txt    N    P    A1    A2    rsID    Hannum_EUR
${SUMSTATS_DIR}/IEAA_EUR_summary_statistics.txt    N    P    A1    A2    rsID    IEAA_EUR
${SUMSTATS_DIR}/PAI1_EUR_summary_statistics.txt    N    P    A1    A2    rsID    PAI1_EUR
${SUMSTATS_DIR}/PhenoAge_EUR_summary_statistics.txt    N    P    A1    A2    rsID    PhenoAge_EUR

#-------------------------------
# Define directories
#-------------------------------
BASE_DIR="/campaign/VB-FM5HPC-001/Vilte/Projects/BrainAgeGWAS_2024"
SUMSTATS_DIR="${BASE_DIR}/genetic-corr/sumstats"
OUTPUT_DIR="${BASE_DIR}/genetic-corr/output"
REF_LD="eur_w_ld_chr"

#-------------------------------
# Munge brain age factor
#-------------------------------
./munge_sumstats.py \
  --sumstats "${SUMSTATS_DIR}/brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz" \
  --p Pval_Estimate \
  --a1 A1 \
  --a2 A2 \
  --signed-sumstats Z_Estimate,0 \
  --snp SNP \
  --out "${SUMSTATS_DIR}/brainageFactor" \
  --merge-alleles w_hm3.snplist

#-------------------------------
# Function to munge traits from TSV
#-------------------------------
munge_traits() {
  local tsv_file="$1"
  tail -n +2 "$tsv_file" | while IFS=$'\t' read -r infile N_col pval a1 a2 snp outname; do
    echo "[MUNGE] $outname..."
    args=(--sumstats "$infile" --p "$pval" --a1 "$a1" --a2 "$a2" --snp "$snp" \
          --out "${SUMSTATS_DIR}/${outname}" --merge-alleles w_hm3.snplist)

    [[ "$N_col" != "-" && -n "$N_col" ]] && args+=(--N "$N_col")
    ./munge_sumstats.py "${args[@]}"
  done
}

# munge traits
munge_traits traits.tsv

#-------------------------------
# Function to run LDSC rg
#-------------------------------
run_ldsc() {
  local tsv_file="$1"
  tail -n +2 "$tsv_file" | while IFS=$'\t' read -r _ _ _ _ _ _ outname; do
    echo "[LDSC] $outname..."
    ./ldsc.py \
      --rg "${SUMSTATS_DIR}/brainageFactor.sumstats.gz,${SUMSTATS_DIR}/${outname}.sumstats.gz" \
      --ref-ld-chr "${REF_LD}/" \
      --w-ld-chr "${REF_LD}/" \
      --out "${OUTPUT_DIR}/brainageFactor-cor-${outname}"
  done
}

# run LDSC on traits
run_ldsc traits.tsv


#-------------------------------
# Additional traits (PD, MDD)
#-------------------------------

# Parkinson's disease
./munge_sumstats.py \
  --sumstats GP2_PD_with_Kaviar_IDs_unique_N.gz \
  --p p_value \
  --N-cas-col ncases \
  --N-con-col ncontrols \
  --a1 effect_allele \
  --a2 other_allele \
  --signed-sumstats beta,0 \
  --snp ID \
  --out "${SUMSTATS_DIR}/Parkinsons" \
  --merge-alleles w_hm3.snplist

./ldsc.py \
  --rg "${SUMSTATS_DIR}/brainageFactor.sumstats.gz,${SUMSTATS_DIR}/Parkinsons.sumstats.gz" \
  --ref-ld-chr "${REF_LD}/" \
  --w-ld-chr "${REF_LD}/" \
  --out "${OUTPUT_DIR}/brainageFactor-cor-PD"

# MDD
./munge_sumstats.py \
  --sumstats pgc-mdd2025_no23andMe_eur_v3-49-24-11_formatted.tsv.gz \
  --p p_value \
  --N-cas-col ncases \
  --N-con-col ncontrols \
  --a1 effect_allele \
  --a2 other_allele \
  --signed-sumstats beta,0 \
  --snp rsid \
  --out "${SUMSTATS_DIR}/pgc-mdd2025_no23andMe_eur" \
  --merge-alleles w_hm3.snplist

./ldsc.py \
  --rg "${SUMSTATS_DIR}/brainageFactor.sumstats.gz,${SUMSTATS_DIR}/pgc-mdd2025_no23andMe_eur.sumstats.gz" \
  --ref-ld-chr "${REF_LD}/" \
  --w-ld-chr "${REF_LD}/" \
  --out "${OUTPUT_DIR}/brainageFactor-cor-mdd2025"
