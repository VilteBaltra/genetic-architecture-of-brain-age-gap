# https://github.com/bulik/ldsc

# this script first munges each of the 6 individual BAG GWASs and obtains heritability estimates
# then runs LDSC regression among the 6 BAG GWASs
#git clone https://github.com/bulik/ldsc.git
cd ldsc

#conda env create --file environment.yml
conda activate ldsc

# specify dirs
wd='path/to/working/dir'
sumstats=${wd}/summary-stats-brainage


# enigma METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22.txt.gz
./munge_sumstats.py \
--sumstats ${sumstats}/METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22.txt.gz \
--snp MarkerName \
--out ${wd}/general-LDSC/results/enigma \
--merge-alleles w_hm3.snplist

# heritability
./ldsc.py \
--h2 ${wd}/general-LDSC/results/enigma.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ${wd}/general-LDSC/results/enigma_h2

# pyment leonardsen.txt.gz
./munge_sumstats.py \
--sumstats ${sumstats}/leonardsen.txt.gz \
--out ${wd}/general-LDSC/results/pyment \
--merge-alleles w_hm3.snplist

# heritability
./ldsc.py \
--h2 ${wd}/general-LDSC/results/pyment.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ${wd}/general-LDSC/results/pyment_h2

#Kaufmann
./munge_sumstats.py \
--sumstats ${sumstats}/Kaufman_sumstats.gz \
--N 20170 \
--out ${wd}/general-LDSC/results/kaufman \
--merge-alleles w_hm3.snplist

# heritability
./ldsc.py \
--h2 ${wd}/general-LDSC/results/kaufman.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ${wd}/general-LDSC/results/kaufman_h2

#Wen brainBAG 
./munge_sumstats.py \
--sumstats ${sumstats}/oag_pheno_normalized_residualized.Brain_age_gap.glm.linear \
--N 30108 \
--snp ID \
--out ${wd}/general-LDSC/results/brainBAG \
--merge-alleles w_hm3.snplist

# heritability
./ldsc.py \
--h2 ${wd}/general-LDSC/results/brainBAG.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ${wd}/general-LDSC/results/brainBAG_h2

# Jawinski
./munge_sumstats.py \
--sumstats ${sumstats}/brainage2025.full.eur.excl2k.gwm.gz \
--snp ID \
--out ${wd}/general-LDSC/results/jawinski \
--merge-alleles w_hm3.snplist

# heritability
./ldsc.py \
--h2 ${wd}/general-LDSC/results/jawinski.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ${wd}/general-LDSC/results/jawinski_h2

# Smith
./munge_sumstats.py \
--sumstats ${sumstats}/Smith_V0140_sumstats_A2swap.txt \
--snp rsid \
--ignore Z \
--N 10612 \
--out ${wd}/general-LDSC/results/smith_V0140 \
--merge-alleles w_hm3.snplist

# heritability
./ldsc.py \
--h2 ${wd}/general-LDSC/results/smith_V0140.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out ${wd}/general-LDSC/results/smith_V0140_h2


### GENETIC CORRELATIONS AMONG 6 BRAIN AGE GAPS

# Define the list of sumstats (no need for .sumstats.gz here)
bash
#!/bin/bash

# Specify dirs
wd='/Users/vb506/Documents/Projects/BrainAge-GWAS/updated-workflow/scripts-2025/SCRIPTS/4.Genetic-correlations/2025'
sumstats="${wd}/summary-stats-brainage"

# Specify traits (no commas!)
traits=("enigma" "pyment" "kaufman" "brainBAG" "jawinski" "smith_V0140")

# Loop through all pairwise combinations
for ((i=0; i<${#traits[@]}; i++)); do
  for ((j=i; j<${#traits[@]}; j++)); do
    trait1=${traits[i]}
    trait2=${traits[j]}
    
    echo "Running rg for ${trait1} and ${trait2}..."

    ./ldsc.py \
      --rg "${wd}/general-LDSC/results/${trait1}.sumstats.gz","${wd}/general-LDSC/results/${trait2}.sumstats.gz" \
      --ref-ld-chr eur_w_ld_chr/ \
      --w-ld-chr eur_w_ld_chr/ \
      --out "${wd}/general-LDSC/results/${trait1}-${trait2}_rg"
  done
done


#### 2. Add all the output to a single file (in bash)
# my edited attempt
cd $wd/general-LDSC/results

# Cleaned output file
OUTFILE="all.rg.cleaned.tsv"

# Clear output file if it exists
> "$OUTFILE"

# Header flag
HEADER_WRITTEN=false

# Loop through all .log files
for FILE in *_rg.log; do
    # Extract the 3 lines of interest
    BLOCK=$(grep -A 2 "Summary of Genetic Correlation Results" "$FILE")

    # Extract the header only once, replace multiple spaces with a single tab
    if [ "$HEADER_WRITTEN" = false ]; then
        # Extract 2nd line and convert spaces to tabs, trim edges
        HEADER_LINE=$(echo "$BLOCK" | sed -n '2p' | tr -s ' ' | sed 's/^ //;s/ $//' | tr ' ' '\t')
        echo -e "$HEADER_LINE" >> "$OUTFILE"
        HEADER_WRITTEN=true
    fi

    # Extract and clean the data line
    DATA=$(echo "$BLOCK" | sed -n '3p' | \
        sed -E 's@.*/@@g' | \
        tr -s ' ' | sed 's/^ //;s/ $//' | tr ' ' '\t')

    echo -e "${FILE}\t${DATA}" >> "$OUTFILE"
done


