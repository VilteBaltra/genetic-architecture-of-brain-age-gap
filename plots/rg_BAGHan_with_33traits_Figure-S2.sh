
### Genetic correlation plot of BAGHan and 33 traits 
##(adding Parkinson's disease)

source activate /shared/home/vb506/.conda/envs/ldsc/
cd /campaign/VB-FM5HPC-001/Vilte/ldsc

# METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22.txt.gz
# munge ENIGMA (noGenR)
./munge_sumstats.py \
  --sumstats METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-07-22.txt.gz \
  --p p_value \
  --a1 Allele1 \
  --a2 Allele2 \
  --ignore N \
  --signed-sumstats Zscore,0 \
  --snp MarkerName \
  --out enigma-noGenr \
  --merge-alleles w_hm3.snplist


# rg of brainage and PD
./ldsc.py \
    --rg enigma-noGenr.sumstats.gz,Parkinsons.sumstats.gz \
    --ref-ld-chr eur_w_ld_chr/ \
    --w-ld-chr eur_w_ld_chr/ \
    --out enigma-noGenr-cor-PD


# new MDD (pgc-mdd2025_no23andMe_eur_v3-49-24-11_formatted.tsv.gz)
# rg of enigma and MDD
./ldsc.py \
    --rg enigma-noGenr.sumstats.gz,pgc-mdd2025_no23andMe_eur.sumstats.gz \
    --ref-ld-chr eur_w_ld_chr/ \
    --w-ld-chr eur_w_ld_chr/ \
    --out enigma-noGenr.sumstats-cor-mdd2025


### Calculate FDR and PLOT ###
## FDR
import os
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

# load data
#cwd = os.getcwd()
#file_path = os.path.join(cwd, 'LDSC/results.tsv')
cwd="/Users/vb506/Documents/rg-additional-PD-MD/enigma/"
file_path = os.path.join(cwd, 'enigma-rg.tsv')
df = pd.read_csv(file_path, sep='\t')

# apply standard BH-FDR correction (Benjamini-Hochberg)
rejected, pvals_corrected = fdrcorrection(df['Pvalue'], alpha=0.05, method='indep')

# add results to DataFrame
df['FDR_BH'] = pvals_corrected
df['Significant_FDR_05'] = rejected

# sort by GeneticCorrelation
df = df.sort_values("Genetic correlation")

# save to new .tsv
#output_path = os.path.join(cwd, 'LDSC/results_with_fdr.tsv')
output_path = os.path.join(cwd, 'results_enigma_with_fdr_PD.tsv')
df.to_csv(output_path, sep='\t', index=False)



## PLOT
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

# === Define your 10 categorical colors (your custom palette) ===
colors_hex_cat = [
    '#fbe6c5', '#f7c9a7', '#f3aa91', '#ee8a82', '#e2797a',
    '#d56973', '#c8586c', '#ab4762', '#8d3757', '#70284a'
]
cmap_cat = colors_hex_cat

#file_path = cwd+'/LDSC/results_with_fdr.tsv'
file_path = cwd+'/results_enigma_with_fdr_PD.tsv'
df = pd.read_csv(file_path, sep='\t')
#df = df.sort_values("GeneticCorrelation")
print(df)

#fontsize_labels = 16
#fontsize_ticks = 16
fontsize_ticks = 14
fontsize_labels = 14
cap_length = 0.05

plt.figure(figsize=(13, 10))

y_pos = range(len(df))

# colours: significant in red, nonsig in light gray
sig_color = cmap_cat[6] # pick color from 0 to 9
nonsig_color = '#aaaaaa'  # light gray

for i, (gc, se, pval, fdr_sig) in enumerate(zip(df["Genetic correlation"], df["SE"], df["Pvalue"], df['Significant_FDR_05'])):
    color = sig_color if fdr_sig else nonsig_color
    ci95 = 1.96 * se 
    
    # plot the circle with error bar
    plt.errorbar(gc, i, xerr=ci95, fmt='o', color=color,
                 ecolor=color, elinewidth=2, capsize=4,
                 markersize=10, markeredgewidth=1.5,
                 markeredgecolor='white', alpha=0.9)

# Apply updated labels to y-axis
plt.yticks(y_pos, df['Trait'], fontsize=fontsize_ticks)
plt.xticks(fontsize=fontsize_ticks)
plt.xlabel("Genetic correlation (rg)", fontsize=fontsize_labels)
#plt.title("Genetic Correlations of Brain Age with Other Traits", fontsize=fontsize_title, weight='bold')
plt.axvline(0, color='gray', linestyle='--', linewidth=1) #vertical line on the 0

# show only vertical grid lines, subtle and light
plt.grid(axis='x', linestyle='--', color='gray', alpha=0.15, linewidth=0.7)
plt.grid(axis='y', visible=False)  # hide horizontal grid lines

sig_patch = mpatches.Patch(color=sig_color, label='FDR-adjusted p < 0.05')
nonsig_patch = mpatches.Patch(color=nonsig_color, label='FDR-adjusted p ≥ 0.05')

plt.legend(
    handles=[sig_patch, nonsig_patch],
    fontsize=fontsize_ticks,
    frameon=False,
    loc='upper left'
)

#plt.gca().set_facecolor('#f9f9f9')
plt.tight_layout()

# Save the plot to a file
plt.savefig('genetic_correlation_enigma_fdr_withPD.png', dpi=600, bbox_inches='tight')

plt.show()
