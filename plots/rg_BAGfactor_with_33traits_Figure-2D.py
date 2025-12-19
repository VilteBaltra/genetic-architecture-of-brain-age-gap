### This script first applies FDR correction to LDSC output and then plots the results
## Adding FDR correction to genetic correlations
import os
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# load data
#cwd = os.getcwd()
#file_path = os.path.join(cwd, 'LDSC/results.tsv')
cwd="/Volumes/files/Psychology/ResearchProjects/EWalton/BrainHealth/BrainAge_GWAS/output/all_plots/forPlotting_Shir_final/"
file_path = os.path.join(cwd, 'LDSC/results.tsv')
df = pd.read_csv(file_path, sep='\t')

# apply standard BH-FDR correction (Benjamini-Hochberg)
rejected, pvals_corrected = fdrcorrection(df['Pvalue'], alpha=0.05, method='indep')

# add results to DataFrame
df['FDR_BH'] = pvals_corrected
df['Significant_FDR_05'] = rejected

# sort by GeneticCorrelation
df = df.sort_values("GeneticCorrelation")

# save to new .tsv
#output_path = os.path.join(cwd, 'LDSC/results_with_fdr.tsv')
output_path = os.path.join(cwd, 'LDSC/results_with_fdr_PD.tsv')
df.to_csv(output_path, sep='\t', index=False)


## Genetic correlation plot of BAG factor with 33 traits

# === Define your 10 categorical colors (your custom palette) ===
colors_hex_cat = [
    '#fbe6c5', '#f7c9a7', '#f3aa91', '#ee8a82', '#e2797a',
    '#d56973', '#c8586c', '#ab4762', '#8d3757', '#70284a'
]
cmap_cat = colors_hex_cat

file_path = cwd+'/LDSC/results_with_fdr_PD.tsv'
df = pd.read_csv(file_path, sep='\t')
#df = df.sort_values("GeneticCorrelation")
print(df)

#fontsize_labels = 16
#fontsize_ticks = 16
fontsize_ticks = 14
fontsize_labels = 14
cap_length = 0.05

plt.figure(figsize=(10, 8))

y_pos = range(len(df))

# colours: significant in red, nonsig in light gray
sig_color = cmap_cat[6] # pick color from 0 to 9
nonsig_color = '#aaaaaa'  # light gray

custom_labels = [
    "Longevity (90th percentile)",
    "Parkinson's disease",
    "Renal age gap",
    "Plasminogen activator inhibitor-1",
    "Coronary artery disease",
    "Pulse pressure",
    "Cardiovascular age gap",
    "Resting heart rate",
    "Intrinsic epigenetic age acceleration",
    "Autism spectrum disorder",
    "Loneliness / social isolation",
    "Sleep duration",
    "Major depressive disorder",
    "Body mass index",
    "Immune age gap",
    "Smoking initiation",
    "Glycated haemoglobin (HbA1c)",
    "Schizophrenia",
    "PhenoAge clock",
    "Systolic blood pressure",
    "Musculoskeletal age gap",
    "Pulmonary age gap",
    "ADHD",
    "C-reactive protein",
    "Type 2 diabetes",
    "Alzheimer's disease",
    "Diastolic blood pressure",
    "Bipolar disorder",
    "GrimAge clock",
    "Metabolic age gap",
    "Hannum clock",
    "Cigarettes per day",
    "Granulocyte proportion"
]


df['plot_labels'] = custom_labels

for i, (gc, se, pval, fdr_sig) in enumerate(zip(df["GeneticCorrelation"], df["SE"], df["Pvalue"], df['Significant_FDR_05'])):
    color = sig_color if fdr_sig else nonsig_color
    ci95 = 1.96 * se 
    
    # plot the circle with error bar
    plt.errorbar(gc, i, xerr=ci95, fmt='o', color=color,
                 ecolor=color, elinewidth=2, capsize=4,
                 markersize=10, markeredgewidth=1.5,
                 markeredgecolor='white', alpha=0.9)



# Apply updated labels to y-axis
plt.yticks(y_pos, df['plot_labels'], fontsize=fontsize_ticks)
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
plt.savefig('LDSC/genetic_correlation_fdr_withPD.png', dpi=600, bbox_inches='tight')

plt.show()
