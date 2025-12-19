import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

cwd = os.getcwd()
file_path = cwd + '/MRlap-results-all-updated.xlsx'

# Read both sheets
df1 = pd.read_excel(file_path, sheet_name=0)
df2 = pd.read_excel(file_path, sheet_name=1)

z = 1.96  # 95% CI

# Set up side-by-side panels
fig, axes = plt.subplots(1, 2, figsize=(13, 6))
plt.rcParams.update({'font.size': 14})

# --- Panel A ---
traits = df1["Trait"]
y_pos = np.arange(len(traits))

axes[0].errorbar(df1["observed_effect"], y_pos + 0.1, 
                 xerr=z * df1["observed_effect_se"], fmt='o', color='red', capsize=4, label='IVW effect')
axes[0].errorbar(df1["corrected_effect"], y_pos - 0.1, 
                 xerr=z * df1["corrected_effect_se"], fmt='o', color='darkred', capsize=4, label='MRlap effect')

for i in y_pos:
    axes[0].axhline(i, color='lightgray', linestyle='--', linewidth=0.5)

axes[0].set_yticks(y_pos)
axes[0].set_yticklabels(traits)
axes[0].axvline(0, color='black', linestyle='--', linewidth=1)
axes[0].invert_yaxis()
axes[0].set_xlabel("Effect size (95% CI)")
axes[0].text(-0.1, 1.05, 'A', transform=axes[0].transAxes, fontsize=16, fontweight='bold')

# --- Panel B ---
traits = df2["Trait"]
y_pos = np.arange(len(traits))

axes[1].errorbar(df2["observed_effect"], y_pos + 0.1, 
                 xerr=z * df2["observed_effect_se"], fmt='o', color='red', capsize=4)
axes[1].errorbar(df2["corrected_effect"], y_pos - 0.1, 
                 xerr=z * df2["corrected_effect_se"], fmt='o', color='darkred', capsize=4)

for i in y_pos:
    axes[1].axhline(i, color='lightgray', linestyle='--', linewidth=0.5)

axes[1].set_yticks(y_pos)
axes[1].set_yticklabels(traits)
axes[1].axvline(0, color='black', linestyle='--', linewidth=1)
axes[1].invert_yaxis()
axes[1].set_xlabel("Effect size (95% CI)")
axes[1].text(-0.1, 1.05, 'B', transform=axes[1].transAxes, fontsize=16, fontweight='bold')

# Shared legend
fig.legend(loc='upper center', ncol=2, fontsize=14)

plt.tight_layout(rect=[0, 0, 1, 0.95])  # leave space for legend
plt.savefig("MRlap_IVW_combined.png", dpi=300)
plt.show()





