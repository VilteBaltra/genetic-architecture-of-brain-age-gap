import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =====================================================
# GenR data (panel A)
# =====================================================
cwd = "/Volumes/Files/Psychology/ResearchProjects/EWalton/BrainHealth/BrainAge_GWAS/output/all_plots/forPlotting_Shir_final/"
df2 = pd.read_csv(cwd + "/trio/SM.table_renamed2.csv", delimiter=',') # original version
#df2 = pd.read_csv(cwd + "/trio/SM_table_renamed_20251209.csv", delimiter='\t') # new version (not integrated yet)

genr = df2[df2["Term"] == "child genetic effect (unadjusted)"].copy()
genr["CI_lower"] = genr["Estimate"] - 1.96 * genr["Std. Error"]
genr["CI_upper"] = genr["Estimate"] + 1.96 * genr["Std. Error"]
genr = genr[genr["model"].isin(["BAGHan", "BAGLeonardsen", "BAGCole"])]

label_map = {
    "BAGHan": r"BAG$_{\mathrm{Han}}$",
    "BAGLeonardsen": r"BAG$_{\mathrm{Leonardsen}}$",
    "BAGCole": r"BAG$_{\mathrm{Cole}}$"
}
genr["label"] = genr["model"].map(label_map)

# =====================================================
# UKBB Data (panel B)
# =====================================================
#df = pd.read_excel("/Users/vb506/Documents/UKBB_validation_subset.xlsx")
df = pd.read_excel("/Users/vb506/Documents/Documents-Nav-PC/UKBB_validation_subset.xlsx")

def se_r2(R2, n):
    """Compute SE of R²."""
    return np.sqrt((4 * R2 * (1 - R2)**2) / (n - 2))

df['SE_Han'] = se_r2(df['R2_Han'], df['n_Han'])
df['SE_Jaw'] = se_r2(df['R2_Jawinski'], df['n_Jawinski'])
df = df.sort_values("Ancestry", ascending=False).reset_index(drop=True)

# =====================================================
# Plotting
# =====================================================
plt.rcParams.update({'font.size': 14})

fig, ax = plt.subplots(
    1, 2, figsize=(16, 6),
    gridspec_kw={'width_ratios': [1, 1.7]}  # panel A = 1x, panel B = 1.7x seemed to work best
)

# ---------------- Panel A: GenR (reminder: plot here R2 eventually!) ----------------
panel_a_width = 0.55
ax[0].bar(genr["label"], genr["Estimate"], width=panel_a_width,
          color="#C55D60", edgecolor="black", alpha=0.9,
          yerr=genr["Std. Error"], capsize=5)
ax[0].set_ylabel("Estimate")
ax[0].set_ylim(0, max(genr["Estimate"] + genr["Std. Error"]) * 1.25)

## annotate GenR bars
#for i, row in genr.iterrows():
#    p_label = "p<0.001" if row["Pr(>|t|)"] < 0.001 else f"p={row['Pr(>|t|)']:.3f}"
#    ax[0].text(i, row["Estimate (b)"] + row["Std. Error"] + 0.005,
#               p_label, ha="center", fontsize=11)
# annotate GenR bars v2
for i, row in genr.iterrows():
    p_label = "p<0.001" if row["Pr(>|t|)"] < 0.001 else f"p={row['Pr(>|t|)']:.3f}"
    ax[0].text(i, row["Estimate"] + row["Std. Error"] + 0.005,
               p_label, ha="center", fontsize=11)


# ---------------- Panel B: UKBB ----------------
labels = df["Ancestry"]
n = len(labels)
bar_width = 1.3
group_positions = np.arange(n) * 3  # spacing

# Bars
ax[1].bar(group_positions - bar_width/2, df["R2_Jawinski"], width=bar_width,
          color="#EBBAB9", edgecolor="black", yerr=df["SE_Jaw"], capsize=5,
          label=r"BAG$_{\mathrm{Jawinski}}$")
ax[1].bar(group_positions + bar_width/2, df["R2_Han"], width=bar_width,
          color="#CC2936", edgecolor="black", alpha=0.7,
          yerr=df["SE_Han"], capsize=5,
          label=r"BAG$_{\mathrm{Han}}$")

# y-axis limit
ymax = max(
    (df["R2_Han"] + df["SE_Han"]).max(),
    (df["R2_Jawinski"] + df["SE_Jaw"]).max()
) + 0.01
ax[1].set_ylim(0, ymax * 1.15)

# annotating bars
for i, row in df.iterrows():
    # Jawinski
    p_label_jaw = "p<0.001" if row["p_Jawinski"] < 0.001 else f"p={row['p_Jawinski']:.3f}"
    ax[1].text(group_positions[i] - bar_width/2,
               row["R2_Jawinski"] + row["SE_Jaw"] + 0.007,
               p_label_jaw, ha="center", fontsize=10)
    ax[1].text(group_positions[i] - bar_width/2,
               row["R2_Jawinski"] + row["SE_Jaw"] + 0.003,
               f"(n={row['n_Jawinski']})", ha="center", fontsize=9, color="grey")

    # Han
    p_label_han = "p<0.001" if row["p_Han"] < 0.001 else f"p={row['p_Han']:.3f}"
    ax[1].text(group_positions[i] + bar_width/2,
               row["R2_Han"] + row["SE_Han"] + 0.007,
               p_label_han, ha="center", fontsize=10)
    ax[1].text(group_positions[i] + bar_width/2,
               row["R2_Han"] + row["SE_Han"] + 0.003,
               f"(n={row['n_Han']})", ha="center", fontsize=9, color="grey")

# X-axis formatting
ax[1].set_xticks(group_positions)
ax[1].set_xticklabels(labels)
ax[1].set_ylabel("Partial R²")
ax[1].legend()

plt.tight_layout()
#plt.savefig("validation_UKBB_GenR_panels_notoverlapping.png", dpi=600)
plt.show()
