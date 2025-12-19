
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import matplotlib.patches as mpatches
import os
from matplotlib.colors import ListedColormap
from plotly.colors import sample_colorscale
from statsmodels.stats.multitest import multipletests


# === Define your 10 categorical colors (your custom palette) ===
colors_hex_cat = [
    '#fbe6c5', '#f7c9a7', '#f3aa91', '#ee8a82', '#e2797a',
    '#d56973', '#c8586c', '#ab4762', '#8d3757', '#70284a'
]
cmap_cat = colors_hex_cat


cap_length = 0.05

font='Arial'


# GenR PGS - BAG
# (figure 3A)

# =====================================================
# GenR data (panel A)
# =====================================================
df2 = pd.read_csv(cwd + "/trio/SM_table_renamed_20251209.csv", delimiter='\t')

# keep only BAG results for child-unadjusted
# remove genetic transmission results
# remove child and parents adjusted effects
mask = (
    df2["Term"].str.contains("unadjusted", case=False, na=False) &
    df2["model"].str.contains("BAG", case=False, na=False)
)

df2 = df2[mask]


def se_r2(R2, n):
    """Compute SE of R²."""
    return np.sqrt((4 * R2 * (1 - R2)**2) / (n - 2))

df2['se_r2'] = se_r2(df2['partial.r2'], df2['n'])

genr = df2.sort_values("partial.r2", ascending=False).reset_index(drop=True)

label_map = {
    "BAGHan": r"BAG$_{\mathrm{Han}}$",
    "BAGLeonardsen": r"BAG$_{\mathrm{Leonardsen}}$",
    "BAGCole": r"BAG$_{\mathrm{Cole}}$"
}
genr["label"] = genr["model"].map(label_map)

# =====================================================
# UKBB Data (panel B)
# =====================================================
df = pd.read_excel(cwd + "/trio/UKBB_validation_subset.xlsx")

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

# ---------------- Panel A: GenR  ----------------
panel_a_width = 0.55
ax[0].bar(genr["label"], genr["partial.r2"], width=panel_a_width,
          color="#C55D60", edgecolor="black", alpha=0.9,
          yerr=genr["se_r2"], capsize=5)
ax[0].set_ylabel(r"Partial R$^2$")

# use when differnt scale to UKBB
#ax[0].set_ylim(0, max(genr["partial.r2"] + genr["se_r2"]) * 1.1)

# use when same scale as UKBB
ax[0].set_ylim(0, 0.11)


# annotate GenR bars v2
for i, row in genr.iterrows():
    p_label = "p<0.001" if row["Pr(>|t|)"] < 0.001 else f"p={row['Pr(>|t|)']:.3f}"
    
    # use when different scale to UKBB
    #ax[0].text(i,
    #           row["partial.r2"] + row["se_r2"] + 0.001,
    #           p_label, ha="center", fontsize=10)
    #ax[0].text(i,
    #           row["partial.r2"] + row["se_r2"] + 0.0003,
    #           f"(n={row['n']})", ha="center", fontsize=9, color="grey")
    
    # use when same scale as UKBB
    ax[0].text(i,
               row["partial.r2"] + row["se_r2"] + 0.01,
               p_label, ha="center", fontsize=10)
    ax[0].text(i,
               row["partial.r2"] + row["se_r2"] + 0.006,
               f"(n={row['n']})", ha="center", fontsize=9, color="grey")
    

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
plt.savefig(cwd + "/trio/validation_UKBB_GenR_panels_notoverlapping_20251210.png", dpi=600)
plt.show()

# GenR PGS - health traits
# (unadjusted for parents)

df2 = pd.read_csv(cwd + "/trio/SM_table_renamed_20251209.csv", delimiter='\t')

fontsize_title = 10
fontsize_ticks = 10


def add_subscript(label, subscript_part):
    return f"BAG$_{{\\mathrm{{{subscript_part}}}}}$"

# Apply only to rows in the "model" column
df2["model"] = df2["model"].apply(
    lambda x: add_subscript("BAG", str(x)[3:]) if str(x).startswith("BAG") else str(x)
)

print(df2.head())

# rename:
replacements = {
    'child (parent-unadjusted)': 'effect size'
}

df2['Term'] = df2['Term'].replace(replacements)


# remove BAG results
# remove genetic transmission results
# remove child and parents adjusted effects
unwanted = ['child', 'mother','father', 'BAG']

mask = df2['Term'].str.contains('|'.join(unwanted), case=False, na=False) | \
       df2['model'].str.contains('|'.join(unwanted), case=False, na=False)

df2 = df2[~mask]

is_brain_age = df2['model'].str.contains('BAG', case=False, na=False)
pvals_group1 = df2.loc[~is_brain_age, "Pr(>|t|)"]


# Perform FDR correction 
reject1, pvals_corrected1, _, _ = multipletests(pvals_group1, alpha=0.05, method='fdr_bh')

# Create a new column with NaNs initially
df2['Pvalue.corrected'] = np.nan

# Fill adjusted p-values back
df2.loc[~is_brain_age, 'Pvalue.corrected'] = pvals_corrected1

print(df2.head())

# ---- SETTINGS ----
sig_color = cmap_cat[6]
nonsig_color = '#aaaaaa'

fontsize_ticks = 14
fontsize_labels = 12
fontsize_title = 14

# ---- SORT ROWS BY EFFECT SIZE ----
# Sort by Estimate from low → high
df2_sorted = df2.sort_values("Estimate").reset_index(drop=True)
df2_sorted["y_index"] = range(len(df2_sorted))

# Determine colors
if "Pvalue.corrected" in df2_sorted.columns:
    df2_sorted["color"] = [
        sig_color if p <= 0.05 else nonsig_color
        for p in df2_sorted["Pvalue.corrected"]
    ]
else:
    df2_sorted["color"] = "black"

# Prepare figure
fig, ax = plt.subplots(figsize=(6, 6))

# ---- Plot error bars ----
for _, row in df2_sorted.iterrows():
    est = row["Estimate"]
    se = row["Std. Error"]
    ci95 = 1.96 * se
    yi = row["y_index"]
    c = row["color"]

    ax.errorbar(
        x=est,
        y=yi,
        xerr=ci95,
        fmt='o',
        color=c,
        ecolor=c,
        elinewidth=2,
        capsize=4,
        markersize=8,
        markeredgewidth=1.5,
        markeredgecolor='white',
        alpha=0.9
    )

# ---- Add zero line ----
ax.axvline(0, linestyle='--', color='gray', linewidth=1)

# ---- Left y-axis labels: MODEL names ----
ax.set_yticks(df2_sorted["y_index"])
ax.set_yticklabels(df2_sorted["model"], fontsize=fontsize_ticks)

# ---- Remove right y-axis entirely ----
ax2 = ax.twinx()
ax2.set_yticks([])
ax2.set_ylabel("")
ax2.set_yticklabels([])

# ---- Style / labels ----
ax.set_xlabel("Estimate", fontsize=fontsize_title)
ax.tick_params(axis='x', labelsize=14)
#ax.set_ylabel("Model", fontsize=fontsize_title)

ax.set_facecolor('#f9f9f9')
ax.grid(axis='x', linestyle='--', color='gray', alpha=0.15, linewidth=0.7)
ax.grid(axis='y', visible=False)

# Frame styling
spine_color = '#cccccc'
for spine in ax.spines.values():
    spine.set_color(spine_color)
    spine.set_linewidth(1)

# ---- Legend ----
sig_patch = mpatches.Patch(color=sig_color, label='FDR-adjusted p < 0.05')
nonsig_patch = mpatches.Patch(color='#aaaaaa', label='FDR-adjusted p ≥ 0.05')

ax.legend(handles=[sig_patch, nonsig_patch],
          fontsize=8,
          frameon=False,
          loc='lower right')

plt.tight_layout()
plt.savefig(cwd + "/trio/PGS.child_health_20251209.png", dpi=300, bbox_inches='tight')
plt.show()

# trio plot

df2 = pd.read_csv(cwd + "/trio/SM_table_renamed_20251209.csv", delimiter='\t')

fontsize_title = 10
fontsize_ticks = 10


def add_subscript(label, subscript_part):
    return f"BAG$_{{\\mathrm{{{subscript_part}}}}}$"

# Apply only to rows in the "model" column
df2["model"] = df2["model"].apply(
    lambda x: add_subscript("BAG", str(x)[3:]) if str(x).startswith("BAG") else str(x)
)

print(df2.head())

# remove BAG results
# remove genetic transmission results
# remove child unadjusted effects
unwanted = ['parent-unadjusted']

mask = df2['Term'].str.contains('|'.join(unwanted), case=False, na=False) 

df2 = df2[~mask]


from statsmodels.stats.multitest import multipletests
df2['Pvalue.corrected'] = multipletests(df2['Pr(>|t|)'], method='fdr_bh')[1]


df2 = df2.sort_values(by='model', key=lambda x: x.str.contains('BAG', case=False, na=False), ascending=False)


g = sns.FacetGrid(
    df2,
    col="model",
    col_wrap=3,
    sharex=False,
    sharey=True,
    height=3,
    aspect=1.2
)

g.fig.set_size_inches(10, 10)  # width=14 inches, height=6 inches

def plot_estimates(data, color, **kwargs):
    ax = plt.gca()
    data = data.sort_values("Term")
    y_pos = np.arange(len(data))

    sig_color = cmap_cat[6]
    nonsig_color = '#aaaaaa'

    if "Pvalue.corrected" in data.columns:
        colors = [sig_color if val <= 0.05 else nonsig_color for val in data["Pvalue.corrected"]]
    else:
        colors = ['black'] * len(data)

    for i, (est, se, c) in enumerate(zip(data["Estimate"], data["Std. Error"], colors)):
        ci95 = 1.96 * se
        ax.errorbar(
            x=est,
            y=i,
            xerr=ci95,
            fmt='o',
            color=c,
            ecolor=c,
            elinewidth=2,
            capsize=4,
            markersize=8,
            markeredgewidth=1.5,
            markeredgecolor='white',
            alpha=0.9
        )

    ax.axvline(0, linestyle='--', color='gray', linewidth=1)
    ax.set_yticks(y_pos)

    padding = 0.5
    ax.set_ylim(-padding, len(data) - 1 + padding)

    ax.set_yticklabels(data["Term"].to_numpy(dtype=str), fontsize=fontsize_ticks)  # use fontsize_ticks here
    ax.tick_params(axis='y', pad=2, labelsize=fontsize_ticks)  # ensure y tick labels font size
    ax.tick_params(axis='x', labelsize=fontsize_ticks)         # use fontsize_ticks for x tick labels

    # Styling for frame and ticks
    ax.set_facecolor('#f9f9f9')
    ax.grid(axis='x', linestyle='--', color='gray', alpha=0.15, linewidth=0.7)
    ax.grid(axis='y', visible=False)

    spine_color = '#cccccc'  # light gray frame color
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(spine_color)
        spine.set_linewidth(1)

    tick_color = '#555555'  # medium gray tick color
    ax.tick_params(axis='x', colors=tick_color)
    ax.tick_params(axis='y', colors=tick_color)

g.map_dataframe(plot_estimates)

# Set axis labels with fontsize_labels (need to set manually on each axis)
for ax in g.axes.flat:
    ax.set_xlabel("Estimate", fontsize=fontsize_title)
    ax.set_ylabel("", fontsize=fontsize_title)  # empty y-label but fontsize applied

#g.set_titles(col_template="{col_name}", size=fontsize_title)  # use fontsize_title here
g.set_titles(col_template="{col_name}", size=fontsize_labels)  # use fontsize_title here

sig_patch = mpatches.Patch(color=cmap_cat[6], label='FDR-adjusted p<0.05')
nonsig_patch = mpatches.Patch(color='#aaaaaa', label='FDR-adjusted p ≥ 0.05')

g.fig.legend(
    handles=[sig_patch, nonsig_patch],
    fontsize=fontsize_ticks,  # legend font size matches ticks
    frameon=False,
    loc='upper right',
    bbox_to_anchor=(1, -0.001),
    borderaxespad=0.8
)

plt.tight_layout()
plt.savefig(cwd + "/trio/trio_20251209.png", dpi=300, bbox_inches='tight')
plt.show()
