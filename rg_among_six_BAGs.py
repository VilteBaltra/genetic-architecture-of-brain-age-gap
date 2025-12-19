import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches

# === Settings ===
colors_hex_cat = [
    '#fbe6c5', '#f7c9a7', '#f3aa91', '#ee8a82', '#e2797a',
    '#d56973', '#c8586c', '#ab4762', '#8d3757', '#70284a'
]
cmap_cat = colors_hex_cat
fontsize_ticks = 14
fontsize_labels = 14
font = "Arial"

# === Load data ===
file_path = "gCorr_BA/gCorr-gsem-updated_v2.csv"
df2 = pd.read_csv(file_path, sep='\t', index_col=0)

# Ensure index/columns are strings
df2.index = df2.index.astype(str)
df2.columns = df2.columns.astype(str)

# === LaTeX subscript labels ===
def make_label(col):
    if col.startswith("BAG"):
        return f"BAG$_{{\\mathrm{{{col[3:]}}}}}$"
    return col

df2.columns = [make_label(c) for c in df2.columns]
df2.index = [make_label(c) for c in df2.index]

# === Manual order ===
manual_order = [
    "BAG$_{\\mathrm{Jawinski}}$",
    "BAG$_{\\mathrm{Leonardsen}}$",
    "BAG$_{\\mathrm{Smith}}$",
    "BAG$_{\\mathrm{Wen}}$",
    "BAG$_{\\mathrm{Kaufmann}}$",
    "BAG$_{\\mathrm{Han}}$"
]

row_order = [df2.index.get_loc(r) for r in manual_order]
col_order = [df2.columns.get_loc(c) for c in manual_order]

df_ordered = df2.iloc[row_order, col_order]

# === Mask upper triangle including diagonal ===
mask = np.triu(np.ones_like(df_ordered, dtype=bool), k=0)

# === Create annotation array (NaN for masked cells) ===
annot_array = df_ordered.round(2).copy()
annot_array = annot_array.where(~mask)  # NaN for upper triangle & diagonal

# === Plot heatmap ===
plt.figure(figsize=(6, 6))
ax = sns.heatmap(
    df_ordered,
    annot=annot_array,
    cmap=cmap_cat,
    mask=mask,
    square=True,
    vmin=0,
    vmax=1,
    annot_kws={"size": fontsize_ticks, "fontname": font},
    xticklabels=df_ordered.columns,
    yticklabels=df_ordered.index,
    linewidth=0,
    cbar_kws={'shrink': 0.64, 'pad': -0.07}  # smaller and closer colorbar
)

# === Hide tick labels only for diagonal/upper triangle ===
n = df_ordered.shape[0]
# X-axis (columns)
xticks = []
for i in range(n):
    lbl = df_ordered.columns[i]
    if any(~mask[:, i]):  # check if any row is unmasked
        xticks.append(lbl)
    else:
        xticks.append('')
# Y-axis (rows)
yticks = []
for i in range(n):
    lbl = df_ordered.index[i]
    if any(~mask[i, :]):  # check if any column is unmasked
        yticks.append(lbl)
    else:
        yticks.append('')

ax.set_xticklabels(xticks, rotation=45, ha='right', fontsize=fontsize_ticks, fontname=font)
ax.set_yticklabels(yticks, rotation=0, fontsize=fontsize_ticks, fontname=font)

# === Faint borders for lower triangle ===
for i in range(n):
    for j in range(n):
        if i > j:
            rect = patches.Rectangle((j, i), 1, 1, fill=False, edgecolor='lightgray', linewidth=0.5)
            ax.add_patch(rect)

ax.grid(False)

# Remove tick marks
ax.tick_params(axis='x', length=0)
ax.tick_params(axis='y', length=0)

# Axis labels
ax.set_xlabel(" ", fontsize=fontsize_labels, fontname=font)
ax.set_ylabel(" ", fontsize=fontsize_labels, fontname=font)

# Colorbar label
cbar = ax.collections[0].colorbar
cbar.set_label("rg", fontsize=fontsize_labels, fontname=font)

plt.savefig("gCorr_BA/genetic_correlation_plot_clustered_BAGHan_on_X.png", dpi=600, bbox_inches='tight')
plt.show()
