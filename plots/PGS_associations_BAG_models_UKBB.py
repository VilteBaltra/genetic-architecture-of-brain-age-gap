# import packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Read Excel file (skip first two rows)
df = pd.read_excel(
    "~/Downloads/Supplementary-Tables-updated_1.xlsx",
    sheet_name=26,  # 0-based index (R sheet=27)
    skiprows=2
)

print(df.head())
print(df.info())

# 2. Filter to FDR rows
df_FDR = df[df["Measure"] == "FDR"].copy()

# 3. Get totals per category
df_total = (
    df[df["Measure"] == "Total"]
    .loc[:, ["Category", "BAGFactor", "BAGHan", "BAGJawinski",
             "BAGKaufmann", "BAGLeonardsen", "BAGSmith", "BAGWen"]]
    .rename(columns=lambda x: x + "_total" if x != "Category" else x)
)

# 4. Combine FDR and Total
df_combined = df_FDR.merge(df_total, on="Category", how="left")

# 5. Compute proportions for each BAG model
bag_models = ["BAGFactor", "BAGHan", "BAGJawinski",
               "BAGKaufmann", "BAGLeonardsen", "BAGSmith", "BAGWen"]

for bag in bag_models:
    df_combined[f"{bag}_prop"] = df_combined[bag] / df_combined[f"{bag}_total"]

# 6. Identify top category per BAG model
top_categories = {
    bag: df_combined.loc[df_combined[f"{bag}_prop"].idxmax(), "Category"]
    for bag in bag_models
}

print("Top categories per BAG model:")
for bag, cat in top_categories.items():
    print(f"{bag}: {cat}")

# 7. Reshape for plotting (wide → long)
df_long = (
    df_combined
    .loc[:, ["Category"] + [f"{bag}_prop" for bag in bag_models]]
    .melt(id_vars="Category", var_name="BAG_model", value_name="Proportion")
)

df_long["BAG_model"] = df_long["BAG_model"].str.replace("_prop", "", regex=False)
df_long["Category"] = pd.Categorical(df_long["Category"], categories=df_combined["Category"], ordered=True)

# # =====================================================
# # Styling
# # =====================================================
# plt.rcParams.update({'font.size': 16})  # larger, consistent font

# fig, ax = plt.subplots(1, 1, figsize=(18, 8))  # single panel

# panel_width = 1  # bar width
colors = sns.color_palette("Set2", n_colors=7)

# Map BAG_model to colors
color_map = dict(zip(bag_models, colors))

# 8. Create mapping for legend labels with *non-italic* subscripts
pretty_labels = {
    "BAGFactor":     r"BAG$_{\mathrm{Factor}}$",
    "BAGHan":        r"BAG$_{\mathrm{Han}}$",
    "BAGJawinski":   r"BAG$_{\mathrm{Jawinski}}$",
    "BAGKaufmann":   r"BAG$_{\mathrm{Kaufmann}}$",
    "BAGLeonardsen": r"BAG$_{\mathrm{Leonardsen}}$",
    "BAGSmith":      r"BAG$_{\mathrm{Smith}}$",
    "BAGWen":        r"BAG$_{\mathrm{Wen}}$"
}


# =====================================================
# Panel: FDR Proportions with slim bars and spacing
# =====================================================
plt.rcParams.update({'font.size': 20})
fig, ax = plt.subplots(1, 1, figsize=(16, 10))

num_categories = len(df_combined["Category"])
spacing = 1.2       # increase spacing between categories
panel_width = 1   # total width for all bars in a group

# Compute x positions for category groups
x_positions = np.arange(num_categories) * spacing

# Plot bars
for i, bag in enumerate(bag_models):
    # offsets for dodged bars within each category
    offset = (i - (len(bag_models) - 1)/2) * (panel_width / len(bag_models))
    ax.bar(
        x=x_positions + offset,
        height=df_combined[f"{bag}_prop"],
        width=panel_width / len(bag_models) * 0.9,  # slightly narrower than full slot
        color=color_map[bag],
        edgecolor="black",
        alpha=0.9,
        label=pretty_labels[bag]
    )

# X-axis
ax.set_xticks(x_positions)
ax.set_xticklabels(df_combined["Category"], rotation=45, ha="right")

# Y-axis label
ax.set_ylabel("Proportion of FDR associations", fontsize=16)

# Legend above plot
ax.legend(
    title="BAG model",
    title_fontsize=15,
    fontsize=15,
    loc='upper right'
    #bbox_to_anchor=(0.5, 1.12)
    #ncol=len(bag_models)
    #frameon=False
)

# Optional grid
ax.yaxis.grid(True, linestyle='--', alpha=0.5)


# Extend y-axis slightly above tallest bar
max_height = df_combined[[f"{bag}_prop" for bag in bag_models]].max().max()
ax.set_ylim(0, max_height * 1.15)

# Remove extra space on X-axis
ax.set_xlim(x_positions[0] - panel_width, x_positions[-1] + panel_width)


plt.tight_layout()
#plt.savefig("FDR_associations_BAG_models_panel_style_spacing.png", dpi=600)
plt.show()
