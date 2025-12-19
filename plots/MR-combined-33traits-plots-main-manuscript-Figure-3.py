
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# define trait groups
trait_groups = {
    'Mental health': [
        'Schizophrenia', 'Major depressive disorder', "Alzheimer's disease",
        'Bipolar disorder', 'ADHD', 'Autism spectrum disorder*',
        "Parkinson's disease", 'Loneliness / social isolation'
    ],
    'Physical health': [
        'Type 2 diabetes', 'Body mass index', 'C-reactive protein',
        'Systolic blood pressure*', 'Diastolic blood pressure*', 'Pulse pressure',
        'Resting heart rate', 'Coronary artery disease', 'Glycated hemoglobin (HbA1c)',
        'Plasminogen activator inhibitor-1', 'Granulocyte proportion'
    ],
    'Aging': [
        'Metabolic age gap', 'GrimAge clock', 'PhenoAge clock', 
        'Intrinsic epigenetic age acceleration', 'Hannum clock',
        'Cardiovascular age gap', 'Pulmonary age gap', 'Renal age gap',
        'Immune age gap', 'Musculoskeletal age gap', 'Longevity (90th percentile)*'
    ],
    'Lifestyle': [
        'Smoking initiation*', 'Cigarettes per day', 'Sleep duration*'
    ]
}

# === Colors, legend, offsets, styling ===
colors_hex_cat = [
    '#fbe6c5', '#f7c9a7', '#f3aa91', '#ee8a82', '#e2797a',
    '#d56973', '#c8586c', '#ab4762', '#8d3757', '#70284a'
]
method_colors = {
    'Inverse variance weighted': colors_hex_cat[9],
    'Weighted median': colors_hex_cat[6],
    'MR Egger': colors_hex_cat[2],
}
legend_handles = [mpatches.Patch(color=color, label=method) for method, color in method_colors.items()]
fontsize_ticks = 20
fontsize_labels = 20
cap_length = 0.05
method_offset = {'Inverse variance weighted': -0.25, 'Weighted median': 0, 'MR Egger': 0.25}

# === Load data ===
df_a = pd.read_excel("main-MR/forward-and-reverse-MR-33traits-results.xlsx", sheet_name=0) # first sheet
df_b = pd.read_excel("main-MR/forward-and-reverse-MR-33traits-results.xlsx", sheet_name=1) # second sheet
# df_a = df_a[df_a['outcome'] == "BAG Han"]
# df_b = df_b[df_b['exposure'] == "BAG Han"]
# df_a = df_a[df_a['outcome'] == "BAG Jawinski"]
# df_b = df_b[df_b['exposure'] == "BAG Jawinski"]
# df_a = df_a[df_a['outcome'] == "BAG Wen"]
# df_b = df_b[df_b['exposure'] == "BAG Wen"]
# df_a = df_a[df_a['outcome'] == "BAG Smith"]
# df_b = df_b[df_b['exposure'] == "BAG Smith"]
# df_a = df_a[df_a['outcome'] == "BAG Kaufmann"]
# df_b = df_b[df_b['exposure'] == "BAG Kaufmann"]
# df_a = df_a[df_a['outcome'] == "BAG Leonardsen"]
# df_b = df_b[df_b['exposure'] == "BAG Leonardsen"]
df_a = df_a[df_a['outcome'] == "BAG factor"]
df_b = df_b[df_b['exposure'] == "BAG factor"]

###
# Step 1: invert mapping
trait_to_category = {
    trait: category
    for category, traits in trait_groups.items()
    for trait in traits
}

# Step 2: map
df_a["group"] = df_a["exposure"].map(trait_to_category)
df_b["group"] = df_b["outcome"].map(trait_to_category)

# Step 3 (optional): show unmatched traits
print("Unmatched traits:", df_a[df_a["group"].isna()]["exposure"].unique())
print("Unmatched traits:", df_b[df_b["group"].isna()]["outcome"].unique())

###

# store four trait groups names: 'Aging', 'Lifestyle', 'Mental health', 'Physical health'
trait_groups = sorted(set(df_a['group'].dropna().unique()).union(
                 set(df_b['group'].dropna().unique())))


# === Loop per trait group ===

# Define custom x-axis limits for Panel A (forward MR)
x_limits_a_dict = {
    'Physical health': (-0.25, 0.25),
    'Aging': (-0.3, 0.3),
    'Mental health': (-0.7, 0.7),
    'Lifestyle': (-0.4, 0.4)  # optional
}

# Define custom x-axis limits for Panel B (reverse MR)
x_limits_b_dict = {
    'Physical health': (-2, 2),
    'Aging': (-1.9, 1.9),
    'Mental health': (-3.6, 3.6)
    #'Lifestyle': (-1.6, 1.6)
}

for group_name in trait_groups:

    df_a_group = df_a[df_a['group'] == group_name].copy()
    df_b_group = df_b[df_b['group'] == group_name].copy()

# replaced below with above
# for group_name, traits in trait_groups.items():
    
#     df_a_group = df_a[df_a['exposure'].isin(traits)].copy()
#     df_b_group = df_b[df_b['outcome'].isin(traits)].copy()
    
    # Recalculate y positions for this group
    exposures_a_group = df_a_group['exposure'].unique()[::-1]
    exposure_positions_a_group = {exposure: i for i, exposure in enumerate(exposures_a_group)}
    df_a_group['y'] = df_a_group.apply(lambda row: exposure_positions_a_group[row['exposure']] + method_offset[row['method']], axis=1)

    exposures_b_group = df_b_group['outcome'].unique()[::-1]
    exposure_positions_b_group = {exposure: i for i, exposure in enumerate(exposures_b_group)}
    df_b_group['y'] = df_b_group.apply(lambda row: exposure_positions_b_group[row['outcome']] + method_offset[row['method']], axis=1)

    # --- Create figure ---
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(18, 10), sharex=False)
    ax1.grid(False)
    ax2.grid(False)

    # ---------------------------------------------------------------------
    # PANEL A — FORWARD MR
    # ---------------------------------------------------------------------
    # group-specific trimming
    if group_name in x_limits_a_dict:
        x_min_a, x_max_a = x_limits_a_dict[group_name]
    else:
        max_abs_a = max(abs(df_a_group['lCI'].min()), abs(df_a_group['hCI'].max()))
        x_min_a = -max_abs_a * 1.05
        x_max_a = max_abs_a * 1.05
    
    clip_buffer_a = 0.02 * (x_max_a - x_min_a)

    x_min_buf_a = x_min_a + clip_buffer_a
    x_max_buf_a = x_max_a - clip_buffer_a
    arrow_gap_a = 0.01


    for _, row in df_a_group.iterrows():
        y = row['y']
        color = method_colors.get(row['method'], '#aaaaaa')

        ci_l, ci_u = row['lCI'], row['hCI']
        clipped_l = max(ci_l, x_min_buf_a)
        clipped_u = min(ci_u, x_max_buf_a)
        is_clip_l = ci_l < x_min_buf_a
        is_clip_r = ci_u > x_max_buf_a

        ci_start = clipped_l
        ci_end = clipped_u


        ax1.plot(row['b'], y, 'o', color=color,
                 markersize=12, markeredgewidth=2, markeredgecolor='white')
        ax1.plot([ci_start, ci_end], [y, y], color=color, linewidth=3)


        arrow_length = 0.03 * (x_max_a - x_min_a)  # 3% of x-axis width

        if is_clip_l:
            ax1.annotate('',
                         xy=(clipped_l, y),
                         xytext=(clipped_l - arrow_length, y),
                         arrowprops=dict(arrowstyle='<|-', color=color, lw=2, mutation_scale=15))

        if is_clip_r:
            ax1.annotate('',
                         xy=(clipped_u, y),
                         xytext=(clipped_u + arrow_length, y),
                         arrowprops=dict(arrowstyle='<|-', color=color, lw=2, mutation_scale=15))
        

    ax1.axvline(0, color='gray', linestyle='--')
    ax1.set_yticks(range(len(exposures_a_group)))
    ax1.set_yticklabels(exposures_a_group, fontsize=24)
    #ax1.set_xlim(x_min_a, x_max_a)
    x_padding_a = 0.1 * (x_max_a - x_min_a)
    ax1.set_xlim(x_min_a - x_padding_a, x_max_a + x_padding_a)
    
    ax1.set_xlabel("Beta (95% CI)", fontsize=20)
    ax1.tick_params(axis='x', labelsize=20)
    #ax1.set_title("A", loc='left', fontsize=18, weight='bold')
    ax1.invert_yaxis()

    # ---------------------------------------------------------------------
    # PANEL B — REVERSE MR
    # ---------------------------------------------------------------------

    # group-specific trimming
    if group_name in x_limits_b_dict:
        x_min_b, x_max_b = x_limits_b_dict[group_name]
    else:
        max_abs_b = max(abs(df_b_group['lCI'].min()), abs(df_b_group['hCI'].max()))
        x_min_b = -max_abs_b * 1.05
        x_max_b = max_abs_b * 1.05

    clip_buffer_b = 0.02 * (x_max_b - x_min_b)
    x_min_buf_b = x_min_b + clip_buffer_b
    x_max_buf_b = x_max_b - clip_buffer_b
    arrow_gap_b = 0.03 * (x_max_b - x_min_b)

    for _, row in df_b_group.iterrows():
        y = row['y']
        color = method_colors.get(row['method'], '#aaaaaa')

        ci_l, ci_u = row['lCI'], row['hCI']
        clipped_l = max(ci_l, x_min_buf_b)
        clipped_u = min(ci_u, x_max_buf_b)
        is_clip_l = ci_l < x_min_buf_b
        is_clip_r = ci_u > x_max_buf_b

        #ci_start = clipped_l + (arrow_gap_b if is_clip_l else 0)
        #ci_end   = clipped_u - (arrow_gap_b if is_clip_r else 0)
        ci_start = clipped_l
        ci_end = clipped_u


        ax2.plot(row['b'], y, 'o', color=color,
                 markersize=12, markeredgewidth=2, markeredgecolor='white')
        ax2.plot([ci_start, ci_end], [y, y], color=color, linewidth=3)
                         
        arrow_length = 0.03 * (x_max_b - x_min_b)  # 3% of x-axis width

        if is_clip_l:
            ax2.annotate('',
                         xy=(clipped_l, y),
                         xytext=(clipped_l - arrow_length, y),
                         arrowprops=dict(arrowstyle='<|-', color=color, lw=2, mutation_scale=15))

        if is_clip_r:
            ax2.annotate('',
                         xy=(clipped_u, y),
                         xytext=(clipped_u + arrow_length, y),
                         arrowprops=dict(arrowstyle='<|-', color=color, lw=2, mutation_scale=15))

        



    ax2.axvline(0, color='gray', linestyle='--')
    ax2.set_yticks(range(len(exposures_b_group)))
    ax2.set_yticklabels([], fontsize=20)
    
    #ax2.set_xlim(x_min_b, x_max_b)
    x_padding_b = 0.1 * (x_max_b - x_min_b)
    ax2.set_xlim(x_min_b - x_padding_b, x_max_b + x_padding_b)
    
    ax2.set_xlabel("Beta (95% CI)", fontsize=20)
    ax2.tick_params(axis='x', labelsize=22)
    #ax2.set_title("B", loc='left', fontsize=18, weight='bold')
    ax2.invert_yaxis()

    # Shared legend
    fig.legend(handles=legend_handles, title="MR Method", fontsize=20,
               title_fontsize=18,
               loc='lower center', bbox_to_anchor=(0.5, -0.05),
               ncol=3, frameon=False)

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    plt.savefig(f"mr_panels_{group_name.replace('/', '_')}.png", dpi=600, bbox_inches='tight')
    plt.show()
