# Install packages if not already installed
packages <- c("openxlsx", "ggplot2", "dplyr")
installed <- rownames(installed.packages())
if (!all(packages %in% installed)) {
  install.packages(setdiff(packages, installed))
}

# Load libraries
library(openxlsx)
library(ggplot2)
library(dplyr)

# Read the Excel file
#file_path <- "/Users/naveennainwal/Documents/Vilte-2025-NaveenPC/2025/3.MR/scripts/input-GWASs/results/forward/MR-output-forward/forward-MR-results/combined/MR_results_trait_to_BAG_2025-07-16.xlsx"
#file_path <- "/Users/naveennainwal/Documents/MR_results_trait_to_BAG.std_2025-07-30.xlsx" # outcome std
#file_path <- "/Users/naveennainwal/Documents/MR-all-std-updated/MR_results_trait_to_BAG.std_2025-07-30.xlsx" # all std
file_path <- "/Users/vb506/Documents/MR-analysis-local/main-MR/forward-and-reverse-MR-33traits-results.xlsx" # all std 33 traits
df <- read.xlsx(file_path, sheet = 1)

# Filter for IVW method and compute confidence intervals
df_ivw <- df %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    ci_lower = b - 1.96 * se,
    ci_upper = b + 1.96 * se
  )

# Create desired outcome order: brainageFactor first
custom_order <- c("brainageFactor", setdiff(unique(df_ivw$outcome), "brainageFactor"))
df_ivw$outcome <- factor(df_ivw$outcome, levels = custom_order)

# Step 1: Add subscript labels
df_ivw <- df_ivw %>%
  mutate(
    significant = (ci_lower > 0) | (ci_upper < 0),
    label = gsub("BAG (.+)", "BAG['\\1']", outcome)  # subscript-compatible
  )

# # Step 2: Plot using the new 'label' column and parse text
# p <- ggplot(df_ivw, aes(x = b, y = label)) +
#   geom_point(aes(color = significant)) +
#   geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = significant), height = 0.2) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
#   facet_wrap(~ exposure, scales = "free_y", ncol = 2, nrow = 6) +
#   scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
#   scale_y_discrete(labels = function(x) parse(text = x)) +  # <- This is key!
#   labs(
#     x = "Beta (95% CI)",
#     y = "",
#     color = "Significant"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(legend.position = "none")

# Step 2: Plot using the new 'label' column and parse text
p <- ggplot(df_ivw, aes(x = b, y = label)) +
  geom_point(aes(color = significant)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = significant), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ exposure, scales = "free_y", ncol = 3, nrow = 11) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
  scale_y_discrete(labels = function(x) parse(text = x)) +  # <- This is key!
  labs(
    x = "Beta (95% CI)",
    y = "",
    color = "Significant"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

# Save the plot
setwd("/Users/vb506/Documents/MR-analysis-local/main-MR/")
#ggsave("MR_results_trait_to_BAG.std_Figure_S8.png", plot = p, width = 10, height = 12, dpi = 600)
ggsave("MR_results_trait_to_BAG.std_Figure_S8_33traits.png", plot = p, width = 12, height = 12, dpi = 600)



### new


# install packages 
packages <- c("openxlsx", "ggplot2", "dplyr", "patchwork")
installed <- rownames(installed.packages())
if (!all(packages %in% installed)) {
  install.packages(setdiff(packages, installed))
}

# load libraries
library(openxlsx)
library(ggplot2)
library(dplyr)
library(patchwork)

# read the Excel file
file_path <- "/Users/vb506/Documents/MR-analysis-local/main-MR/forward-and-reverse-MR-33traits-results.xlsx"
df <- read.xlsx(file_path, sheet = 1)  # reverse MR sheet

# Filter for IVW and compute confidence intervals
df_ivw <- df %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    ci_lower = b - 1.96 * se,
    ci_upper = b + 1.96 * se
  )

# order exposures (brainageFactor first)
custom_order <- c("brainageFactor", setdiff(unique(df_ivw$outcome), "brainageFactor"))
df_ivw$outcome <- factor(df_ivw$outcome, levels = custom_order)

# subscript labels for BAG(...)
df_ivw <- df_ivw %>%
  mutate(
    significant = (ci_lower > 0) | (ci_upper < 0),
    label = gsub("BAG (.+)", "BAG['\\1']", outcome)
  )

# ---> new: split outcomes into 3 columns for shared-Y-axis design <---
unique_outcomes <- unique(df_ivw$exposure)
df_ivw$col_id <- rep(1:3, length.out = length(unique_outcomes))[match(df_ivw$exposure, unique_outcomes)]

# Function to build a single column plot
make_plot <- function(dat, show_y_labels = TRUE) {
  ggplot(dat, aes(x = b, y = label)) +
    geom_point(aes(color = significant)) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = significant),
                   height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~ exposure, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
    scale_y_discrete(labels = if (show_y_labels) function(x) parse(text = x) else NULL) +
    labs(
      x = NULL,                       # REMOVE individual x-axis title
      y = if (show_y_labels) "" else NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      axis.title.y = if (!show_y_labels) element_blank() else NULL,
      axis.text.y  = if (!show_y_labels) element_blank() else NULL
    ) +
    coord_cartesian(xlim = c(-0.7, 0.7))
}


# create 3 aligned columns
p1 <- make_plot(df_ivw %>% filter(col_id == 1), show_y_labels = TRUE)
p2 <- make_plot(df_ivw %>% filter(col_id == 2), show_y_labels = FALSE)
p3 <- make_plot(df_ivw %>% filter(col_id == 3), show_y_labels = FALSE)

# combine plots with shared y-axis (only left column has labels)
final_plot <- p1 | p2 | p3

#save the final figure
setwd("/Users/vb506/Documents/MR-analysis-local/main-MR/")
ggsave("MR_results_trait_to_BAG.std_shared_yaxis_33traits.png",
       final_plot, width = 11, height = 16, dpi = 600)





