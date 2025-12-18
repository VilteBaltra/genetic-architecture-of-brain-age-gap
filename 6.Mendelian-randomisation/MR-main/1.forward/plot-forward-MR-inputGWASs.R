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
file_path <- "MR-all-std-updated/MR_results_trait_to_BAG.std_2025-07-30.xlsx" # all std
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

# Step 2: Plot using the new 'label' column and parse text
p <- ggplot(df_ivw, aes(x = b, y = label)) +
  geom_point(aes(color = significant)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = significant), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ exposure, scales = "free_y", ncol = 2, nrow = 6) +
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
ggsave("MR_results_trait_to_BAG.std_Figure_S8.png", plot = p, width = 10, height = 12, dpi = 600)


