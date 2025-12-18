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
file_path <- "MR-all-std-updated/MR_results_BAG_to_trait.std_2025-07-31.xlsx"
df <- read.xlsx(file_path, sheet = 1)

# Filter for IVW method and compute confidence intervals
df_ivw <- df %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    ci_lower = b - 1.96 * se,
    ci_upper = b + 1.96 * se
  )

# Create desired exposure order: brainageFactor first
custom_order <- c("brainageFactor", setdiff(unique(df_ivw$exposure), "brainageFactor"))
df_ivw$exposure <- factor(df_ivw$exposure, levels = custom_order)

# Define significance: CI does NOT include zero
df_ivw <- df_ivw %>%
  mutate(significant = (ci_lower > 0) | (ci_upper < 0))

# Set CI plotting bounds
x_min_bound <- -0.8  # lower x limit
x_max_bound <- 0.6 # optional upper bound

# Create new columns for plotting truncated CIs
df_ivw <- df_ivw %>%
  mutate(
    ci_lower_plot = pmax(ci_lower, x_min_bound),
    ci_upper_plot = pmin(ci_upper, x_max_bound),
    arrow_left = ci_lower < x_min_bound,
    arrow_right = ci_upper > x_max_bound
  )

# Add subscript-compatible labels
df_ivw <- df_ivw %>%
  mutate(
    label = gsub("BAG (.+)", "BAG['\\1']", exposure)
  )

p <- ggplot(df_ivw, aes(x = b, y = label)) +
  geom_point(aes(color = significant)) +
  
  # CI bars cropped within bounds
  geom_errorbarh(
    aes(xmin = ci_lower_plot, xmax = ci_upper_plot, color = significant),
    height = 0.2
  ) +
  
  # Arrows for lower bound
  geom_segment(
    data = df_ivw %>% filter(arrow_left),
    aes(x = x_min_bound, xend = x_min_bound - 0.05, y = label, yend = label, color = significant),
    arrow = arrow(length = unit(0.15, "cm"), ends = "last", type = "closed"),
    inherit.aes = FALSE
  ) +
  
  # Arrows for upper bound
  geom_segment(
    data = df_ivw %>% filter(arrow_right),
    aes(x = x_max_bound, xend = x_max_bound + 0.05, y = label, yend = label, color = significant),
    arrow = arrow(length = unit(0.15, "cm"), ends = "last", type = "closed"),
    inherit.aes = FALSE
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ outcome, scales = "free_y", ncol = 2, nrow = 6) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  labs(
    x = "Beta (95% CI)",
    y = "",
    #title = "Inverse Variance Weighted MR Estimates by Exposure",
    color = "Significant"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +
  xlim(x_min_bound - 0.1, x_max_bound + 0.1)

ggsave("MR_results_BAG_to_trait.std_Figure_S9.png", plot = p, width = 10, height = 12, dpi = 600)
