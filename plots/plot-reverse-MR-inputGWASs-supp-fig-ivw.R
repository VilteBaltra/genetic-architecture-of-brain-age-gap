
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
df <- read.xlsx(file_path, sheet = 2)  # reverse MR sheet

# Filter for IVW and compute confidence intervals
df_ivw <- df %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    ci_lower = b - 1.96 * se,
    ci_upper = b + 1.96 * se
  )

# order exposures (brainageFactor first)
custom_order <- c("brainageFactor", setdiff(unique(df_ivw$exposure), "brainageFactor"))
df_ivw$exposure <- factor(df_ivw$exposure, levels = custom_order)

# subscript labels for BAG(...)
df_ivw <- df_ivw %>%
  mutate(
    significant = (ci_lower > 0) | (ci_upper < 0),
    label = gsub("BAG (.+)", "BAG['\\1']", exposure)
  )

# ---> new: split outcomes into 3 columns for shared-Y-axis design <---
unique_outcomes <- unique(df_ivw$outcome)
df_ivw$col_id <- rep(1:3, length.out = length(unique_outcomes))[match(df_ivw$outcome, unique_outcomes)]

# Function to build a single column plot
make_plot <- function(dat, show_y_labels = TRUE) {
  ggplot(dat, aes(x = b, y = label)) +
    geom_point(aes(color = significant)) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = significant),
                   height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~ outcome, scales = "free_y", ncol = 1) +
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
ggsave("MR_results_BAG_to_trait.std_shared_yaxis_33traits.png",
       final_plot, width = 11, height = 16, dpi = 600)


