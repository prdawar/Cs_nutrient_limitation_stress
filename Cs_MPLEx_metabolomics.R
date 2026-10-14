library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

setwd("working directory")

# Read results and calculate BH adjustment before display filtering.
read_contrast <- function(file, contrast_name) {
  x <- read.csv(file, row.names = "X")
  colnames(x) <- c("FoldChange", "log2FC", "pvalue", "padjusted")
  
  valid <- is.finite(x$pvalue) & x$pvalue >= 0 & x$pvalue <= 1
  stopifnot(all(is.na(x$pvalue) | valid))
  
  x$padjusted <- NA_real_
  x$padjusted[valid] <- p.adjust(x$pvalue[valid], method = "BH")
  
  x %>%
    rownames_to_column("ID") %>%
    mutate(
      contrast = contrast_name,
      direction = if_else(log2FC > 0, "Up", "Down"))
}

meta_metab <- read.csv("./metabolomics_meta.csv") %>%
  mutate(ID = as.character(ID)) %>%
  filter(!grepl("Unknown", Name))

all_results <- bind_rows(
  read_contrast(
    "./onetenthvsfullPDA/fold_change_onetenthvsfullPDA.csv",
    "1/10 vs Full"),
  read_contrast(
    "./onetenthvshalfPDA/fold_change_onetenthvshalfPDA.csv",
    "1/10 vs Half"),
  read_contrast(
    "./halfvsfullPDA/fold_change_halfvsfullPDA.csv",
    "Half vs Full")) %>%
  inner_join(meta_metab, by = "ID") %>%
  mutate(
    contrast = factor(
      contrast,
      levels = c("1/10 vs Full", "1/10 vs Half", "Half vs Full")))

class_order <- c(
  "Organic acids", "Carbohydrates", "Nucleic acids",
  "Glycerophospholipids", "Fatty Acyls")

# Figure 5C:
# Filled circles = BH-adjusted p < 0.05.
# Open triangles = unadjusted p < 0.05, but BH-adjusted p >= 0.05.
# Both categories require |log2FC| >= 1.
scatter_df <- all_results %>%
  filter(
    is.finite(log2FC),
    !is.na(pvalue),
    !is.na(padjusted),
    pvalue < 0.05,
    abs(log2FC) >= 1,
    Super.class %in% class_order) %>%
  # Retain the original rule for repeated metabolite names.
  group_by(Name, contrast, direction) %>%
  slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Super.class = factor(Super.class, levels = class_order),
    significance = factor(
      if_else(
        padjusted < 0.05,
        "BH-adjusted p < 0.05",
        "Unadjusted p < 0.05 only"),
      levels = c(
        "BH-adjusted p < 0.05",
        "Unadjusted p < 0.05 only")))

figure5C <- ggplot(
  scatter_df,
  aes(
    x = Super.class,
    y = log2FC,
    color = direction,
    shape = significance)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50") +
  geom_jitter(width = 0.18, height = 0, size = 3, alpha = 0.7) +
  facet_wrap(~contrast, ncol = 1) +
  scale_color_manual(
    values = c("Up" = "#D55E00", "Down" = "#0072B2")) +
  scale_shape_manual(
    values = c(
      "BH-adjusted p < 0.05" = 16,
      "Unadjusted p < 0.05 only" = 17),
    drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = seq(-10, 10, by = 2),
    limits = c(-10, 10)) +
  labs(
    x = NULL,
    y = expression(log[2]~"fold change"),
    color = "Direction",
    shape = "Statistical evidence",
    title = "Metabolite Abundance Changes Under Nutrient Limitation",
    caption = paste(
      "All displayed results have |log2FC| ≥ 1.",
      "Open triangles denote nominal results retained for exploratory interpretation.")) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.box = "vertical",
    plot.caption = element_text(hjust = 0, size = 10))

print(figure5C)

ggsave(
  "Figure5C_revised.png",
  figure5C,
  width = 10,
  height = 12,
  dpi = 300,
  bg = "white")


write.csv(
  scatter_df,
  "Figure5C_plotted_results.csv",
  row.names = FALSE)
