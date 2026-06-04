library(cowplot)
library(data.table)
library(DEP)
library(devtools)
library(dplyr)
library(DreamAI)
library(edgeR)
library(EnhancedVolcano)
library(ggplot2)
library(glmGamPoi)
library(GO.db)
library(gprofiler2)
library(limma)
library(Matrix)
library(msImpute)
library(MSnbase)
library(MSnSet.utils)
library(PCAtools)
library(pheatmap)
library(proDA)
library(readxl)
library(stringr)
library(sva)
library(tidyr)
library(tidyverse)
library(uwot)

getwd()
setwd("../metabolomics/")

meta_metab <- read.csv("./metabolomics_meta.csv") %>%
  filter(!grepl("Unknown", Name))

onetenth_vs_fullPDA <- read.csv("./onetenthvsfullPDA/fold_change_onetenthvsfullPDA.csv", row.names = "X")
colnames(onetenth_vs_fullPDA) <- c("FoldChange", "log2FC", "pvalue", "padjusted")
half_vs_fullPDA <- read.csv("./halfvsfullPDA/fold_change_halfvsfullPDA.csv", row.names = "X")
colnames(half_vs_fullPDA) <- c("FoldChange", "log2FC", "pvalue", "padjusted")
onetenth_vs_halfPDA <- read.csv("./onetenthvshalfPDA/fold_change_onetenthvshalfPDA.csv", row.names = "X")
colnames(onetenth_vs_halfPDA) <- c("FoldChange", "log2FC", "pvalue", "padjusted")

onetenth_vs_fullPDA_sig <- onetenth_vs_fullPDA %>%
  rownames_to_column(var = "ID") %>%
  filter(pvalue < 0.05 & abs(log2FC) >= 1) %>%
  inner_join(meta_metab, by = "ID") %>%
  ungroup() %>%
  mutate(
    contrast = "1/10 vs Full",
    direction = if_else(log2FC > 0, "Up", "Down"))

onetenth_vs_halfPDA_sig <- onetenth_vs_halfPDA %>%
  rownames_to_column(var = "ID") %>%
  filter(pvalue < 0.05 & abs(log2FC) >= 1) %>%
  inner_join(meta_metab, by = "ID") %>%
  ungroup() %>%
  mutate(
    contrast = "1/10 vs Half",
    direction = if_else(log2FC > 0, "Up", "Down"))

half_vs_fullPDA_sig <- half_vs_fullPDA %>%
  rownames_to_column(var = "ID") %>%
  filter(pvalue < 0.05 & abs(log2FC) >= 1) %>%
  inner_join(meta_metab, by = "ID") %>%
  ungroup() %>%
  mutate(
    contrast = "Half vs Full",
    direction = if_else(log2FC > 0, "Up", "Down"))

all_sig <- bind_rows(
  onetenth_vs_fullPDA_sig,
  onetenth_vs_halfPDA_sig,
  half_vs_fullPDA_sig)

plot_df <- all_sig %>%
  filter(Super.class %in% c("Carbohydrates", "Fatty Acyls", "Glycerophospholipids", "Organic acids", "Nucleic acids")) %>%
  group_by(Name, contrast, direction) %>%
  slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    subclass = if_else(is.na(Sub.class) | Sub.class == "", "Unassigned", Sub.class), 
    Super.class = factor(
      Super.class,
      levels = c("Organic acids", "Nucleic acids", "Carbohydrates", "Glycerophospholipids", "Fatty Acyls"))) %>%
  count(contrast, Super.class, direction, name = "n") %>%
  mutate(n_plot = if_else(direction == "Down", -n, n)) %>%
  group_by(contrast, Super.class) %>%
  mutate(total = sum(abs(n_plot))) %>%
  ungroup()

ggplot(plot_df, aes(x = Super.class, y = n_plot, fill = direction)) +
  geom_col(width = 0.75) +
  coord_flip() +
  facet_grid(~contrast) +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2")) +
  labs(
    x = "",
    y = "Number of significant metabolites",
    fill = "Direction",
    title = "Up- and down-regulated metabolites by subclass") +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.major.y = element_blank())

scatter_df <- all_sig %>%
  filter(Super.class %in% c("Organic acids", "Carbohydrates", 
                            "Glycerophospholipids", "Fatty Acyls", "Nucleic acids")) %>%
  group_by(Name, contrast, direction) %>%
  slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Super.class = factor(
      Super.class,
      levels = c("Organic acids", "Carbohydrates", "Nucleic acids", 
                 "Glycerophospholipids", "Fatty Acyls")),
    direction = if_else(log2FC > 0, "Up", "Down"))

ggplot(scatter_df, aes(x = Super.class, y = log2FC, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_jitter(width = 0.18, height = 0, size = 3, alpha = 0.5) +
  facet_wrap(~ contrast, ncol = 1) +
  scale_color_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2")) +
  scale_y_continuous(breaks = seq(-8, 10, by = 2)) +
  coord_cartesian(ylim = c(-8, 10)) +
  labs(
    x = "",
    y = "log2 fold change",
    color = "Direction",
    title = "Log2FC distribution of significant metabolites across select classes"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#ggsave("metabolomics_subset_scatterplot.png", width = 8, height = 10, bg = "white")
