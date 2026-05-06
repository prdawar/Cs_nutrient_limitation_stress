library(tidyverse)
library(tidyr)
library(dplyr)
library("devtools")
#require("remotes")
#install_github("WangLab-MSSM/DreamAI/Code")
library(DreamAI)
library(PCAtools)
library(sva)
library(proDA)
library(ggplot2)
library(uwot)
library(stringr)
library(msImpute)
#devtools::install_version("Matrix", version = "1.6.1.1")
library(Matrix)
library(cowplot)
library(data.table)
library(edgeR)
library(glmGamPoi)
library(gprofiler2)
library(pheatmap)
library(MSnSet.utils)
library(EnhancedVolcano)
library(DEP)
library(MSnbase)
library(GO.db)
library(limma)
library(readxl)

getwd()
setwd("../../Users/dawa726/OneDrive - PNNL/Desktop/BRAVE_PD/Cs_MPLEx_peptides/")

meta_metab <- read.csv("./metabolomics_meta.csv") %>%
  filter(!grepl("Unknown", Name))

onetenth_vs_fullPDA <- read.csv("./onetenthvsfullPDA/fold_change_onetenthvsfullPDA.csv", row.names = "X")
half_vs_fullPDA <- read.csv("./halfvsfullPDA/fold_change_halfvsfullPDA.csv", row.names = "X")
onetenth_vs_halfPDA <- read.csv("./onetenthvshalfPDA/fold_change_onetenthvshalfPDA.csv", row.names = "X")

#EnhancedVolcano(onetenth_vs_fullPDA,
#                lab = NA,
#                x = 'log2.FC.',
#                y = 'pvalue',
#                #selectLab = c(""),
#                pCutoff = 0.0093728, 
#                FCcutoff = 1, 
#                labSize = 0,
#                legendLabels=c('Not sig.',
#                               'Log (base 2) FC',
#                               'p-adj < 0.05',
#                               'p-adj < 0.05 & Log (base 2) FC'),
#                legendPosition = 'right',
#                legendLabSize = 12,
#                legendIconSize = 4,
#                cutoffLineType = 'twodash',
#                cutoffLineWidth = 1,
#                pointSize = 4, 
#                title = "1/10 PDA vs PDA; metabolomics", 
#                subtitle = NULL) +
#  theme_minimal(base_size = 20)

#EnhancedVolcano(half_vs_fullPDA,
#                lab = NA,
#                x = 'log2.FC.',
#               y = 'pvalue',
#                #selectLab = c(""),
#                pCutoff = 0.020863, 
#                FCcutoff = 1, 
#                labSize = 0,
#                legendLabels=c('Not sig.',
#                               'Log (base 2) FC',
#                               'p-adj < 0.05',
#                               'p-adj < 0.05 & Log (base 2) FC'),
#                legendPosition = 'right',
#                legendLabSize = 12,
#                legendIconSize = 4,
#                cutoffLineType = 'twodash',
#                cutoffLineWidth = 1,
#                pointSize = 4, 
#                title = "1/2 PDA vs PDA; metabolomics", 
#                subtitle = NULL) +
#  theme_minimal(base_size = 20)

#EnhancedVolcano(onetenth_vs_halfPDA,
#                lab = NA,
#                x = 'log2.FC.',
#                y = 'pvalue',
#                #selectLab = c(""),
#                pCutoff = 0.0088382, 
#                FCcutoff = 1, 
#                labSize = 0,
#                legendLabels=c('Not sig.',
#                               'Log (base 2) FC',
#                               'p-adj < 0.05',
#                               'p-adj < 0.05 & Log (base 2) FC'),
#                legendPosition = 'right',
#                legendLabSize = 12,
#                legendIconSize = 4,
#                cutoffLineType = 'twodash',
#                cutoffLineWidth = 1,
#                pointSize = 4, 
#                title = "1/2 PDA vs PDA; metabolomics", 
#                subtitle = NULL) +
#  theme_minimal(base_size = 20)

#onetenth_vs_fullPDA_sig <- onetenth_vs_fullPDA %>%
#  rownames_to_column(var = "ID") %>%
#  filter(padjusted < 0.05) %>%
#  inner_join(meta_metab, by = "ID") %>%
#  group_by(Formula) %>%
#  slice_min(padjusted, n = 1, with_ties = FALSE) %>%
#  ungroup()
#onetenth_vs_halfPDA_sig <- onetenth_vs_halfPDA %>%
#  rownames_to_column(var = "ID") %>%
#  filter(padjusted < 0.05) %>%
#  inner_join(meta_metab, by = "ID") %>%
#  group_by(Formula) %>%
#  slice_min(padjusted, n = 1, with_ties = FALSE) %>%
#  ungroup()
#half_vs_fullPDA_sig <- half_vs_fullPDA %>%
#  rownames_to_column(var = "ID") %>%
#  filter(padjust < 0.05) %>%
#  inner_join(meta_metab, by = "ID") %>%
#  group_by(Formula) %>%
#  slice_min(padjust, n = 1, with_ties = FALSE) %>%
#  ungroup()

meta_metab <- read.csv("./metabolites_annotations.csv") %>%
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
    direction = if_else(log2FC > 0, "Up", "Down")
  )

onetenth_vs_halfPDA_sig <- onetenth_vs_halfPDA %>%
  rownames_to_column(var = "ID") %>%
  filter(pvalue < 0.05 & abs(log2FC) >= 1) %>%
  inner_join(meta_metab, by = "ID") %>%
  ungroup() %>%
  mutate(
    contrast = "1/10 vs Half",
    direction = if_else(log2FC > 0, "Up", "Down")
  )

half_vs_fullPDA_sig <- half_vs_fullPDA %>%
  rownames_to_column(var = "ID") %>%
  filter(pvalue < 0.05 & abs(log2FC) >= 1) %>%
  inner_join(meta_metab, by = "ID") %>%
  ungroup() %>%
  mutate(
    contrast = "Half vs Full",
    direction = if_else(log2FC > 0, "Up", "Down")
  )

all_sig <- bind_rows(
  onetenth_vs_fullPDA_sig,
  onetenth_vs_halfPDA_sig,
  half_vs_fullPDA_sig
)

plot_df <- all_sig %>%
  filter(Super_class %in% c("Carbohydrates", "Fatty Acyls", "Glycerophospholipids", "Organic acids", "Nucleic acids")) %>%
  group_by(Name, contrast, direction) %>%
  slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    subclass = if_else(is.na(Sub_class) | Sub_class == "", "Unassigned", Sub_class), 
    Super_class = factor(
      Super_class,
      levels = c("Organic acids", "Nucleic acids", "Carbohydrates", "Glycerophospholipids", "Fatty Acyls")
    )
  ) %>%
  count(contrast, Super_class, direction, name = "n") %>%
  mutate(n_plot = if_else(direction == "Down", -n, n)) %>%
  group_by(contrast, Super_class) %>%
  mutate(total = sum(abs(n_plot))) %>%
  ungroup()

ggplot(plot_df, aes(x = Super_class, y = n_plot, fill = direction)) +
  geom_col(width = 0.75) +
  coord_flip() +
  facet_grid(~contrast) +
  scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c("Up" = "#D55E00", "Down" = "#0072B2")) +
  labs(
    x = "",
    y = "Number of significant metabolites",
    fill = "Direction",
    title = "Up- and down-regulated metabolites by subclass"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.major.y = element_blank()
  )

library(tidyverse)

scatter_df <- all_sig %>%
  filter(Super_class %in% c("Organic acids", "Carbohydrates", "Glycerophospholipids", "Fatty Acyls", "Nucleic acids")) %>%
  group_by(Name, contrast, direction) %>%
  slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Super_class = factor(
      Super_class,
      levels = c("Organic acids", "Carbohydrates", "Nucleic acids", "Glycerophospholipids", "Fatty Acyls")
    ),
    direction = if_else(log2FC > 0, "Up", "Down")
  )

ggplot(scatter_df, aes(x = Super_class, y = log2FC, color = direction)) +
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
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
