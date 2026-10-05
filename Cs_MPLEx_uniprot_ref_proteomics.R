library(dplyr)
library(tidyverse)
library(proDA)
library(DreamAI)
library(sva)
library(DEP)
library(gprofiler2)
library(ComplexHeatmap)
library(PCAtools)
library(EnhancedVolcano)
library(scales)
library(factoextra)
library(RColorBrewer)
library(biomaRt)
library(rlang)
library(aLFQ)
library(limma)
library(clusterProfiler)
library(eulerr)
library(UpSetR)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)


###playing with peptide data and rolling up the intensities to protein level.
setwd("C:/Users/dawa726/OneDrive - PNNL/Desktop/BRAVE_PD/Cs_MPLEx_peptides/Cs_MPLEx_UniProt_reference/cs_Mplex_MBR_5/")
mod_pep <- read_tsv("combined_modified_peptide.tsv") %>%
  setNames(., gsub(" Intensity", "", names(.))) %>%
  setNames(., gsub(" ", ".", names(.)))

meta <- mod_pep[,1:15]

meta2 <- meta %>%
  dplyr::select(`Modified.Sequence`, Protein)

peptide <- mod_pep %>%
  filter(grepl("COLSU", Protein)) %>%
  dplyr::select(Modified.Sequence, contains("_")) %>%
  dplyr::select(-contains("Spectral"),
                -contains("MaxLFQ"),
                -contains("Match"),
                -contains("Localization")) %>%
  column_to_rownames(var = "Modified.Sequence") 

peptide_long <- peptide %>%
  rownames_to_column(var = "Modified.Sequence") %>%
  pivot_longer(!Modified.Sequence, names_to = "SampleID", values_to = "Intensity") %>%
  filter(Intensity > 0) %>%
  mutate(
    num_4 = str_count(Modified.Sequence, "\\[57\\.0215\\]"),
    num_C = str_count(Modified.Sequence, "C"),
    at_N = str_detect(Modified.Sequence, "^n\\[57\\.0215\\]"),
    at_C = str_detect(Modified.Sequence, "C\\[57\\.0215\\]"),
    at_K = str_detect(Modified.Sequence, "K\\[57\\.0215\\]")
  ) %>%
  mutate(Label = case_when(num_4 == 0 & num_C == 0 & at_N == F & at_K == 0 ~ "Cannot be alkylated", 
                           num_4 == num_C & at_C == T & at_N == F & at_K == F ~ "alkylated",
                           num_C > num_4 & at_N == F & at_K == F ~ "under-alkylated",
                           TRUE ~ "over-alkylated"
  )) %>%
  dplyr::select(-num_4,-num_C,-at_N,-at_C,-at_K) %>%
  filter(Label != "over-alkylated")

peptide_long %>%
  mutate(Intensity = na_if(Intensity, 0)) %>%
  mutate(log2Intensity = log2(Intensity)) %>% 
  mutate(condition = as.character(case_when(grepl("full", SampleID) ~ "PDA", 
                                            grepl("half", SampleID) ~ "half-PDA", 
                                            grepl("one", SampleID) ~ "OneTenth-PDA"))) %>% 
  filter(!is.na(log2Intensity)) %>% 
  group_by(`Modified.Sequence`, condition) %>%
  add_count(name = "n") %>% 
  distinct(SampleID, n, condition, log2Intensity) %>%
  group_by(`Modified.Sequence`) %>%
  summarize(median = median(log2Intensity), sd=sd(log2Intensity)) %>%
  ggplot() + 
  aes(x = median) +
  geom_histogram(fill = "#87CEFA", binwidth = 0.3, alpha = 0.6) +
  labs(title = "Median Intensity - peptide",
       x = "median log2intensity",
       y = "Number of Proteins") +
  scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 1000)) +
  #scale_x_continuous(limits = c(17.5, 37.5), breaks = seq(17.5, 37.5, 2.5)) +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("median_peptide_intensities.png", width = 9, height = 6, bg = "white")

protein <- peptide_long %>%
  inner_join(., meta2) %>%
  distinct(Protein, `Modified.Sequence`, SampleID, Intensity) %>%
  filter(!is.na(Intensity)) %>%
  mutate(Intensity = log2(Intensity)) %>%
  dplyr::rename(protein_list = Protein,
                sample_list = SampleID,
                id = `Modified.Sequence`,
                quant = Intensity) %>%
  as.list() %>%
  iq::fast_MaxLFQ(.) %>%
  .[[1]]
protein <- as.data.frame(protein)

#calculates the number of proteins identified in each column.
protein_ids <- nrow(protein)-colSums(is.na(protein))

protein_long <- protein %>%
  rownames_to_column(var = "Protein") %>%
  pivot_longer(-Protein, names_to = "SampleID",
               values_to = "Intensity") %>%
  mutate(condition = as.character(case_when(grepl("full", SampleID) ~ "PDA",
                                            grepl("half", SampleID) ~ "half-PDA",
                                            grepl("one", SampleID) ~ "OneTenth-PDA"))) %>%
  filter(!is.na(Intensity))

protein_long %>%
  group_by(Protein, condition) %>%
  add_count(name = "n") %>%
  distinct(SampleID, n, condition, Intensity) %>%
  filter(n > 2) %>%
  group_by(Protein) %>%
  summarize(median = median(Intensity), sd=sd(Intensity)) %>%
  ggplot() + 
  aes(x = median) +
  geom_histogram(fill = "#4682B4", binwidth = 0.3, alpha = 0.6) +
  labs(title = "Median Intensity - Protein",
       x = "median log2intensity",
       y = "Number of Proteins") +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, 100)) +
  #scale_x_continuous(limits = c(17.5, 37.5), breaks = seq(17.5, 37.5, 2.5)) +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("median_protein_intensities.png", width = 9, height = 6, bg = "white")

protein_long_filtered <- protein_long %>%
  group_by(Protein, condition) %>%
  add_count(name = "n") %>% 
  distinct(SampleID, n, condition, Intensity) %>%
  filter(n > 2)

length(unique(protein_long_filtered$Protein))
#4590

protein_input <- protein_long_filtered %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  #filter(!is.na(n)) %>%
  dplyr::select(SampleID, condition, n) %>% 
  distinct() %>%
  ungroup() %>%
  group_by(condition) %>%
  summarize(
    mean_n = mean(n), 
    sd_n = sd(n))

peptide_input <- peptide_long %>%
  mutate(condition = as.character(case_when(grepl("full", SampleID) ~ "PDA",
                                            grepl("half", SampleID) ~ "half-PDA",
                                            grepl("one", SampleID) ~ "OneTenth-PDA"))) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  #filter(!is.na(n)) %>%
  dplyr::select(SampleID, condition, n) %>% 
  distinct() %>%
  ungroup() %>%
  group_by(condition) %>%
  summarize(
    mean_n = mean(n),  
    sd_n = sd(n))

combined_stats <- bind_rows(
  protein_input %>% mutate(type = "Protein"),
  peptide_input %>% mutate(type = "Peptide")) %>% 
  mutate(condition = factor(condition, levels = c("PDA", "half-PDA", "OneTenth-PDA")))

ggplot(combined_stats, aes(x = condition, y = mean_n, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  geom_errorbar(aes(ymin = mean_n - sd_n, ymax = mean_n + sd_n), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 32000), breaks = seq(0, 32000, 4000))+
  scale_fill_manual(values = c("Protein" = "#4682B4", "Peptide" = "#87CEFA")) +
  ggtitle("Protein and Peptide Identifications") +
  xlab("Condition") +
  ylab("mean identifications") +
  theme_classic(base_size = 20)
#ggsave("Protein_and_Peptide_Identifications.png", width = 8, height = 5, bg = "white")

protein_wide_filtered <- protein_long_filtered %>% 
  ungroup() %>% 
  dplyr::select(Protein, SampleID, Intensity) %>% 
  pivot_wider(names_from = SampleID,  
              values_from = Intensity) %>% 
  column_to_rownames(var = "Protein")

protein_wide_filtered <- as.data.frame(median_normalization(as.matrix(protein_wide_filtered)))

###Using "Colletotrichum higginsianum; chig" annotations - Blast transfer
###Blast
#system('"C:/Program Files/NCBI/blast-2.17.0+/bin/makeblastdb.exe" -in "../../../../../../../references/Colletotrichum_sublineola/UP000092177_759273.fasta/UP000092177_759273.fasta" -dbtype prot -out uniprot_chig_db')
#system('"C:/Program Files/NCBI/blast-2.17.0+/bin/blastp.exe" -query "./Cs_UniProt.fasta" -db uniprot_chig_db -out blastp_results.txt -outfmt 6')

blast1 <- read.table("blastp_results.txt", header = FALSE, sep = "\t")
colnames(blast1) <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)

#blast1 %>%
#  ggplot(aes(x = bitscore, y = length)) +
#  geom_point(alpha = 0.6) +
#scale_y_log10() +
#  labs(title = "Relationship between Bit Score and E-value (Top Hits)",
#       x = "Bit Score",
#       y = "length") +
#  theme_minimal(base_size = 20)

#blast1 %>%
#  ggplot(aes(x = bitscore, y = evalue)) +
#  geom_point(alpha = 0.6) +
#  scale_y_log10() +
#  labs(title = "Relationship between Bit Score and E-value (Top Hits)",
#       x = "Bit Score",
#       y = "log10(evalue)") +
#  theme_minimal(base_size = 20)

top_hits_blast <- blast1 %>%
  group_by(qseqid) %>%
  filter(evalue < 0.001) %>%
  filter(bitscore >= 50) %>%
  slice_max(order_by = pident , n = 1, with_ties = FALSE) %>%
  ungroup()

annota_final <- top_hits_blast %>%
  dplyr::select(qseqid, sseqid) %>%
  distinct()
colnames(annota_final) <- c("COLSU", "COLHI")  

annota_final <- annota_final %>%
  mutate(uniprot_colsu = sub(".*\\|(.*)\\|.*", "\\1", COLSU))


combined_proteins <- protein_wide_filtered %>%
  rownames_to_column(var = "Protein") %>%
  pivot_longer(!Protein, names_to = "SampleID", values_to = "Intensity") %>%
  filter(!is.na(Intensity)) %>%
  mutate(condition = as.character(case_when(grepl("full", SampleID) ~ "PDA", 
                                            grepl("half", SampleID) ~ "half-PDA", 
                                            grepl("one", SampleID) ~ "OneTenth-PDA"))) %>%
  group_by(Protein, condition) %>%
  add_count(name = "n") %>%
  filter(n > 2) %>%
  dplyr::select(Protein, condition) %>%
  dplyr::rename(COLSU = Protein) %>%
  distinct()
#write.csv(combined_proteins, "combined_proteins.csv")

sets <- combined_proteins %>%
  filter(!is.na(COLSU), !is.na(condition)) %>%
  group_by(condition) %>%
  summarise(members = list(unique(COLSU)), .groups = "drop") %>%
  tibble::deframe()
sets <- sets[c("PDA", "half-PDA", "OneTenth-PDA")]
fit <- euler(sets)
plot(
  fit,
  quantities = TRUE,
  labels = FALSE,               
  legend = list(side = "right"), 
  fills = list(
    fill  = c("dodgerblue2", "#FF6A6A", "#90EE90"),
    alpha = 0.5
  ),
  edges  = list(lwd = 1)
)

protein_counts <- combined_proteins %>%
  group_by(COLSU) %>%
  summarise(n_conditions = n_distinct(condition), .groups = "drop")

unique_proteins <- combined_proteins %>%
  full_join(protein_counts, by = "COLSU") %>%
  filter(n_conditions == 1) %>%
  dplyr::select(-n_conditions)
unique_condition_proteins <- unique_proteins %>%
  left_join(annota_final, by = "COLSU") %>%
  group_split(condition) %>%
  set_names(map_chr(., ~ unique(.x$condition)))

OneTenth_PDA_pres <- unique_proteins %>%
  filter(condition == "OneTenth-PDA") %>%
  inner_join(annota_final, by = "COLSU") %>%
  filter(!is.na(COLHI)) %>%
  mutate(COLHI = sub(".*\\|(.*)\\|.*", "\\1", COLHI))

top_10_OneTenth_PDA_pres <- protein_wide_filtered %>%
  rownames_to_column(var = "Protein") %>%
  pivot_longer(!Protein, names_to = "SampleID", values_to = "Intensity") %>%
  filter(!is.na(Intensity)) %>%
  mutate(Intensity = 2^Intensity) %>%
  mutate(condition = as.character(case_when(grepl("full", SampleID) ~ "PDA", 
                                            grepl("half", SampleID) ~ "half-PDA", 
                                            grepl("one", SampleID) ~ "OneTenth-PDA"))) %>%
  group_by(Protein, condition) %>%
  summarise(log2_mean_Intensity = log2(mean(Intensity, na.rm = TRUE)),
            .groups = "drop") %>%
  group_by(condition) %>%
  mutate(rank_condition = dense_rank(desc(log2_mean_Intensity))) %>%
  ungroup() %>%
  filter(condition == "OneTenth-PDA") %>% 
  dplyr::rename(COLSU = Protein) %>%
  semi_join(OneTenth_PDA_pres, by = "COLSU") %>% 
  arrange(rank_condition) %>% 
  slice_head(n = 10) %>%
  mutate(COLSU_2 = sub(".*\\|(.*)\\|.*", "\\1", COLSU))

OneTenth_PDA_pres_ranks <- protein_wide_filtered %>%
  rownames_to_column(var = "Protein") %>%
  pivot_longer(!Protein, names_to = "SampleID", values_to = "Intensity") %>%
  filter(!is.na(Intensity)) %>%
  mutate(Intensity = 2^Intensity) %>%
  mutate(condition = as.character(case_when(grepl("full", SampleID) ~ "PDA", 
                                            grepl("half", SampleID) ~ "half-PDA", 
                                            grepl("one", SampleID) ~ "OneTenth-PDA"))) %>%
  group_by(Protein, condition) %>%
  summarise(log2_mean_Intensity = log2(mean(Intensity, na.rm = TRUE)),
            .groups = "drop") %>%
  group_by(condition) %>%
  mutate(rank_condition = dense_rank(desc(log2_mean_Intensity))) %>%
  ungroup() %>%
  filter(condition == "OneTenth-PDA") %>% 
  dplyr::rename(COLSU = Protein) %>%
  semi_join(OneTenth_PDA_pres, by = "COLSU") %>% 
  arrange(rank_condition) %>% 
  mutate(COLSU_2 = sub(".*\\|(.*)\\|.*", "\\1", COLSU))

protein_wide_filtered %>%
  rownames_to_column(var = "Protein") %>%
  pivot_longer(!Protein, names_to = "SampleID", values_to = "Intensity") %>%
  filter(!is.na(Intensity)) %>%
  mutate(Intensity = 2^Intensity) %>%
  mutate(condition = as.character(case_when(grepl("full", SampleID) ~ "PDA", 
                                            grepl("half", SampleID) ~ "half-PDA", 
                                            grepl("one", SampleID) ~ "OneTenth-PDA"))) %>%
  group_by(Protein, condition) %>%
  summarise(log2_mean_Intensity = log2(mean(Intensity, na.rm = TRUE)),
            .groups = "drop") %>%
  group_by(condition) %>%
  mutate(rank_condition = dense_rank(desc(log2_mean_Intensity))) %>%
  ungroup() %>%
  ggplot(aes(x = rank_condition, y = log2_mean_Intensity, color = condition)) +
  geom_point(size = 3) +
  geom_point(
    data = . %>% 
      filter(Protein %in% OneTenth_PDA_pres$COLSU),
    aes(x = rank_condition, y = log2_mean_Intensity),   # keep mapping explicit
    inherit.aes = FALSE,
    color = "black", size = 2, alpha = 0.2
  ) +
  labs(x = "Rank within condition",
       y = "log2(mean Intensity)",
       title = "Protein abundance rank vs mean log-intensity") +
  coord_cartesian(ylim = c(20, 35)) +
  scale_color_manual(values = c(
    "PDA" = "dodgerblue2",
    "half-PDA" = "#FF6A6A",
    "OneTenth-PDA" = "#90EE90"
  )) +
  geom_label_repel(
    data = top_10_OneTenth_PDA_pres,
    aes(x = rank_condition, y = log2_mean_Intensity, label = top_10_OneTenth_PDA_pres$COLSU_2),
    inherit.aes = FALSE,
    color = "black",
    fill = "white",
    size = 3,
    max.overlaps = 25,
    point.padding = 0.1,
    label.padding = 0.1,
    nudge_x = 5,
    nudge_y = -5
  ) +
  theme_minimal(base_size = 20) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  guides(color = guide_legend(title = "Condition"))

meta4 <- data.frame(colnames(protein_wide_filtered)) %>%
  mutate(SampleID = `colnames.protein_wide_filtered.`) %>%
  mutate(condition = as.character(case_when(grepl("full", SampleID) ~ "PDA",
                                            grepl("half", SampleID) ~ "half-PDA",
                                            grepl("one", SampleID) ~ "OneTenth-PDA"))) %>%
  dplyr::select(SampleID, condition)

meta4$condition <- factor(
  meta4$condition,
  levels = c("PDA", "half-PDA", "OneTenth-PDA")
)

protein_imputed <- DreamAI(protein_wide_filtered, 
                           k = 5, maxiter_MF = 10, ntree = 100, 
                           maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2), 
                           gamma_ADMIN = NA, gamma = 50, CV = FALSE, 
                           fillmethod = "row_mean", maxiter_RegImpute = 10, 
                           conv_nrmse = 1e-06, iter_SpectroFM = 40, 
                           method = "KNN",out = c("KNN"))

protein_imputed <- as.data.frame(protein_imputed$KNN)
protein_imputed <- protein_imputed[, meta4$SampleID]

set.seed(101)
pca1 <- PCAtools::pca(protein_imputed, scale = T, center = T)
pca1$metadata <- meta4

eigencorplot(pca1,
             components = getComponents(pca1),
             metavars = colnames(pca1$metadata),
             scale = FALSE,
             col = c("darkgreen", "white", "purple"), cexCorval = "black")
PCAtools::biplot(pca1, x = "PC1", y =  "PC2", 
                 lab = meta4$SampleID,
                 showLoadings = F, 
                 sizeLoadingsNames = F,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = 'right', 
                 colby = 'condition', colkey = c(
                   "PDA" = "dodgerblue2",
                   "half-PDA" = "#FF6A6A",
                   "OneTenth-PDA" = "#90EE90"),
                 encircle = F,
                 ellipse = TRUE,
                 ellipseLevel = 0.7,
                 xlim = c(-50, 100),
                 ylim = c(-50, 75),
                 pointSize = 5, ellipseAlpha = 0.2)

meta5 <- data.frame(colnames(protein_wide_filtered)) %>%
  mutate(SampleID = `colnames.protein_wide_filtered.`) %>%
  mutate(condition = as.character(case_when(grepl("full", SampleID) ~ "PDA",
                                            grepl("half", SampleID) ~ "half_PDA",
                                            grepl("one", SampleID) ~ "OneTenth_PDA"))) %>%
  dplyr::select(SampleID, condition)

lev <- c("PDA", "half_PDA", "OneTenth_PDA")
f <- factor(meta5$condition, levels=lev)
design <- model.matrix(~0+f)
design
colnames(design) <- lev
fit <- lmFit(protein_wide_filtered, design)
contrast.matrix <- makeContrasts(half_PDA-PDA, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
half_PDA_vs_PDA <- topTable(fit2, adjust.method = "BH", sort.by = "P", number = 7000)
EnhancedVolcano(half_PDA_vs_PDA,
                lab = NA,
                x = 'logFC',
                y = 'P.Value',
                #selectLab = c(""),
                pCutoff = 1.331183e-03, 
                FCcutoff = 1, 
                labSize = 0,
                legendLabels=c('Not sig.',
                               'Log (base 2) FC',
                               'p-adj < 0.05',
                               'p-adj < 0.05 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 1,
                pointSize = 4, 
                #xlim = c(-10, 10), 
                ylim = c(0, 10), 
                title = "half-PDA vs PDA", 
                subtitle = NULL,
                #labCol = 'black',
                #labFace = 'bold',
                #boxedLabels = TRUE,
                #colAlpha = 4/5,
                #drawConnectors = TRUE,
                #widthConnectors = 0.75,
                #colConnectors = 'black', #max.overlaps = 50
) +
  theme_minimal(base_size = 20)
#ggsave("half_PDA_vs_PDA.png", width = 15, height = 10, bg = "white")

contrast.matrix2 <- makeContrasts(OneTenth_PDA-PDA, levels = design)
fit3 <- contrasts.fit(fit, contrast.matrix2)
fit3 <- eBayes(fit3)
OneTenth_PDA_vs_PDA <- topTable(fit3, adjust.method = "BH", sort.by = "P", number = 7000)
EnhancedVolcano(OneTenth_PDA_vs_PDA,
                lab = NA,
                x = 'logFC',
                y = 'P.Value',
                #selectLab = c(""),
                pCutoff = 0.014054955, 
                FCcutoff = 1, 
                labSize = 0,
                legendLabels=c('Not sig.',
                               'Log (base 2) FC',
                               'p-adj < 0.05',
                               'p-adj < 0.05 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 1,
                pointSize = 4, 
                #xlim = c(-10, 10), 
                ylim = c(0, 10), 
                title = "Onetenth-PDA vs PDA", 
                subtitle = NULL,
                #labCol = 'black',
                #labFace = 'bold',
                #boxedLabels = TRUE,
                #colAlpha = 4/5,
                #drawConnectors = TRUE,
                #widthConnectors = 0.75,
                #colConnectors = 'black', #max.overlaps = 50
) +
  theme_minimal(base_size = 20)
#ggsave("OneTenth_PDA_vs_PDA.png", width = 15, height = 10, bg = "white")

contrast.matrix3 <- makeContrasts(OneTenth_PDA-half_PDA, levels = design)
fit4 <- contrasts.fit(fit, contrast.matrix3)
fit4 <- eBayes(fit4)
OneTenth_PDA_vs_half_PDA <- topTable(fit4, 
                                     adjust.method = "BH", 
                                     sort.by = "P", 
                                     number = 7000)
EnhancedVolcano(OneTenth_PDA_vs_half_PDA,
                lab = NA,
                x = 'logFC',
                y = 'P.Value',
                #selectLab = c(""),
                pCutoff = 0.01534861, 
                FCcutoff = 1, 
                labSize = 0,
                legendLabels=c('Not sig.',
                               'Log (base 2) FC',
                               'p-adj < 0.05',
                               'p-adj < 0.05 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 1,
                pointSize = 4, 
                #xlim = c(-10, 10), 
                ylim = c(0, 10), 
                title = "Onetenth-PDA vs half-PDA", 
                subtitle = NULL,
                #labCol = 'black',
                #labFace = 'bold',
                #boxedLabels = TRUE,
                #colAlpha = 4/5,
                #drawConnectors = TRUE,
                #widthConnectors = 0.75,
                #colConnectors = 'black', #max.overlaps = 50
) +
  theme_minimal(base_size = 20)
#ggsave("OneTenth_PDA_vs_half_PDA.png", width = 15, height = 10, bg = "white")
#write.csv(OneTenth_PDA_vs_half_PDA, file = "OneTenth_PDA_vs_half_PDA.csv")
#write.csv(OneTenth_PDA_vs_PDA, file = "OneTenth_PDA_vs_PDA.csv")
#write.csv(half_PDA_vs_PDA, file = "half_PDA_vs_PDA.csv")

protein_imputed %>%
  rownames_to_column(var = "ProteinID") %>%
  pivot_longer(!ProteinID, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(
    condition = case_when(
      grepl("full", SampleID) ~ "PDA",
      grepl("half", SampleID) ~ "half-PDA",
      grepl("one", SampleID) ~ "OneTenth-PDA"
    ),
    condition = factor(condition, levels = c("PDA", "half-PDA", "OneTenth-PDA")),
    Intensity = 2^Intensity
  ) %>%
  filter(!is.na(Intensity), !is.na(condition)) %>%
  group_by(ProteinID, condition) %>%
  summarise(
    Avg = mean(Intensity, na.rm = TRUE),
    CV  = sd(Intensity, na.rm = TRUE) / Avg,
    .groups = "drop"
  ) %>%
  filter(!is.na(CV), is.finite(CV)) %>% 
  ggplot(aes(x = condition, y = CV, fill = condition)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = median, geom = "point", color = "black") +
  scale_fill_manual(values = c(
    "PDA" = "dodgerblue2",
    "half-PDA" = "#FF6A6A",
    "OneTenth-PDA" = "#90EE90"
  )) +
  ylab("Coefficient of Variation") +
  xlab("") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 0.5)) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")
#ggsave("CVs_per_condition.png", bg = "white")

### Step 1. Filter significant proteins for each comparison; Find shared and unique significant proteins
sig_OneTenth_PDA_vs_half_PDA <- OneTenth_PDA_vs_half_PDA %>%
  rownames_to_column(var = "ProteinID") %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)

sig_OneTenth_PDA_vs_PDA <- OneTenth_PDA_vs_PDA %>%
  rownames_to_column(var = "ProteinID") %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)

sig_half_PDA_vs_PDA <- half_PDA_vs_PDA %>%
  rownames_to_column(var = "ProteinID") %>%
  filter(adj.P.Val < 0.05, abs(logFC) >= 1)

sets_DAP <- list(
  `OneTenth vs Half` = unique(as.character(sig_OneTenth_PDA_vs_half_PDA$ProteinID)),
  `OneTenth vs PDA`  = unique(as.character(sig_OneTenth_PDA_vs_PDA$ProteinID)),
  `half vs PDA`  = unique(as.character(sig_half_PDA_vs_PDA$ProteinID))
)

#fit_DAP <- euler(sets_DAP)
#plot(
#  fit_DAP,
#  quantities = TRUE,
#  labels = FALSE,               
#  legend = list(side = "right"), 
#  fills = list(
#    fill  = c("#CDB5CD", "slategray2", "grey"
#  ),
#  edges  = list(lwd = 1)
#))

#shared_sig_proteins %>%
#  dplyr::filter(
#    !is.na(logFC_OneTenth_half),
#    !is.na(logFC_OneTenth_PDA),
#    logFC_OneTenth_half != 0,
#    logFC_OneTenth_PDA != 0,
#    sign(logFC_OneTenth_half) != sign(logFC_OneTenth_PDA)
#  ) %>%
#  dplyr::pull(ProteinID)

upset(fromList(sets_DAP),
      order.by = "freq",
      sets.bar.color = "grey40",
      mainbar.y.label = "Intersection size",
      sets.x.label = "Set size", text.scale = 1.5)

shared_sig_proteins <- inner_join(sig_OneTenth_PDA_vs_half_PDA,
                                  sig_OneTenth_PDA_vs_PDA, 
                                  by = "ProteinID",
                                  suffix = c("_OneTenth_half", "_OneTenth_PDA"))


unique_sig_onetenth_vs_half <- sig_OneTenth_PDA_vs_half_PDA %>%
  filter(!ProteinID %in% shared_sig_proteins$ProteinID)

unique_sig_onetenth_vs_PDA <- sig_OneTenth_PDA_vs_PDA %>%
  filter(!ProteinID %in% shared_sig_proteins$ProteinID)

chig_rankedProteins_OneTenth_PDA_vs_PDA_sig <- OneTenth_PDA_vs_PDA %>%
  rownames_to_column(var = "COLSU") %>%
  filter(adj.P.Val < 0.05) %>%
  filter(!is.na(COLSU)) %>%
  mutate(ranking_metric = -log10(P.Value) * sign(logFC)) %>%
  group_by(COLSU) %>%
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>%
  inner_join(annota_final, by = "COLSU") %>%
  dplyr::select(COLHI, ranking_metric) %>%
  filter(!is.na(COLHI)) %>%
  mutate(COLHI = sub(".*\\|(.*)\\|.*", "\\1", COLHI)) %>%
  group_by(COLHI) %>% 
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>%  # ensures uniqueness
  arrange(desc(ranking_metric)) %>%
  deframe()

anyDuplicated(names(chig_rankedProteins_OneTenth_PDA_vs_PDA_sig))

head(chig_rankedProteins_OneTenth_PDA_vs_PDA_sig)

chig_gsea_OneTenth_PDA_vs_PDA <- clusterProfiler::gseKEGG(
  geneList = chig_rankedProteins_OneTenth_PDA_vs_PDA_sig,
  organism = "chig",        # change this to your species' KEGG code
  keyType = "uniprot",     # because you have UniProt IDs
  eps = 0.0,
  minGSSize = 0,
  maxGSSize = 1000,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  verbose = FALSE
)

chig_gsea_OneTenth_PDA_vs_PDA@result %>% 
  filter(pvalue < 0.05) %>%
  ggplot(aes(x = NES,
             y = reorder(Description, NES),
             size = setSize,
             color = pvalue)) +
  geom_point(alpha = 0.65) +
  #scale_size(limits = c(2, 16), breaks = seq(2, 16, by = 2)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "OneTenth_PDA_vs_PDA",
       x = "NES",
       y = "Pathway Description",
       color = "pvalue",
       size = "Gene Count") +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1)) +
  theme_minimal(base_size = 20)

chig_rankedProteins_OneTenth_PDA_vs_half_PDA_sig <- OneTenth_PDA_vs_half_PDA %>%
  rownames_to_column(var = "COLSU") %>%
  filter(adj.P.Val < 0.05) %>%
  filter(!is.na(COLSU)) %>%
  mutate(ranking_metric = -log10(P.Value) * sign(logFC)) %>%
  group_by(COLSU) %>%
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>%
  inner_join(annota_final, by = "COLSU") %>%
  dplyr::select(COLHI, ranking_metric) %>%
  filter(!is.na(COLHI)) %>%
  mutate(COLHI = sub(".*\\|(.*)\\|.*", "\\1", COLHI)) %>%
  group_by(COLHI) %>% 
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>%  # ensures uniqueness
  arrange(desc(ranking_metric)) %>%
  deframe()

anyDuplicated(names(chig_rankedProteins_OneTenth_PDA_vs_half_PDA_sig))

head(chig_rankedProteins_OneTenth_PDA_vs_half_PDA_sig)

chig_gsea_OneTenth_PDA_vs_half_PDA <- clusterProfiler::gseKEGG(
  geneList = chig_rankedProteins_OneTenth_PDA_vs_half_PDA_sig,
  organism = "chig",        # change this to your species' KEGG code
  keyType = "uniprot",     # because you have UniProt IDs
  eps = 0.0,
  minGSSize = 0,
  maxGSSize = 1000,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  verbose = FALSE
)

chig_gsea_OneTenth_PDA_vs_half_PDA@result %>%
  filter(pvalue < 0.05) %>%
  ggplot(aes(x = NES,
             y = reorder(Description, NES),
             size = setSize,
             color = pvalue)) +
  geom_point(alpha = 0.65) +
  #scale_size(limits = c(2, 16), breaks = seq(2, 16, by = 2)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "OneTenth_PDA_vs_half_PDA",
       x = "NES",
       y = "Pathway Description",
       color = "pvalue",
       size = "Gene Count") +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1)) +
  theme_minimal(base_size = 20)

chig_rankedProteins_half_PDA_vs_PDA_sig <- half_PDA_vs_PDA %>%
  rownames_to_column(var = "COLSU") %>%
  filter(adj.P.Val < 0.05) %>%
  filter(!is.na(COLSU)) %>%
  mutate(ranking_metric = -log10(P.Value) * sign(logFC)) %>%
  group_by(COLSU) %>%
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>%
  inner_join(annota_final, by = "COLSU") %>%
  dplyr::select(COLHI, ranking_metric) %>%
  filter(!is.na(COLHI)) %>%
  mutate(COLHI = sub(".*\\|(.*)\\|.*", "\\1", COLHI)) %>%
  group_by(COLHI) %>% 
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>%  # ensures uniqueness
  arrange(desc(ranking_metric)) %>%
  deframe()

anyDuplicated(names(chig_rankedProteins_half_PDA_vs_PDA_sig))

tail(chig_rankedProteins_half_PDA_vs_PDA_sig)

chig_gsea_half_PDA_vs_PDA <- clusterProfiler::gseKEGG(
  geneList = chig_rankedProteins_half_PDA_vs_PDA_sig,
  organism = "chig",        # change this to your species' KEGG code
  keyType = "uniprot",     # because you have UniProt IDs
  eps = 0.0,
  minGSSize = 0,
  maxGSSize = 1000,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  verbose = FALSE
)

chig_gsea_half_PDA_vs_PDA@result %>%
  filter(pvalue < 0.05) %>%
  ggplot(aes(x = NES,
             y = reorder(Description, NES),
             size = setSize,
             color = pvalue)) +
  geom_point(alpha = 0.65) +
  #scale_size(limits = c(2, 16), breaks = seq(2, 16, by = 2)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "half_PDA_vs_PDA",
       x = "NES",
       y = "Pathway Description",
       color = "pvalue",
       size = "Gene Count") +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1)) +
  theme_minimal(base_size = 20)

res_PDA <- chig_gsea_OneTenth_PDA_vs_PDA@result %>%
  select(Description, NES, pvalue, p.adjust, setSize) %>%
  rename(
    NES_PDA    = NES,
    pvalue_PDA = pvalue,
    padj_PDA   = p.adjust,
    setSize_PDA= setSize
  )

res_half <- chig_gsea_OneTenth_PDA_vs_half_PDA@result %>%
  select(Description, NES, pvalue, p.adjust, setSize) %>%
  rename(
    NES_half    = NES,
    pvalue_half = pvalue,
    padj_half   = p.adjust,
    setSize_half= setSize
  )

gsea_both <- full_join(res_PDA, res_half, by = "Description")

gsea_both %>% 
  pivot_longer(
    cols = -Description,
    names_to = c(".value", "contrast"),
    names_pattern = "(NES|pvalue|padj|setSize)_(PDA|half)"
  ) %>%
  #filter(!is.na(NES)) %>%
  filter(pvalue < 0.05) %>%
  ggplot(aes(x = NES,
             y = fct_reorder(Description, NES, .fun = mean, na.rm = TRUE),
             size = setSize,
             color = pvalue,
             shape = contrast)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 1)) +
  facet_wrap(~ contrast) +
  labs(title = "KEGG GSEA (OneTenth PDA comparisons)",
       x = "NES",
       y = "Pathway Description",
       color = "pvalue",
       size = "setSize",
       shape = "contrast") +
  theme_minimal(base_size = 20)

###
#OneTenth_PDA_pres2 <- unique(OneTenth_PDA_pres$COLHI)

#bg <- annota_final %>% 
#  dplyr::select(COLHI) %>%
#  filter(!is.na(COLHI)) %>%
#  mutate(COLHI = sub(".*\\|(.*)\\|.*", "\\1", COLHI)) %>% pull(COLHI)

#enrich_result_OneTenth_PDA_pres <- enrichKEGG(
#  gene = OneTenth_PDA_pres2,
#  organism = "chig",
#  keyType = "uniprot", 
#  universe = bg)

#enrich_result_OneTenth_PDA_pres@result %>% view()

#shared_sig_proteins_new <- shared_sig_proteins %>%
#  dplyr::select(ProteinID) %>%
#  dplyr::rename(COLSU = ProteinID) %>%
#  inner_join(annota_final, by = "COLSU") %>%
#  dplyr::select(COLHI) %>%
#  mutate(COLHI = sub(".*\\|(.*)\\|.*", "\\1", COLHI)) %>% pull(COLHI)

#enrich_result_shared <- enrichKEGG(
#  gene = shared_sig_proteins_new,
#  organism = "chig",
#  keyType = "uniprot", 
#  universe = bg)

#enrich_result_shared@result %>% dplyr::filter(pvalue < 0.05) %>% view()

###Calculating Protein ranked abundance in each condition
#mod_pep2 <- read.delim("combined_modified_peptide.tsv", 
#                       check.names = FALSE) 
#sample_columns <- colnames(mod_pep2)[!grepl("PDA", colnames(mod_pep2))]
#meta3 <- mod_pep2[, sample_columns] 
#mod_pep2[mod_pep2 == 0] <- NA
#mod_pep2 <- mod_pep2 %>%
#  dplyr::select(-ends_with("Count")) %>%
#  dplyr::select(-ends_with("Type")) %>%
#  dplyr::select(c(ends_with("maxLFQ Intensity"),  colnames(meta3)))
#filtered_peptides <- peptide_long %>% 
#  dplyr::select(Modified.Sequence) %>% 
#  unique()

#mod_pep2 <- mod_pep2 %>%
#  filter(`Modified Sequence` %in% #filtered_peptides$Modified.Sequence) %>%
#  pivot_longer(cols = -colnames(meta3), 
#               names_to = "run_id", 
#               values_to = "peptide_intensity") %>%
#  mutate(concentration = "?") %>%
#  filter(peptide_intensity > 0)

#ibaqfd <- mod_pep2 %>%
#  mutate(peptide_sequence = `Peptide Sequence`) %>%
#  mutate(protein_id = Protein) %>%
#  mutate(precursor_charge = Charges) %>%
#  mutate(peptide_id = `Modified Sequence`) %>%
#  dplyr::select(protein_id, peptide_sequence, run_id, peptide_id, peptide_intensity, concentration, precursor_charge)

#iBAQ <- ProteinInference(ibaqfd, 
#                         peptide_method = "iBAQ",
#                         consensus_peptides = FALSE,
#                         consensus_proteins = FALSE,
#                         fasta = "../../../../../../references/Colletotrichum_sublineola/Cs_UniProt.fasta")

#rIBAQ_PDA <- iBAQ %>%
#  #filter(!grepl("contam_sp", protein_id)) %>%
#  mutate(condition = as.character(case_when(grepl("full", run_id) ~ "PDA",
#                                            grepl("half", run_id) ~ "half-PDA",
#                                            grepl("one", run_id) ~ "OneTenth-PDA"))) %>%
#  filter(condition == "PDA") %>%
#  group_by(protein_id) %>%
#  summarize(mean1 = mean(response)) %>%
#  mutate(Rank = rank(-mean1)) %>%
#  filter(protein_id %in% unique(protein_long_filtered$Protein))

#rIBAQ_half <- iBAQ %>%
  #filter(!grepl("contam_sp", protein_id)) %>%
#  mutate(condition = as.character(case_when(grepl("full", run_id) ~ "PDA",
#                                            grepl("half", run_id) ~ "half-PDA",
#                                            grepl("one", run_id) ~ "OneTenth-PDA"))) %>%
#  filter(condition == "half-PDA") %>%
#  group_by(protein_id) %>%
#  summarize(mean1 = mean(response)) %>%
#  mutate(Rank = rank(-mean1)) %>%
#  filter(protein_id %in% unique(protein_long_filtered$Protein))

#rIBAQ_One_tenth <- iBAQ %>%
  #filter(!grepl("contam_sp", protein_id)) %>%
#  mutate(condition = as.character(case_when(grepl("full", run_id) ~ "PDA",
#                                            grepl("half", run_id) ~ "half-PDA",
#                                            grepl("one", run_id) ~ "OneTenth-PDA"))) %>%
#  filter(condition == "OneTenth-PDA") %>%
#  group_by(protein_id) %>%
#  summarize(mean1 = mean(response)) %>%
#  mutate(Rank = rank(-mean1)) %>%
#  filter(protein_id %in% unique(protein_long_filtered$Protein))

#combined_rIBAQ <- bind_rows(
#  rIBAQ_PDA %>% mutate(condition = "PDA"),
#  rIBAQ_half %>% mutate(condition = "half-PDA"),
#  rIBAQ_One_tenth %>% mutate(condition = "OneTenth-PDA")
) %>%
#  mutate(condition = factor(condition,
#                            levels = c("PDA", "half-PDA", "OneTenth-PDA")))

#top_labels <- bind_rows(rIBAQ_PDA %>% mutate(condition = "PDA") %>% filter(Rank <= 5), rIBAQ_half %>% mutate(condition = "half-PDA") %>% filter(Rank <= 5), rIBAQ_One_tenth %>% mutate(condition = "OneTenth-PDA") %>% filter(Rank <= 5))

#annota_final3_new <- read.csv(file = "./annota_final3_new.csv", header = TRUE)

#unique_one_tenth <- combined_rIBAQ %>%
#  filter(condition == "OneTenth-PDA") %>%
#  distinct(protein_id) %>%
#  anti_join(
#    combined_rIBAQ %>%
#      filter(condition != "OneTenth-PDA") %>%
#      distinct(protein_id),
#    by = "protein_id")

#unique_one_tenth_top10 <- combined_rIBAQ %>% 
#  filter(condition == "OneTenth-PDA") %>% 
#  semi_join(unique_one_tenth, by = "protein_id") %>% 
#  arrange(Rank) %>% 
#  slice_head(n = 10) %>%
#  mutate(protein_id2 = str_extract(protein_id, "(?<=\\|)[^|]+(?=\\|)"))

#ggplot(combined_rIBAQ, aes(x = Rank, y = log2(mean1), color = condition)) +
#  geom_point(size = 3) +
#  geom_point(
#    data = combined_rIBAQ %>% 
#      filter(protein_id %in% unique_one_tenth$protein_id),
#    aes(x = Rank, y = log2(mean1)),   # keep mapping explicit
#    inherit.aes = FALSE,
#    color = "black", size = 3
#  ) +
#  geom_label_repel(
#    data = unique_one_tenth_top10,
#    aes(x = Rank, y = log2(mean1), label = protein_id2),
#    inherit.aes = FALSE,
#    color = "black",
#    fill = "white",
#    size = 3,
#    max.overlaps = 25,
#    point.padding = 0.1,
#    label.padding = 0.1,
#    nudge_x = 5,
#    nudge_y = -10
#  ) +
#  ylab("log2(mean protein abundance)") +
#  coord_cartesian(ylim = c(10, 40)) +
#  scale_color_manual(values = c(
#    "PDA" = "dodgerblue2",
#    "half-PDA" = "#FF6A6A",
#    "OneTenth-PDA" = "#90EE90"
#  )) +
#  theme_minimal(base_size = 20) +
#  theme(panel.grid.minor = element_blank(),
#        panel.grid.major = element_blank(),
#        panel.border = element_blank(),
#        axis.line = element_line()) +
#  guides(color = guide_legend(title = "Condition"))
#ggsave("ranked_abundance.png", width = 9, height = 6, bg = "white")


###sample-sample correlation
correlation_matrix <- cor(as.data.frame(lapply(protein_wide_filtered, as.numeric)), method = "pearson", use = "complete.obs")

hc <- hclust(dist(correlation_matrix))
reordered_correlation_matrix <- correlation_matrix[hc$order, hc$order]
#pheatmap(reordered_correlation_matrix_new2, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", color = colorRampPalette(c("white", "darkblue"))(100), main = "C10 Correlation heatmap", fontsize_row = 15, fontsize_col = 15, cellwidth = 17, cellheight = 17, fontsize = 15, border_color = NA, filename = "./pheatmap_SPINS_vs_DS.png")

png("sample_correlation.png", width = 2500, height = 2500, bg = "white", res = 300)
corrplot::corrplot(correlation_matrix, method = 'shade', diag = TRUE, tl.cex = 1, order = 'alphabet', col.lim = c(0,1)) %>% corrplot::corrRect(c(1,6,11,15), lwd = 2, col = "red")
dev.off()

full_PDA_Protein_matrix <- correlation_matrix[grep("full", rownames(correlation_matrix)), grep("full", colnames(correlation_matrix))]
average_correlation_full <- mean(full_PDA_Protein_matrix[upper.tri(full_PDA_Protein_matrix) | lower.tri(full_PDA_Protein_matrix)])
print(average_correlation_full)

half_PDA_Protein_matrix <- correlation_matrix[grep("half", rownames(correlation_matrix)), grep("half", colnames(correlation_matrix))]
average_correlation_half <- mean(half_PDA_Protein_matrix[upper.tri(half_PDA_Protein_matrix) | lower.tri(half_PDA_Protein_matrix)])
print(average_correlation_half)

onetenth_PDA_Protein_matrix <- correlation_matrix[grep("one", rownames(correlation_matrix)), grep("one", colnames(correlation_matrix))]
average_correlation_onetenth <- mean(onetenth_PDA_Protein_matrix[upper.tri(onetenth_PDA_Protein_matrix) | lower.tri(onetenth_PDA_Protein_matrix)])
print(average_correlation_onetenth)

###mclust and protein-protein covariance
#Credit - https://github.com/Cajun-data/nanoPOTS_Arabidopsis/tree/main/Functions
#Fulcher, J.M., Dawar, P., Balasubramanian, V.K. et al. Single-cell proteomics of Arabidopsis leaf mesophyll reveals dynamic protein responses to water-deficit stress. Genome Biol 27, 32 (2026). https://doi.org/10.1186/s13059-025-03919-6

library(mclust)
library(tidyverse)
library(stats)
library(cluster)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

### Example usage:
###  test1 <- fast_cor(
###   data,
###   use = "pairwise.complete.obs",
###   method = "spearman"
###     )
### "data" is class "matrix" with samples as rows and observations as columns


fast_cor <- function(data, use = "pairwise.complete.obs", method = "spearman") 
{
  {
    n = NULL
    n <- psych::pairwiseCount(data) - 2 ## Create matrix of degrees freedom
  }
  if (method == "spearman") {
    data <- base::apply(X = data, MARGIN = 2, data.table::frankv)
  }
  r = NULL
  p1 = NULL
  p2 = NULL
  r <- coop::pcor(x = data, use = use)
  {
    
    t <- sqrt(n)*r/sqrt(1 - r^2) ## T-test
    p1 <- stats::pt(t, n)  ## One side p-value
    p2 <- stats::pt(t, n, lower.tail = FALSE) ## Unfortunately run twice...
  }
  flt.Corr.Matrix <- function(cormat, pmat1 = NULL, pmat2 = NULL,
                              df = NULL) {
    ut <- base::upper.tri(cormat)
    flt_data <- data.frame(Var1 = base::rownames(cormat)[base::row(cormat)[ut]], 
                           Var2 = base::rownames(cormat)[base::col(cormat)[ut]], 
                           cor = cormat[ut])
    if (!is.null(pmat1)) 
      flt_data$p1 <- pmat1[ut]
    if (!is.null(pmat2)) 
      flt_data$p2 <- pmat2[ut]
    if (!is.null(df)) 
      flt_data$df <- df[ut]
    return(flt_data)
  }
  {
    result <- flt.Corr.Matrix(cormat = r,
                              df = n,
                              pmat1 = p1,
                              pmat2 = p2)
    result <-  base::transform(result, p = base::pmin(p1, p2)*2)
    result$FDR <- stats::p.adjust(result$p, method = "BH")
    result <- base::subset(result, select = -c(p1,p2))
    result <- list(result, r)
  }
  return(result)
}

###Cluster comparison function
Clust_compare <-function(data, 
                         clusterNumbers= c(15, 30),
                         nameAlgorithm = c('mclust'),
                         models = c('VEI','VII', 'EEI',
                                    'EII','EVI'),
                         shrinkage_values = c(0.01,0)) {
  
  dfList <- list()
  print('Initializing...')
  data <- t(scale(t(data)))
  for (algorithm in nameAlgorithm) {
    for (clusters in clusterNumbers) {
      ##MCLUST (REMEMBER SPECIFIC TYPES)
      if (algorithm == 'mclust'){
        print('Running mclust')
        for(model in models) {
          set.seed(420)
          mclust::mclust.options(subset = 4500)
          default_shrinkage <- shrinkage_values[1]
          tryCatch({
            result <- mclust::Mclust(data, 
                                     modelNames = model, 
                                     G = clusters,
                                     prior = priorControl(shrinkage = default_shrinkage))
            
          }, error = function(e) {
            cat("Error with shrinkage=", default_shrinkage,"\n")
            result <- NULL
          })
          if(is.null(result)) {
            for(shrink_val in shrinkage_values[-1]){
              tryCatch({
                result <- mclust::Mclust(data, 
                                         modelNames = model, 
                                         G = clusters,
                                         prior = priorControl(shrinkage =shrink_val))
              }, error = function(e) {
                cat("Error with shrinkage=", shrink_val,"\n")
              })
            }
          }
          
          df <- data.frame(
            Gene = names(result$classification),
            Cluster = result$classification,
            uncertainty = result$uncertainty)
          df_name <- paste( model, "_clusters_", clusters, sep = "")
          dfList[[df_name]] <- df
        }
      }
    }
  }
  dfList <- dfList[sapply(dfList, nrow)>0]
  return(dfList)
}

all_cor <- fast_cor(t(protein_imputed), method = "pearson")

X <- t(scale(t(protein_imputed)))     # same scaling you used
set.seed(444)
mclust.options(subset = 4000)

models <- c("VEI","VII","EEI","EII","EVI")
G_grid  <- 1:40                        # change upper bound as needed

fit <- Mclust(
  data = X,
  G = G_grid,
  modelNames = models,
  prior = mclust::priorControl(shrinkage = 0.01)
)

fit$G          # <-- selected number of clusters
fit$modelName  # <-- selected covariance model

# Visualize the BIC surface across G and model
plot(fit, what = "BIC")

# Cluster membership + uncertainty
table(fit$classification)
summary(fit$uncertainty)

##27 clusters seems apt
clust_summary <- Clust_compare(protein_imputed,
                               clusterNumbers= c(27),
                               nameAlgorithm = c("mclust"
                               ),
                               models = "EVI",
                               shrinkage_values = c(0.01,0))

# Inputs:
# all_cor <- fast_cor(t(protein_imputed), method = "pearson")
# clust_summary <- Clust_compare(protein_imputed, clusterNumbers=c(30), models="EVI", shrinkage_values=c(0.01,0))

All_cors <- all_cor[[1]]
cor_mat  <- all_cor[[2]]

###Pull clustering table
clust_df1 <- if (!is.null(names(clust_summary)) && "EVI_clusters_27" %in% names(clust_summary)) {
  clust_summary[["EVI_clusters_27"]]
} else {
  clust_summary[[1]]
}

clust_df1 <- clust_df1 %>%
  transmute(
    Gene = as.character(Gene),
    Cluster = as.character(Cluster),
    uncertainty = as.numeric(uncertainty)
  )

###Keep only HIGH-certainty genes (low uncertainty)
u_cut <- 0.1   # adjust if you want stricter (e.g., 0.15) or looser
clust_df <- clust_df1 %>% filter(uncertainty <= u_cut)

###Keep only genes present in correlation matrix
genes_in_mat <- intersect(clust_df$Gene, rownames(cor_mat))
clust_df <- clust_df %>% filter(Gene %in% genes_in_mat)
cor_mat  <- cor_mat[genes_in_mat, genes_in_mat, drop = FALSE]

###ONLY use FDR<0.05 for ordering clusters (everything else same)
All_cors_fdr <- All_cors %>% dplyr::filter(FDR < 0.05)

pairs_med <- All_cors_fdr %>%
  filter(Var1 %in% genes_in_mat, Var2 %in% genes_in_mat) %>%
  inner_join(clust_df %>% select(Gene, Cluster), by = c("Var1" = "Gene")) %>%
  rename(Cluster1 = Cluster) %>%
  inner_join(clust_df %>% select(Gene, Cluster), by = c("Var2" = "Gene")) %>%
  rename(Cluster2 = Cluster) %>%
  filter(Cluster1 != Cluster2) %>%
  group_by(Cluster1, Cluster2) %>%
  summarise(MedCor = mean(cor, na.rm = TRUE), .groups = "drop")

clusters <- sort(unique(clust_df$Cluster))
cc_mat <- matrix(0, nrow = length(clusters), ncol = length(clusters),
                 dimnames = list(clusters, clusters))
diag(cc_mat) <- 1

if (nrow(pairs_med) > 0) {
  cc_mat[cbind(pairs_med$Cluster1, pairs_med$Cluster2)] <- pairs_med$MedCor
  cc_mat[cbind(pairs_med$Cluster2, pairs_med$Cluster1)] <- pairs_med$MedCor
  
  cc_mat <- pmax(pmin(cc_mat, 1), -1)
  hc <- hclust(as.dist(1 - cc_mat), method = "complete")
  cluster_levels <- hc$labels[hc$order]
} else {
  # fallback if FDR<0.05 leaves no usable inter-cluster edges
  cluster_levels <- clusters
}

###5) Order genes by cluster -> uncertainty -> gene
clust_order <- clust_df %>%
  mutate(Cluster = factor(Cluster, levels = cluster_levels)) %>%
  arrange(Cluster, uncertainty, Gene)

gene_order  <- clust_order$Gene
cor_mat_ord <- cor_mat[gene_order, gene_order, drop = FALSE]

###6) Annotations: cluster only (since we filtered to high-certainty already)
set.seed(290)
cluster_cols <- structure(circlize::rand_color(length(cluster_levels)), names = cluster_levels)

top_ha <- HeatmapAnnotation(
  Module = clust_order$Cluster,
  col = list(Module = cluster_cols),
  show_annotation_name = FALSE,
  border = FALSE
)

left_ha <- rowAnnotation(
  Module = clust_order$Cluster,
  col = list(Module = cluster_cols),
  show_annotation_name = FALSE,
  border = FALSE
)

###Heatmap (full correlations among kept genes; NOT masked by FDR)
ht <- Heatmap(
  cor_mat_ord,
  name = "Pearson r",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = top_ha,
  left_annotation = left_ha,
  use_raster = TRUE
)

draw(ht, heatmap_legend_side = "right")

all_cluster_proteins <- All_cors_fdr %>%
  filter(Var1 %in% genes_in_mat, Var2 %in% genes_in_mat) %>%
  inner_join(clust_df %>% select(Gene, Cluster), 
             by = c("Var1" = "Gene")) %>%
  rename(Cluster1 = Cluster) %>%
  inner_join(clust_df %>% select(Gene, Cluster), 
             by = c("Var2" = "Gene")) %>%
  rename(Cluster2 = Cluster)

###Use the same correlation matrix
cor_mat <- if (is.list(all_cor)) all_cor[[2]] else all_cor

###define your hit sets
hits8  <- c("tr|A0A066X3R0|A0A066X3R0_COLSU", "tr|A0A066Y192|A0A066Y192_COLSU")
hits18 <- c("tr|A0A066XJE7|A0A066XJE7_COLSU")

###merged protein universe: union of both "keep" sets
res8 <- all_cluster_proteins %>% dplyr::filter(Cluster1 == 8, Cluster2 == 8) %>% mutate(Cluster1 = sub(".*\\|(.*)\\|.*", "\\1", Cluster1)) %>% mutate(Cluster1 = sub(".*\\|(.*)\\|.*", "\\1", Cluster1))
neighbors <- res8 %>% dplyr::filter(Var1 %in% hits8 | Var2 %in% hits8) %>% dplyr::transmute(p = ifelse(Var1 %in% hits8, Var2, Var1)) %>% dplyr::pull(p) %>% unique() 
keep <- sort(unique(c(hits8, neighbors)))

res18 <- all_cluster_proteins %>% dplyr::filter(Cluster1 == 18, Cluster2 == 18)
neighbors_18 <- res18 %>% dplyr::filter(Var1 %in% hits18 | Var2 %in% hits18) %>% dplyr::transmute(p = ifelse(Var1 %in% hits18, Var2, Var1)) %>% dplyr::pull(p) %>% unique() 
keep18 <- sort(unique(c(hits18, neighbors_18)))

keep_all <- sort(unique(c(keep, keep18)))

###merged matrix
#mat_merged <- cor_mat[keep_all, keep_all, drop = FALSE]
#ann_merged <- data.frame(
#  Cluster = dplyr::case_when(
#    keep_all %in% keep  & keep_all %in% keep18 ~ "Both",
#    keep_all %in% keep  ~ "8",
#    keep_all %in% keep18 ~ "18",
#    TRUE ~ "Other"
#  ),
#  Hit = ifelse(keep_all %in% c(hits8, hits18), "Hit", "Other"),
#  row.names = keep_all)

#ann_colors <- list(
#  Cluster = c(`8` = "skyblue", `18` = "orange", Other = "grey80"),
#  Hit     = c(Hit = "red", Other = "grey95"))

#ord <- order(factor(ann_merged$Cluster, levels = c("8", "18", "Other")), keep_all)
#mat_merged <- mat_merged[ord, ord, drop = FALSE]
#ann_merged <- ann_merged[ord, , drop = FALSE]
###add a gap between 8 and 18 blocks (optional)
#gap_pos <- sum(ann_merged$Cluster == "8")
#ph_merged <- pheatmap::pheatmap(
#  mat_merged,
#  main = "Protein-Protein Covariation",
#  breaks = seq(-1, 1, length.out = 101),
#  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
#  cluster_rows = TRUE,
#  cluster_cols = TRUE,
#  show_rownames = TRUE,
#  show_colnames = TRUE,
#  border_color = NA,
#  annotation_row = ann_merged,
#  annotation_col = ann_merged,
#  annotation_colors = ann_colors,
#  gaps_row = gap_pos,
#  gaps_col = gap_pos, fontsize = 7)

#png("Cluster8_18_hits_neighbors_heatmap2.png", width = 6000, height = 5000, res = 300)
#grid::grid.newpage()
#grid::grid.draw(ph_merged$gtable)
#dev.off()

#res8 <- all_cluster_proteins %>% 
#  dplyr::filter(Cluster1 == 8, Cluster2 == 8) %>% 
#  mutate(Cluster1 = sub(".*\\|(.*)\\|.*", "\\1", Cluster1)) %>%
#  mutate(Cluster1 = sub(".*\\|(.*)\\|.*", "\\1", Cluster1))

#prots <- sort(unique(c(res8$Var1, res8$Var2)))

#cor_mat <- if (is.list(all_cor)) all_cor[[2]] else all_cor
#mat <- cor_mat[prots, prots]
#ph_all <- pheatmap::pheatmap(
#  mat,
#  cluster_rows = TRUE,
#  cluster_cols = TRUE,
#  show_rownames = FALSE,
#  show_colnames = FALSE, 
#  breaks = seq(-1, 1, length.out = 101),
#  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
#  main = "Cluster 8 correlation heatmap (subset of all_cor)")

#hits <- c("tr|A0A066X3R0|A0A066X3R0_COLSU","tr|A0A066Y192|A0A066Y192_COLSU")
#neighbors <- res8 %>%
#  dplyr::filter(Var1 %in% hits | Var2 %in% hits) %>%
#  dplyr::transmute(p = ifelse(Var1 %in% hits, Var2, Var1)) %>%
#  dplyr::pull(p) %>%
#  unique()
#keep <- sort(unique(c(hits, neighbors)))
#mat_small <- mat[keep, keep]

#ann <- data.frame(
#  Hit = ifelse(rownames(mat_small) %in% hits, "Hit", "Other"),
#  row.names = rownames(mat_small))

#ann_colors <- list(
#  Hit = c(Hit = "red", Other = "grey80"))

#ph <- pheatmap::pheatmap(mat_small, 
#                         main = "Hits + neighbors (Cluster 8)",
#                         breaks = seq(-1, 1, length.out = 101), 
#                         legend_breaks = c(-1, -0.5, 0, 0.5, 1),
#                         cluster_rows = TRUE,
#                         cluster_cols = TRUE,
#                         fontsize_row = 5,
#                         fontsize_col = 5,
#                         border_color = NA,
#                         annotation_row = ann,
#                         annotation_colors = ann_colors
#)

#png("Cluster8_hits_neighbors_heatmap.png", width = 3000, height = 3000, res = 300)
#grid::grid.newpage()
#grid::grid.draw(ph$gtable)
#dev.off()

#res18 <- all_cluster_proteins %>%
#  dplyr::filter(Cluster1 == 18, Cluster2 == 18)

#prots_18 <- sort(unique(c(res18$Var1, res18$Var2)))

#cor_mat_18 <- if (is.list(all_cor)) all_cor[[2]] else all_cor
#mat18 <- cor_mat_18[prots_18, prots_18]

#ph_all_18 <- pheatmap::pheatmap(
#  mat18,
#  cluster_rows = TRUE,
#  cluster_cols = TRUE,
#  show_rownames = FALSE,
#  show_colnames = FALSE,
#  breaks = seq(-1, 1, length.out = 101),
#  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
#  main = "Cluster 18 correlation heatmap (subset of all_cor)")

#hits_18 <- c("tr|A0A066XJE7|A0A066XJE7_COLSU")

#neighbors_18 <- res18 %>%
#  dplyr::filter(Var1 %in% hits_18 | Var2 %in% hits_18) %>%
#  dplyr::transmute(p = ifelse(Var1 %in% hits_18, Var2, Var1)) %>%
#  dplyr::pull(p) %>%
#  unique()

#keep18 <- sort(unique(c(hits_18, neighbors_18)))
#mat_small18 <- mat18[keep18, keep18, drop = FALSE]

#ann18 <- data.frame(
#  Hit = ifelse(rownames(mat_small18) %in% hits_18, "Hit", "Other"),
#  row.names = rownames(mat_small18))

#ann_colors18 <- list(Hit = c(Hit = "red", Other = "grey80"))

#ph18 <- pheatmap::pheatmap(
#  mat_small18,
#  main = "Cluster 18: hits + neighbors",
#  breaks = seq(-1, 1, length.out = 101),
#  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
#  cluster_rows = TRUE,
#  cluster_cols = TRUE,
#  fontsize_row = 5,
#  fontsize_col = 5,
#  border_color = NA,
#  annotation_row = ann18,
#  annotation_colors = ann_colors18)

#png("Cluster18_hits_neighbors_heatmap.png", width = 2000, height = 1000, res = 300)
#grid::grid.newpage()
#grid::grid.draw(ph18$gtable)
#dev.off()

###Get all unique proteins involved in Cluster 8
prots8 <- sort(unique(c(res8$Var1, res8$Var2)))
###Subset your intensity matrix (rows = ProteinID, cols = samples)
mat8 <- protein_imputed[prots8, , drop = FALSE]
###Convert to long format
df8_long <- mat8 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ProteinID") %>%
  tidyr::pivot_longer(-ProteinID, names_to = "Sample", values_to = "log2Intensity") %>%
  dplyr::filter(!is.na(log2Intensity))

###Define condition for each sample
df8_long <- df8_long %>%
  dplyr::mutate(
    Condition = dplyr::case_when(
      str_detect(Sample, regex("one_tenth", ignore_case = TRUE)) ~ "OneTenth",
      str_detect(Sample, regex("half", ignore_case = TRUE)) ~ "Half",
      str_detect(Sample, regex("full", ignore_case = TRUE)) ~ "Full",
      TRUE ~ "Other"
    )
  )

med_df <- df8_long %>%
  group_by(Condition) %>%
  summarise(med = median(log2Intensity, na.rm = TRUE), .groups = "drop")

ggplot(df8_long, aes(x = log2Intensity, fill = Condition)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.35,
                 aes(y = after_stat(density))) +
  geom_vline(data = med_df, aes(xintercept = med, color = Condition),
             linewidth = 0.8, linetype = "dashed", show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  scale_x_continuous(limits = c(20, 35), breaks = seq(20, 35, 5)) +
  labs(
    title = "Cluster 8 proteins: log2 intensity distributions with medians",
    x = "log2 intensity",
    y = "Density",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 16)

###Get all unique proteins involved in Cluster 18
prots18 <- sort(unique(c(res18$Var1, res18$Var2)))
###Subset your intensity matrix (rows = ProteinID, cols = samples)
mat18 <- protein_imputed[prots18, , drop = FALSE]
###Convert to long format
df18_long <- mat18 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ProteinID") %>%
  tidyr::pivot_longer(-ProteinID, names_to = "Sample", values_to = "log2Intensity") %>%
  dplyr::filter(!is.na(log2Intensity))

###Define condition for each sample
df18_long <- df18_long %>%
  dplyr::mutate(
    Condition = dplyr::case_when(
      str_detect(Sample, regex("one_tenth", ignore_case = TRUE)) ~ "OneTenth",
      str_detect(Sample, regex("half", ignore_case = TRUE)) ~ "Half",
      str_detect(Sample, regex("full", ignore_case = TRUE)) ~ "Full",
      TRUE ~ "Other"
    )
  )

med_df_18 <- df18_long %>%
  group_by(Condition) %>%
  summarise(med = median(log2Intensity, na.rm = TRUE), .groups = "drop")

ggplot(df18_long, aes(x = log2Intensity, fill = Condition)) +
  geom_histogram(bins = 25, position = "identity", alpha = 0.35,
                 aes(y = after_stat(density))) +
  #geom_vline(data = med_df_18, aes(xintercept = med, color = Condition), linewidth = 0.18, linetype = "dashed", show.legend = FALSE) +
  #scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  #scale_x_continuous(limits = c(20, 35), breaks = seq(20, 35, 5)) +
  labs(
    title = "Cluster 18 proteins: log2 intensity distributions with medians",
    x = "log2 intensity",
    y = "Density",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 16)

###Get all unique proteins involved in Cluster 13
res13 <- all_cluster_proteins %>% 
  dplyr::filter(Cluster1 == 13, Cluster2 == 13)

prots13 <- sort(unique(c(res13$Var1, res13$Var2)))

mat13 <- protein_imputed[prots13, , drop = FALSE]

df13_long <- mat13 %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ProteinID") %>%
  tidyr::pivot_longer(-ProteinID, names_to = "Sample", values_to = "log2Intensity") %>%
  dplyr::filter(!is.na(log2Intensity))

df13_long <- df13_long %>%
  dplyr::mutate(
    Condition = dplyr::case_when(
      str_detect(Sample, regex("one_tenth", ignore_case = TRUE)) ~ "OneTenth",
      str_detect(Sample, regex("half", ignore_case = TRUE)) ~ "Half",
      str_detect(Sample, regex("full", ignore_case = TRUE)) ~ "Full",
      TRUE ~ "Other"
    )
  )
med_df_13 <- df13_long %>%
  group_by(Condition) %>%
  summarise(med = median(log2Intensity, na.rm = TRUE), .groups = "drop")

ggplot(df13_long, aes(x = log2Intensity, fill = Condition)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.35,
                 aes(y = after_stat(density))) +
  #geom_vline(data = med_df_13, aes(xintercept = med, color = Condition),
  #           linewidth = 0.8, linetype = "dashed", show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  scale_x_continuous(limits = c(20, 35), breaks = seq(20, 35, 5)) +
  labs(
    title = "Cluster 13 proteins: log2 intensity distributions with medians",
    x = "log2 intensity",
    y = "Density",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 16)

cluster8_prots <- unique(df8_long$ProteinID)

# helper to standardize each topTable
prep_tt <- function(tt, comparison_name) {
  tt %>%
    tibble::rownames_to_column("ProteinID") %>%
    transmute(
      ProteinID,
      logFC = logFC,
      adj.P.Val = adj.P.Val,
      Comparison = comparison_name
    )
}

###combine + subset to Cluster 8 proteins
volc_df <- bind_rows(
  prep_tt(half_PDA_vs_PDA,          "Half vs Full"),
  prep_tt(OneTenth_PDA_vs_PDA,      "OneTenth vs Full"),
  prep_tt(OneTenth_PDA_vs_half_PDA, "OneTenth vs Half")) %>%
  filter(ProteinID %in% cluster8_prots) %>%
  mutate(adj.P.Val = pmax(adj.P.Val, 1e-300),
         negLog10FDR = -log10(adj.P.Val),
         Significant = (adj.P.Val < 0.05 & abs(logFC) >= 1))
###Volcano (all comparisons overlaid; color-coded by comparison)
ggplot(volc_df, aes(x = logFC, y = negLog10FDR, color = Comparison)) +
  geom_point(aes(alpha = Significant), size = 3) +
  scale_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.3), guide = "none") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.75) +
  labs(
    title = "Volcano plot (Cluster 8 proteins) across comparisons",
    x = "log2 fold-change",
    y = "-log10(FDR)",
    color = "Comparison"
  ) +
  theme_minimal(base_size = 16)

cluster18_prots <- unique(df18_long$ProteinID)
###combine + subset to Cluster 8 proteins
volc_df_18 <- bind_rows(
  prep_tt(half_PDA_vs_PDA,          "Half vs Full"),
  prep_tt(OneTenth_PDA_vs_PDA,      "OneTenth vs Full"),
  prep_tt(OneTenth_PDA_vs_half_PDA, "OneTenth vs Half")
) %>%
  filter(ProteinID %in% cluster18_prots) %>%
  mutate(
    adj.P.Val = pmax(adj.P.Val, 1e-300),
    negLog10FDR = -log10(adj.P.Val),
    Significant = (adj.P.Val < 0.05 & abs(logFC) >= 1)
  )
###Volcano (all comparisons overlaid; color-coded by comparison)
ggplot(volc_df_18, aes(x = logFC, y = negLog10FDR, color = Comparison)) +
  geom_point(aes(alpha = Significant), size = 3) +
  scale_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.3), guide = "none") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.75) +
  labs(
    title = "Volcano plot (Cluster 18 proteins) across comparisons",
    x = "log2 fold-change",
    y = "-log10(FDR)",
    color = "Comparison"
  ) +
  theme_minimal(base_size = 16)

cluster13_prots <- unique(df13_long$ProteinID)
###combine + subset to Cluster 8 proteins
volc_df_13 <- bind_rows(
  prep_tt(half_PDA_vs_PDA,          "Half vs Full"),
  prep_tt(OneTenth_PDA_vs_PDA,      "OneTenth vs Full"),
  prep_tt(OneTenth_PDA_vs_half_PDA, "OneTenth vs Half")
) %>%
  filter(ProteinID %in% cluster13_prots) %>%
  mutate(
    adj.P.Val = pmax(adj.P.Val, 1e-300),
    negLog10FDR = -log10(adj.P.Val),
    Significant = (adj.P.Val < 0.05 & abs(logFC) >= 1)
  )
###Volcano (all comparisons overlaid; color-coded by comparison)
ggplot(volc_df_13, aes(x = logFC, y = negLog10FDR, color = Comparison)) +
  geom_point(aes(alpha = Significant), size = 3) +
  scale_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.3), guide = "none") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.75) +
  labs(
    title = "Volcano plot (Cluster 18 proteins) across comparisons",
    x = "log2 fold-change",
    y = "-log10(FDR)",
    color = "Comparison"
  ) +
  theme_minimal(base_size = 16)

#volc_df_test <- volc_df %>%
#  dplyr::mutate(COLSU = ProteinID) %>%
#  filter(Significant == TRUE) %>%
#  filter(!is.na(COLSU)) %>%
#  mutate(ranking_metric = -log10(adj.P.Val) * sign(logFC)) %>% 
#  group_by(COLSU) %>%
#  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>%
#  inner_join(annota_final, by = "COLSU") %>%
#  dplyr::select(COLHI, ranking_metric) %>%
#  filter(!is.na(COLHI)) %>%
#  mutate(COLHI = sub(".*\\|(.*)\\|.*", "\\1", COLHI)) %>%
#  group_by(COLHI) %>% 
#  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>%  # ensures uniqueness
#  arrange(desc(ranking_metric)) %>%
#  deframe()

#anyDuplicated(names(volc_df_test))

#tail(volc_df_test)

#chig_gsea_volc_df_test <- clusterProfiler::gseKEGG(
#  geneList = volc_df_test,
#  organism = "chig",        # change this to your species' KEGG code
#  keyType = "uniprot",     # because you have UniProt IDs
#  eps = 0.0,
#  minGSSize = 0,
#  maxGSSize = 1000,
#  pAdjustMethod = "BH",
#  pvalueCutoff = 1,
#  verbose = FALSE)

##06152026
##Calculating Protein ranked abundance in each condition
mod_pep2 <- mod_pep 
sample_columns <- colnames(mod_pep2)[!grepl("PDA", colnames(mod_pep2))]
meta3 <- mod_pep2[, sample_columns] 
mod_pep2[mod_pep2 == 0] <- NA
mod_pep2 <- mod_pep2 %>%
  dplyr::select(-ends_with("Count")) %>%
  dplyr::select(-ends_with("Type")) %>%
  dplyr::select(c(ends_with("maxLFQ"),  colnames(meta3)))
filtered_peptides <- peptide_long %>% 
  dplyr::select(Modified.Sequence) %>% 
  unique()

mod_pep2 <- mod_pep2 %>%
  filter(`Modified.Sequence` %in% filtered_peptides$Modified.Sequence) %>%
  pivot_longer(
    cols = -all_of(colnames(meta3)),
    names_to = "run_id",
    values_to = "peptide_intensity"
  ) %>%
  mutate(
    peptide_intensity = as.numeric(peptide_intensity),
    concentration = "?"
  ) %>%
  filter(!is.na(peptide_intensity), peptide_intensity > 0)

ibaqfd <- mod_pep2 %>%
  mutate(peptide_sequence = `Peptide.Sequence`) %>%
  mutate(protein_id = Protein) %>%
  mutate(precursor_charge = Charges) %>%
  mutate(peptide_id = `Modified.Sequence`) %>%
  dplyr::select(protein_id, peptide_sequence, run_id, peptide_id, peptide_intensity, concentration, precursor_charge)

iBAQ <- ProteinInference(ibaqfd, 
                         peptide_method = "iBAQ",
                         consensus_peptides = FALSE,
                         consensus_proteins = FALSE,
                         fasta = "../../references/Colletotrichum_sublineola/Cs_UniProt.fasta")
rIBAQ <- iBAQ %>%
  filter(!grepl("contam_sp", protein_id)) %>%
  mutate(
    condition = case_when(
      grepl("full", run_id) ~ "PDA",
      grepl("half", run_id) ~ "half-PDA",
      grepl("one", run_id)  ~ "OneTenth-PDA",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(condition)) %>%
  group_by(protein_id, condition) %>%
  summarize(
    mean1 = mean(response, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(condition) %>%
  mutate(
    Rank = rank(-mean1, ties.method = "min")
  ) %>%
  ungroup() %>%
  filter(protein_id %in% unique(protein_long_filtered$Protein))