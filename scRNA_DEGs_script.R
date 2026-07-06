library(scales)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(gplots)
library(writexl)
library(EnhancedVolcano)
library(ggrepel)
library(UpSetR)

#allgenes,not just significant ones, to create volcano plots
DE_ALL <- X10d_DESeq2_pseudobulk_shrinkage_0_3res_Salmonella_vs_Control

#Ensure the necessary columns (log2FC and p-value) are numeric (double)
DE_ALL$log2FoldChange <- as.numeric(DE_ALL$log2FoldChange)
DE_ALL$pvalue <- as.numeric(DE_ALL$pvalue)
DE_ALL$padj <- as.numeric(DE_ALL$padj)
DE_ALL$cluster <- as.character(DE_ALL$cluster)
DE_ALL$cluster <- factor(DE_ALL$cluster, levels = cluster_order)
rownames(DE_ALL) <- DE_ALL$gene

##change cluster order,first set cluster order identities
cluster_order <- c("0", "4", "5", "9", "12", "16", "3", "10", "17", "1", 
                   "7", "8", "6", "14", "2", "13", "11", "18", "15")
cluster_labels <- c('Mono_1', 'Mono_2', 'Mono_3', 'Mono_4',
                    'pDC', 'cDC', 'Bcell_1', 'Bcell_2', 'ASC', 
                    'abT_CD4_1', 'abT_CD4_2', 'abT_NK_1', 'abT_NK_2',
                    'abT_prolif', 'gdT_CD2n', 'gdT_CD2p', 
                    'Mixed_Tcells_1', 'Mixed_Tcells_2', 'Mixed_T_Mono')
DE_ALL$cluster <- as.character(DE_ALL$cluster)
DE_ALL$cluster <- factor(DE_ALL$cluster, levels = cluster_order)
cluster_labels_map <- setNames(cluster_labels, cluster_order)

##repeat volcanoplots from Sharu ##SupplementalFigure5
ggplot(DE_ALL, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50"),
                     name = "Significant DE") +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       title = "Differentially Expressed Genes by Cluster (0.3res)",
       subtitle = "Red = p_adj < 0.05 & |log2FC| > 1") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 35), oob = squish) +
  scale_x_continuous(limits = c(-6, 6), oob = squish) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(
    ~ cluster,
    ncol = 6,
    labeller = labeller(cluster = as_labeller(cluster_labels_map))
  )

##isolate specific cell type data to graph independently
Mono1 <- DE_ALL %>%
  filter(cluster == "0")
Mono2 <- DE_ALL %>%
  filter(cluster == "4")
Mono3 <- DE_ALL %>%
  filter(cluster == "5")
Mono4 <- DE_ALL %>%
  filter(cluster == "9")
##change row names to gene names
rownames(Mono1) <- Mono1$gene
rownames(Mono2) <- Mono2$gene
rownames(Mono3) <- Mono3$gene
rownames(Mono4) <- Mono4$gene
cDC <- DE_ALL %>%
  filter(cluster == "16")
rownames(cDC) <- cDC$gene
abT_NK_2 <- DE_ALL %>%
  filter(cluster == "6")
rownames(abT_NK_2) <- abT_NK_2$gene

#plots all clusters from DE_ALL
EnhancedVolcano(toptable = DE_ALL,
                title = "All Clusters - Control vs Salmonella",
                lab = rownames(DE_ALL),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoffCol = 'padj',
                pCutoff = .05,
                FCcutoff = 1.0,
                cutoffLineWidth = 0,
                col = c('lightgrey', 'grey40', 'grey40', 'red3'),
                pointSize = 1.5,
                labSize = 6.0,
                colAlpha = 1,
                legendLabels=c('Not sig.','Log2FC','padj',
                               'padj & Log2FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0,
                ylim = c(0, 40),
                xlim = c(-6, 6)) + 
                facet_wrap(
                ~ cluster,
                ncol = 6, labeller = labeller(cluster = as_labeller(cluster_labels_map)))

#create a list of genes for Mono1 cluster
Mono1genes <- c('IDO1', 'S100A12', 'TCEA3', 'GBP2', 'LCN2', 'EDH1',
                'C15orf48', 'S100A9', 'CD163', 'C1QC', 'SMAD3', 'SMAD6',
                'CD1.1', 'ICOSLC', 'SLO6A6', 'CXCR4', 'S100A8', 'IDO1', 'SELL',
                'LRG1', 'IL6', 'SLC6A6','JDP1', 'CD274', 'C5AR2', 'WNT10B',
                'TGFB1', 'CXCL10', 'CCL21', 'SOX5', 'MS4A7')
#create volcano plot of just Mono1 cluster using the Mono1genes list
EnhancedVolcano(toptable = Mono1,
                title = "Mono_1 - Control vs Salmonella",
                lab = rownames(Mono1),
                selectLab = Mono1genes,
                drawConnectors = TRUE,
                max.overlaps = Inf,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoffCol = 'padj',
                ylab = bquote(~-log[10]~"p-adj"),
                xlab = bquote(~log[2]~"FoldChange"),
                pCutoff = .05,
                FCcutoff = 1.0,
                cutoffLineWidth = 0,
                col = c('lightgrey', 'grey40', 'grey40', 'red3'),
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 4/5,
                caption = '555 DEGs',
                legendLabels=c('Not sig.','Log2FC','p-adj',
                               'p-adj & Log2FC'),
                legendPosition = 'bottom',
                legendLabSize = 16,
                legendIconSize = 3,
                ylim = c(0, 35),
                xlim = c(-5.8, 5.8))

##gene list
Mono2genes <- c('HPS5', 'MSA4A7', 'S100A12', 'TCEA3', 'GBP2', 'LCN2', 'EDH1', 'MS4A8',
                'C15orf48', 'S100A9', 'CD163', 'CCR1', 'SMAD3', 'SMAD6', 'GB2',
                'ICOSLC', 'SLO6A6', 'CXCR4', 'S100A8', 'IDO1', 'MS4A8', 'LAP3', 'SELL',
                'LRG1', 'ID3', 'SLC03A1','HDAC9', 'BLK', 'CD274', 'IL21R',
                'S1PR1', 'VEGFA', 'UBP1', 'FLT3', 'SIGLEC1', 'RETN', 'EHD1', 'UPP1', 'CCR7', 'EMP1', 'CLEC2B', 'GBP1')

EnhancedVolcano(toptable = Mono2,
                title = "Mono_2 - Control vs Salmonella",
                lab = rownames(Mono2),
                selectLab = Mono2genes,
                drawConnectors = TRUE,
                max.overlaps = Inf,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoffCol = 'padj',
                ylab = bquote(~-log[10]~"p-adj"),
                xlab = bquote(~log[2]~"FoldChange"),
                pCutoff = .05,
                FCcutoff = 1.0,
                cutoffLineWidth = 0,
                col = c('lightgrey', 'grey40', 'grey40', 'red3'),
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 4/5,
                caption = '709 DEGs',
                legendLabels=c('Not sig.','Log2FC','p-adj',
                               'p-adj & Log2FC'),
                legendPosition = 'bottom',
                legendLabSize = 16,
                legendIconSize = 3,
                ylim = c(0, 40),
                xlim = c(-5, 5))

Mono3genes <- c('HPS5', 'MSA4A7', 'S100A12', 'TCEA3', 'GBP2', 'LCN2', 'EDH1', 'MS4A8',
                'S100A9', 'SMAD3', 'SMAD6', 'GB2', 'IL6', 'LTF', 'SERPING1',
                'ICOSLC', 'SLO6A6', 'CXCR4', 'S100A8', 'MS4A8','SELL',
                'ID3', 'SLC03A1','HDAC9', 'BLK', 'CD274', 'ID2', 'KLF9',
                'VEGFA', 'FLT3', 'EHD1', 'UPP1', 'CCR7', 'EMP1', 'CLEC2B', 'GBP1')
EnhancedVolcano(toptable = Mono3,
                title = "Mono_3 - Control vs Salmonella",
                lab = rownames(Mono3),
                selectLab = Mono3genes,
                drawConnectors = TRUE,
                max.overlaps = Inf,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoffCol = 'padj',
                ylab = bquote(~-log[10]~"p-adj"),
                xlab = bquote(~log[2]~"FoldChange"),
                pCutoff = .05,
                FCcutoff = 1.0,
                cutoffLineWidth = 0,
                col = c('lightgrey', 'grey40', 'grey40', 'red3'),
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 4/5,
                caption = '423 DEGs',
                legendLabels=c('Not sig.','Log2FC','p-adj',
                               'p-adj & Log2FC'),
                legendPosition = 'bottom',
                legendLabSize = 16,
                legendIconSize = 3,
                ylim = c(0, 40),
                xlim = c(-4, 6))

Mono4genes <- c('MSA4A7', 'S100A12', 'GBP2', 'MS4A8', 'S100A8', 'S100A9', 'CD163', 
                'CCR1', 'SMAD3', 'SMAD6', 'CXCR4', 'IDO1', 'SELL', 'IL18', 'IDO1',
                'CD274', 'WNT10B', 'GBP2', 'CYP1B1', 'VEGFA', 'PPARG', 'LSAMP', 'IL15',
                'LCN2', 'C4BPA', 'TGFB1', 'KLF9')
EnhancedVolcano(toptable = Mono4,
                title = "Mono_4 - Control vs Salmonella",
                lab = rownames(Mono4),
                selectLab = Mono4genes,
                drawConnectors = TRUE,
                max.overlaps = Inf,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoffCol = 'padj',
                ylab = bquote(~-log[10]~"p-adj"),
                xlab = bquote(~log[2]~"FoldChange"),
                pCutoff = .05,
                FCcutoff = 1.0,
                cutoffLineWidth = 0,
                col = c('lightgrey', 'grey40', 'grey40', 'red3'),
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 4/5,
                caption = '293 DEGs',
                legendLabels=c('Not sig.','Log2FC','p-adj',
                               'p-adj & Log2FC'),
                legendPosition = 'bottom',
                legendLabSize = 16,
                legendIconSize = 3,
                ylim = c(0, 25),
                xlim = c(-6, 6))

cDCgenes <- c('MSA4A7', 'GBP2', 'MS4A8', 'S100A8', 'CD163', 
              'CCR1', 'CCR2', 'IDO1', 'LAP3','NOD1', 'SERPING1', 'SLC11A1', 'LYZ',
              'CD274', 'WNT10B', 'GBP1', 'CLEC7A', 'LSAMP', 'IL15', 'TRIB3', 'RAMP1',
              'LCN2', 'C4BPA', 'TLR4', 'C1QB', 'CCL11', 'C5AR1', 'ACTN1', 'BLK')
EnhancedVolcano(toptable = cDC,
                title = "cDC - Control vs Salmonella",
                lab = rownames(cDC),
                selectLab = cDCgenes,
                drawConnectors = TRUE,
                max.overlaps = Inf,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoffCol = 'padj',
                ylab = bquote(~-log[10]~"p-adj"),
                xlab = bquote(~log[2]~"FoldChange"),
                pCutoff = .05,
                FCcutoff = 1.0,
                cutoffLineWidth = 0,
                col = c('lightgrey', 'grey40', 'grey40', 'red3'),
                pointSize = 1.5,
                labSize = 4,
                colAlpha = 4/5,
                caption = '97 DEGs',
                legendLabels=c('Not sig.','Log2FC','p-adj',
                               'p-adj & Log2FC'),
                legendPosition = 'bottom',
                legendLabSize = 16,
                legendIconSize = 3,
                ylim = c(0, 15),
                xlim = c(-3, 7))
