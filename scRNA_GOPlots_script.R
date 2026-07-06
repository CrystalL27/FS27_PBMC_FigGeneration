library(scales)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(EnhancedVolcano)
library(gplots)
library(writexl)
library(EnhancedVolcano)
library(ggrepel)
library(clusterProfiler)
library(scCustomize)


All_GO <- X07_GO_enrichment_DESeq_res0_3_cpm_shrink
#plotting GOTerms from scClusters, which are still numbered 
#change cluster to a Factor
##if have to make a factor to reorder clusters

#use to change cluster order and names
All_GO$cluster <- factor(All_GO$cluster, levels = as.character(0:18))

cluster_order <- c("0", "4", "5", "9", "12", "16", "3", "10", "17", "1", 
                   "7", "8", "6", "14", "2", "13", "11", "18", "15")

All_GO$cluster <- factor(All_GO$cluster, levels = cluster_order)
All_GO$FoldEnrichment <- as.numeric(All_GO$FoldEnrichment)
All_GO$p.adjust <- as.numeric(All_GO$p.adjust)
All_GO$Count <- as.numeric(All_GO$Count)


cluster_labels <- c('Mono_1', 'Mono_2', 'Mono_3', 'Mono_4',
                    'pDC', 'cDC', 'Bcell_1', 'Bcell_2', 'ASC', 
                    'abT_CD4_1', 'abT_CD4_2', 'abT_NK_1', 'abT_NK_2',
                    'abT_prolif', 'gdT_CD2n', 'gdT_CD2p', 
                    'Mixed_Tcells_1', 'Mixed_Tcells_2', 'Mixed_T_Mono')
cluster_labels_map <- setNames(cluster_labels, cluster_order)

MonoClusters_GO <- All_GO %>%
  filter(cluster %in% c("0", "4", "5", "9"))
MonoClusters_GO <- MonoClusters_GO %>%
  group_by(cluster) %>%
  slice_min(p.adjust, n = 16)

ggplot(data=MonoClusters_GO, aes(x=factor(direction, levels = c("up", "down")), 
                                 y=fct_reorder(Term_GO, direction), 
                                 fill=p.adjust, 
                                 size=FoldEnrichment)) +
  scale_fill_gradient(low = "tomato3",
                      high = "yellow2") + 
  geom_point(shape=21) +
  theme_bw() +
  facet_wrap(~cluster, ncol=4, labeller = labeller(cluster = as_labeller(cluster_labels_map))) +
  labs(y="GO Term", 
       x="DEG Direction", 
       fill="AdjPvalue", 
       size="Fold Enrichment") + RotatedAxis()
