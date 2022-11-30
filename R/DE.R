rm(list = ls())

library(openxlsx)
library(DESeq2)
library(dplyr)
library(magrittr)
library(edgeR)
library(ggplot2)
library(pheatmap)

# Analyses:
#   1. At each timepoint, stimulation vs control
#     - Three volcano plots per timepoint
#     - Barplot with nr DEG per timepoint per stimulation: 3 stimulations over 4 timepoints
#     - Venn Diagrams with overlap at each timepoint for the three stimulations
#
#   2. Check the cytokines from Bushra
#     - IL1B, IL6, TNFa, etc.
#
#       For each cytokine:
#         - x=time, y = avg_log2FC, geom_bar with geom_errorbar for avg_log2FC standard error from deseq2. Fill bars according to significance compared to RPMI control
#         - facet_grid with stimulation.
#


setwd('/vol/projects/CIIM/COMBI/bulk_seq/analysis/')
load('../output/data.RData')

# Paste the timepoint and stimulation together for easier comparison
meta$time_stim <- paste0(meta$timepoint, '_', meta$stimulation)
meta$age_class <- as.factor(ifelse(meta$age < 40, 'young', 'old'))
meta$sex <- factor(meta$sex)

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = meta, 
                              design = ~ age_class + sex + time_stim)

dds <- DESeq(dds)

# Timepoint 0: baseline
res1 <- results(dds, contrast = c('time_stim', 't0_Influenza', 't0_RPMI')) %>% as.data.frame() %>% arrange(padj) %>% mutate(gene = rownames(.))
res2 <- results(dds, contrast = c('time_stim', 't0_R848', 't0_RPMI'))%>% as.data.frame() %>% arrange(padj) %>% mutate(gene = rownames(.))
res3 <- results(dds, contrast = c('time_stim', 't0_SARSCOV2', 't0_RPMI'))%>% as.data.frame() %>% arrange(padj) %>% mutate(gene = rownames(.))

write.xlsx(res1, file = '../output/DE/t0_influenza_RPMI.xlsx')
write.xlsx(res2, file = '../output/DE/t0_R848_RPMI.xlsx')
write.xlsx(res3, file = '../output/DE/t0_SARSCOV2_RPMI.xlsx')



# Timepoint t1
res1 <- results(dds, contrast = c('time_stim', 't1_Influenza', 't1_RPMI'))%>% as.data.frame() %>% arrange(padj) %>% mutate(gene = rownames(.))
res2 <- results(dds, contrast = c('time_stim', 't1_R848', 't1_RPMI'))%>% as.data.frame() %>% arrange(padj)%>% mutate(gene = rownames(.))
res3 <- results(dds, contrast = c('time_stim', 't1_SARSCOV2', 't1_RPMI'))%>% as.data.frame() %>% arrange(padj)%>% mutate(gene = rownames(.))

write.xlsx(res1, file = '../output/DE/t1_influenza_RPMI.xlsx')
write.xlsx(res2, file = '../output/DE/t1_R848_RPMI.xlsx')
write.xlsx(res3, file = '../output/DE/t1_SARSCOV2_RPMI.xlsx')



# Timepoint t2
res1 <- results(dds, contrast = c('time_stim', 't2_Influenza', 't2_RPMI'))%>% as.data.frame() %>% arrange(padj) %>% mutate(gene = rownames(.))
res2 <- results(dds, contrast = c('time_stim', 't2_R848', 't2_RPMI'))%>% as.data.frame() %>% arrange(padj) %>% mutate(gene = rownames(.))
res3 <- results(dds, contrast = c('time_stim', 't2_SARSCOV2', 't2_RPMI'))%>% as.data.frame() %>% arrange(padj) %>% mutate(gene = rownames(.))

write.xlsx(res1, file = '../output/DE/t2_influenza_RPMI.xlsx')
write.xlsx(res2, file = '../output/DE/t2_R848_RPMI.xlsx')
write.xlsx(res3, file = '../output/DE/t2_SARSCOV2_RPMI.xlsx')



# Timepoint t3
res1 <- results(dds, contrast = c('time_stim', 't3_Influenza', 't3_RPMI'))%>% as.data.frame() %>% arrange(padj) %>% mutate(gene = rownames(.))
res2 <- results(dds, contrast = c('time_stim', 't3_R848', 't3_RPMI'))%>% as.data.frame() %>% arrange(padj) %>% mutate(gene = rownames(.))
res3 <- results(dds, contrast = c('time_stim', 't3_SARSCOV2', 't3_RPMI'))%>% as.data.frame() %>% arrange(padj) %>% mutate(gene = rownames(.))

write.xlsx(res1, file = '../output/DE/t3_influenza_RPMI.xlsx')
write.xlsx(res2, file = '../output/DE/t3_R848_RPMI.xlsx')
write.xlsx(res3, file = '../output/DE/t3_SARSCOV2_RPMI.xlsx')












