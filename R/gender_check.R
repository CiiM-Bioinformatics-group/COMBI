# Check genders
rm(list = ls())

library(openxlsx)
library(DESeq2)
library(dplyr)
library(magrittr)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(ggsci)

setwd('/vol/projects/CIIM/COMBI/bulk_seq/')
load('output/data.RData')

meta$time_stim <- paste0(meta$timepoint, '_', meta$stimulation)
meta$age_class <- as.factor(ifelse(meta$age < 40, 'young', 'old'))
meta$sex <- factor(meta$sex)

# Filter those low genes
cutoff <- 20
table(rowSums(counts) > cutoff)
counts <- counts[rowSums(counts) > cutoff, ] # 22K genes remaining

genes.to.check <- c('XIST', 'RPS4Y1')

# PCA color by gender
counts <- counts[genes.to.check, ]

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = meta, 
                              design = ~ timepoint)

vsd <- varianceStabilizingTransformation(dds)


df <- cbind(t(counts), meta)

table(df %>% filter(XIST > 100) %>% pull(sex))
table(df %>% filter(RPS4Y1 > 100) %>% pull(sex))

pdf('output/gender_check.pdf', width = 4, height = 4)
ggplot(df) + 
  geom_point(aes(XIST, RPS4Y1, color = sex)) +
  theme_bw() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = .5)) +
  labs(x = 'XIST', y = 'RPS4Y1', title = 'Gender check')
dev.off()