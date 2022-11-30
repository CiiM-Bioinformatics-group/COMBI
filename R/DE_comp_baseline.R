rm(list = ls())

library(openxlsx)
library(DESeq2)
library(dplyr)
library(magrittr)
library(edgeR)
library(ggplot2)
library(pheatmap)

setwd('/vol/projects/CIIM/COMBI/bulk_seq/')
load('output/data.RData')

meta$time_stim <- paste0(meta$timepoint, '_', meta$stimulation)
meta$age_class <- as.factor(ifelse(meta$age < 40, 'young', 'old'))
meta$sex <- factor(meta$sex)

runDEseq2 <- function(counts, meta, stimul, time, outfile) {

    meta_sub <- meta %>%
      filter(stimulation == stimul | stimulation == 'RPMI') %>%
      filter(timepoint == time | timepoint == 't0')

    meta_sub$stim_time <- paste0(meta_sub$stimulation, '_', meta_sub$timepoint)
    counts_sub <- counts[, rownames(meta_sub)]

    # filter those low genes
    cutoff <- 20
    table(rowSums(counts_sub) > cutoff)
    counts_sub <- counts_sub[rowSums(counts_sub) > cutoff, ]

    stopifnot(all(colnames(counts_sub) == rownames(meta_sub)))

    print(table(meta_sub$stimulation, meta_sub$timepoint))


    dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                                  colData = meta_sub,
                                  design = ~ age_class + sex + stim_time)

    dds <- DESeq(dds)

    res <- results(dds, contrast = c('stim_time', paste0(stimul, '_', time), 'RPMI_t0')) %>%
      as.data.frame() %>%
      # arrange(padj) %>%
      mutate(gene = rownames(.))

    # merge with the amount of reads per group. Median
    samples <- meta_sub %>% filter(stimulation == stimul) %>% rownames()
    median1 <- counts_sub[, samples] %>% as.matrix() %>% rowMedians()
    sum1 <- counts_sub[, samples] %>% as.matrix() %>% rowSums()

    samples <- meta_sub %>% filter(stimulation == 'RPMI') %>% rownames()
    median2 <- counts_sub[, samples] %>% as.matrix() %>% rowMedians()
    sum2 <- counts_sub[, samples] %>% as.matrix() %>% rowSums()

    res <- cbind(res,
                 data.frame(
                   'median1' = median1,
                   'median2' = median2,
                   'sum1' = sum1,
                   'sum2' = sum2
                 )) %>%
      arrange(padj)

    write.xlsx(res, outfile, overwrite = T)
}



comps <- list(
  c(stim = 'Influenza', time = 't0', outfile = 'output/DE/t0_influenza_vs_baseline.xlsx'),
  c(stim = 'SARSCOV2', time = 't0', outfile = 'output/DE/t0_SARSCOV2_vs_baseline.xlsx'),
  c(stim = 'R848', time = 't0', outfile = 'output/DE/t0_R848_vs_baseline.xlsx'),

  c(stim = 'Influenza', time = 't1', outfile = 'output/DE/t1_influenza_vs_baseline.xlsx'),
  c(stim = 'SARSCOV2', time = 't1', outfile = 'output/DE/t1_SARSCOV2_vs_baseline.xlsx'),
  c(stim = 'R848', time = 't1', outfile = 'output/DE/t1_R848_vs_baseline.xlsx'),

  c(stim = 'Influenza', time = 't2', outfile = 'output/DE/t2_influenza_vs_baseline.xlsx'),
  c(stim = 'SARSCOV2', time = 't2', outfile = 'output/DE/t2_SARSCOV2_vs_baseline.xlsx'),
  c(stim = 'R848', time = 't2', outfile = 'output/DE/t2_R848_vs_baseline.xlsx'),

  c(stim = 'Influenza', time = 't3', outfile = 'output/DE/t3_influenza_vs_baseline.xlsx'),
  c(stim = 'SARSCOV2', time = 't3', outfile = 'output/DE/t3_SARSCOV2_vs_baseline.xlsx'),
  c(stim = 'R848', time = 't3', outfile = 'output/DE/t3_R848_vs_baseline.xlsx')
)

for (comp in comps) {
  runDEseq2(counts = counts,
            meta = meta,
            stimul = comp[['stim']],
            time = comp[['time']],
            outfile = comp[['outfile']])
}
