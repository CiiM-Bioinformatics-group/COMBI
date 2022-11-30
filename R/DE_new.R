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

runDEseq2 <- function(counts, meta, stimul, control, time, outfile) {

    meta_sub <- meta %>%
      filter(stimulation == stimul | stimulation == control) %>%
      filter(timepoint == time)

    counts_sub <- counts[, rownames(meta_sub)]

    # filter those low genes
    cutoff <- 20
    table(rowSums(counts_sub) > cutoff)
    counts_sub <- counts_sub[rowSums(counts_sub) > cutoff, ]

    stopifnot(all(colnames(counts_sub) == rownames(meta_sub)))

    print(table(meta_sub$stimulation, meta_sub$timepoint))


    dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                                  colData = meta_sub,
                                  design = ~ age_class + sex + stimulation)

    dds <- DESeq(dds)

    res <- results(dds, contrast = c('stimulation', stimul, control)) %>% 
      as.data.frame() %>%
      # arrange(padj) %>%
      mutate(gene = rownames(.))

    # merge with the amount of reads per group. Median
    stimul.samples <- meta_sub %>% filter(stimulation == stimul) %>% rownames()
    median.stimul <- counts_sub[, stimul.samples] %>% as.matrix() %>% rowMedians()
    sum.stimuli <- counts_sub[, stimul.samples] %>% as.matrix() %>% rowSums()

    control.samples <- meta_sub %>% filter(stimulation == control) %>% rownames()
    median.controls <- counts_sub[, control.samples] %>% as.matrix() %>% rowMedians()
    sum.controls <- counts_sub[, control.samples] %>% as.matrix() %>% rowSums()

    res <- cbind(res,
                 data.frame(
                   'median_reads_stimulation' = median.stimul,
                   'median_reads_control' = median.controls,
                   'sum_reads_stimulation' = sum.stimuli,
                   'sum_reads_control' = sum.controls
                 )) %>%
      arrange(padj)

    write.xlsx(res, outfile, overwrite = T)
}



comps <- list(
  c(stim = 'Influenza', control = 'RPMI', time = 't0', outfile = 'output/DE/t0_influenza_RPMI.xlsx'),
  c(stim = 'SARSCOV2', control = 'RPMI', time = 't0', outfile = 'output/DE/t0_SARSCOV2_RPMI.xlsx'),
  c(stim = 'R848', control = 'RPMI', time = 't0', outfile = 'output/DE/t0_R848_RPMI.xlsx'),

  c(stim = 'Influenza', control = 'RPMI', time = 't1', outfile = 'output/DE/t1_influenza_RPMI.xlsx'),
  c(stim = 'SARSCOV2', control = 'RPMI', time = 't1', outfile = 'output/DE/t1_SARSCOV2_RPMI.xlsx'),
  c(stim = 'R848', control = 'RPMI', time = 't1', outfile = 'output/DE/t1_R848_RPMI.xlsx'),

  c(stim = 'Influenza', control = 'RPMI', time = 't2', outfile = 'output/DE/t2_influenza_RPMI.xlsx'),
  c(stim = 'SARSCOV2', control = 'RPMI', time = 't2', outfile = 'output/DE/t2_SARSCOV2_RPMI.xlsx'),
  c(stim = 'R848', control = 'RPMI', time = 't2', outfile = 'output/DE/t2_R848_RPMI.xlsx'),

  c(stim = 'Influenza', control = 'RPMI', time = 't3', outfile = 'output/DE/t3_influenza_RPMI.xlsx'),
  c(stim = 'SARSCOV2', control = 'RPMI', time = 't3', outfile = 'output/DE/t3_SARSCOV2_RPMI.xlsx'),
  c(stim = 'R848', control = 'RPMI', time = 't3', outfile = 'output/DE/t3_R848_RPMI.xlsx')
)

for (comp in comps) {
  runDEseq2(counts = counts,
            meta = meta,
            stim = comp[['stim']],
            control = comp[['control']],
            time = comp[['time']],
            outfile = comp[['outfile']])
}
