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
# meta$age_class <- as.factor(ifelse(meta$age < 40, 'young', 'old'))
# meta$sex <- factor(meta$sex)



runDEseq2 <- function(counts, meta, stimul, t, t_base, outfile) {

    meta_sub <- meta %>%
      dplyr::filter(stimulation == stimul) %>%
      dplyr::filter(timepoint == t | timepoint == t_base)

    table(meta_sub$stimulation)
    table(meta_sub$timepoint)
    table(meta_sub$studyID)
    counts_sub <- counts[, rownames(meta_sub)]

    print(table(meta_sub$timepoint, meta_sub$stimulation))

    # filter those low genes
    cutoff <- 20
    table(rowSums(counts_sub) > cutoff)
    counts_sub <- counts_sub[rowSums(counts_sub) > cutoff, ]

    stopifnot(all(colnames(counts_sub) == rownames(meta_sub)))
    dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                                  colData = meta_sub,
                                  design = ~ studyID + timepoint)
                                  # design = ~ sex + age_class + timepoint + studyID)

    dds <- DESeq(dds)

    res <- results(dds, contrast = c('timepoint', t, t_base)) %>%
      as.data.frame() %>%
      # arrange(padj) %>%
      mutate(gene = rownames(.))

    # merge with the amount of reads per group. Median
    samples <- meta_sub %>% filter(timepoint == t) %>% rownames()
    median1 <- counts_sub[, samples] %>% as.matrix() %>% rowMedians()
    sum1 <- counts_sub[, samples] %>% as.matrix() %>% rowSums()

    samples <- meta_sub %>% filter(timepoint == t_base) %>% rownames()
    median2 <- counts_sub[, samples] %>% as.matrix() %>% rowMedians()
    sum2 <- counts_sub[, samples] %>% as.matrix() %>% rowSums()

    res <- cbind(res,
                 data.frame(
                   'median1' = median1,
                   'median2' = median2,
                   'sum1' = sum1,
                   'sum2' = sum2
                 )) %>%
      arrange(padj) %>% 
      mutate(significance = ifelse(padj < 0.05, T, F))

    write.xlsx(res, outfile, overwrite = T)
}

comps <- list(
  c(stim = 'RPMI', t = 't1', t_base = 't0', outfile = 'output/DE/RPMI_timepoints_t1_t0.xlsx'),
  c(stim = 'RPMI', t = 't2', t_base = 't0', outfile = 'output/DE/RPMI_timepoints_t2_t0.xlsx'),
  c(stim = 'RPMI', t = 't3', t_base = 't0', outfile = 'output/DE/RPMI_timepoints_t3_t0.xlsx'),

  c(stim = 'Influenza', t = 't1', t_base = 't0', outfile = 'output/DE/influenza_timepoints_t1_t0.xlsx'),
  c(stim = 'Influenza', t = 't2', t_base = 't0', outfile = 'output/DE/influenza_timepoints_t2_t0.xlsx'),
  c(stim = 'Influenza', t = 't3', t_base = 't0', outfile = 'output/DE/influenza_timepoints_t3_t0.xlsx'),

  c(stim = 'R848', t = 't1', t_base = 't0', outfile = 'output/DE/R848_timepoints_t1_t0.xlsx'),
  c(stim = 'R848', t = 't2', t_base = 't0', outfile = 'output/DE/R848_timepoints_t2_t0.xlsx'),
  c(stim = 'R848', t = 't3', t_base = 't0', outfile = 'output/DE/R848_timepoints_t3_t0.xlsx'),

  c(stim = 'SARSCOV2', t = 't1', t_base = 't0', outfile = 'output/DE/sars_timepoints_t1_t0.xlsx'),
  c(stim = 'SARSCOV2', t = 't2', t_base = 't0', outfile = 'output/DE/sars_timepoints_t2_t0.xlsx'),
  c(stim = 'SARSCOV2', t = 't3', t_base = 't0', outfile = 'output/DE/sars_timepoints_t3_t0.xlsx')
)

for (comp in comps) {
  runDEseq2(counts = counts,
            meta = meta,
            stimul = comp[['stim']],
            t = comp[['t']],
            t_base = comp[['t_base']],
            outfile = comp[['outfile']])
}
