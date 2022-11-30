library(openxlsx)
library(ggplot2)
library(dplyr)
library(ggsci)
library(magrittr)

rm(list = ls())
setwd('/vol/projects/CIIM/COMBI/bulk_seq/')

pval.thresh <- 0.05
logfc.thresh <- 0.5

isSig <- function(x, pt = pval.thresh, lt = logfc.thresh) {
  x$significance <- ifelse(x$padj < pt & abs(x$log2FoldChange) > lt, TRUE, FALSE)
  x$direction <- ifelse(x$log2FoldChange < 0, 'Downregulated', 'Upregulated')
  return(x)
}

doVolcano <- function(result, outfile, plottitle, xlim, ycap, nlabel, nudge_x, nudge_y, pt = pval.thresh, lt = logfc.thresh) {

  result$log10padj <- -log10(result$padj)
  result$log10padj <- ifelse(result$log10padj > ycap, ycap, result$log10padj)

  png(outfile, width = 5, height = 5, res = 600, units = 'cm')
  print(
    ggplot() +
      geom_point(data = result %>% filter(significance == FALSE), aes(x = log2FoldChange, y = log10padj), color = 'lightgray', size = 0.5) +
      geom_point(data = result %>% filter(significance == TRUE), aes(x = log2FoldChange, y = log10padj, color = direction), size = 0.5) +
      geom_hline(yintercept = -log10(pt), lty= 2) +
      geom_vline(xintercept = lt, lty=2) +
      geom_vline(xintercept = -lt, lty=2) +
      theme_classic() +
      labs(x = expression(paste(log[2], ' fold-change')),
           y = expression(paste(-log[10], '(adj. p-value)')),
           title = plottitle) +
      theme(plot.title = element_blank(),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6),
            aspect.ratio = 1,
            legend.position = 'none')  +
      scale_color_manual(values = rev(pal_nejm()(2))) +
      xlim(xlim) +
      ylim(c(0, ycap)) +
      ggrepel::geom_label_repel(data = result %>% filter(significance) %>% head(label, n = nlabel),
                                aes(x = log2FoldChange, y = log10padj, label = gene),
                                label.padding = 0.1, box.padding=0.1,
                                label.size = NA,
                                show.legend = F, nudge_x = nudge_x, nudge_y = nudge_y,
                                size = 1.5, max.overlaps = 50,
                                min.segment.length = unit(2, 'mm')) )

  dev.off()

}

results <- list(
  'T0 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t0_influenza_RPMI.xlsx'),
  'T0 R848' =       read.xlsx(xlsxFile = 'output/DE/t0_R848_RPMI.xlsx'),
  'T0 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t0_SARSCOV2_RPMI.xlsx'),

  'T1 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t1_influenza_RPMI.xlsx'),
  'T1 R848' =       read.xlsx(xlsxFile = 'output/DE/t1_R848_RPMI.xlsx'),
  'T1 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t1_SARSCOV2_RPMI.xlsx'),

  'T2 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t2_influenza_RPMI.xlsx'),
  'T2 R848' =       read.xlsx(xlsxFile = 'output/DE/t2_R848_RPMI.xlsx'),
  'T2 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t2_SARSCOV2_RPMI.xlsx'),

  'T3 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t3_influenza_RPMI.xlsx'),
  'T3 R848' =       read.xlsx(xlsxFile = 'output/DE/t3_R848_RPMI.xlsx'),
  'T3 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t3_SARSCOV2_RPMI.xlsx')
)

results <- lapply(X = results, FUN = isSig)

# Volcanos
# T0
doVolcano(results[['T0 Influenza']],  outfile = 'output/volcano_t0_influenza_RPMI.png', plottitle = 'T0 Influenza',  xlim = c(-10, 10), ycap = 50, nlabel = 5)
doVolcano(results[['T0 R848']],       outfile = 'output/volcano_t0_R848_RPMI.png',      plottitle = 'T0 R848',       xlim = c(-25, 25), ycap = 60, nlabel = 10)
doVolcano(results[['T0 SARS-CoV-2']], outfile = 'output/volcano_t0_SARSCOV2_RPMI.png',  plottitle = 'T0 SARS-CoV-2', xlim = c(-6, 6),   ycap = 10, nlabel = 10)

# T1
doVolcano(results[['T1 Influenza']],  outfile = 'output/volcano_t1_influenza_RPMI.png', plottitle = 'T1 Influenza',  xlim = c(-10, 10), ycap = 50, nlabel = 6)
doVolcano(results[['T1 R848']],       outfile = 'output/volcano_t1_R848_RPMI.png',      plottitle = 'T1 R848',       xlim = c(-25, 25), ycap = 60, nlabel = 10)
doVolcano(results[['T1 SARS-CoV-2']], outfile = 'output/volcano_t1_SARSCOV2_RPMI.png',  plottitle = 'T1 SARS-CoV-2', xlim = c(-6, 6),   ycap = 10, nlabel = 10)

# T2
doVolcano(results[['T2 Influenza']],  outfile = 'output/volcano_t2_influenza_RPMI.png', plottitle = 'T2 Influenza',  xlim = c(-10, 10), ycap = 50, nlabel = 10)
doVolcano(results[['T2 R848']],       outfile = 'output/volcano_t2_R848_RPMI.png',      plottitle = 'T2 R848',       xlim = c(-25, 25), ycap = 60, nlabel = 10)
doVolcano(results[['T2 SARS-CoV-2']], outfile = 'output/volcano_t2_SARSCOV2_RPMI.png',  plottitle = 'T2 SARS-CoV-2', xlim = c(-6, 6),   ycap = 10, nlabel = 10)

# T3
doVolcano(results[['T3 Influenza']],  outfile = 'output/volcano_t3_influenza_RPMI.png', plottitle = 'T3 Influenza',  xlim = c(-10, 10), ycap = 50, nlabel = 5)
doVolcano(results[['T3 R848']],       outfile = 'output/volcano_t3_R848_RPMI.png',      plottitle = 'T3 R848',       xlim = c(-25, 25), ycap = 60, nlabel = 10)
doVolcano(results[['T3 SARS-CoV-2']], outfile = 'output/volcano_t3_SARSCOV2_RPMI.png',  plottitle = 'T3 SARS-CoV-2', xlim = c(-6, 6),   ycap = 10, nlabel = 10)


# Get nr of up and down
up <- lapply(results, FUN = function(x) { x %>% filter(significance == TRUE) %>% filter(direction == 'Upregulated') %>%  nrow() })
down <- lapply(results, FUN = function(x) {x %>% filter(significance == TRUE) %>% filter(direction == 'Downregulated') %>% nrow() })

up <- do.call(cbind.data.frame, up)
down <- do.call(cbind.data.frame, down)

up <- as.data.frame( t(up) ) %>% mutate('Direction' =  'Upregulated') %>% dplyr::rename('Nr' = V1)
down <- as.data.frame( t(down) ) %>% mutate('Direction' =  'Downregulated') %>% dplyr::rename('Nr' = V1)

up$comparison <- rownames(up)
down$comparison <- rownames(down)
rownames(up) <- NULL
rownames(down) <- NULL

stopifnot(all(rownames(up) == rownames(down)))
df <- rbind(up, down)

df <- cbind(df,
            reshape2::colsplit(df$comparison, pattern = ' ', names = c('time', 'stimulation')))
df$time <- ifelse(df$time == 'T0', 'Baseline', df$time)
df %>% group_by(time, stimulation) %>% summarise(tot = sum(Nr)) -> an

pdf('output/overview_nr_DEG.pdf', width = 6, height = 3)
ggplot(df) +
  geom_bar(aes(x = time, y = Nr, fill = Direction), stat = 'identity') +
  facet_wrap("stimulation", ncol = 4, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Timepoint',
       y = 'Nr. genes',
       title = 'Nr. DEG after stimulation compared to RPMI') +
  theme(strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = rev(pal_nejm()(2)))
dev.off()



# nr DEG for old vs young

results2 <- list(
  'T0 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t0_influenza_old_young.xlsx'),
  'T0 R848' =       read.xlsx(xlsxFile = 'output/DE/t0_R848_old_young.xlsx'),
  'T0 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t0_SARSCOV2_old_young.xlsx'),
  'T0 RPMI' =       read.xlsx(xlsxFile = 'output/DE/t0_RPMI_old_young.xlsx'),

  'T1 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t1_influenza_old_young.xlsx'),
  'T1 R848' =       read.xlsx(xlsxFile = 'output/DE/t1_R848_old_young.xlsx'),
  'T1 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t1_SARSCOV2_old_young.xlsx'),
  'T1 RPMI' =       read.xlsx(xlsxFile = 'output/DE/t1_RPMI_old_young.xlsx'),

  'T2 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t2_influenza_old_young.xlsx'),
  'T2 R848' =       read.xlsx(xlsxFile = 'output/DE/t2_R848_old_young.xlsx'),
  'T2 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t2_SARSCOV2_old_young.xlsx'),
  'T2 RPMI' =       read.xlsx(xlsxFile = 'output/DE/t2_RPMI_old_young.xlsx'),

  'T3 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t3_influenza_old_young.xlsx'),
  'T3 R848' =       read.xlsx(xlsxFile = 'output/DE/t3_R848_old_young.xlsx'),
  'T3 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t3_SARSCOV2_old_young.xlsx'),
  'T3 RPMI' =       read.xlsx(xlsxFile = 'output/DE/t3_RPMI_old_young.xlsx')
)

results2 <- lapply(X = results2, FUN = isSig)
up <- do.call(cbind.data.frame,
              lapply(results2, FUN = function(x) { x %>% filter(significance == TRUE) %>%
                  filter(direction == 'Upregulated') %>%  nrow() }) )

down <- do.call(cbind.data.frame,
                lapply(results2, FUN = function(x) {x %>% filter(significance == TRUE) %>%
                    filter(direction == 'Downregulated') %>% nrow() }) )

up <- as.data.frame( t(up) ) %>% mutate('Direction' =  'Upregulated') %>% dplyr::rename('Nr' = V1)
down <- as.data.frame( t(down) ) %>% mutate('Direction' =  'Downregulated') %>% dplyr::rename('Nr' = V1)

up$comparison <- rownames(up)
down$comparison <- rownames(down)
rownames(up) <- NULL
rownames(down) <- NULL

stopifnot(all(rownames(up) == rownames(down)))
df <- rbind(up, down)

df <- cbind(df,
            reshape2::colsplit(df$comparison, pattern = ' ', names = c('time', 'stimulation')))

df$stimulation <- factor(df$stimulation, levels = c('RPMI', 'Influenza', 'R848', 'SARS-CoV-2'))

png('output/overview_nr_DEG_old_young.png', width = 6, height = 3)
ggplot(df) +
  geom_bar(aes(x = time, y = Nr, fill = Direction), stat = 'identity') +
  facet_wrap("stimulation", ncol = 4, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Timepoint',
       y = 'Nr. genes',
       title = 'Nr. DEG old vs young') +
  theme(strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = rev(pal_nejm()(2)))
dev.off()



# nr DEG compared to baseline

results3 <- list(
  'T0 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t0_influenza_vs_baseline.xlsx'),
  'T0 R848' =       read.xlsx(xlsxFile = 'output/DE/t0_R848_vs_baseline.xlsx'),
  'T0 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t0_SARSCOV2_vs_baseline.xlsx'),

  'T1 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t1_influenza_vs_baseline.xlsx'),
  'T1 R848' =       read.xlsx(xlsxFile = 'output/DE/t1_R848_vs_baseline.xlsx'),
  'T1 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t1_SARSCOV2_vs_baseline.xlsx'),

  'T2 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t2_influenza_vs_baseline.xlsx'),
  'T2 R848' =       read.xlsx(xlsxFile = 'output/DE/t2_R848_vs_baseline.xlsx'),
  'T2 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t2_SARSCOV2_vs_baseline.xlsx'),

  'T3 Influenza' =  read.xlsx(xlsxFile = 'output/DE/t3_influenza_vs_baseline.xlsx'),
  'T3 R848' =       read.xlsx(xlsxFile = 'output/DE/t3_R848_vs_baseline.xlsx'),
  'T3 SARS-CoV-2' = read.xlsx(xlsxFile = 'output/DE/t3_SARSCOV2_vs_baseline.xlsx')
)


results3 <- lapply(X = results3, FUN = isSig)

up <- do.call(cbind.data.frame,
              lapply(results3, FUN = function(x) { x %>% filter(significance == TRUE) %>%
                  filter(direction == 'Upregulated') %>%  nrow() }) )

down <- do.call(cbind.data.frame,
                lapply(results3, FUN = function(x) {x %>% filter(significance == TRUE) %>%
                    filter(direction == 'Downregulated') %>% nrow() }) )

up <- as.data.frame( t(up) ) %>% mutate('Direction' =  'Upregulated') %>% dplyr::rename('Nr' = V1)
down <- as.data.frame( t(down) ) %>% mutate('Direction' =  'Downregulated') %>% dplyr::rename('Nr' = V1)

up$comparison <- rownames(up)
down$comparison <- rownames(down)
rownames(up) <- NULL
rownames(down) <- NULL

stopifnot(all(rownames(up) == rownames(down)))
df <- rbind(up, down)

df <- cbind(df,
            reshape2::colsplit(df$comparison, pattern = ' ', names = c('time', 'stimulation')))


png('output/overview_nr_DEG_baseline.png', width = 6, height = 3)
ggplot(df) +
  geom_bar(aes(x = time, y = Nr, fill = Direction), stat = 'identity') +
  facet_wrap("stimulation", ncol = 4, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Timepoint',
       y = 'Nr. genes',
       title = 'Nr. DEG compared to baseline') +
  theme(strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = rev(pal_nejm()(2)))
dev.off()




# Enrichment modules

library(rjson)
library(fgsea)
library(ComplexHeatmap)

data <- fromJSON(file = "/vol/projects/mzoodsma/BTM/BTM/test.json")

pathways <- list()
for (name in names(data)) { pathways[[data[[name]][['name']]]] <- data[[name]][['genes']] }

deg <- results$`T0 Influenza` %>% filter(significance ==TRUE) %>% filter(direction == 'Upregulated')

# Input FGSEA:
  # Names should be the genes
  # Values should be the ranking (scoreType = pos) or t-test statistic (scoreType=  std)

genes <- deg %>%  pull(gene)
pvals <- deg %>% pull(stat)

names(pvals) <- genes

res <- fgsea(pathways, pvals, scoreType='pos')
res <- rbind(head(res, 10), tail(res, 10))


mat <- matrix(data = res$NES, dimnames = list(res$pathway, 'NES'))
Heatmap(mat)

enr <- function(results, pathways) {
  deg <- results %>% filter(significance ==TRUE) #%>% filter(direction == 'Upregulated')

  genes <- deg %>%  pull(gene)
  pvals <- deg %>% pull(stat)

  names(pvals) <- genes

  res <- fgsea(pathways, pvals, scoreType='std')
  res %<>% arrange(desc(NES))

  res <- rbind(head(res, 5), tail(res, 5))
  mat <- matrix(data = res$NES, dimnames = list(res$pathway, 'NES'))
  mat <- data.frame(mat)
  return(mat)
}


enr.all <- lapply(results, enr, pathways = pathways)
enr.all

# Influenza
res <- list('T0 Influenza' = enr.all$`T0 Influenza`,
            "T1 Influenza" = enr.all$`T1 Influenza`,
            "T2 Influenza" = enr.all$`T2 Influenza`,
            "T3 Influenza" =  enr.all$`T3 Influenza`)

res2 <- lapply(res, function(x) {x$rn <- rownames(x); return(x)})
for (name in names(res2)) { res2[[name]][, name] <- res2[[name]]$NES; res2[[name]] %<>% select(-NES) }

new <- Reduce(function(x, y) {full_join(x, y, by = 'rn')}, res2)
rownames(new) <- new$rn
new %<>% select(-rn)

png('output/heatmap_influenza_BTM_enr.png', width = 6, height = 10)
Heatmap(as.matrix(new), cluster_rows = F, cluster_columns = F,
        na_col = 'lightgray', column_names_side = 'top', column_names_rot = 45,
        row_names_side = 'left', rect_gp = gpar(col = 'white', lwd = 2))
dev.off()
