rm(list = ls())

library(ggVennDiagram)
library(ggplot2)
library(dplyr)
library(magrittr)
library(openxlsx)

setwd('/vol/projects/CIIM/COMBI/bulk_seq/')

pval.thresh <- 0.05
logfc.thresh <- 0.5

isSig <- function(x, pt = pval.thresh, lt = logfc.thresh) {
  x$significance <- ifelse(x$padj < pt & abs(x$log2FoldChange) > lt, TRUE, FALSE)
  x$direction <- ifelse(x$log2FoldChange < 0, 'Downregulated', 'Upregulated')
  return(x)
}
getAll <- function(x) {return (x %>% filter(significance == TRUE) %>% pull(gene))}
getUp <- function(x) {
  return( x %>% filter(significance == TRUE) %>% filter(log2FoldChange > 0) %>% pull(gene))
}
getDown <- function(x) {
  return( x %>% filter(significance == TRUE) %>% filter(log2FoldChange < 0) %>% pull(gene))
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

results <- lapply(results, isSig)
deg.up <- lapply(results, getUp)
deg.down <- lapply(results, getDown)

overlaps <- list()

# Per timepoint, overlap between stimulations
for (t in c('T0', 'T1', 'T2', 'T3')) {
  names <- paste0(t, c(' Influenza', ' R848', ' SARS-CoV-2'))
  
  sub.up <- deg.up[names]
  sub.down <- deg.down[names]  

  names(sub.up) <- c('Influenza', 'R848', 'SARS-CoV-2')  
  names(sub.down) <- c('Influenza', 'R848', 'SARS-CoV-2')  
  
  venn <- Venn(sub.up)
  data <- process_data(venn)  
  
  overlaps[[paste0(t, 'up')]] <- venn_region(data) %>% as.data.frame() %>% select(name, item)
  
  pdf(paste0('output/overlap_DEG_', t, '_up.pdf'), width = 3, height = 3)
  print(ggplot() +
    geom_sf(aes(fill = count), data = venn_region(data), show.legend = F) +
    geom_sf(data = venn_setedge(data)) +
    geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
    geom_sf_label(aes(label = count), data = venn_region(data), label.size = 0, alpha = 0) +
    theme_void() +
    scale_fill_gradient(low = 'white', high = 'red') +
    labs(title = paste0('Upregulated\n', t)) +
    theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  )
  dev.off()
  
  venn <- Venn(sub.down)
  data <- process_data(venn)  
  overlaps[[paste0(t, 'down')]] <- venn_region(data) %>% as.data.frame() %>% select(name, item)
  
    pdf(paste0('output/overlap_DEG_', t, '_down.pdf'), width = 3, height = 3)
    print(ggplot() +
      geom_sf(aes(fill = count), data = venn_region(data), show.legend = F) +
      geom_sf(data = venn_setedge(data)) +
      geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
      geom_sf_label(aes(label = count), data = venn_region(data), label.size = 0, alpha = 0) +
      theme_void() +
      scale_fill_gradient(low = 'white', high = 'red') +
      labs(title = paste0('Downregulated\n', t)) +
      theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  )
  dev.off()
}

write.xlsx(x = overlaps, file = 'output/overlap_DEG_timepoints.xlsx')

overlaps <- list()

# Per stimulation, overlap between timepoints
for (s in c('Influenza', 'R848', 'SARS-CoV-2')) {
  names <- paste0(c('T0 ', 'T1 ', 'T2 ', 'T3 '), s)
  
  sub.up <- deg.up[names]
  sub.down <- deg.down[names]  
  
  names(sub.down) <- c('T0', 'T1', 'T2', 'T3')
  names(sub.up) <- c('T0', 'T1', 'T2', 'T3')
  
  venn <- Venn(sub.up)
  data <- process_data(venn)  
  overlaps[[paste0(s, 'down')]] <- venn_region(data) %>% as.data.frame() %>% select(name, item)
  
  pdf(paste0('output/overlap_DEG_stim_', s, '_up.pdf'), width = 3, height = 3)
  print(ggplot() +
          geom_sf(aes(fill = count), data = venn_region(data), show.legend = F) +
          geom_sf(data = venn_setedge(data)) +
          geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
          geom_sf_label(aes(label = count), data = venn_region(data), label.size = 0, alpha = 0) +
          theme_void() +
          scale_fill_gradient(low = 'white', high = 'red') +
          labs(title = paste0('Upregulated\n', s)) +
          theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  )
  dev.off()
  
  venn <- Venn(sub.down)
  data <- process_data(venn)  
  overlaps[[paste0(s, 'up')]] <- venn_region(data) %>% as.data.frame() %>% select(name, item)
  
  pdf(paste0('output/overlap_DEG_stim_', s, '_down.pdf'), width = 3, height = 3)
  print(ggplot() +
          geom_sf(aes(fill = count), data = venn_region(data), show.legend = F) +
          geom_sf(data = venn_setedge(data)) +
          geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
          geom_sf_label(aes(label = count), data = venn_region(data), label.size = 0, alpha = 0) +
          theme_void() +
          scale_fill_gradient(low = 'white', high = 'red') +
          labs(title = paste0('Downregulated\n', s)) +
          theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  )
  dev.off()
}

write.xlsx(x = overlaps, file = 'output/overlap_DEG_stimulations.xlsx')

