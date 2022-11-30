# Manual PCA on the DEG


rm(list = ls())
dev.off()
setwd('/vol/projects/CIIM/COMBI/bulk_seq/')

load('output/data.RData')
# Any gene that is DEG between any of the stimulations

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

getSig <- function(x, pt = 0.05, lt = 0.5) {
  y <- x %>% filter(padj < pt & abs(x$log2FoldChange) > lt) %>% pull(gene)
  return(y)
}

deg <- lapply(results, getSig)

# For each stimulation compared to RPMI

for (s in c('Influenza', 'R848', 'SARS-CoV-2')) {
  
  
  meta_sub <- meta %>% filter(stim == 'RMPI' | stim == s)
  print(table(meta_sub$timepoint, meta_sub$stim))
  # Get the DEG
  n <- names(deg)[grepl(s, names(deg))]
  de <- deg[n] %>% unlist() %>% unname() %>% unique()

  counts_sub <- counts[de, rownames(meta_sub)]  

  prcomp(t(counts_sub), center=T, scale.=T) -> pca
  
  pdf(paste0('output/PCA_DEG_', s, '.pdf') ,width = 4, height = 4)
  print(ggplot(data = data.frame(pca$x)) +
    geom_point(aes(PC1, PC2, color = meta_sub$stim, shape = meta_sub$age_class)) +
    labs(x = 'PC1', y = 'PC2', color = 'Stim') +
    theme_classic() +
    theme(aspect.ratio = 1))
  dev.off()  
  
}

