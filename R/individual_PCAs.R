library(openxlsx)
library(ggplot2)
library(dplyr)
library(DESeq2)

rm(list = ls())
setwd('/vol/projects/CIIM/COMBI/bulk_seq/')
load('output/data.RData')

meta$age_class <- as.factor(ifelse(meta$age < 40, 'young', 'old'))
pls <- list()

for (s in c('R848', 'Influenza', 'SARSCOV2')) {
  
  for (t in c('t0', 't1', 't2', 't3')) {

    meta_sub <- meta %>% 
      filter(timepoint == t) %>% 
      filter(stimulation == s)
    
    counts_sub <- counts[, rownames(meta_sub)]
    
    dds <- DESeqDataSetFromMatrix(countData = counts_sub, 
                                  colData = meta_sub, 
                                  design = ~ 1)
    
    vsd <- varianceStabilizingTransformation(dds)
    p <- plotPCA(vsd, intgroup = 'sampleID')

    pca <- plotPCA(vsd, intgroup = colnames(meta_sub), returnData=T)
    
    ylabel <- p$labels$y
    xlabel <- p$labels$x
    
    pl <- ggplot(pca) +
      geom_point(aes(PC1, PC2, color = age_class)) +
      theme_classic() +
      labs(x = xlabel, y = ylabel, title = paste0(s, '\n', t)) +
      theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1, 
            legend.position = 'none') +
      scale_color_tron()
      
    name <- paste0(t, ' ',  s)
    pls[[name]] <- ggplotGrob(pl)

  }
}

library(cowplot)
leg <- get_legend(ggplot(pca) + geom_point(aes(PC1, PC2, color = meta_sub$age_class)) + 
                    scale_color_tron() + labs(color = 'Age') + theme_classic())

pdf('output/individual_PCAs.pdf', width = 10, height = 8)
plot_grid(
  plot_grid(plotlist = pls, ncol = 4, nrow = 3), 
  leg, rel_widths = c(9, 1))
dev.off()

# pdf('output/individual_PCAs.pdf', width = 5, height = 4)
ggplot(df) +
  geom_point(aes(PC1, PC2, color = ifelse(age < 40, 'Young', 'Old'))) +
  theme_bw() +
  facet_grid(stimulation ~ timepoint) +
  # labs(x = xlabel, y = ylabel, color = 'Age') +
  # labs(x = 'PC1', y = 'PC2', color = 'Age') +
  theme(strip.background = element_rect(fill = 'white'), 
        strip.text = element_text(face = 'bold'), 
        aspect.ratio = 1) +
  scale_x_continuous(breaks = seq(-20, 20, 20)) +
  scale_y_continuous(breaks = seq(-20, 20, 20)) +
  scale_color_manual(values = rev(ggsci::pal_nejm()(2)))
dev.off()

