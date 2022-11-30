# All the PCAs
# All together, individual timepoints, stimulation and so on

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



##### PCA by stimulation and timepoints. Everything together
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = meta, 
                              design = ~ stim)
vsd <- varianceStabilizingTransformation(dds)

pca <- plotPCA(vsd, intgroup = c('stimulation', 'timepoint'), returnData=T)

stim <- ggplot(pca) +
  geom_point(aes(PC1, PC2, color = stimulation)) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  scale_color_nejm() +
  labs(x = 'PC1 [38%]', y = 'PC2 [22%]', color = 'Stimulation', title = 'PCA by stimulation')

time <- ggplot(pca) +
  geom_point(aes(PC1, PC2, color = timepoint)) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  scale_color_npg() +
  labs(x = 'PC1 [38%]', y = 'PC2 [22%]', color = 'Time', title = 'PCA by time')

pdf('output/PCA_time_stim.pdf', width = 8, height = 4)
ggpubr::ggarrange(time, stim, ncol = 2, align = 'hv')
dev.off()


###### PCA per stimulation. Do we see timepoint difference?
for (s in c('Influenza', 'SARSCOV2', 'R848', 'RPMI')) {
  
  meta_sub <- meta %>% filter(stimulation == s)
  counts_sub <- counts[, rownames(meta_sub)]
  
  stopifnot(colnames(counts_sub) == rownames(meta_sub))
  
  dds <- DESeqDataSetFromMatrix(countData = counts_sub, 
                                colData = meta_sub, 
                                design = ~ timepoint)
  
  vsd <- varianceStabilizingTransformation(dds)
  
  pca <- plotPCA(vsd, intgroup = 'timepoint', returnData=T)
  p <- ggpubr::ggscatterhist(data = pca, x = 'PC1', y = 'PC2', 
                             color = 'timepoint', fill = 'timepoint', palette = 'npg', 
                             alpha = 1, margin.plot = 'density', legend = 'right', title = s)
  p$sp <- p$sp + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  
  pdf(paste0('output/PCA_alone_', s, '.pdf'), width = 4, height = 4, onefile = F)
  print(p)
  dev.off()
}

###### PCA per timepoint. Do we see difference between the stimulations over the timepoints?
load('output/data.RData')
t <- 't3'
for (t in c('t0', 't1', 't2', 't3')) {
  
  meta_sub <- meta %>% filter(timepoint == t)
  counts_sub <- counts[, rownames(meta_sub)]
  
  stopifnot(colnames(counts_sub) == rownames(meta_sub))
  
  dds <- DESeqDataSetFromMatrix(countData = counts_sub, 
                                colData = meta_sub, 
                                design = ~ stimulation)
  
  vsd <- varianceStabilizingTransformation(dds)  
  
  plotPCA(vsd, 'stimulation')
  
  pca <- plotPCA(vsd, intgroup = 'stimulation', returnData=T)
  
  pdf(paste0('output/PCA_timepoint_', t, '.pdf'), width = 4, height = 4, onefile = F)
  print(ggplot(pca) +
    geom_point(aes(PC1, PC2, color = stimulation)) +
    scale_color_nejm() +
    theme_classic() +
    theme(aspect.ratio = 1, 
          plot.title = element_text(hjust = 0.5)) +
    labs(x = 'PC1 [39%]', y = 'PC2 [25%]', color = 'Stimulation', title = paste0('Time: ', t)))
  dev.off()
}

?plotPCA

##### UMAP per stimulation separately. Do we see timepoint difference?
library(umap)

for (s in c('Influenza', 'SARSCOV2', 'R848')) {
  
  meta_sub <- meta %>% filter(stimulation == s)
  counts_sub <- counts[, rownames(meta_sub)]
  
  um <- umap(t(counts_sub))
  stopifnot(all(rownames(um$layout) == rownames(meta_sub)))

  colnames(um$layout) <- c('UMAP_1', 'UMAP_2')
  df <- cbind(meta_sub, um$layout)

  df$timepoint <- as.factor(df$timepoint)
  
  p <- ggpubr::ggscatterhist(data = df, x = 'UMAP_1', y = 'UMAP_2', 
                        color = 'timepoint', margin.plot = 'density', 
                        fill = "timepoint", 
                        palette = 'npg', title = paste0('UMAP ', s))  
  
  p$sp <- p$sp + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  
  pdf(paste0('output/UMAP_', s, '.pdf'), width = 4, height = 4, onefile=F)
  print(p)
  dev.off()
}



##### Heatmap of correlations to PCs

for (s in c('Influenza', 'SARSCOV2', 'R848')) {
  
  meta_sub <- meta %>% filter(stimulation == s)
  counts_sub <- counts[, rownames(meta_sub)]
  print(unique(meta_sub$stimulation))
  counts_sub <- counts_sub[!apply(X = counts_sub, MARGIN = 1, FUN = var) == 0, ]
  
  dds <- DESeqDataSetFromMatrix(counts_sub, meta_sub, design = ~1)
  dds <- vst( dds, blind = T )
  rv <- rowVars(assay(dds))
  select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  pca <- prcomp(t(assay(dds)[select,]))
  
  variables <- meta_sub %>% select(age, sex, timepoint)
  
  variables$sex <- as.numeric(as.factor(variables$sex))
  variables$timepoint <- as.numeric(as.factor(variables$timepoint))

  n.pcs <- 10
  pcs <- pca$x[, 1:n.pcs]
  
  corr.pvalues <- matrix(nrow = ncol(variables), 
                         ncol = ncol(pcs), 
                         dimnames = list(colnames(variables), colnames(pcs)))
  
  for (pc in colnames(pcs)) {
    for (vars in colnames(variables)) {
      res <- cor.test(pcs[, pc], variables[, vars], method = 'spearman')

      if (res$estimate < 0) {
        pval <- -1 * log10(res$p.value)
      } else {
        pval <- 1 * log10(res$p.value)
      }

      corr.pvalues[vars, pc] <- pval
      
      
    }
  }
  
  threshold <- log10(0.05)
  corr.pvalues[abs(corr.pvalues) < abs(threshold)] <- NA
  corr.pvalues <- reshape2::melt(corr.pvalues)
  
  corr.pvalues$value <- abs(corr.pvalues$value)
  
  pdf(paste0('output/heatmap_corr_', s, '.pdf'), width = 6, height = 3)
  print(ggplot(data = corr.pvalues) + 
    geom_tile(aes(Var2, Var1, fill = value)) +
    coord_equal() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), 
          axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_fill_gradient2(low = 'white', na.value = 'white', high = 'red') +
    labs(title = s, fill = expression(-Log[10](p-value)))
  )
  dev.off()
}

