
rm(list = ls())
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(ggsci)

# To make QC plots:
#   Mapped / unmapped reads ratio barplot
#   featureCounts proportion barplot
#   DESeq2 PCA for all stimulations
#   Sample-sample heatmap
#   Distributions of male / female over timepoints and stimulations
#     Is there any confounded factor? If not, we can correct for this in the model

setwd('/vol/projects/CIIM/COMBI/bulk_seq/analysis/')
load('../output/data.RData')


# DESeq2 PCAs
# Complete for all samples
# Per timepoint
# Per stimulation

# Complete PCA
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = meta, 
                              design = ~ stim)
vsd <- varianceStabilizingTransformation(dds)

pca <- plotPCA(vsd, intgroup = c('stimulation', 'timepoint'), returnData=T)

pdf('../output/PCA_stimulation.pdf', width = 4, height = 4)
stim <- ggplot(pca) +
  geom_point(aes(PC1, PC2, color = stimulation)) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  scale_color_nejm() +
  labs(x = 'PC1 [38%]', y = 'PC2 [22%]', color = 'Stimulation', title = 'PCA by stimulation')
stim
dev.off()

pdf('../output/PCA_time.pdf', width = 4, height = 4)
time <- ggplot(pca) +
  geom_point(aes(PC1, PC2, color = timepoint)) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  scale_color_npg() +
  labs(x = 'PC1 [38%]', y = 'PC2 [22%]', color = 'Time', title = 'PCA by time')
time
dev.off()

pdf('../output/PCA_time_stim.pdf', width = 8, height = 4)
ggpubr::ggarrange(time, stim, ncol = 2, align = 'hv')
dev.off()


###### PCAs for the timepoints and stimulations separately
x <- 'R848'
class = 'stimulation'

meta_sub <- meta %>% filter(!!rlang::sym(class) == x)
counts_sub <- counts[, colnames(counts) %in% rownames(meta_sub)]

stopifnot(colnames(counts_sub) == rownames(meta_sub))

dds <- DESeqDataSetFromMatrix(countData = counts_sub, 
                              colData = meta_sub, 
                              design = ~ timepoint)

vsd <- varianceStabilizingTransformation(dds)

pca <- plotPCA(vsd, intgroup = 'timepoint', returnData=T)

pdf(paste0('../output/PCA_time_', x, '.pdf'), width = 4, height = 4)
ggplot(pca) +
  geom_point(aes(PC1, PC2, color = timepoint)) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust =0.5)) +
  labs(title = paste0('Stimulation: ', x))
dev.off()



# Sample-sample heatmaps
library(ComplexHeatmap)
library(ggsci)

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = meta, 
                              design = ~ stimulation)

vsd <- varianceStabilizingTransformation(dds)
sampleDists<- cor(assay(vsd))
sampleDistMatrix <- as.matrix(sampleDists)

#rownames(sampleDistMatrix) <- paste(vsd$timepoint, vsd$stimulation, sep="-")
#colnames(sampleDistMatrix) <- NULL
#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

cols <- pal_nejm()(4)
stims <- c('RPMI' = cols[[1]], 'Influenza' = cols[[2]], 'R848' = cols[[3]], 'SARSCOV2' = cols[[4]])
cols <- pal_npg()(4)
times <- c('t0' = cols[[1]], 't1' = cols[[2]], 't2' = cols[[3]], 't3' = cols[[4]])

row_ha = rowAnnotation(Timepoint = as.factor(vsd$timepoint), 
                       Stimulation = as.character(vsd$stimulation), 
                       annotation_name_side ='top', annotation_name_rot = 45, 
                       col = list(Stimulation = stims, timepoint = times))

pdf('../output/heatmap_samples.pdf', width = 7, height = 5)
Heatmap(matrix = sampleDistMatrix, 
        right_annotation = row_ha, 
        name = 'Correlation', show_row_names = F, 
        show_column_names = F, clustering_distance_rows = 'pearson', clustering_distance_columns = 'pearson')
dev.off()



###### MultiQC plots to indicate which samples we leave out
rm(list = ls())

df <- read.csv('../results/multiqc/star/multiqc_data/multiqc_featureCounts.txt', sep = '\t', row.names=1) %>% 
  mutate(sampleID = rownames(.)) %>% 
  select(Total, sampleID) %>% 
  reshape2::melt(.)

# pdf('../output/barplot_reads.pdf', width = 6, height = 6)
nreads <- ggplot(df) +
  geom_bar(aes(x = value, y = sampleID), position = 'stack', stat = 'identity') +
  theme_classic() +
  scale_fill_npg() +
  theme(axis.text.y = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Nr. reads', y = 'Samples', title = 'Total nr reads per sample')

df2 <- read.csv('../results/multiqc/star/multiqc_data/mqc_star_alignment_plot_1.txt', sep = '\t', row.names=1) %>% 
  mutate(sampleID = rownames(.)) %>% 
  set_colnames(c('Uniquely mapped', 'Mapped to multiple loci', 'Mapped to too many loci', 'Unmapped (too short)', 'Unmapped (other)', 'sampleID')) %>% 
  reshape2::melt(.)


df2$variable <- factor(df2$variable, 
                       levels = rev(c('Uniquely mapped', 'Mapped to multiple loci', 'Mapped to too many loci', 'Unmapped (too short)', 'Unmapped (other)')))

status <- ggplot(df2) +
  geom_bar(aes(x = value, y = sampleID, fill = variable), position = 'stack', stat = 'identity') +
  theme_classic() +
  # scale_fill_nejm() +
  scale_fill_manual(values = rev(pal_nejm()(5))) +
  theme(axis.text.y = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Nr. reads', y = 'Samples', title = 'STAR mapping per sample', fill= 'Status')

pdf('../output/QC_barplots.pdf', width = 12, height = 6)
ggpubr::ggarrange(nreads, status, align = 'h')

dev.off()


# median nr of reads
df <- read.csv('../results/multiqc/star/multiqc_data/multiqc_featureCounts.txt', sep = '\t', row.names=1)
head(df)
df %>% filter(Total > 2e7) -> filt

median(filt$Total) / 1e6





###### PCA for each stimulation comapred to RPMI
rm(list = ls())
load('../output/data.RData')

# Complete PCA
total <- list()
cols <- list(
  'Influenza' = "#BC3C29FF", 
  'R848' = "#0072B5FF", 
  'SARSCOV2' = "#20854EFF"
)
stim <- 'R848'

meta_sub <- meta[which(meta$stimulation %in% c(stim, 'RPMI')), ]
counts_sub <- counts[, colnames(counts) %in% rownames(meta_sub)]

dds <- DESeqDataSetFromMatrix(countData = counts_sub, 
                              colData = meta_sub, 
                              design = ~ stim)
vsd <- varianceStabilizingTransformation(dds)

plotPCA(vsd, intgroup = c('stimulation'))
pca <- plotPCA(vsd, intgroup = 'stimulation', returnData=T)

total[[stim]] <- ggplot(pca) +
  geom_point(aes(PC1, PC2, color = stimulation)) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = rev(c(cols[[stim]], '#E18727FF'))) +
  labs(x = 'PC1 [17%]', y = 'PC2 [19%]', color = 'Stimulation', title = 'SARS-CoV-2')

pdf('../output/PCA_stim_separate.pdf', width = 10, height = 3)
ggpubr::ggarrange(plotlist = total, nrow = 1, ncol = 3, common.legend = F)
dev.off()
