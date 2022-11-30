rm(list = ls())

library(dplyr)
library(magrittr)

setwd('/vol/projects/CIIM/COMBI/bulk_seq/analysis/')

counts <- read.csv('../results/star/featurecounts.merged.counts.tsv', header=T, row.names=1, sep = '\t')
counts %<>% select(-gene_name) 

mapfile <- read.csv('../design.csv')

# Remove the two resequenced samples that are duplicated in the mapfile
rem <- paste0(c('2-V4-R8A2', '6-V2-IA2'), collapse = '|')
mapfile <- mapfile[grep(mapfile$fastq_1, perl = T, pattern = rem, invert = T), ]
rownames(mapfile) <- paste0(mapfile$group, '_R', mapfile$replicate)

mapfile <- mapfile[match(colnames(counts), rownames(mapfile)), ]

stopifnot(all(colnames(counts) == rownames(mapfile)))

meta <- read.csv('../metadata_bulkRNA_COMBI.csv', header=T, row.names=1)
meta$sampleID <- gsub(meta$sampleID, pattern = '_', replacement = '-')

meta %<>% arrange(match(sampleID, mapfile$sampleID))
stopifnot(all(meta$sampleID == mapfile$sampleID))


meta <- cbind(mapfile, meta %>% select(-sampleID, -timepoint))


# Remove the unwanted samples
read.csv('../results/multiqc/star/multiqc_data/multiqc_featureCounts.txt', sep = '\t') -> featurecounts
featurecounts %>% filter(Assigned < 0.5*1e7) %>% pull(Sample) -> rem

counts <- counts[, !colnames(counts) %in% rem]
meta <- meta[!rownames(meta) %in% rem, ]

# # Add the batch information frmo Christine
# batches <- read.csv('../batches.csv', header=T, row.names=1)
# batches$Sample.name <- gsub(pattern = '_', replacement = '-', x = batches$Sample.name)
# colnames(batches) <- c('sampleID', 'sample_nr_BGI', 'batch')
# 
# meta <- left_join(x = meta, y = batches, by = 'sampleID')
# rownames(meta) <- paste0(meta$group, '_R', meta$replicate)

rm(rem, mapfile, featurecounts)
stopifnot(all(rownames(meta) == colnames(counts)))

save.image('../output/data.RData')
