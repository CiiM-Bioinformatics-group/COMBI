try(dev.off())
rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(fgsea)
library(reshape2)
library(ggsci)
library(ComplexHeatmap)

setwd('/vol/projects/CIIM/COMBI/bulk_seq')

degs <- list(
  'R848_T1vsT0' = read.xlsx('output/DE/R848_timepoints_t1_t0.xlsx'), 
  'R848_T2vsT0' = read.xlsx('output/DE/R848_timepoints_t2_t0.xlsx'), 
  'R848_T3vsT0' = read.xlsx('output/DE/R848_timepoints_t3_t0.xlsx'), 
  
  'influenza_T1vsT0' = read.xlsx('output/DE/influenza_timepoints_t1_t0.xlsx'), 
  'influenza_T2vsT0' = read.xlsx('output/DE/influenza_timepoints_t2_t0.xlsx'), 
  'influenza_T3vsT0' = read.xlsx('output/DE/influenza_timepoints_t3_t0.xlsx'), 
  
  'sars_T1vsT0' = read.xlsx('output/DE/sars_timepoints_t1_t0.xlsx'), 
  'sars_T2vsT0' = read.xlsx('output/DE/sars_timepoints_t2_t0.xlsx'), 
  'sars_T3vsT0' = read.xlsx('output/DE/sars_timepoints_t3_t0.xlsx'), 
  
  'rpmi_T1vsT0' = read.xlsx('output/DE/RPMI_timepoints_t1_t0.xlsx'), 
  'rpmi_T2vsT0' = read.xlsx('output/DE/RPMI_timepoints_t2_t0.xlsx'), 
  'rpmi_T3vsT0' = read.xlsx('output/DE/RPMI_timepoints_t3_t0.xlsx')
)

runGSEA <- function(res, pathways) {
  
  # Names are the genes, values are t-test statistic
  ranked_list <- res$stat
  names(ranked_list) <- res$gene
  
  # Remove the duplicates values
  # ranked_list <- ranked_list[!duplicated(names(ranked_list))]
  sort(ranked_list, decreasing = T) -> ranked_list
  
  fgseaRes <- fgsea(pathways    = pathways, 
                    stats       = ranked_list, 
                    nPermSimple = 100)
  
  return (fgseaRes)
}

btm <- gmtPathways('/vol/projects/CIIM/resources/BTM/BTM/datasets/BTM_for_GSEA_20131008.gmt')
# Remove the 'TBA'
btm <- btm[!grepl('TBA*', names(btm))]

res <- lapply(degs, runGSEA, pathways = btm)

for (x in names(res)) { res[[x]]$grp <- x }

do.call(rbind, res) -> res
res <- as.data.frame(res)

unique_pways <- res %>% filter(padj < 0.0001) %>% pull(pathway)
# saveRDS(object = res, file = '../output/gsea_btms_overtime.rds')
dcast(pathway ~ grp, data=res %>% filter(pathway %in% unique_pways), value.var = 'NES') %>% as.data.frame() %>% tibble::column_to_rownames('pathway')-> mat
Heatmap(mat) -> hm
ro <- row_order(hm)

res %<>% filter(padj < 0.0001)
dcast(pathway ~ grp, data=res, value.var = 'NES') %>% as.data.frame() %>% tibble::column_to_rownames('pathway')-> mat

head(mat)
mat$influenza_T1vsT0 <- NA
mat %<>% select(c(rpmi_T1vsT0, rpmi_T2vsT0, rpmi_T3vsT0, 
                  influenza_T1vsT0, influenza_T2vsT0, influenza_T3vsT0, 
                 R848_T1vsT0, R848_T2vsT0, R848_T3vsT0, 
                 sars_T1vsT0, sars_T2vsT0, sars_T3vsT0))

pdf('output/btms_overtime.pdf', width = 5, height = 4)
Heatmap(mat, cluster_columns = F, cluster_rows = F, row_order = ro, name = 'NES',
        na_col = 'lightgray', column_names_side = 'bottom', column_names_rot = 45, 
        row_names_gp = gpar(fontsize=6), column_names_gp = gpar(fontsize = 6),  
        column_split = factor(c('RPMI', 'RPMI', 'RPMI', 'Influenza', 'Influenza', 'Influenza', 'R848', 'R848', 'R848', 'SARS', 'SARS', 'SARS'), levels = c('RPMI', 'Influenza', 'R848', 'SARS')), column_gap = unit(2, 'mm'),
        column_labels = c(rep(c('t1', 't2', 't3'), 4)), column_title_gp = gpar(fontsize=8),
        row_names_side = 'left', rect_gp = gpar(col = 'white', lwd = 2))
dev.off()

?Heatmap
