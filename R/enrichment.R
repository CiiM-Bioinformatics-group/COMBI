library(openxlsx)
library(ggplot2)
library(dplyr)
library(clusterProfiler)

setwd('/vol/projects/CIIM/COMBI/bulk_seq/')

pval.thresh <- 0.05
logfc.thresh <- 0.5

enrich <- function(genes, thresh.incl=5, n.show=5) {
  
  if (length(up) == 0 | length(down) == 0) {
    print('No up/downregulated genes. Skipping enrichment')
    return()
  }
  
  enr.GO <- enrichGO(gene = genes,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
  return(enr.GO)
  # 
  # res = list('GO Up' = enr.GO.up,
  #            'GO Down' = enr.GO.down)
  # merged <- merge_result(res)
  # 
  # return(merged)
  # 
  # # Construct the plot manually
  # df <- merged@compareClusterResult %>% filter(Count >= thresh.incl)
  # 
  # # If nothing left to plot, return and continue
  # if (nrow(df) == 0) { return() }
  # 
  # df$old_cluster <- df$Cluster
  # df <- cbind(df,
  #             reshape2::colsplit(df$Cluster, pattern = ' ', names = c('enr', 'direction')))
  # 
  # # Mutate the df to the final df we use for plotting
  # df %<>%
  #   mutate(GeneRatio2 = sapply(df$GeneRatio, function(x) eval(parse(text=x)))) %>%
  #   mutate(Description = stringr::str_wrap(.$Description, width = 50)) %>%
  #   group_by(old_cluster) %>%
  #   dplyr::slice(1:n.show)
  # 
  # # order the df as we want and set the description factor last
  # # to ensure proper plotting order
  # df$enr <- factor(df$enr, levels = c('GO', 'KEGG'))
  # df$direction <- factor(df$direction, levels = c('Down', 'Up'))
  # df %<>% arrange(enr, direction)
  # df$Description <- factor(df$Description, levels = df$Description)
  # 
  # # Hacky stuff
  # ratio <- nrow(df) * 0.5
  # 
  # pdf(outpath_pdf, width = 10, height = nrow(df) * 0.5)
  # print(ggplot(data = df) +
  #         geom_point(aes(x = direction, y = Description, size = GeneRatio2, color = p.adjust)) +
  #         DOSE::theme_dose() +
  #         theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
  #               strip.background = element_rect(fill = 'white', color = 'black'),
  #               axis.title.y = element_blank(),
  #               axis.title.x = element_blank(), aspect.ratio = ratio) +
  #         scale_colour_gradient(limits=c(0, 0.05), low="red", high = 'blue') +
  #         scale_x_discrete(drop=F) +
  #         labs(color = 'Adj. P', size = 'Gene ratio') +
  #         facet_wrap('enr'))
  # dev.off()
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

getUp <- function(x, pt = pval.thresh, lt = logfc.thresh) {
  x$significance <- ifelse(x$padj < pt & abs(x$log2FoldChange) > lt, TRUE, FALSE)
  return( x %>% filter(significance == TRUE) %>% filter(log2FoldChange > 0) %>% pull(gene))
}
getDown <- function(x, pt = pval.thresh, lt = logfc.thresh) {
  x$significance <- ifelse(x$padj < pt & abs(x$log2FoldChange) > lt, TRUE, FALSE)
  return( x %>% filter(significance == TRUE) %>% filter(log2FoldChange < 0) %>% pull(gene))
}
conv <- function(x) {
  return(bitr(geneID = x, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb="org.Hs.eg.db") %>% pull(ENTREZID))
}

up <- lapply(results, getUp)
down <- lapply(results, getDown)

up <- lapply(up, conv)
down <- lapply(down, conv)

up.enr <- lapply(up, enrich)
down.enr <- lapply(down, enrich)

up.enr <- merge_result(up.enr) 
down.enr <- merge_result(down.enr)

# Up / Down
n.show <- 10
n.thresh <- 5
df <- up.enr
df <- data.frame(df@compareClusterResult)

df <- cbind(df, 
            reshape2::colsplit(df$Cluster, pattern = ' ', names = c('Time', 'Stimulation'))) %>% 
        mutate(GeneRatio2 = sapply(df$GeneRatio, function(x) eval(parse(text=x)))) %>%
        filter(Count >= n.thresh) %>% 
        # mutate(Description = stringr::str_wrap(.$Description, width = 50)) %>% 
        group_by(Cluster) %>%
        dplyr::slice(1:n.show)

# Order everything in the way that we want
df$Time <- ifelse(df$Time == 'T0', 'Baseline', df$Time)
df$Time <- factor(df$Time, levels = c('Baseline', 'T1', 'T2', 'T3'))
df$Stimulation <- factor(df$Stimulation, levels = c('Influenza', 'R848', 'SARS-CoV-2'))
df$Description <- factor(df$Description, levels = unique(df$Description))

# Hacky stuff
ratio <- nrow(df) * .05


pdf('output/enr_down_stim_vs_RPMI.pdf', width = 12, height = 15)
print(ggplot(data = df) +
        geom_point(aes(x = Time, y = Description, size = GeneRatio2, color = p.adjust)) +
        DOSE::theme_dose() +
        theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
              strip.background = element_rect(fill = 'white', color = 'black'),
              strip.text = element_text(face = 'bold', size = 15),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(), 
              plot.title = element_text(hjust = 0.5)) +
              # aspect.ratio = ratio) +
        scale_colour_gradient(limits=c(0, 0.05), low="darkred", high = 'darkblue') +
        scale_x_discrete(drop=F) +
        labs(color = 'Adj. P', size = 'Gene ratio', title = 'Downregulated processes') +
        facet_wrap('Stimulation', ncol=3, nrow = 1)
)
dev.off()
