try(dev.off())
rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)
library(fgsea)

setwd('/vol/projects/CIIM/COMBI/bulk_seq/analysis')

pulendran <- read.xlsx('../sigBTMs_Pulendran.xlsx')
data <- readRDS('../output/gsea_BTM_RPMIovertime.rds')

# T1 against day 21 vs BL
own <- data %>% filter(grp == 't1')
pul <- pulendran %>% filter(comparison == 'BL_vs_Day21')

own %<>% filter(pathway %in% pul$pathway) %>% arrange(match(pathway, pul$pathway))
own$sig <- ifelse(own$padj < 0.05, T, F)
stopifnot(all(own$pathway == pul$pathway))

label <- data.frame(own$pathway, own$NES, pul$NES)
label$incl <- ifelse(sign(label$own.NES) == sign(label$pul.NES), T, F)

png('../output/compNES_arunachalamday21.png', width = 5, height = 5, res=300, units='in')
ggplot(mapping = aes(x = own$NES, y = pul$NES)) +
  geom_point(aes(color = own$sig)) +
  theme_bw() +
  labs(x = 'NES (this study)', 
       y = 'NES (Arunachalam et al., 2021)', 
       color = 'Sig in COMBI?',
       title = 'Comparison vs Arunchalam et al., \nDay21 vs BL') +
  theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = .5)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'lightblue')) +
  ggrepel::geom_label_repel(data = label %<>% filter(incl), 
                            aes(x = own.NES, y = pul.NES, label = own.pathway), 
                            label.size = 0, size = 1.5, max.overlaps = 100)
dev.off()  





# T1 against day 22 vs BL
own <- data %>% filter(grp == 't1')
pul <- pulendran %>% filter(comparison == 'BL_vs_Day22')

own %<>% filter(pathway %in% pul$pathway) %>% arrange(match(pathway, pul$pathway))
own$sig <- ifelse(own$padj < 0.05, T, F)
stopifnot(all(own$pathway == pul$pathway))

label <- data.frame(own$pathway, own$NES, pul$NES)
label$incl <- ifelse(sign(label$own.NES) == sign(label$pul.NES), T, F)

png('../output/compNES_arunachalamday22.png', width = 5, height = 5, res=300, units='in')
ggplot(mapping = aes(x = own$NES, y = pul$NES)) +
  geom_point(aes(color = own$sig)) +
  theme_bw() +
  labs(x = 'NES (this study)', 
       y = 'NES (Arunachalam et al., 2021)',        
       color = 'Sig in COMBI?',
       title = 'Comparison vs Arunchalam et al., \nDay22 vs BL') +
  theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = .5)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'lightblue')) +
  ggrepel::geom_label_repel(data = label %<>% filter(incl), 
                            aes(x = own.NES, y = pul.NES, label = own.pathway), 
                            label.size = 0, size = 1.5, max.overlaps = 100)
dev.off()  




# T1 against day 23 vs BL
own <- data %>% filter(grp == 't1')
pul <- pulendran %>% filter(comparison == 'BL_vs_Day23')

own %<>% filter(pathway %in% pul$pathway) %>% arrange(match(pathway, pul$pathway))
own$sig <- ifelse(own$padj < 0.05, T, F)
stopifnot(all(own$pathway == pul$pathway))

label <- data.frame(own$pathway, own$NES, pul$NES)
label$incl <- ifelse(sign(label$own.NES) == sign(label$pul.NES), T, F)

png('../output/compNES_arunachalamday23.png', width = 5, height = 5, res=300, units='in')
ggplot(mapping = aes(x = own$NES, y = pul$NES)) +
  geom_point(aes(color = own$sig)) +
  theme_bw() +
  labs(x = 'NES (this study)', 
       y = 'NES (Arunachalam et al., 2021)', 
       color = 'Sig in COMBI?',
       title = 'Comparison vs Arunchalam et al., \nDay23 vs BL') +
  theme(aspect.ratio = 1, panel.grid = element_blank(), plot.title = element_text(hjust = .5)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c('TRUE' = 'darkred', 'FALSE' = 'lightblue')) +
  ggrepel::geom_label_repel(data = label %<>% filter(incl), 
                            aes(x = own.NES, y = pul.NES, label = own.pathway), 
                            label.size = 0, size = 1.5, max.overlaps = 100)
dev.off()  
  