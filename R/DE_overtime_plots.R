try(dev.off())
rm(list = ls())

library(dplyr)
library(openxlsx)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggsci)

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

lapply(degs, function(x) {x %>% filter(significance) %>% nrow()}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('id') %>% 
  cbind(., colsplit(.$id, '_', c('id2', 'time'))) -> df
head(df)
df$time <- ifelse(df$time == 'T1vsT0', 'T1\nvs\nBaseline', df$time)
df$time <- ifelse(df$time == 'T2vsT0', 'T2\nvs\nBaseline', df$time)
df$time <- ifelse(df$time == 'T3vsT0', 'T3\nvs\nBaseline', df$time)


pdf('output/lineplot_nrdeg.pdf', width = 4, height = 3)
ggplot(df) +
  geom_point(aes(x = time, y = V1, color = id2)) +
  geom_line(aes(x = time, y = V1, group = id2, color = id2)) +
  theme_bw() +
  theme(aspect.ratio = 1, 
        axis.title.x = element_blank(),
        panel.grid = element_blank()) +
  labs(x = 'Time', y = 'nr. DEG', color = 'Stimulation', 
       title = 'Nr. DEG over time within stimulations') +
  ggsci::scale_color_futurama()
dev.off()
