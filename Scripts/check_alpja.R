#checando con la data

library(tidyverse)

tabla_div<- read.delim("Data/qs.csv", sep = ",")

ojiss<- tabla_div %>% filter(ids=="AJIS") %>% inner_join(metadata %>% rownames_to_column(var="SampleID"))

library(ggpubr)

ojiss %>% filter(order=="q0") %>%  ggboxplot(x="Poligono", y="val", fill="Poligono")+
  theme_bw()

#volviendolo a correr con hillr

library(hillR)
table_ojiss<- read.delim("Data/filt_table_ojis/table_ojis_filt.txt", check.names = F, row.names = 1)
filt_table_ojis<- read_qza("Data/filt_table_ajis.qza")$data

q0_ojiss<- hill_taxa(table_ojiss, MARGIN = 2) %>% as.data.frame() %>% dplyr::rename(val=".") %>%
  rownames_to_column(var="SampleID") %>% 
  inner_join(metadata %>% rownames_to_column(var="SampleID"))

q0_ojiss_qza<- hill_taxa(filt_table_ojis, MARGIN = 2)

q0_ojiss  %>%  ggboxplot(x="Poligono", y="val", fill="Poligono")+
  theme_bw()

