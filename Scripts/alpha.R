library(qiime2R)
library(tidyverse)

table_oc<-read_qza("Data/table_oc.qza")$data
table_oim<- read_qza("Data/table_oim.qza")$data %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "^\\w..[^_]")) %>% 
  column_to_rownames(var = "ids") %>% t()
table_oism<- read_qza("Data/table_oism.qza")$data  %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "[^_]+(?=_................\\w[^_]*$)")) %>% 
  column_to_rownames(var = "ids") %>% t()
table_ojis<- read_qza("Data/table_ojis.qza")$data %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "[^_]+(?=_............\\w[^_]*$)")) %>% 
  column_to_rownames(var = "ids") %>% t()
table_oiss<- read_qza("Data/table_oiss.qza")$data%>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "[^_]+(?=................\\w[^_]*$)")) %>% 
  column_to_rownames(var = "ids") %>% t()
table_ac<- read_qza("Data/table_ac.qza")$data
table_aism<- read_qza("Data/table_aism.qza")$data
table_aisc<- read_qza("Data/table_aisc.qza")$data
table_ajis<- read_qza("Data/table_ajis.qza")$data

filt_table_oc<-read_qza("Data/filt_table_oc.qza")$data
filt_table_oim<- read_qza("Data/filt_table_oim.qza")$data %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "^\\w..[^_]")) %>% 
  column_to_rownames(var = "ids") %>% t()
filt_table_oism<- read_qza("Data/filt_table_oism.qza")$data %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "[^_]+(?=_................\\w[^_]*$)")) %>% 
  column_to_rownames(var = "ids") %>% t()
filt_table_ojis<- read_qza("Data/filt_table_ojis.qza")$data %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "[^_]+(?=_............\\w[^_]*$)")) %>% 
  column_to_rownames(var = "ids") %>% t()
filt_table_oiss<- read_qza("Data/filt_table_oiss.qza")$data%>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "[^_]+(?=................\\w[^_]*$)")) %>% 
  column_to_rownames(var = "ids") %>% t()
filt_table_ac<- read_qza("Data/filt_table_ac.qza")$data
filt_table_aism<- read_qza("Data/filt_table_aism.qza")$data
filt_table_aisc<- read_qza("Data/filt_table_aisc.qza")$data
filt_table_ajis<- read_qza("Data/filt_table_ajis.qza")$data



table_oim<- table_oim[,match(colnames(table_oc), colnames(table_oim))]
table_oism<- table_oism[,match(colnames(table_oc), colnames(table_oism))]
table_ojis<- table_ojis[,match(colnames(table_oc), colnames(table_ojis))]
table_oiss<- table_oiss[,match(colnames(table_oc), colnames(table_oiss))]

table_ac<- table_ac[,match(colnames(table_oc), colnames(table_ac))]
table_aism <- table_aism[,match(colnames(table_oc), colnames(table_aism))]
table_aisc <- table_aisc[,match(colnames(table_oc), colnames(table_aisc))]
table_ajis <- table_ajis[,match(colnames(table_oc), colnames(table_ajis))]

filt_table_oc <- filt_table_oc[,match(colnames(table_oc), colnames(filt_table_oc))]
filt_table_oim <- filt_table_oim[,match(colnames(table_oc), colnames(filt_table_oim))]
filt_table_oism <- filt_table_oism[,match(colnames(table_oc), colnames(filt_table_oism))]
filt_table_ojis <- filt_table_ojis[,match(colnames(table_oc), colnames(filt_table_ojis))]
filt_table_oiss <- filt_table_oiss[,match(colnames(table_oc), colnames(filt_table_oiss))]

filt_table_ac <- filt_table_ac[,match(colnames(table_oc), colnames(filt_table_ac))]
filt_table_aism <- filt_table_aism[,match(colnames(table_oc), colnames(filt_table_aism))]
filt_table_aisc <- filt_table_aisc[,match(colnames(table_oc), colnames(filt_table_aisc))]
filt_table_ajis <- filt_table_ajis[,match(colnames(table_oc), colnames(filt_table_ajis))]



list_table<- list(table_oc, table_oim, table_oism, table_ojis, table_oiss,
                  table_ac, table_aism, table_aisc, table_ajis)
names(list_table)<- c("OC", "OIM", "OISM", "OJIS", "OISS",
                      "AC", "AISM", "AISC", "AJIS")
list_single<- list(filt_table_oc, filt_table_oim, filt_table_oism, filt_table_ojis, filt_table_oiss,
                   filt_table_ac, filt_table_aism, filt_table_aisc, filt_table_ajis)
names(list_single)<- c("OCS", "OIMS", "OISMS", "OJISS", "OISSS",
                       "ACS", "AISMS", "AISCS", "AJISS")


all<- c(list_table, list_single)

meta<- read.delim("Data/its_map.txt")

library(hillR)

hill_fun0<- function(x){hill_taxa(comm = x, q = 0, MARGIN = 2)}
hill_fun1<- function(x){hill_taxa(comm = x, q = 1, MARGIN = 2)}
hill_fun2<- function(x){hill_taxa(comm = x, q = 2, MARGIN = 2)}

q0s<- lapply(all, hill_fun0)
q1s<- lapply(all, hill_fun1)
q2s<- lapply(all, hill_fun2)

q0<- as.data.frame(q0s) %>% #colMeans() %>% 
  round() %>% as.data.frame() %>%# dplyr::rename(q=".") %>% 
  mutate(order="q0") %>% rownames_to_column(var = "SampleID")
q1<- as.data.frame(q1s) %>% #colMeans() %>% 
  round() %>% as.data.frame() %>%# dplyr::rename(q=".") %>% 
  mutate(order="q1")%>% rownames_to_column(var = "SampleID")
q2<- as.data.frame(q2s) %>% #colMeans() %>% 
  round() %>% as.data.frame() %>% #dplyr::rename(q=".") %>% 
  mutate(order="q2")%>% rownames_to_column(var = "SampleID")

qs<- q0 %>% full_join(q1) %>% full_join(q2) %>% pivot_longer(
  cols = c(-SampleID,-order), values_to = "val", names_to = "ids")


qs$ids<- factor(qs$ids,
  levels =c("OC" ,   "OIM",   "OISM",  "OJIS",  "OISS" ,
  "OCS",   "OIMS"  ,"OISMS", "OJISS", "OISSS",
  "AC"  ,  "AISM" , "AISC",  "AJIS", 
 "ACS"   ,"AISMS", "AISCS" ,"AJISS" ))

#library(ggpubr)

#fig <- plot_ly(ggplot2::diamonds, x = ~cut, y = ~price, color = ~clarity, type = "box")
#fig <- fig %>% layout(boxmode = "group")

#fig

  
qs %>% ggboxplot(x = "ids", y = "val", fill="ids", facet.by = "order")+
  facet_wrap(vars(order), scales = "free")+theme_bw()+
  theme(axis.text.x = element_blank() ,
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ylab("Effective number of features")

qs %>% inner_join(meta)%>% filter(order=="q0") %>% 
  ggbarplot(x = "Poligono", y = "val", fill="Poligono", 
            facet.by = "ids", add="mean_sd")+
  facet_wrap(vars(ids), ncol= 5,scales = "free")+theme_bw()+
  theme(#axis.text.x = element_blank() ,
        axis.title.x = element_blank())+
        #axis.ticks.x = element_blank())+
  ylab("Effective number of features")+theme(legend.position = "none")#+stat_compare_means()

qs %>% inner_join(meta)%>% filter(order=="q0") %>% 
  ggbarplot(x = "Poligono", y = "val", fill="Poligono", 
            facet.by = "ids", add="mean_sd")+
  facet_wrap(vars(ids), ncol= 5,scales = "free")+theme_bw()+
  theme(#axis.text.x = element_blank() ,
        axis.title.x = element_blank())+
        #axis.ticks.x = element_blank())+
  ylab("Effective number of features")+theme(legend.position = "none")#+stat_compare_means()
qs %>% inner_join(meta)%>% filter(order=="q1") %>% 
  ggbarplot(x = "Poligono", y = "val", fill="Poligono", 
            facet.by = "ids", add="mean_sd")+
  facet_wrap(vars(ids), ncol= 5,scales = "free")+theme_bw()+
  theme(#axis.text.x = element_blank() ,
        axis.title.x = element_blank())+
        #axis.ticks.x = element_blank())+
  ylab("Effective number of features")+theme(legend.position = "none")#+stat_compare_means()

write.csv(qs, "Data/qs.csv", row.names = F)
