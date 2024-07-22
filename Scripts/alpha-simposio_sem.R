library(qiime2R)
library(tidyverse)

filt_table_ojis<- read_qza("Data/filt_table_ojis.qza")$data %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "[^_]+(?=_............\\w[^_]*$)")) %>% 
  column_to_rownames(var = "ids") %>% t()

meta<- read.delim("Data/its_map.txt")

library(hillR)

hill0<- hill_taxa(comm = filt_table_ojis, q = 0, MARGIN = 2)
hill1<- hill_taxa(comm = filt_table_ojis, q = 1, MARGIN = 2)
hill2<- hill_taxa(comm = filt_table_ojis, q = 2, MARGIN = 2)


q0<- as.data.frame(hill0) %>% #colMeans() %>% 
  round() %>% as.data.frame() %>% dplyr::rename(q="hill0") %>% 
  mutate(order="q0") %>% rownames_to_column(var = "SampleID")
q1<- as.data.frame(hill1) %>% #colMeans() %>% 
  round() %>% as.data.frame()%>% dplyr::rename(q="hill1") %>% 
  mutate(order="q1")%>% rownames_to_column(var = "SampleID")
q2<- as.data.frame(hill2) %>% #colMeans() %>% 
  round() %>% as.data.frame() %>% dplyr::rename(q="hill2") %>% 
  mutate(order="q2")%>% rownames_to_column(var = "SampleID")

qs<- q0 %>% full_join(q1) %>% full_join(q2) %>% pivot_longer(
  cols = c(-SampleID,-order), values_to = "val", names_to = "ids") %>% inner_join(meta)



#library(ggpubr)

#fig <- plot_ly(ggplot2::diamonds, x = ~cut, y = ~price, color = ~clarity, type = "box")
#fig <- fig %>% layout(boxmode = "group")

#fig
library(ggpubr)

qs %>%filter(order=="q2") %>%  ggboxplot(x = "Season", y = "val", fill="Season", facet.by = "Poligono")+
  facet_wrap( vars(Poligono), scales = "free")+theme_linedraw()+
  theme(axis.text.x = element_blank() ,
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ylab("Número efectivo de especies")+stat_compare_means(label = "p.signif")+
  scale_fill_manual(values = c("#FE922AFF","#4490FEFF"),labels=c("Secas", "Lluvias"))+
  labs(fill="Temporada")+
  theme(strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16))

a
qs %>% ggboxplot(x = "Season", y = "val", fill="Season", facet.by = "Poligono")+
  facet_grid(vars(order), vars(Poligono), scales = "free")+theme_bw()+
  theme(axis.text.x = element_blank() ,
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  ylab("Effective number of features")+stat_compare_means()

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

ggsave('Plot/alpha.png', width = 9, height = 7, dpi = 300, plot =a)
