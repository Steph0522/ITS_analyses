#beta
library(tidyverse)
library(betapart)
library(qiime2R)
library(vegan)
library(colorRamps)

table<- read_qza("Data/table_SMOQ_rar.qza")$data

metadata<- read.delim("Data/its_map.txt")


#convrt data
presabs<- function(x){x[x>0]=1 
return(t(x))}
tables_presabs<- presabs(table)
table_presabs<- tables_presabs %>% as.data.frame() %>% rownames_to_column(var = "SampleID") %>% 
  arrange(SampleID) %>% column_to_rownames(var = "SampleID")

#betapart function with jaccard 
x<- beta.pair(table_presabs, index.family = "jaccard")

#extract each part
jac<- x$beta.jac
jtu<- x$beta.jtu
jne<- x$beta.jne

jac<- x$beta.sor
jtu<- x$beta.sim
jne<- x$beta.sne

env1<- table_presabs%>% as.data.frame() %>% rownames_to_column(var="SampleID") %>% inner_join(metadata)


jacs<- betadisper(jac,factor(env1$Poligono))
jtus<- betadisper(jtu,factor(env1$Poligono))
jnes<- betadisper(jne,factor(env1$Poligono))

 

library(tidyverse)
library(ggordiplots)

function_plot_beta<- function(x,env){
  y <- gg_ordiplot(x, groups = env$Poligono, hull = FALSE, spiders = TRUE, 
                   ellipse = FALSE, plot = FALSE, label = TRUE)
  xlabs <- y$plot$labels$x
  ylabs <- y$plot$labels$y
  z<-ggplot()+ geom_point(data = y$df_ord %>% rownames_to_column(var="SampleID") %>% 
                            inner_join(env),
                          aes(x = x, y = y, color = Group, shape=Season), size = 3) + 
    xlab(xlabs) +  ylab(ylabs)+
    geom_segment(data = y$df_spiders, aes(x = cntr.x, xend = x, y = cntr.y, yend = y, color = Group), 
                 show.legend = FALSE)
  
  a<-z+geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group))+guides(
      color=guide_legend(title="Polygon"))+theme_linedraw() +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo", name="Polygon")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(aspect.ratio =3/10)
  return(a)
}
 


#jaccard
plot_jac<-  function_plot_beta(jacs, env1)
plot_jac=plot_jac+ guides(colour = guide_legend(nrow = 1, title = "Polygon"),
                          shape = guide_legend(nrow = 1) )+theme(legend.position = "top")

#turnover
plot_turn<- function_plot_beta(jtus, env1)

#nestedness

plot_nes<-function_plot_beta(jnes, env1)


library(cowplot)
leg<- get_legend(plot_jac)
a<-plot_grid(plot_jac+theme(legend.position = "none")+ylab("DIM2")+xlab("")+theme(aspect.ratio =8/10)+ggtitle("Jaccard dissimilarity"),  
             plot_turn+theme(legend.position = "none")+ylab("DIM2")+xlab("")+theme(aspect.ratio =8/10)+ggtitle("Turnover component of Jaccard"),  
             plot_nes+theme(legend.position = "none")+ylab("DIM2")+xlab("DIM1")+theme(aspect.ratio =8/10)+ggtitle("Nestedness component of Jaccard"),  
             ncol = 1, align = "hv", labels = c("A", "B", "C"), label_x = 0.3)
b= plot_grid(leg,a, ncol = 1, rel_heights = c(0.2,1))


b
ggsave("Plots/gao_beta_jaccard.jpg",width = 5, height =10, dpi = 300, plot = b, device = "jpg")


#betapart season

nam<- colnames(table)
nam2<- str_extract(nam, pattern = "^\\d+")
nam3<- unique(nam2)

drys<- paste0(nam3,"D")
rainys<- paste0(nam3,"R")

tables_dry<- table_presabs %>% as.data.frame() %>%rownames_to_column(var="SampleID") %>% 
  filter(SampleID %in% drys) %>% mutate(SampleID = str_extract(SampleID, "^\\d+")) %>% 
  column_to_rownames(var = "SampleID") 

tables_rainy<- table_presabs  %>% as.data.frame() %>%rownames_to_column(var="SampleID") %>% 
  filter(SampleID %in% rainys) %>% mutate(SampleID = str_extract(SampleID, "^\\d+")) %>% 
  column_to_rownames(var = "SampleID") 


beta.t <-beta.temp(tables_dry, tables_rainy,index.family="jaccard")


nams<-paste0(rownames(beta.t), "D")


table_jac<-rbind(beta.t[[3]])
colnames(table_jac)<-nams

table_turn<-rbind(beta.t[[1]]) 
colnames(table_turn)<-nams

table_nes<-rbind(beta.t[[2]])
colnames(table_nes)<-nams


table_jacs<- table_jac %>%t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% inner_join(metadata) %>% 
  dplyr::rename(beta=V1)

table_turns<- table_turn %>%t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% inner_join(metadata) %>% 
  dplyr::rename(beta=V1)

table_ness<- table_nes %>%t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% inner_join(metadata) %>% 
  dplyr::rename(beta=V1)

model<-aov(beta~Poligono, data=table_turns)
out <- agricolae::HSD.test(model,"Poligono", group=TRUE,console=FALSE)
g<- data.frame(out$groups) %>% rownames_to_column(var = "Poligono") %>% dplyr::select(-beta)
data<- data.frame(Poligono=c("P1", "P2", "P3", "P4", "P5", "P6"),
                  beta=rep(0.85,6)) %>% inner_join(g, by = "Poligono")
library(ggpubr)
t<-ggboxplot(data = table_turns, x = "Poligono", y="beta", fill="Poligono")+
  scale_fill_viridis_d(option ="turbo", name="Polygon")+#color of points 
  theme_grey() +theme(axis.text = element_text(colour = "black", size = 14),
                     axis.title = element_text(colour = "black", size = 16),
                     legend.text = element_text(size = 10),
                     legend.title = element_text(size = 12), 
                     legend.position = "right", 
                     legend.box = "vertical",
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=1,))+
  ylab("Temporal turnonver of Jaccard dissimilarity")+
  xlab("Polygon")+
  stat_compare_means(method = "anova", label.y = 0.86, size=5)+
  geom_label(data = data,aes(Poligono, beta, label=groups), size=5)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=4, color="black", fill="white")
  
t
#ggsave('Plots/temporal_turnover_jaccard.png', width =8, height = 6, dpi = 300, plot =t)





#convrt data
presabs<- function(x){x[x>0]=1 
return(t(x))}
tables_presabs<- presabs(table)
table_presabs<- tables_presabs %>% as.data.frame() %>% rownames_to_column(var = "SampleID") %>% 
  arrange(SampleID) %>% column_to_rownames(var = "SampleID")
```


#betapart season

nam<- colnames(table)
nam2<- str_extract(nam, pattern = "^\\d+")
nam3<- unique(nam2)

drys<- paste0(nam3,"D")
rainys<- paste0(nam3,"R")

tables_dry<- table_presabs %>% as.data.frame() %>%rownames_to_column(var="SampleID") %>% 
  filter(SampleID %in% drys) %>% mutate(SampleID = str_extract(SampleID, "^\\d+")) %>% 
  column_to_rownames(var = "SampleID") 

tables_rainy<- table_presabs  %>% as.data.frame() %>%rownames_to_column(var="SampleID") %>% 
  filter(SampleID %in% rainys) %>% mutate(SampleID = str_extract(SampleID, "^\\d+")) %>% 
  column_to_rownames(var = "SampleID") 


beta.t <-beta.temp(tables_dry, tables_rainy,index.family="jaccard")


nams<-paste0(rownames(beta.t), "D")


table_jac<-rbind(beta.t[[3]])
colnames(table_jac)<-nams

table_turn<-rbind(beta.t[[1]]) 
colnames(table_turn)<-nams

table_nes<-rbind(beta.t[[2]])
colnames(table_nes)<-nams





#write.csv(table_turns, "Data/table_turns.csv")
#write.csv(table_jacs, "Data/table_jacs.csv")
#write.csv(table_nes, "Data/table_nes.csv")

pairdis<- pair_dis(table, qvalue = 0)

library(reshape2)
tab<-reshape2::melt(as.matrix(pairdis$L1_SqN), varnames = c("row", "col")) %>% drop_na()%>%
  dplyr::rename("SampleID"="row") %>% inner_join(metadata) %>% dplyr::select(
    row=SampleID, SampleID=col, value, Pol1=Poligono ) %>% inner_join(metadata) %>% 
  dplyr::select(row, col=SampleID, value, Pol1, Pol2=Poligono) %>% 
  unite("compar", c("Pol1", "Pol2"), remove = F, sep = "_vs_") %>% filter(
    compar %in% c("P2_vs_P1","P3_vs_P2","P4_vs_P3","P5_vs_P4","P6_vs_P5",
                  "P3_vs_P1" ,"P4_vs_P2","P5_vs_P3","P6_vs_P4",
                  "P4_vs_P1" ,"P5_vs_P2","P6_vs_P3",
                  "P5_vs_P1", "P6_vs_P2", "P6_vs_P1"))


tab %>% ggpubr::ggboxplot(x="compar", y="value", fill="compar")