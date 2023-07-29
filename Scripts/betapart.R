#beta
library(tidyverse)
library(betapart)
library(qiime2R)
library(vegan)
library(colorRamps)

table_oc<-read_qza("Data/filt_table_oc.qza")$data
table_oim<- read_qza("Data/filt_table_oim.qza")$data %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "^\\w..[^_]")) %>% 
  column_to_rownames(var = "ids") %>% t()
table_oism<- read_qza("Data/filt_table_oism.qza")$data  %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "[^_]+(?=_................\\w[^_]*$)")) %>% 
  column_to_rownames(var = "ids") %>% t()
table_ojis<- read_qza("Data/filt_table_ojis.qza")$data %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "[^_]+(?=_............\\w[^_]*$)")) %>% 
  column_to_rownames(var = "ids") %>% t()
table_oiss<- read_qza("Data/filt_table_oiss.qza")$data%>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = "ids") %>% 
  mutate(ids=str_extract(ids, "[^_]+(?=................\\w[^_]*$)")) %>% 
  column_to_rownames(var = "ids") %>% t()
table_ac<- read_qza("Data/filt_table_ac.qza")$data
table_aism<- read_qza("Data/filt_table_aism.qza")$data
table_aisc<- read_qza("Data/filt_table_aisc.qza")$data
table_ajis<- read_qza("Data/filt_table_ajis.qza")$data



table_oim<- table_oim[,match(colnames(table_oc), colnames(table_oim))]
table_oism<- table_oism[,match(colnames(table_oc), colnames(table_oism))]
table_ojis<- table_ojis[,match(colnames(table_oc), colnames(table_ojis))]
table_oiss<- table_oiss[,match(colnames(table_oc), colnames(table_oiss))]

table_ac<- table_ac[,match(colnames(table_oc), colnames(table_ac))]
table_aism <- table_aism[,match(colnames(table_oc), colnames(table_aism))]
table_aisc <- table_aisc[,match(colnames(table_oc), colnames(table_aisc))]
table_ajis <- table_ajis[,match(colnames(table_oc), colnames(table_ajis))]

list_table<- list(table_oc, table_oim, table_oism, table_ojis, table_oiss,
                  table_ac, table_aism, table_aisc, table_ajis)
names(list_table)<- c("OC", "OIM", "OISM", "OJIS", "OISS",
                      "AC", "AISM", "AISC", "AJIS")
metadata<- read.delim("Data/its_map.txt")

#convrt data
presabs<- function(x){x[x>0]=1 
return(t(x))}
tables_presabs<- lapply(list_table, presabs)

#betapart function with jaccard 
fds<- lapply(tables_presabs, beta.pair,index.family = "jaccard")

#extract each part
jac<- function(x){x$beta.jac}
jacs<- lapply(fds, jac)

jtu<- function(x){x$beta.jtu}
jtus<- lapply(fds, jtu)

jne<- function(x){x$beta.jne}
jnes<- lapply(fds, jne)

env1<- tables_presabs[[1]] %>% as.data.frame() %>% rownames_to_column(var="SampleID") %>% inner_join(metadata)
env2<- tables_presabs[[2]] %>% as.data.frame() %>% rownames_to_column(var="SampleID")%>% inner_join(metadata)
env3<- tables_presabs[[3]] %>% as.data.frame() %>% rownames_to_column(var="SampleID")%>% inner_join(metadata)
env4<- tables_presabs[[4]] %>% as.data.frame() %>% rownames_to_column(var="SampleID")%>% inner_join(metadata)
env5<- tables_presabs[[5]] %>% as.data.frame() %>% rownames_to_column(var="SampleID")%>% inner_join(metadata)
env6<- tables_presabs[[6]] %>% as.data.frame() %>% rownames_to_column(var="SampleID")%>% inner_join(metadata)
env7<- tables_presabs[[7]] %>% as.data.frame() %>% rownames_to_column(var="SampleID")%>% inner_join(metadata)
env8<- tables_presabs[[8]] %>% as.data.frame() %>% rownames_to_column(var="SampleID")%>% inner_join(metadata)
env9<- tables_presabs[[9]] %>% as.data.frame() %>% rownames_to_column(var="SampleID")%>% inner_join(metadata)


envs<- list(env1, env2,env3,env4,env5,env6,env7,env8,env9)

jd<- function(x,y){
  betadisper(x,factor(y$Poligono)) 
  }
jds<- mapply(jd, jacs, envs, SIMPLIFY = F)
tud<- mapply(jd, jtus, envs,SIMPLIFY = F)
ned<- mapply(jd, jnes, envs, SIMPLIFY = F)

  
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
plot_jac<-  mapply(function_plot_beta, jds, envs, SIMPLIFY = F)

#turnover
plot_turn<- mapply(function_plot_beta, tud, envs, SIMPLIFY = F)

#nestedness

plot_nes<- mapply(function_plot_beta, ned, envs, SIMPLIFY = F)

saveRDS(plot_jac, "Data/plot_jac.RDS")
saveRDS(plot_turn, "Data/plot_turn.RDS")
saveRDS(plot_nes, "Data/plot_nes.RDS")


#library(cowplot)
#leg<- get_legend(plot_jac[[1]])
#a<-plot_grid(plot_jac[[1]]+theme(legend.position = "none"), plot_turn[[1]]+theme(legend.position = "none"), plot_nes[[1]]+theme(legend.position = "none"), ncol = 3, align = "hv")
#b<-plot_grid(plot_jac[[2]]+theme(legend.position = "none"), plot_turn[[2]]+theme(legend.position = "none"), plot_nes[[2]]+theme(legend.position = "none"), ncol = 3, align = "hv")
#c<-plot_grid(plot_jac[[3]]+theme(legend.position = "none"), plot_turn[[3]]+theme(legend.position = "none"), plot_nes[[3]]+theme(legend.position = "none"), ncol = 3, align = "hv")
#d<-plot_grid(plot_jac[[4]]+theme(legend.position = "none"), plot_turn[[4]]+theme(legend.position = "none"), plot_nes[[4]]+theme(legend.position = "none"), ncol = 3, align = "hv")
#e<-plot_grid(plot_jac[[5]]+theme(legend.position = "none"), plot_turn[[5]]+theme(legend.position = "none"), plot_nes[[5]]+theme(legend.position = "none"), ncol = 3, align = "hv")
#f<-plot_grid(plot_jac[[6]]+theme(legend.position = "none"), plot_turn[[6]]+theme(legend.position = "none"), plot_nes[[6]]+theme(legend.position = "none"), ncol = 3, align = "hv")
#g<-plot_grid(plot_jac[[7]]+theme(legend.position = "none"), plot_turn[[7]]+theme(legend.position = "none"), plot_nes[[7]]+theme(legend.position = "none"), ncol = 3, align = "hv")
#h<-plot_grid(plot_jac[[8]]+theme(legend.position = "none"), plot_turn[[8]]+theme(legend.position = "none"), plot_nes[[8]]+theme(legend.position = "none"), ncol = 3, align = "hv")
#i<-plot_grid(plot_jac[[9]]+theme(legend.position = "none"), plot_turn[[9]]+theme(legend.position = "none"), plot_nes[[9]]+theme(legend.position = "none"), ncol = 3, align = "hv")

#plot_grid(a,b,c,d,e,f,g,h,i, nrow = 9)

library(plotly)
subplot(plot_jac[[1]]+theme(legend.position = "none"), 
        plot_turn[[1]]+theme(legend.position = "none"),
        plot_nes[[1]]+theme(legend.position = "none"))

subplot(plot_jac[[2]]+theme(legend.position = "none"), 
        plot_turn[[2]]+theme(legend.position = "none"),
        plot_nes[[2]]+theme(legend.position = "none"))

subplot(plot_jac[[3]]+theme(legend.position = "none"), 
        plot_turn[[3]]+theme(legend.position = "none"),
        plot_nes[[3]]+theme(legend.position = "none"))

subplot(plot_jac[[4]]+theme(legend.position = "none"), 
        plot_turn[[4]]+theme(legend.position = "none"),
        plot_nes[[4]]+theme(legend.position = "none"))

subplot(plot_jac[[5]]+theme(legend.position = "none"), 
        plot_turn[[5]]+theme(legend.position = "none"),
        plot_nes[[5]]+theme(legend.position = "none"))

subplot(plot_jac[[6]]+theme(legend.position = "none"), 
        plot_turn[[6]]+theme(legend.position = "none"),
        plot_nes[[6]]+theme(legend.position = "none"))

subplot(plot_jac[[7]]+theme(legend.position = "none"), 
        plot_turn[[7]]+theme(legend.position = "none"),
        plot_nes[[7]]+theme(legend.position = "none"))

subplot(plot_jac[[8]]+theme(legend.position = "none"), 
        plot_turn[[8]]+theme(legend.position = "none"),
        plot_nes[[8]]+theme(legend.position = "none"))

subplot(plot_jac[[9]]+theme(legend.position = "none"), 
        plot_turn[[9]]+theme(legend.position = "none"),
        plot_nes[[9]]+theme(legend.position = "none"))

#ggsave("gao_beta_jaccard.pdf",width = 16, height = 10, dpi = 300, plot = e, device = "pdf")


#betapart season

nam<- colnames(table_ac)
nam2<- str_extract(nam, pattern = "^\\d+")
nam3<- unique(nam2)

drys<- paste0(nam3,"D")
rainys<- paste0(nam3,"R")

dry_func<- function(x){y<- x  %>% as.data.frame() %>%rownames_to_column(var="SampleID") %>% 
  filter(SampleID %in% drys) %>% mutate(SampleID = str_extract(SampleID, "^\\d+")) %>% 
  column_to_rownames(var = "SampleID") 
return(y)}

rainy_func<- function(x){y<- x  %>% as.data.frame() %>%rownames_to_column(var="SampleID") %>% 
  filter(SampleID %in% rainys) %>% mutate(SampleID = str_extract(SampleID, "^\\d+")) %>% 
  column_to_rownames(var = "SampleID") 
return(y)}

tables_dry<- lapply(tables_presabs, dry_func)
tables_rainy<- lapply(tables_presabs, rainy_func)

beta.t.func <- function(x,y){beta.temp(x, y,index.family="jaccard")}

beta.t<- mapply(beta.t.func, tables_dry, tables_rainy, SIMPLIFY = F)

nams<-paste0(rownames(beta.t[[1]]), "D")

table_jac<-do.call("rbind", lapply(beta.t, "[[", 3)) 
colnames(table_jac)<-nams

table_turn<-do.call("rbind", lapply(beta.t, "[[", 1)) 
colnames(table_turn)<-nams

table_nes<-do.call("rbind", lapply(beta.t, "[[", 2)) 
colnames(table_nes)<-nams

table_jacs<- table_jac %>%t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% inner_join(metadata) %>% 
  dplyr::select(-overhang:-ReversePrimerSequence) %>% 
  pivot_longer(cols = OC:AJIS,names_to = "Method", values_to = "beta")

table_turns<- table_turn %>%t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% inner_join(metadata) %>% 
  dplyr::select(-overhang:-ReversePrimerSequence) %>% 
  pivot_longer(cols = OC:AJIS,names_to = "Method", values_to = "beta")

table_jacs<- table_jac %>%t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% inner_join(metadata) %>% 
  dplyr::select(-overhang:-ReversePrimerSequence) %>% 
  pivot_longer(cols = OC:AJIS,names_to = "Method", values_to = "beta")

table_ness<- table_nes %>%t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% inner_join(metadata) %>% 
  dplyr::select(-overhang:-ReversePrimerSequence) %>% 
  pivot_longer(cols = OC:AJIS,names_to = "Method", values_to = "beta")


fig <- plot_ly(table_jacs , x = ~Method, y = ~beta, color = ~Poligono, type = "box")
fig <- fig %>% layout(boxmode = "group")

fig

write.csv(table_turns, "Data/table_turns.csv")
write.csv(table_jacs, "Data/table_jacs.csv")
write.csv(table_nes, "Data/table_nes.csv")



