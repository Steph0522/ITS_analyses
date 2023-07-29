library(qiime2R)
library(tidyverse)
set.seed(124)

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

#taxonomys
taxa_oc<- read_qza("Data/taxonomy_oc.qza")$data
taxa_oim<- read_qza("Data/taxonomy_oim.qza")$data
taxa_oism<- read_qza("Data/taxonomy_oism.qza")$data
taxa_ojis<- read_qza("Data/taxonomy_ojis.qza")$data
taxa_oiss<- read_qza("Data/taxonomy_oiss.qza")$data
taxa_ac<- read_qza("Data/taxonomy_sklearn_ac.qza")$data
taxa_aism<- read_qza("Data/taxonomy_sklearn_aism.qza")$data
taxa_aisc<- read_qza("Data/taxonomy_sklearn_aisc.qza")$data
taxa_ajis<- read_qza("Data/taxonomy_ajis.qza")$data


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
list_taxa<- list(taxa_oc, taxa_oim, taxa_oism, taxa_ojis,
                 taxa_oiss, taxa_ac, taxa_aism, taxa_aisc, taxa_ajis)


names(list_taxa)<- c("OC", "OIM", "OISM", "OJIS", "OISS",
                     "AC", "AISM", "AISC", "AJIS")

metadata<- read.delim("Data/its_map.txt")


#functions

transform_clr<- function(x){
  library(ALDEx2)
  set.seed(123)
  aldex.clr.transform <- aldex.clr(x, mc.samples = 128, denom="all",
                                   verbose = FALSE, useMC=FALSE)
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1) )
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1) )
  return(aldex.clr.transform.data)}

pca_compositional<- function(table_transformed){
    otu_pca<- prcomp(table_transformed)}

PC1.f <- function(pca_compositional){paste("PC1 : ", round(pca_compositional$sdev[1]^2/sum(pca_compositional$sdev^2),3)*100, "%",sep="")}
PC2.f <- function(pca_compositional){paste("PC2 : ", round(pca_compositional$sdev[2]^2/sum(pca_compositional$sdev^2),3)*100, "%",sep="")}

pca_plot<- function(pca_compositional, scales, taxonomys, feature){ggplot() +
    geom_segment(data=data.frame(pca_compositional$rotation) %>%   #arrows
                   rownames_to_column(var = "Feature.ID")%>%  
                   mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                   top_n(5, a) %>% #keep 10 furthest away points
                   mutate(PC1=PC1*scales, PC2=PC2*scales),
                 aes(x=0, xend=PC1, y=0, yend=PC2),
                 arrow = arrow(length = unit(0.3,"cm")))+
    geom_point(data=data.frame(pca_compositional$x) %>% #individuals
                 rownames_to_column(var = "SampleID")%>%
                 left_join(metadata, by = "SampleID"),
               aes(x=PC1, y=PC2, fill=Poligono,shape=Season, color=Poligono), size=4) +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo" )+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical") +
    geom_polygon(data=data.frame(pca_compositional$x) %>% #individuals
                   rownames_to_column(var = "SampleID")%>%
                   left_join(metadata, by = "SampleID")%>%
                   drop_na() %>%
                   group_by(Poligono) %>% 
                   slice(chull(PC1, PC2)),
                 aes(x=PC1, y=PC2, fill=Poligono, color=Poligono),
                 alpha = 0.3,
                 show.legend = FALSE)+
    ggrepel::geom_label_repel(data=data.frame(pca_compositional$rotation) %>%   #arrows
                                rownames_to_column(var = "Feature.ID")%>%  
                                mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                                top_n(5, a) %>% #keep 10 furthest away points
                                mutate(PC1=PC1*scales, PC2=PC2*scales)%>%
                                left_join(
                                  taxonomys)%>% dplyr::select(
                                    Taxon, PC1, PC2, Feature.ID)%>%
                                mutate_at(
                                  c("Taxon"), ~str_replace(.,";s__unidentified", "")) %>% mutate(
                                    tax= str_extract(Taxon, "[^;s__]\\w+$")) %>%
                                mutate_at(c("tax"), funs(tax = case_when(
                                  tax=="Fungi" ~ "Unidentified",
                                  tax=="sajor_caju" ~ "Lentinus",
                                  TRUE~as.character(tax)))),
                              aes(x=PC1, y=PC2, label= tax),
                              segment.colour = NA, col = 'black', fill= "#EEEEEE",
                              fontface="italic",  box.padding = 0.2, size=4)}

pca_new<-function(pca, scales, taxonomys, feature){
  metadata1<- as.data.frame(pca$x) %>% rownames_to_column(var = "SampleID") %>% 
    inner_join(metadata)
  y<-ggordiplots::gg_ordiplot(pca, metadata1$Poligono, hull = FALSE, 
                              spiders = TRUE,  ellipse = FALSE,   pt.size = 4,
                              plot =FALSE, label = FALSE)
  
  # Basic ordination plot:
  xlab <- y$plot$labels$x
  ylab <- y$plot$labels$y
  z<-ggplot()+ geom_point(data = y$df_ord %>% rownames_to_column(var="SampleID") %>% 
                            inner_join(metadata1),
                          aes(x = x, y = y, color = Group, shape=Season), size = 3) + xlab(xlab) + 
    ylab(ylab)+
    
    # Plot spiders:
    geom_segment(data = y$df_spiders, aes(x = cntr.x, xend = x, y = cntr.y, yend = y, color = Group), 
                 show.legend = FALSE)+
      #geom_label(
    #data = y$df_mean.ord,
    #aes(x = x, y = y, label=Group), 
  #  label.padding = unit(0.15, "lines"),label.size = 0.4  )+
  guides(
    color=guide_legend(title="Sites"))+theme_linedraw() +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo", name="Poligono")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    geom_text(data=data.frame(pca$rotation) %>%   #arrows
                                rownames_to_column(var = "Feature.ID")%>%  
                                mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                                top_n(5, a) %>% #keep 10 furthest away points
                                mutate(PC1=PC1*scales, PC2=PC2*scales)%>% left_join(
                                  taxonomys)%>% dplyr::select(
                                    Taxon, PC1, PC2, Feature.ID)%>%
                                mutate_at(
                                  c("Taxon"), ~str_replace(
                                    .,";s__unidentified", ""))%>%
                                mutate_at(
                                  c("Taxon"), ~str_replace(
                                    .,";g__unidentified", "")) %>% mutate(
                                      tax= str_extract(Taxon, "[^;s__]\\w+$")) %>%
                                mutate_at(c("tax"), funs(tax = case_when(
                                  tax=="Fungi" ~ "Unidentified",
                                  tax=="pseudograminearum"~"Fusarium",
                                  tax=="oryzae"~ "Aspergillus oryzae",
                                  tax=="oreades"~ "",
                                  tax=="solani"~"Rhizoctonia solani",
                                  tax=="romaleae"~"Encephalitozoon romaleae",
                                  tax=="Pseudogymnoascus verrucosus"~"",
                                  tax=="sajor_caju" ~ "Lentinus",
                                  TRUE~as.character(tax)))),
                              aes(x=PC1, y=PC2, label= tax),
                             # segment.colour = NA,
               col = 'black', fill= "#EEEEEE",
                              fontface="italic", 
               label.r = unit(0.1, "cm"),
               size=4, 
               nudge_y = 0.25,nudge_x = 0.25,
               label.padding = unit(0.05, "cm"),
    )
  print(z)}


permanova_compo<- function(table_transformed){
  tab<-table_transformed %>% as.data.frame() %>% rownames_to_column(
    var = "SampleID") %>% inner_join(metadata)
  library(vegan)
  perm<- how(nperm = 999)
  
  perm<-adonis2(table_transformed~Poligono, data = tab, method = 
                  "euclidian", permutations =perm)
  print(perm)}

permdisp_compo<- function(table_transformed){
  tab<-table_transformed %>% as.data.frame() %>% rownames_to_column(
    var = "SampleID") %>% inner_join(metadata)
  library(vegan)
  perm<- how(nperm = 999)
  mat<- dist(table_transformed, method = "euclidean")
  permdisp<-betadisper(mat, tab$Poligono)
  permdisp2<- permutest( permdisp, permutations = 999)
  print(permdisp2)
}

#applying
table_transformed<- lapply(list_table, transform_clr)
pcas<- lapply(table_transformed, pca_compositional)
pc1<- lapply(pcas, PC1.f)
pc2<- lapply(pcas, PC2.f)
perma<- lapply(table_transformed, permanova_compo)
permd<- lapply(table_transformed, permdisp_compo)  

a<-pca_new(pca =  pcas[[1]], taxonomys = taxa_oc, scales =1200)+theme(legend.position = "none")
b<-pca_new(pca = pcas[[2]], taxonomys = taxa_oim, scales = 900)+theme(legend.position = "none")
c<-pca_new(pca = pcas[[3]], taxonomys = taxa_oism, scales = 1000)+theme(legend.position = "none")
d<-pca_new(pca = pcas[[4]], taxonomys = taxa_ojis, scales = 1800)+theme(legend.position = "none")
e<-pca_new(pca = pcas[[5]], taxonomys = taxa_oiss, scales = 800)+theme(legend.position = "none")
f<-pca_new(pca = pcas[[6]], taxonomys = taxa_ac, scales = 500)+theme(legend.position = "none")
g<-pca_new(pca = pcas[[7]], taxonomys = taxa_aism, scales = 800)+theme(legend.position = "none")
h<-pca_new(pca = pcas[[8]], taxonomys = taxa_aisc, scales = 800)+theme(legend.position = "none")
i<-pca_new(pca = pcas[[9]], taxonomys = taxa_ajis, scales = 800)+theme(legend.position = "none")




library(plotly)
ggplotly(a)

subplot(a,b, c)

tab<-do.call("rbind", lapply(permd, "[[", 1)) %>% 
  drop_na() %>% rownames_to_column(var = "Method") %>% 
  mutate(Method=str_extract(Method,"^\\w+")) %>% 
  dplyr::select(Method,F,`Pr(>F)`) %>% mutate_at("F", ~round(.,2))

tab2<-do.call("rbind", lapply(perma, "[", 3:5)) %>% 
  rownames_to_column(var="Method")%>%
  filter(!str_detect(Method, 'Residual|Total')) %>% 
  mutate(Method=str_extract(Method,"^\\w+")) %>% 
  dplyr::select(Method,R2,F,`Pr(>F)`) %>% mutate_at("F", ~round(.,2))%>% 
  mutate_at("R2", ~round(.,2)) %>% as.data.frame()

library(ggpubr)
library(plotly)
ggtexttable(tab2, theme = ttheme("blank"), rows = NULL) %>% 
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(10), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(.), column.side = "left", from.row = 2, linetype = 2) %>%
  table_cell_font(row = 2:tab_nrow(.), column = 4, face = "bold")

ggtexttable(tab, theme = ttheme("blank"), rows = NULL) %>% 
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(10), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(.), column.side = "left", from.row = 2, linetype = 2) %>%
  table_cell_font(row = c(7,8,9), column = 3, face = "bold")
  
saveRDS(pcas, "Data/pcas.RDS")
write.csv(tab, "tab.csv")
write.csv(tab2, "tab2.csv")
