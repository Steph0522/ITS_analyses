# Functions of beta diversity


relabunda<- function(x){(as.data.frame(t(t(x)/colSums(x)))*100)}


log_norm <- function(otu) {
  log(otu + 1)
}
beta_div_dist_bray <- function(otu, method = "bray") {
  require(phyloseq)
  physeq <- phyloseq(otu_table(otu, taxa_are_rows = TRUE))
  as.matrix(phyloseq::distance(physeq, method))
} 

beta_div_dist_hill<- function(otu, q=NULL){
  require(hillR)
  require(hilldiv)
otu<-otu  
  if (q==0) {
b<- hilldiv::pair_dis(otu, qvalue = 0, metric = "U")
return(b[["L1_UqN"]])} 
if (q==1) {
b<-1-hill_taxa_parti_pairwise(t(otu), q=1,
                              output = "matrix", 
                              pairs = "full",
                              .progress = FALSE, 
                              show_warning = FALSE)$local_similarity
  return(b)}
if (q==2) {
  b<- hilldiv::pair_dis(otu, qvalue = 2, metric = "C")
return(b[["L1_CqN"]])

      }
  
}

pcoa_all <- function(dist) {
  require(ape)
  require(vegan)
  map <- filter(map, SampleID %in% rownames(dist))
  dist <- dist[match(map$SampleID, rownames(dist)), match(map$SampleID, colnames(dist))]
  pcoa <- wcmdscale(as.dist(dist), eig=TRUE)
  return(pcoa)
}

pcoa_axes <- function(dist) {
  require(ape)
  map <- filter(map, SampleID %in% rownames(dist))
  dist <- dist[match(map$SampleID, rownames(dist)), match(map$SampleID, colnames(dist))]
  pcoa <- pcoa(as.dist(dist))
  as.data.frame(pcoa$vectors) %>%
    mutate(SampleID = rownames(.)) %>%
    inner_join(map, by = "SampleID")
}
pcoa_eigval <- function (dist) {
  require(ape)
  map <- filter(map, SampleID %in% rownames(dist))
  dist <- dist[match(map$SampleID, rownames(dist)), match(map$SampleID, colnames(dist))]
  pcoa <- pcoa(as.dist(dist))
  eigval <- round(pcoa$values$Relative_eig * 100, digits = 2)
  data.frame( PC = 1:length(eigval), Eigval = eigval, CumEigval = cumsum(eigval))
}

dist.tidy.filter.bray <-function(dist){ 
  require(reshape2)
  dist[upper.tri(dist)] <- NA 
  
  dist.tidy<-dist %>% melt(., varnames = c(
    "SampleID.x", "SampleID.y")) %>% 
    inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
    inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
    filter(!is.na(value)) %>% dplyr::rename(Distance=value)
  
    dist.filt <- dist.tidy %>% 
    filter(Distance > 0)  %>% 
    full_join(distance_dm) %>% dplyr::rename(SpatialDistance = value)  %>% 
      full_join(env.dist.tidy) %>% 
      full_join(env.dist.tidy2) %>% 
    mutate(Similarity = 1 - Distance)
  return(dist.filt)
}


dist.tidy.filter.hill <-function(dist){ 
  require(reshape2)
  bc.dist[upper.tri(dist)] <- NA 
    dist.tidy<-dist %>% melt(., varnames = c(
    "SampleID.x", "SampleID.y")) %>% 
    inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
    inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
    filter(!is.na(value)) %>% dplyr::rename(Distance=value)
   
    dist.filt <-dist.tidy %>% 
    filter(Distance > 0)  %>% 
    full_join(distance_dm) %>% dplyr::rename(SpatialDistance = value) %>% 
    full_join(env.dist.tidy) %>% 
    full_join(env.dist.tidy2) %>% 
    mutate(Similarity = 1 - Distance)
  return(dist.filt)
}



pcoa_plot<-function(pca){
  require(rcartocolor)
  y<-ggordiplots::gg_ordiplot(pca, map$Sites, hull = FALSE, 
                              spiders = TRUE,  ellipse = FALSE,   pt.size = 4,
                              plot =FALSE, label = FALSE)
  z <- y$plot
  a<-z+geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group), colour="white", fill="black",
    label.padding = unit(0.1, "lines"),label.size = 0.3,fontface = "bold",
  )+
    guides(
    color=guide_legend(title="Sites"))+theme_linedraw() +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    rcartocolor::scale_color_carto_d(palette = "Safe") +
  #  scale_color_viridis_d(option ="turbo" )+#color of points +
    #scale_y_continuous(limits = c(-0.2,0.2))+
    #scale_x_continuous(limits = c(-0.2,0.2))+
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
          #  plot.margin = margin(t = 0,  # Top margin
         #                        r = 0,  # Right margin
           #                      b = 0,  # Bottom margin
            #                       l = 0))# Left margin #+
  #geom_text(data=data.frame(pca$rotation) %>%   #arrows
  #           rownames_to_column(var = "Feature.ID")%>%  
  #          mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
  #         mutate(PC1=PC1*scales, PC2=PC2*scales)%>% left_join(
  #          taxonomy)%>% dplyr::select(
  #           Taxon, PC1, PC2, Feature.ID,a)%>%
  #      mutate_at(
  #       c("Taxon"), ~str_replace(.,";s__unidentified", "")) %>% filter(
  #        Feature.ID %in% vars) %>% separate(Taxon, c(letters[1:7]),sep = separ),  #keep 10 furthest away points
  # aes(x=PC1, y=PC2, label= f),
  #fontface="italic",  box.padding = 0.5, size=4)
}


otu.match <- function(x){x[,match(map$SampleID, colnames(x))]}
otu.single <- function(otu){otu[rowSums(otu>0)>1,]}
otu.norm <- function(otu){otu %>% relabunda() %>% log_norm()}

permanova_beta<- function(matrix, metadata){
tab<-matrix%>% as.data.frame() %>% rownames_to_column(
    var = "SampleID") %>% inner_join(metadata)
  library(vegan)
  perm<- how(nperm = 999)
  
  perm<-adonis2(matrix~Sites, data = tab, permutations =perm)
  }

permdisp_beta<- function(matrix, metadata){
  tab<-matrix%>% as.data.frame() %>% rownames_to_column(
    var = "SampleID") %>% inner_join(metadata)
  library(vegan)
  perm<- how(nperm = 999)
  permdisp<-suppressWarnings(betadisper(as.dist(matrix), tab$Sites))
  permdisp2<- suppressWarnings(permutest( permdisp, permutations = 999))
  

}


transform_clr<- function(x){
  library(ALDEx2)
  set.seed(123)
  aldex.clr.transform <- aldex.clr(x, mc.samples = 999, denom="all",
                                   verbose = FALSE, useMC=FALSE)
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1) )
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1) )
  return(aldex.clr.transform.data)}


pca_compo<- function(x){
  otu_pca<- prcomp(x)}

pca_compositional<- function(x){
  library(ALDEx2)
  set.seed(123)
  aldex.clr.transform <- aldex.clr(x, mc.samples = 999, denom="all",
                                   verbose = FALSE, useMC=FALSE)
  aldex.clr.transform.data<-  t(getMonteCarloSample(aldex.clr.transform,1) )
  otu_pca<- prcomp(aldex.clr.transform.data)}

PC1.f <- function(pcx){paste("PC1 : ", round(pcx$sdev[1]^2/sum(pcx$sdev^2),3)*100, "%",sep="")}
PC2.f <- function(pcx){paste("PC2 : ", round(pcx$sdev[2]^2/sum(pcx$sdev^2),3)*100, "%",sep="")}

pca_plot<- function(tab, scales, taxonomys, feature){ggplot() +
    geom_segment(data=data.frame(tab$rotation) %>%   #arrows
                   rownames_to_column(var = "Feature.ID")%>%  
                   mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                   top_n(5, a) %>% #keep 10 furthest away points
                   mutate(PC1=PC1*scales, PC2=PC2*scales),
                 aes(x=0, xend=PC1, y=0, yend=PC2),
                 arrow = arrow(length = unit(0.3,"cm")))+
    geom_point(data=data.frame(tab$x) %>% #individuals
                 rownames_to_column(var = "SampleID")%>%
                 left_join(metadata, by = "SampleID"),
               aes(x=PC1, y=PC2, fill=Sites),shape=21, size=4) +
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
    geom_polygon(data=data.frame(tab$x) %>% #individuals
                   rownames_to_column(var = "SampleID")%>%
                   left_join(metadata, by = "SampleID")%>%
                   drop_na() %>%
                   group_by(Sites) %>% 
                   slice(chull(PC1, PC2)),
                 aes(x=PC1, y=PC2, fill=Sites, color=Sites),
                 alpha = 0.3,
                 show.legend = FALSE)+
    ggrepel::geom_label_repel(data=data.frame(tab$rotation) %>%   #arrows
                                rownames_to_column(var = "Feature.ID")%>%  
                                mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
                                top_n(5, a) %>% #keep 10 furthest away points
                                mutate(PC1=PC1*scales, PC2=PC2*scales)%>% left_join(
                                  taxonomys)%>% dplyr::select(
                                    Taxon, PC1, PC2, Feature.ID)%>%
                                mutate_at(
                                  c("Taxon"), ~str_replace(.,";s__unidentified", "")) %>% mutate(
                                    tax= str_extract(Taxon, "[^__]+$")) %>%
                                mutate_at(c("tax"), funs(tax = case_when(
                                  tax=="Fungi" ~ "Unidentified",
                                  tax=="sajor_caju" ~ "Lentinus",
                                  TRUE~as.character(tax)))),
                              aes(x=PC1, y=PC2, label= tax),
                              segment.colour = NA, col = 'black', fill= "#EEEEEE",
                              fontface="italic",  box.padding = 0.2, size=4)}



permanova_compo<- function(table_transformed, metadata){
  tab<-table_transformed %>% as.data.frame() %>% rownames_to_column(
    var = "SampleID") %>% inner_join(metadata)
  library(vegan)
  perm<- how(nperm = 999)
  
  perm<-adonis2(table_transformed~Site, data = tab, method = 
                  "euclidian", permutations =perm)
  print(perm)}

permdisp_compo<- function(table_transformed, metadata){
  tab<-table_transformed %>% as.data.frame() %>% rownames_to_column(
    var = "SampleID") %>% inner_join(metadata)
  library(vegan)
  perm<- how(nperm = 999)
  mat<- dist(table_transformed, method = "euclidean")
  permdisp<-betadisper(mat, tab$Site)
  permdisp2<- permutest( permdisp, permutations = 999)
  
  print(permdisp2)
}


pca_compositional_sites<-function(pca){
  metadata1<- as.data.frame(pca$x) %>% rownames_to_column(var = "SampleID") %>% 
    inner_join(metadata)
  y<-ggordiplots::gg_ordiplot(pca, metadata1$Sites, hull = FALSE, 
                              spiders = TRUE,  ellipse = FALSE,   pt.size = 4,
                              plot =FALSE, label = FALSE)
  z <- y$plot
  a<-z+geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group), 
    label.padding = unit(0.15, "lines"),label.size = 0.4,
  )+guides(
    color=guide_legend(title="Sites"))+theme_linedraw() +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo", name="Sites")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())#+
  #geom_text(data=data.frame(pca$rotation) %>%   #arrows
  #           rownames_to_column(var = "Feature.ID")%>%  
  #          mutate(a=sqrt(PC1^2+PC2^2)) %>% # calculate the distance from the origin
  #         mutate(PC1=PC1*scales, PC2=PC2*scales)%>% left_join(
  #          taxonomy)%>% dplyr::select(
  #           Taxon, PC1, PC2, Feature.ID,a)%>%
  #      mutate_at(
  #       c("Taxon"), ~str_replace(.,";s__unidentified", "")) %>% filter(
  #        Feature.ID %in% vars) %>% separate(Taxon, c(letters[1:7]),sep = separ),  #keep 10 furthest away points
  # aes(x=PC1, y=PC2, label= f),
  #fontface="italic",  box.padding = 0.5, size=4)
  print(a)}


pca_new<-function(pca, scales, taxonomys, feature){
  metadata1<- as.data.frame(pca$x) %>% rownames_to_column(var = "SampleID") %>% 
    inner_join(metadata)
  y<-ggordiplots::gg_ordiplot(pca, metadata1$Sites, hull = FALSE, 
                              spiders = TRUE,  ellipse = FALSE,   pt.size = 4,
                              plot =FALSE, label = FALSE)
  z <- y$plot
  a<-z+geom_label(
    data = y$df_mean.ord,
    aes(x = x, y = y, label=Group), 
    label.padding = unit(0.15, "lines"),label.size = 0.4,
  )+guides(
    color=guide_legend(title="Sites"))+theme_linedraw() +
    geom_vline(xintercept = 0, linetype = 2) +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2) +
    theme_linedraw()+
    scale_fill_viridis_d(option ="turbo", name="Sites")+#color of points 
    scale_color_viridis_d(option ="turbo" )+#color of points 
    theme(axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12), 
          legend.position = "right", 
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    ggrepel::geom_label_repel(data=data.frame(pca$rotation) %>%   #arrows
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
                                      tax= str_extract(Taxon, "[^__]+$")) %>%
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
                              segment.colour = NA, col = 'black', fill= "#EEEEEE",
                              fontface="italic",  box.padding = 0.1, size=4, label.padding = 0.1)
  print(a)}


#bc.dist <-funct beta_div_dist(otu.norm)
cor.b <- function(x){
  cor.test(x$SpatialDistance,x$Similarity, method= "pearson", alternative = "two.sided") %>% tidy()
}
lm.b <- function(x){lm(Similarity ~ SpatialDistance, data = x) %>% tidy() %>% filter(term == "SpatialDistance")}

mantel.b = function(xdis,ydis){
  mantel(xdis, ydis, method="spearman", permutations=999, strata = NULL,
         na.rm = FALSE, parallel = getOption("mc.cores"))
  
}

bc.dist.tidy.filter <-function(bc.dist){ 
  require(reshape2)
  bc.dist[upper.tri(bc.dist)] <- NA 
  bc.dist.tidy<-bc.dist %>% melt(., varnames = c(
    "SampleID.x", "SampleID.y")) %>% 
    inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
    inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
    filter(!is.na(value)) %>% dplyr::rename(Distance=value)
  
  b.dist.filt <- bc.dist.tidy %>% 
    filter(Distance > 0)  %>% 
    full_join(distance_dm) %>% dplyr::rename(SpatialDistance = value) %>% 
    mutate(Similarity = 1 - Distance)
  return(b.dist.filt)
}


bc.dist.tidy.filter.hill <-function(bc.dist){ 
  require(reshape2)
  bc.dist[upper.tri(bc.dist)] <- NA 
  bc.dist.tidy<-bc.dist %>% melt(., varnames = c(
    "SampleID.x", "SampleID.y")) %>% 
    inner_join(map, by = c("SampleID.x" = "SampleID")) %>% 
    inner_join(map, by = c("SampleID.y" = "SampleID")) %>% 
    filter(!is.na(value)) %>% dplyr::rename(Distance=value)
  
  b.dist.filt <- bc.dist.tidy %>% 
    filter(Distance > 0)  %>% 
    full_join(distance_dm) %>% dplyr::rename(SpatialDistance = value) %>% 
    mutate(Similarity = 1 - Distance)
  return(b.dist.filt)
}