#load data

library(tidyverse)

tab <- read.delim("Data/filt_table_ojis/table_ojis_filt.txt", check.names = F) %>% dplyr::rename(Feature.ID=ASV)
taxa<- qiime2R::read_qza("Data/taxonomy_ojis.qza")$data %>% as.data.frame()

tab_tax<- tab %>% inner_join(taxa)

tab_taxa<- tab_tax %>% dplyr::select(OTUID=Feature.ID, `111D`:`623R`, taxonomy=Taxon)

#guilds

guild1<- read.delim("~/Documents/FUNGuild/table_ojis_funguild.taxa.guilds.txt", check.names = F) %>% 
  filter(!trophicMode=="na")
guild<- guild1 %>% group_by(Genus, guild) %>% dplyr::count()
trophic<- guild1 %>% group_by(Genus, trophicMode) %>% dplyr::count()


relabunda<- function(x){(as.data.frame(t(t(x)/colSums(x)))*100)}

metadata<- read.delim("Data/its_map.txt")


vector_order<- c("111D" , "112D", "113D" ,"121D", "122D" ,"123D", "111R" ,"112R", "113R", "121R", "122R" ,"123R",
                 "211D" ,"212D" ,"213D" ,"221D" ,"222D" ,"223D", "211R" ,"212R" ,"213R", "221R", "222R" ,"223R",
                 "311D" ,"312D" ,"313D" ,"321D", "323D" ,"322D","311R" ,"312R" ,"313R" ,"321R", "322R" ,"323R" ,
                 "411D" ,"412D" ,"413D" ,"421D" ,"422D" ,"423D", "411R", "412R" ,"413R" ,"421R" ,"422R" ,"423R",
                 "511D" ,"512D" ,"513D" ,"521D" ,"522D", "523D", "511R", "512R" ,"513R" ,"521R", "522R" ,"523R",
                 "611D" ,"612D" ,"613D", "621D" ,"622D", "623D",  "611R", "612R","613R", "621R" ,"622R", "623R")

#ordering

library(ggh4x)
taxones_color<- read_csv("~/Documents/2023/corredor_scripts/Fungal_Communities_PNIP_PNML/Data/taxones_color.csv") %>% dplyr::rename(Genus=Taxon)

table_guilds<- tab_taxa %>%  inner_join(guild1, by = c("OTUID"="OTU")) %>%  mutate_if(
      is.character, ~replace_na(., "Unassigned")) %>%
  filter(!trophicMode=="na") %>% 
  group_by(guild) %>% summarise_if(is.numeric, sum) %>% 
  column_to_rownames(var = "guild") %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% 
  inner_join(metadata) %>% group_by(Poligono, Season) %>% 
  summarise_if(is.numeric, sum) %>% 
  dplyr::select(-LIB) %>% 
  unite("ids", Poligono:Season, sep = ".") %>% column_to_rownames(var = "ids") %>% 
  t() %>% as.data.frame() %>% mutate(
            all= rowSums(.)) %>% dplyr::arrange(
              -all) %>% relabunda(.) %>%as.data.frame( ) %>%  rownames_to_column(
                var = "Guild")%>% filter(!Guild=="unidentified" ,
                                         !Guild=="Unassigned") %>% slice(
                                           c(1:20))  %>% pivot_longer(
                                             ., cols = -Guild, names_to ="SampleID", 
                                             values_to = "relab" ) %>% filter(
                                               !SampleID=="all")
  

table_trophic<- tab_taxa %>%  inner_join(guild1, by = c("OTUID"="OTU")) %>%  mutate_if(
  is.character, ~replace_na(., "Unassigned")) %>%
  filter(!trophicMode=="na") %>% 
  group_by(trophicMode) %>% summarise_if(is.numeric, sum) %>% 
  column_to_rownames(var = "trophicMode") %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SampleID") %>% 
  inner_join(metadata) %>% group_by(Poligono, Season) %>% 
  summarise_if(is.numeric, sum) %>% 
  dplyr::select(-LIB) %>% 
    unite("ids", Poligono:Season, sep = ".") %>% column_to_rownames(var = "ids") %>% 
  t() %>% as.data.frame() %>% mutate(
    all= rowSums(.)) %>% dplyr::arrange(
      -all) %>% relabunda(.) %>%as.data.frame( ) %>%  rownames_to_column(
        var = "Guild")%>% filter(!Guild=="unidentified" ,
                                 !Guild=="Unassigned") %>% slice(
                                   c(1:30))  %>% pivot_longer(
                                     ., cols = -Guild, names_to ="SampleID", 
                                     values_to = "relab" ) %>% filter(
                                       !SampleID=="all")

barplot_guilds<- table_guilds %>% dplyr::select(
    Guild, SampleID, relab) %>% pivot_wider(
      names_from = Guild, values_from = relab) %>% column_to_rownames(
        var = "SampleID") %>% t() %>% as.data.frame()


barplot_trophic<- table_trophic%>% dplyr::select(
  Guild, SampleID, relab) %>% pivot_wider(
    names_from = Guild, values_from = relab) %>% column_to_rownames(
      var = "SampleID") %>% t() %>% as.data.frame()



data_fungi_order=barplot_guilds%>%
rownames_to_column(var="Guild")


data_fungi_order_t=barplot_trophic%>%
  rownames_to_column(var="Guild")

merge_data<- data_fungi_order %>% mutate(Guild=case_when(
        Guild=="Leptosphaeria maculans species complex" ~ "Leptosphaeria", 
        TRUE ~ as.character(Guild)))%>% column_to_rownames(
          var = "Guild") %>% mutate(proms=rowMeans(.)) %>% dplyr::arrange(-proms) %>% slice(
            1:60)%>% dplyr::select(-proms) %>% mutate_all(., funs(R = case_when(
              . <= 0.001 ~ 0,
              . >  0.001 & .  <= 0.005 ~ 1,
              . >  0.005 & .  <= 0.01 ~ 2,
              . >  0.01 & .  <= 0.10 ~ 3,
              . >  0.10 & .  <= 0.20 ~ 4,
              . >  0.20 & .  <= 1.00 ~ 5,
              . >  1.00 & .  <= 2.00 ~ 6,
              . >  2.00 & .  <= 5.00 ~ 7,
              . >  5.00 & .  <= 10.00 ~ 8,
              . >  10.00 & .  <= 25.00 ~ 9,
              . >  25.00 & .  <= 50.00 ~ 10,
              . >  50.00 & .  <= 75.00 ~ 11,
              . >  75.00 ~ 12))) %>%select_at(
                vars(contains("_R"))) %>% select_all(~str_replace(., "_R", ""))

merge_data2<- data_fungi_order_t %>% mutate(Guild=case_when(
  Guild=="Leptosphaeria maculans species complex" ~ "Leptosphaeria", 
  TRUE ~ as.character(Guild)))%>% column_to_rownames(
    var = "Guild") %>% mutate(proms=rowMeans(.)) %>% dplyr::arrange(-proms) %>% slice(
      1:60)%>% dplyr::select(-proms) %>% mutate_all(., funs(R = case_when(
        . <= 0.001 ~ 0,
        . >  0.001 & .  <= 0.005 ~ 1,
        . >  0.005 & .  <= 0.01 ~ 2,
        . >  0.01 & .  <= 0.10 ~ 3,
        . >  0.10 & .  <= 0.20 ~ 4,
        . >  0.20 & .  <= 1.00 ~ 5,
        . >  1.00 & .  <= 2.00 ~ 6,
        . >  2.00 & .  <= 5.00 ~ 7,
        . >  5.00 & .  <= 10.00 ~ 8,
        . >  10.00 & .  <= 25.00 ~ 9,
        . >  25.00 & .  <= 50.00 ~ 10,
        . >  50.00 & .  <= 75.00 ~ 11,
        . >  75.00 ~ 12))) %>%select_at(
          vars(contains("_R"))) %>% select_all(~str_replace(., "_R", ""))

library(ComplexHeatmap)
library(circlize)
library(viridis)
     

split=c("P1", "P2", "P3", "P4", "P5", "P6")


ha = HeatmapAnnotation("Pol" = anno_block(gp = gpar(
  fill = c("black" ,"black" ,"black", "black", "black", "black")), 
  labels = c("P1", "P2", "P3", "P4","P5",  "P6" ), 
  labels_gp = gpar(col = "white", fontsize = 9, fontface= "bold")))




cols_ho<- list("Season" = c("Dry"="#FE922AFF","Rainy"="#4490FEFF"))

ho = HeatmapAnnotation("Season" = c(rep("Dry", 1), rep("Rainy", "1"),
                                    rep("Dry", 1), rep("Rainy", "1"),
                                    rep("Dry", 1), rep("Rainy", "1"),
                                    rep("Dry", 1), rep("Rainy", "1"),
                                    rep("Dry", 1), rep("Rainy", "1"),
                                    rep("Dry", 1), rep("Rainy", "1")),
                       which = "col", col = cols_ho,
                       annotation_name_gp = gpar(fontsize=8,
                                                 fontface="bold"),
                       show_legend = F, gp = gpar(
                         col = "white", fontize=12), 
                       simple_anno_size = unit(0.25, "cm"),
                       show_annotation_name = T)

my_palette <- viridis::viridis(n = 12, option = "B", direction = -1)



heats<-ComplexHeatmap::Heatmap(
  merge_data,
  col = my_palette,
  row_dend_width = unit(0.4, "cm"),
 width = ncol(merge_data)*unit(6, "mm"), 
 #height = nrow(merge_data)*unit(2.4, "mm"),
  heatmap_legend_param = list(direction = "horizontal",
                              title = "Relative \n abund(%)",
                              grid_height = unit(0.2, "cm"),
                              legend_height = unit(1, "cm"),
                              labels_gp = gpar(fontsize = 7),
                              title_gp = gpar(fontsize = 6, 
                                              fontface="bold"),
                              at = c(0,1,2,3,5,8,10, 50,100),
                              break_dist = 3),
  rect_gp = gpar(col = "white"), 
  cluster_columns = F, cluster_rows = F,
  show_heatmap_legend = T, #top_annotation = c(ha,ho),
  #left_annotation = c(anntro, annguild),
  #column_order = sort(colnames(merge_data)),
 column_split = c(1,1,2,2,3,3,4,4,5,5,6,6) ,
 column_title = NULL,
  show_column_names = F,
  row_names_gp = gpar(fontsize=7),
  column_title_gp = gpar(
    fill = c("#800000" ,"#808000" ,"#008000", "#D35400", "#2E4053" )))

heats

my_palette2 <- viridis::viridis(n = 12, option = "B", direction = -1)


heats2<-ComplexHeatmap::Heatmap(
  merge_data2,
  col = my_palette2,
  row_dend_width = unit(0.4, "cm"),
  width = ncol(merge_data)*unit(6, "mm"), 
  #height = nrow(merge_data)*unit(2.4, "mm"),
  heatmap_legend_param = list(direction = "horizontal",
                              title = "Relative abund(%)",
                              grid_height = unit(0.2, "cm"),
                              legend_height = unit(1, "cm"),
                              labels_gp = gpar(fontsize = 7),
                              title_gp = gpar(fontsize = 6, 
                                              fontface="bold"),
                              at = c(0,1,2,3,5,8,10, 50,100),
                              break_dist = 3),
  rect_gp = gpar(col = "white"), 
  cluster_columns = F, cluster_rows = F,
  show_heatmap_legend = F, top_annotation = c(ha,ho),
  #left_annotation = c(anntro, annguild),
  #column_order = sort(colnames(merge_data)),
  column_split = c(1,1,2,2,3,3,4,4,5,5,6,6) ,
  column_title = NULL,
  show_column_names = F,
  row_names_gp = gpar(fontsize=7),
  column_title_gp = gpar(
    fill = c("#800000" ,"#808000" ,"#008000", "#D35400", "#2E4053" )))

heats2




ht_list = heats2 %v% heats
draw(ht_list)
library(viridis)


lgd1 = Legend(at =  c("Dry", "Rainy"), 
              title = "Season", nrow = 1,
              title_position = "leftcenter",
              legend_gp = gpar(fill = cols_ho$Season),
              labels_gp = gpar(fontsize = 8),
              title_gp = gpar(fontsize = 9, fontface="bold"))


draw(ht_list, heatmap_legend_side = "right",
     annotation_legend_side = "top", 
     merge_legend=F,
     annotation_legend_list = list(lgd1))


heatm<-grid.grabExpr(draw(ht_list, heatmap_legend_side = "top",
                          annotation_legend_side = "top", 
                          merge_legend=F,
                          annotation_legend_list = list(lgd1)))

heatm
ggsave('Plot/heatmap_all_bypol_noroder_guilds.png', width = 11, height = 7, dpi = 300, plot =heatm)

