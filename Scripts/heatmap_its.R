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

table_genus<- tab_taxa %>%  separate(
      taxonomy, c("Kingdom","Phylum","Order","Class","Family","Genus","Species"), sep = ";" ) %>% mutate_at(
        c("Genus"), ~str_replace(., "g__", ""))%>% 
    dplyr::mutate(Genus = stringr::str_trim(Genus, side = "both")) %>% mutate_if(
      is.character, ~replace_na(., "Unassigned")) %>% group_by(
      Genus) %>% summarise_if(is.numeric, sum) %>% column_to_rownames(
          var = "Genus") %>%  mutate(
            all= rowSums(.)) %>% dplyr::arrange(
              -all) %>% relabunda(.) %>%as.data.frame( ) %>%  rownames_to_column(
                var = "Genus")%>% filter(!Genus=="unidentified" ,
                                         !Genus=="Unassigned") %>% 
  filter(!grepl('Incertae_sedis',Genus))%>% slice(
                                           c(1:50))  %>% pivot_longer(
                                             ., cols = -Genus, names_to ="SampleID", 
                                             values_to = "relab" ) %>% filter(
                                               !SampleID=="all")
  cols<- table_genus %>% inner_join(taxones_color) %>% arrange(Genus)
  col <- as.character(cols$color)
  names(col) <- as.character(cols$Genus)

table_phylum <- tab_taxa %>%  separate(
    taxonomy, c("Kingdom","Phylum","Order","Class","Family","Genus","Species"), sep = ";" ) %>% mutate_at(
      c("Phylum"), ~str_replace(., "p__", "")  )%>% mutate_at(
        c("Genus"), ~str_replace(., "g__", ""))%>% group_by(Genus, Phylum) %>% count()
      
   
  
  
barplot_genus<- table_genus %>% inner_join(metadata)%>% dplyr::select(
    Genus, SampleID, relab) %>% pivot_wider(
      names_from = Genus, values_from = relab) %>% column_to_rownames(
        var = "SampleID") %>% t() %>% as.data.frame()



data_fungi_order=barplot_genus[, match(vector_order, colnames(barplot_genus))]%>%
rownames_to_column(var="Genus")


merge_data<- data_fungi_order %>% mutate(Genus=case_when(
        Genus=="Leptosphaeria maculans species complex" ~ "Leptosphaeria", 
        TRUE ~ as.character(Genus)))%>% column_to_rownames(
          var = "Genus") %>% mutate(proms=rowMeans(.)) %>% dplyr::arrange(-proms) %>% slice(
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
col_fun2 = colorRamp2(c(0,  1, 1+1e-5,10,50, 100), 
                      viridis(6, option = "F", direction = -1))
my_palette <- viridis::viridis(n = 12, option = "F", direction = -1)
col_fun = colorRamp2(c(0,  1, 1+1e-5,5,5+1e-5,10,10+1e-5,50,50+1e-5, 75, 75+1e-5,100), 
                     viridis(12, option = "F", direction = -1))

#annotatio seccion y origen (columnas)
annotation_columns<- data.frame(id=colnames(merge_data)) 
rownames(annotation_columns) <- colnames(heatmap)

#set.seed(123)
split = rep(1:6, each = 12)
split= c(rep("1D", 6),rep("1R", 6),
         rep("2D", 6),rep("2R", 6),
         rep("3D", 6),rep("3R", 6),
         rep("4D", 6),rep("4R", 6),
         rep("5D", 6),rep("5R", 6),
         rep("6D", 6),rep("6R", 6)         )

Pol<- c(1,1,1,1,1,1,
          2,2,2,2,2,2,
          3,3,3,3,3,3,
          4,4,4,4,4,4,
          5,5,5,5,5,5,
          6,6,6,6,6,6)

ha = HeatmapAnnotation("Pol" = anno_block(gp = gpar(
  fill = c("black" ,"black" ,"black", "black", "black", "black",
           "black" ,"black" ,"black", "black", "black", "black")), 
  labels = c("P1","P1", "P2","P2", "P3","P3", "P4","P4","P5", "P5", "P6", "P6"), 
  labels_gp = gpar(col = "white", fontsize = 9, fontface= "bold")))


cols_ho<- list("Season" = c("Dry"="#FE922AFF","Rainy"="#4490FEFF"))

ho = HeatmapAnnotation("Season" = c(rep("Dry", 6), rep("Rainy", "6"),
                                    rep("Dry", 6), rep("Rainy", "6"),
                                    rep("Dry", 6), rep("Rainy", "6"),
                                    rep("Dry", 6), rep("Rainy", "6"),
                                    rep("Dry", 6), rep("Rainy", "6"),
                                    rep("Dry", 6), rep("Rainy", "6")),
                       which = "col", col = cols_ho,
                       annotation_name_gp = gpar(fontsize=8,
                                                 fontface="bold"),
                       show_legend = F, gp = gpar(
                         col = "white", fontize=12), 
                       simple_anno_size = unit(0.25, "cm"),
                       show_annotation_name = T)
#annotatio row (filas)
guilds<- guild %>% mutate(Guilds=case_when(
  guild=="Dung Saprotroph-Undefined Saprotroph-Wood Saprotroph"~"Dung Sap.-Undef Sap.-Wood Sap.",
  guild=="Ectomycorrhizal-Fungal Parasite"~"Ectomycorrhizal-Fung.Parasit",
  guild=="Animal Pathogen-Endophyte-Endosymbiont-Epiphyte-Soil Saprotroph-Undefined Saprotroph"~"Animal P-Endoph/symb-Undef. Sap.",
  guild=="Endophyte-Plant Pathogen-Wood Saprotroph"~"Endophyte-Pathogen-Saprotroph",
  guild=="Animal Pathogen-Endophyte-Epiphyte-Fungal Parasite-Plant Pathogen-Wood Saprotroph"~"Animal-Plant P-Endophyte-Saprotroph",
  guild=="Endophyte-Litter Saprotroph-Soil Saprotroph-Undefined Saprotroph"~"Endophyte-Undefined Saprotrop",
  guild=="Animal Endosymbiont-Animal Pathogen-Plant Pathogen-Undefined Saprotroph"~"Animal-Plant endosymb-P-Undef. Sap.",
  guild=="Endophyte-Plant Pathogen-Undefined Saprotroph"~"Endophyte-Plant Pathogen-Undef. Sap.",
  guild=="Ectomycorrhizal-Undefined Saprotroph"~"Ectomycorrhizal-Undef. Sap.",
  guild=="Animal Pathogen-Fungal Parasite-Undefined Saprotroph"~"Animal P-Fungal Parasite-Undef. Sap.",
  guild=="Fungal Parasite-Plant Pathogen-Plant Saprotroph"~"Fungal Parasite-Plant Path-Plant Sap.",
  guild=="Animal Pathogen-Endophyte-Plant Pathogen-Soil Saprotroph-Wood Saprotroph" ~"Animal Pathogen-Endoph-Plant P-Soil/Wood Sap", 
  guild=="Ectomycorrhizal-Orchid Mycorrhizal-Root Associated Biotroph" ~"Ectomycorrhizal-Mycorrhizal-Biotroph",
  guild=="Clavicipitaceous Endophyte-Plant Pathogen" ~"Clavicipitaceous Endophyte-Plant Pat.",
  TRUE~as.character(guild))) 






annotation_rows<- merge_data %>% rownames_to_column(
  var = "Genus") %>% left_join(
    guilds) %>% left_join(table_phylum) %>% dplyr::select_if(
      is.character) %>% replace(
        is.na(.), "Unassigned") %>% left_join(trophic) %>% column_to_rownames(
          var = "Genus") %>% dplyr::rename(Trophic=trophicMode)

cols_phy<-list("Phylum"=c(
  "Basidiomycota" ="#133337",
  "Ascomycota"="#fff68f"  ,
  "Mortierellomycota"="#ffa500",
  "Mucoromycota"="#008080"))

cols_guild <- list('Guild' = c(
  "Unassigned"= "#85929e",
  "Animal Pathogen-Endophyte-Plant Pathogen-Soil Saprotroph-Wood Saprotroph"="#D95F02",
  "Animal Pathogen-Endoph-Plant P-Soil/Wood Sap"="#D95F02",
  "Animal Pathogen-Endophyte-Epiphyte-Fungal Parasite-Plant Pathogen-Wood Saprotroph"="#d658c3",
  "Animal Pathogen-Soil Saprotroph"="#F08080",
  "Animal Endosymbiont-Animal Pathogen-Plant Pathogen-Undefined Saprotroph"="#ca7822",
  "Animal Pathogen"="#ba4a00",
  "Animal Pathogen-Endophyte-Endosymbiont-Epiphyte-Soil Saprotroph-Undefined Saprotroph"="#C70039",
  "Animal P-Endoph/symb-Undef. Sap."="#C70039",
  "Animal Pathogen-Undefined Saprotroph"="#c70051",
  
  "Plant Pathogen"="#E7298A",
  "Fungal Parasite"="#666666", #no sale
  "Animal-Plant P-Endophyte-Saprotroph"="#ff452d",
  "Animal-Plant endosymb-P-Undef. Sap."="#f77c00",
  "Animal Pathogen-Fungal Parasite-Undefined Saprotroph"="#900C3F",
  "Animal P-Fungal Parasite-Undef. Sap."="#ba9194",
  "Fungal Parasite-Plant Pathogen-Plant Saprotroph"="#e22b0f",
  "Fungal Parasite-Plant Path-Plant Sap."="#e22b0f",
  "Plant Pathogen-Undefined Saprotroph"="#b65169",
  
  
  "Arbuscular Mycorrhizal"="#008000",
  "Ectomycorrhizal-Orchid Mycorrhizal-Root Associated Biotroph"="#b4ff68",
  "Ectomycorrhizal-Mycorrhizal-Biotroph"="#b4ff68",
  "Ectomycorrhizal" = "#7FC97F",
  "Ectomycorrhizal-Fungal Parasite"="#58d68d",
  "Ectomycorrhizal-Undefined Saprotroph"="#edf2a3",
  "Ectomycorrhizal-Undef. Sap."="#edf2a3",
  "Endophyte"="#00FFFF",
  "Endophyte-Plant Pathogen"="#117a65",
  "Epiphyte"="#02f0a5",
  "Ericoid Mycorrhizal"="#a2ab16",
  "Lichenized"="#DAF7A6",
  "Ectomycorrhizal-Orchid"="#0000FF", #no sale
  "Endophyte-Plant Pathogen-Undefined Saprotroph" ="#00FF00",
  "Endophyte-Plant Pathogen-Undef. Sap." ="#7dcea0",
  "Endophyte-Pathogen-Saprotroph"="#0b5345",
  "Endophyte-Plant Pathogen-Wood Saprotroph"="#00FF00",
  "Endophyte-Litter Saprotroph-Soil Saprotroph-Undefined Saprotroph"="#73c6b6",
  "Endophyte-Undefined Saprotrop"="#008080",
  "Endophyte-Insect Pathogen"="#00FFFF",
  "Clavicipitaceous Endophyte-Plant Pat."="#5e6326",
  
  "Dung Saprotroph"="#FDC086",
  "Soil Sparotroph"="#851d01",
  "Undefined Saprotroph"="#7d6608",
  "Wood Saprotroph"="#E6AB02",
  "Plant Pathogen-Wood Saprotroph"="#5c2402",
  "Dung Saprotroph-Wood Saprotroph"="#7d5c48",
  "Dung Sap.-Undef Sap.-Wood Sap."="#FDC086",
  "Dung Saprotroph-Undefined Saprotroph-Wood Saprotroph"="#FDC086" ))


cols_tro <- list('Trophic' = c(
  "Pathotroph"=	"#D95F02",
  "Pathotroph-Saprotroph"=	"#581845",
  "Pathotroph-Symbiotroph"="#F08080",
  "Saprotroph"=	"#851d01",
  "Saprotroph-Symbiotroph"=	"#7d5c48",
  "Symbiotroph"=	"#008000",
  "Unassigned"= "#85929e",
  "Pathotroph-Saprotroph-Symbiotroph"=	"#C70039"))

annphyl= HeatmapAnnotation("Phylum" = annotation_rows$Phylum, 
                           which = "row", col = cols_phy,
                           show_legend = T,   
                           show_annotation_name = T,
                           annotation_name_gp =gpar(
                             fontsize = 7, fontface="bold"),
                           annotation_legend_param = list(
                             title_gp = gpar(fontsize = 7, 
                                             fontface="bold"),
                             labels_gp = gpar(fontsize = 7),
                             direction ="horizontal",
                             grid_width = unit(0.3, "cm"),
                             grid_height = unit(0.1, "cm")),
                           
                           simple_anno_size = unit(0.3, "cm"),
                           gp = gpar(col = "white"))


annguild = HeatmapAnnotation("Guilds" = annotation_rows$Guilds, 
                             which = "row", col = cols_guild,
                             show_legend = T,   
                             show_annotation_name = T,
                             annotation_name_gp =gpar(
                               fontsize = 7, fontface="bold"),
                             annotation_legend_param = list(
                               title_gp = gpar(fontsize = 7, 
                                               fontface="bold"),
                               labels_gp = gpar(fontsize = 7),
                               direction ="horizontal",
                               grid_width = unit(0.3, "cm"),
                               grid_height = unit(0.1, "cm")),
                             
                             simple_anno_size = unit(0.3, "cm"),
                             gp = gpar(col = "white"))

anntro = HeatmapAnnotation("Trophic" = annotation_rows$Trophic, 
                           which = "row", col = cols_tro,
                           show_legend = T,   
                           show_annotation_name = T,
                           annotation_name_gp =gpar(
                             fontsize = 7,  fontface="bold"),
                           annotation_legend_param = list(
                             title_gp = gpar(fontsize = 7, 
                                             fontface="bold"),
                             labels_gp = gpar(fontsize =7),
                             direction ="horizontal",
                             grid_width = unit(0.3, "cm"),
                             grid_height = unit(0.1, "cm")),
                           simple_anno_size = unit(0.3, "cm"),
                           gp = gpar(col = "white"))

my_palette <- viridis::viridis(n = 11, option = "B", direction = -1)

heats<-ComplexHeatmap::Heatmap(
  merge_data,
  col = my_palette,
  row_dend_width = unit(0.4, "cm"),
  #width = ncol(merge_data)*unit(2.2, "mm"), 
 #height = nrow(merge_data)*unit(3.4, "mm"),
  heatmap_legend_param = list(direction = "horizontal",
                              title = "Relative abund(%)",
                              grid_height = unit(0.2, "cm"),
                              legend_height = unit(1, "cm"),
                              labels_gp = gpar(fontsize = 7),
                              title_gp = gpar(fontsize = 6, 
                                              fontface="bold"),
                              at = c(0,1,2,3,5,8,10,25, 50, 100),
                              break_dist = 1),
  rect_gp = gpar(col = "white"), 
  cluster_columns = F, cluster_rows = T,
  show_heatmap_legend = T, top_annotation = c(ha,ho),
  right_annotation = c(annguild, anntro, annphyl),
  #column_order = sort(colnames(merge_data)),
  column_split = split, column_title = NULL,
  show_column_names = F,
  row_names_gp = gpar(fontsize=7.5, fontface="italic"),
  column_title_gp = gpar(
    fill = c("#800000" ,"#808000" ,"#008000", "#D35400", "#2E4053" )))

heats
library(viridis)


lgd1 = Legend(at =  c("Dry", "Rainy"), 
              title = "Season", nrow = 1,
              title_position = "leftcenter",
              legend_gp = gpar(fill = cols_ho$Season),
              labels_gp = gpar(fontsize = 8),
              title_gp = gpar(fontsize = 9, fontface="bold"))


draw(heats, heatmap_legend_side = "right",
     annotation_legend_side = "top", 
     merge_legend=F,
     annotation_legend_list = list(lgd1))


heatm<-grid.grabExpr(draw(heats, heatmap_legend_side = "right",
                          annotation_legend_side = "top", 
                          merge_legend=F,
                          annotation_legend_list = list(lgd1)))

heatm
ggsave('Plot/heatmap_all_phyl.png', width = 11, height = 7, dpi = 300, plot =heatm)

