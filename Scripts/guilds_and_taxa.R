
library(tidyverse)

tab <- read.delim("Data/filt_table_ojis/table_ojis_filt.txt", check.names = F) %>% dplyr::rename(Feature.ID=ASV)
taxa<- qiime2R::read_qza("Data/taxonomy_ojis.qza")$data %>% as.data.frame()

tab_tax<- tab %>% inner_join(taxa)

tab_taxa<- tab_tax %>% dplyr::select(OTUID=Feature.ID, `111D`:`623R`, taxonomy=Taxon)

#write_tsv(tab_taxa,"/home/yendi/Documents/FUNGuild/table_ojis_funguild.txt")

#guilds

guild1<- read.delim("~/Documents/FUNGuild/table_ojis_funguild.guilds_matched.txt", check.names = F)

relabunda<- function(x){(as.data.frame(t(t(x)/colSums(x)))*100)}

metadata<- read.delim("Data/its_map.txt")

library(RColorBrewer)

guild_colors<- read.csv("~/Documents/2023/corredor_scripts/Fungal_Communities_PNIP_PNML/Data/colors_guild")
guild_colors$col[5]<-"#aed6f1"

library(tidyverse)
library(ggh4x)
table_guild<- guild1_fil %>%dplyr::select(
  `111D`:`623R`, Guild) %>% group_by(
    Guild) %>% summarise_if(
      is.numeric, sum)  %>%  column_to_rownames(
          var = "Guild") %>%  mutate(
            all= rowSums(.)) %>% dplyr::arrange(-all) %>%
  filter(all>2000) %>% 
  relabunda(.) %>% rownames_to_column(
    var = "Guild") %>% 
  pivot_longer(., cols = -Guild, names_to ="SampleID", 
               values_to = "relab" ) %>% filter(
                 !SampleID=="all") 
cols_guild<- table_guild %>% inner_join(guild_colors) %>% arrange(Guild)
col_guild <- as.character(cols_guild$col)
names(col_guild) <- as.character(cols_guild$Guild)


barplot_guild<- table_guild %>% inner_join(metadata) %>% ggbarplot(
      x = "Season", y = "relab",add = "mean",
      facet.by = "Poligono", fill="Guild", 
      position = position_fill()) +facet_nested(
        .~Poligono, scales = "free_x")+scale_fill_manual(
          name = "Guild",values =guild_colors$col )+
  theme_linedraw()+ylab("Relative abundance")+
  xlab("")+theme(legend.text = element_text(face = "plain"))+
  #  guides(fill = guide_legend(nrow = 30))+
  theme(legend.title = element_text(size = 9),
        #axis.ticks = element_blank(),
        legend.text = element_text(size = 9), 
        axis.text.x = element_text(size = 10),
        legend.key.size = unit(0.6, 'cm'), #change legend key size
        legend.key.height = unit(0.45, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'),
        strip.text.x = element_text(size = 16),
        legend.box = "vertical")
print(barplot_guild)

ggsave("Plot/guilds.png",width = 13, height = 7, dpi = 300, plot = barplot_guild, device = "png")


table_trophic<- guild1_fil %>%dplyr::select(
  `111D`:`623R`, `Trophic Mode`) %>% group_by(
    `Trophic Mode`) %>% summarise_if(
      is.numeric, sum)  %>%  column_to_rownames(
        var = "Trophic Mode") %>%  mutate(
          all= rowSums(.)) %>% dplyr::arrange(-all) %>%
  filter(all>2000) %>% 
  relabunda(.) %>% rownames_to_column(
    var = "Trophic") %>% 
  pivot_longer(., cols = -Trophic, names_to ="SampleID", 
               values_to = "relab" ) %>% filter(
                 !SampleID=="all") 

trophic_colors<- read.delim("~/Documents/2023/corredor_scripts/Fungal_Communities_PNIP_PNML/Data/colors_trophic")

cols_trophic<- table_trophic %>% inner_join(trophic_colors) %>% arrange(Trophic)
col_trophic <- as.character(cols_trophic$col)
names(col_trophic) <- as.character(cols_trophic$Trophic)


barplot_trophic<- table_trophic %>% inner_join(metadata) %>% ggbarplot(
  x = "Season", y = "relab",add = "mean",
  facet.by = "Poligono", fill="Trophic", 
  position = position_fill()) +facet_nested(
    .~Poligono, scales = "free_x")+scale_fill_manual(
      name = "Trophic",values =trophic_colors$col )+
  theme_linedraw()+ylab("Relative abundance")+
  xlab("")+theme(legend.text = element_text(face = "plain"))+
  #  guides(fill = guide_legend(nrow = 30))+
  theme(legend.title = element_text(size = 9),
        #axis.ticks = element_blank(),
        legend.text = element_text(size = 9), 
        axis.text.x = element_text(size = 10),
        legend.key.size = unit(0.6, 'cm'), #change legend key size
        legend.key.height = unit(0.45, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'),
        strip.text.x = element_text(size = 16),
        legend.box = "vertical")
print(barplot_trophic)

ggsave("Plot/trophics.png",width = 13, height = 7, dpi = 300, plot = barplot_trophic, device = "png")

#phyla

table_phyl<- tab_taxa %>%dplyr::select(
  `111D`:`623R`, taxonomy)  %>% 
  separate(taxonomy, into = c("k", "phyl", "c", "o", "f", "gen", "s"), sep = ";")%>%
  group_by(phyl) %>% summarise_if(
      is.numeric, sum)  %>%
  mutate_at(vars(phyl), ~replace_na(., "Unassigned")) %>% 
  column_to_rownames(
        var = "phyl")   %>% 
  relabunda(.) %>% rownames_to_column(
    var = "phyl") %>% 
  pivot_longer(., cols = -phyl, names_to ="SampleID", 
               values_to = "relab" )%>% 
  dplyr::mutate(across('phyl', str_replace, 'p__', ''))

barplot_phyl<- table_phyl %>% inner_join(metadata) %>% ggbarplot(
  x = "Season", y = "relab",add = "mean",
  facet.by = "Poligono", fill="phyl", 
  position = position_fill()) +facet_nested(
    .~Poligono, scales = "free_x")+scale_fill_manual(
      name = "Phyla",values =guild_colors$col )+
  theme_linedraw()+ylab("Relative abundance")+
  xlab("")+theme(legend.text = element_text(face = "plain"))+
  #  guides(fill = guide_legend(nrow = 30))+
  theme(legend.title = element_text(size = 9),
        #axis.ticks = element_blank(),
        legend.text = element_text(size = 9), 
        axis.text.x = element_text(size = 10),
        legend.key.size = unit(0.6, 'cm'), #change legend key size
        legend.key.height = unit(0.45, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'),
        strip.text.x = element_text(size = 16),
        legend.box = "vertical")
print(barplot_phyl)


#gena

table_gen<- tab_taxa %>%dplyr::select(
  `111D`:`623R`, taxonomy)  %>% 
  separate(taxonomy, into = c("k", "p", "c", "o", "f", "gen", "s"), sep = ";")%>%
  group_by(gen) %>% summarise_if(
    is.numeric, sum)  %>%
  mutate_at(vars(gen), ~replace_na(., "Unassigned")) %>% 
  column_to_rownames(
    var = "gen")   %>% 
  relabunda(.) %>%
  mutate(all=rowMeans(.)) %>% dplyr::arrange(-all) %>% 
  filter(all>0.7) %>% 
  dplyr::select(-all) %>% 
  rownames_to_column(
    var = "gen") %>% 
  pivot_longer(., cols = -gen, names_to ="SampleID", 
               values_to = "relab" )%>% 
  dplyr::mutate(across('gen', str_replace, 'g__', ''))

barplot_gen<- table_gen %>% inner_join(metadata) %>% ggbarplot(
  x = "Season", y = "relab",add = "mean",
  facet.by = "Poligono", fill="gen", 
  position = position_fill()) +facet_nested(
    .~Poligono, scales = "free_x")+scale_fill_manual(
      name = "Genera",values =guild_colors$col )+
  theme_linedraw()+ylab("Relative abundance")+
  xlab("")+theme(legend.text = element_text(face = "plain"))+
  #  guides(fill = guide_legend(nrow = 30))+
  theme(legend.title = element_text(size = 9),
        #axis.ticks = element_blank(),
        legend.text = element_text(size = 9), 
        axis.text.x = element_text(size = 10),
        legend.key.size = unit(0.6, 'cm'), #change legend key size
        legend.key.height = unit(0.45, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'),
        strip.text.x = element_text(size = 16),
        legend.box = "vertical")
print(barplot_gen)


ggsave("Plot/gen.png",width = 13, height = 7, dpi = 300, plot = barplot_gen, device = "png")
ggsave("Plot/phyl.png",width = 13, height = 7, dpi = 300, plot = barplot_phyl, device = "png")
