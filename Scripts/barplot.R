library(tidyverse)
table_taxa<- seqtab_nochim %>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "Feature.ID") %>% inner_join(as.data.frame(taxa) %>% rownames_to_column(var="Feature.ID"))
meta<- read.delim("../../../../maps_its/its_map.txt")

table_genus<- table_taxa %>% group_by(Genus) %>% summarise_if(is.numeric, sum) %>%
  replace(is.na(.), "Unassigned") %>% mutate(
    total = rowSums(across(where(is.numeric)))) %>% filter(!total<1500) %>% dplyr::select(-total)

table_genus_pol_season<- table_genus %>% column_to_rownames(var = "Genus") %>% 
  t() %>% as.data.frame() %>% rownames_to_column(var = "SampleID") %>% 
  inner_join(meta) %>% group_by(Poligono,Season) %>%
  summarise_if(is.numeric, sum) %>% unite("ids", Poligono:Season) %>% column_to_rownames(var = "ids") %>% 
  t() %>% as.data.frame()

table_genus_pol_season  %>% rownames_to_column(var="Genus") %>% pivot_longer(
  cols = -Genus, names_to = "SampleID", values_to = "abund") %>% ggplot(
    aes(SampleID, abund, fill=Genus), color="black")+geom_col(position = "fill")
