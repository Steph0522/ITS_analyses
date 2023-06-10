setwd("~/Documents/ITS_corredor/ITS_analysis/")
library(tidyverse)
library(metagMisc)
library("Biostrings")
library(dada2)

#seqs
#load("all_seqs/joined_trimmed/joined_trimmed/resultados/ITS_dada2_results_joined_trimmed.RData")
#load("all_seqs/raw/ITS_dada2_results.RData")
#load("all_seqs/trimmed/resultados/ITS_dada2_results_trimmed.RData")
load("../all_seqs/joined_trimmed_raw/joined_trimmed/ITS_dada2_results_joined_trimmed.RData")

seqtab <- getUniques(seqtab_nochim)
dada_to_fasta(seqtab, out = "seqs_joined_raw.fasta", hash = "md5")


fastaFile <- readDNAStringSet("seqs_joined_raw.fasta")
"Feature ID" = names(fastaFile)
Sequence = paste(fastaFile)


joined_trimmed <- data.frame(`Feature ID`, Sequence, check.names = F)
#raw <- data.frame(`Feature ID`, Sequence, check.names = F)
#trimmed <- data.frame(`Feature ID`, Sequence, check.names = F)
#trimmed_merged <- data.frame(`Feature ID`, Sequence, check.names = F)

remove(dadaFs, dadaRs, errF, errR, mergers, fastaFile, seqtab_nochim, taxa)
remove(seqtab, Sequence, `Feature ID`)

#tables

table_raw<- read.csv("all_seqs/raw/table.csv", row.names = 1)
table_trimmed<- read.csv("all_seqs/trimmed/table.csv", row.names = 1)
table_joined_trimmed<- read.csv("../all_seqs/joined_trimmed_raw/joined_trimmed//table.csv", row.names = 1)
table_trimmed_merged<- read.csv("all_seqs/trimmed_merged/table.csv", row.names = 1)

table_raw_t<- table_raw %>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "Sequence") %>% inner_join(raw) %>% dplyr::select(`Feature ID`, `111D`:`623R`) %>% 
  dplyr::rename("#OTUID"=`Feature ID`)

table_trimmed_t<- table_trimmed%>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "Sequence") %>% inner_join(trimmed) %>% dplyr::select(`Feature ID`, `111D`:`623R`) %>% 
  dplyr::rename("#OTUID"=`Feature ID`)

table_joined_t<- table_joined_trimmed %>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "Sequence") %>% inner_join(joined_trimmed) %>% dplyr::select(`Feature ID`, `111D`:`623R`) %>% 
  dplyr::rename("#OTUID"=`Feature ID`)

table_merged_t<- table_trimmed_merged %>% t() %>% as.data.frame() %>% rownames_to_column(
  var = "Sequence") %>% inner_join(trimmed_merged) %>% dplyr::select(`Feature ID`, `111D`:`623R`) %>% 
  dplyr::rename("#OTUID"=`Feature ID`)

write_tsv(table_merged_t,"table_trimmed_merged.txt")
write_tsv(table_raw_t,"table_raw.txt")
write_tsv(table_trimmed_t,"table_trimmed.txt")
write_tsv(table_joined_t,"table_joined_trimmed_raw.txt")



#taxonomys
taxonomy_raw<- read.csv("all_seqs/raw/taxonomy.csv")
taxonomy_trimmed<- read.csv("all_seqs/trimmed/taxonomy.csv")
taxonomy_joined_trimmed<- read.csv("../all_seqs/joined_trimmed_raw//joined_trimmed/taxonomy.csv")
taxonomy_trimmed_merged<- read.csv("all_seqs/trimmed_merged/taxonomy.csv")

taxonomy_raw_t<- taxonomy_raw %>% unite("Taxon", Kingdom:Species, sep = ";") %>% dplyr::rename(
  Sequence=X) %>% inner_join(raw) %>% dplyr::select(`Feature ID`, Taxon)

taxonomy_trimmed_t<- taxonomy_trimmed %>% unite("Taxon", Kingdom:Species, sep = ";") %>% dplyr::rename(
  Sequence=X) %>% inner_join(trimmed) %>% dplyr::select(`Feature ID`, Taxon)

taxonomy_joined_t<- taxonomy_joined_trimmed%>% unite("Taxon", Kingdom:Species, sep = ";") %>% dplyr::rename(
  Sequence=X) %>% inner_join(joined_trimmed) %>% dplyr::select(`Feature ID`, Taxon)

taxonomy_merged_t<- taxonomy_trimmed_merged%>% unite("Taxon", Kingdom:Species, sep = ";") %>% dplyr::rename(
  Sequence=X) %>% inner_join(trimmed_merged) %>% dplyr::select(`Feature ID`, Taxon)


write_tsv(taxonomy_merged_t,"taxonomy_trimmed_merged.txt")
write_tsv(taxonomy_trimmed_t,"taxonomy_trimmed.txt")
write_tsv(taxonomy_raw_t,"taxonomy_raw.txt")
write_tsv(taxonomy_joined_t,"taxonomy_joined_trimmed_raw.txt")
