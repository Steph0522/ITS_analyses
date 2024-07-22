tab<- read.delim("Data/table_ojis_filt.txt", check.names = F)
tab2<- tab %>%dplyr::rename(OTU=ASV) %>%
  mutate(ASV=paste0("ASV", rownames(.))) %>% dplyr::select(74,2:73) 

tab<-tab[match(names(fasta), tab$ASV),]

nombres<- tab2$ASV

write.csv(tab2, "Data/table_ojis_filt.csv", row.names = F)

names(fasta)<-nombres

meta<- read.delim("Data/its_map.txt", check.names = F)
write.csv(meta, "Data/its_map.csv", row.names = F)

ver<- read.csv("Data/table_ojis_filt.csv", check.names = F)
vers<- read.csv("~/Downloads/example_files/example_Unoise_table.csv", check.names = F)

library(seqinr)
fasta<- read.fasta("Data/filt_seqs_ojis.fasta")

fastas<- read.delim("Data/filt_seqs_ojis.fasta")

write.fasta(sequences = fasta, file.out = "Data/filt_seqs_ojis.fasta", names = names(fasta))
