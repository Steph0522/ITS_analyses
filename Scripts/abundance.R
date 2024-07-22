tab98<- read.delim("/media/yendi/ADATA HD830/METAGENOMAS/CORREDOR_METAGENOMA/metagenoma_corredor/99/99_plus.txt", skip = 1, check.names = F, row.names = 1)
tab98<- read.delim("/media/yendi/ADATA HD830/METAGENOMAS/CORREDOR_METAGENOMA/metagenoma_corredor/98/98.txt", skip = 1, check.names = F, row.names = 1)

tab9_tax<- tab98 %>% separate(taxonomy, into = c("k", "p", "c", "o", "f", "g", "s"), sep = ";")

tab9_bact<- tab9_tax %>% filter(k=="k__Bacteria")
tab9_fungi<- tab9_tax %>% filter(k=="k__Eukaryota") %>% filter(!p==" p__Chordata") %>% 
  filter(!p==" p__Apicomplexa") %>% filter(!p==" p__Ciliophora" ) %>% 
  filter(!p== " p__Bacillariophyta" ) %>% filter(!p==" p__Cercozoa" ) %>% filter(!p==" p__Euglenozoa" ) %>% 
  filter(!p==" p__Parabasalia") %>% filter(!p==" p__Fornicata") %>%
  filter(!p==" p__Evosea") #%>% filter(!p==" p__Microsporidia")
  

tab98_fungi<- tab9_fungi %>% unite("taxonomy", k:s, sep = ";")

tab9_vir<- tab9_tax %>% filter(k=="k__Viruses")
tab9_arq<- tab9_tax %>% filter(k=="k__Archaea")


nrow(tab9_fungi)/nrow(tab9_tax)*100
nrow(tab9_bact)/nrow(tab9_tax)*100
nrow(tab9_vir)/nrow(tab9_tax)*100
nrow(tab9_arq)/nrow(tab9_tax)*100

gen <- tab9_fungi %>% group_by(g) %>% summarise_if(is.numeric, sum) %>% dplyr::arrange(-kraken_fungi_98_report_p_bracken)

f99<- read.delim("/media/yendi/ADATA HD830/METAGENOMAS/CORREDOR_METAGENOMA/metagenoma_corredor/99/99_fungi.txt", skip = 1, check.names = F, row.names = 1)
f98<- read.delim("/media/yendi/ADATA HD830/METAGENOMAS/CORREDOR_METAGENOMA/metagenoma_corredor/98/98_k/98_fungi.txt", skip = 1, check.names = F, row.names = 1)



