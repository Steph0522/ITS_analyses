library(qiime2R)
library(tidyverse)

table_oc<-read_qza("Data/table_oc.qza")$data
table_oim<- read_qza("Data/table_oim.qza")$data
table_oism<- read_qza("Data/table_oism.qza")$data 
table_ojis<- read_qza("Data/table_ojis.qza")$data
table_oiss<- read_qza("Data/table_oiss.qza")$data
table_ac<- read_qza("Data/table_ac.qza")$data
table_aism<- read_qza("Data/table_aism.qza")$data
table_aisc<- read_qza("Data/table_aisc.qza")$data
table_ajis<- read_qza("Data/table_ajis.qza")$data

filt_table_oc<-read_qza("Data/filt_table_oc.qza")$data
filt_table_oim<- read_qza("Data/filt_table_oim.qza")$data
filt_table_oism<- read_qza("Data/filt_table_oism.qza")$data 
filt_table_ojis<- read_qza("Data/filt_table_ojis.qza")$data
filt_table_oiss<- read_qza("Data/filt_table_oiss.qza")$data
filt_table_ac<- read_qza("Data/filt_table_ac.qza")$data
filt_table_aism<- read_qza("Data/filt_table_aism.qza")$data
filt_table_aisc<- read_qza("Data/filt_table_aisc.qza")$data
filt_table_ajis<- read_qza("Data/filt_table_ajis.qza")$data


list_table<- list(table_oc, table_oim, table_oism, table_ojis, table_oiss,
                  table_ac, table_aism, table_aisc, table_ajis)
names(list_table)<- c("OC", "OIM", "OISM", "OJIS", "OISS",
                     "AC", "AISM", "AISC", "AJIS")
list_single<- list(filt_table_oc, filt_table_oim, filt_table_oism, filt_table_ojis, 
                   filt_table_oiss, filt_table_ac, filt_table_aism,
                   filt_table_aisc, filt_table_ajis)
names(list_single)<- c("OC", "OIM", "OISM", "OJIS", "OISS",
                     "AC", "AISM", "AISC", "AJIS")

features<-sapply(list_table, nrow)
features_single <- sapply(list_single, nrow)

reads<- sapply(list_table, sum)
read_single<- sapply(list_single, sum)



multiple.func <- function(x) {
  c(min = min(colSums(x)),
    max = max(colSums(x)),
    mean = mean(colSums(x)),
    median=median(colSums(x)))}

resume_features<- sapply(list_table, multiple.func)
resume_features_single<- sapply(list_single, multiple.func)

library(ggpubr)
table_features_reads<- data.frame(
  features= features, features_single=features_single,
  reads=reads, reads_single=read_single) %>% rownames_to_column(var = "method")

DT::datatable(
  table_features_reads,
  fillContainer = FALSE, options = list(pageLength = 10),
  filter="top")
DT::datatable(
  as.data.frame(t(round(resume_features))) %>% rownames_to_column(var = "method"),
  fillContainer = FALSE, options = list(pageLength = 10),
  filter=list(position = 'top', clear = TRUE, plain = FALSE)
)
DT::datatable(
  as.data.frame(t(round(resume_features_single))) %>% rownames_to_column(var = "method"),
  fillContainer = FALSE, options = list(pageLength = 10),
  filter=list(position = 'top', clear = TRUE, plain = FALSE)
)

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

list_taxa<- list(taxa_oc, taxa_oim, taxa_oism, taxa_ojis,
              taxa_oiss, taxa_ac, taxa_aism, taxa_aisc, taxa_ajis)


names(list_taxa)<- c("OC", "OIM", "OISM", "OJIS", "OISS",
                       "AC", "AISM", "AISC", "AJIS")


count_p<- suppressPackageStartupMessages(suppressWarnings(function(x,y){x %>% as.data.frame() %>% 
  rownames_to_column(var = "Feature.ID") %>% inner_join(
    y) %>% separate(
      "Taxon", c("k", "p", "c", "o", "f", "g", "s"), sep = ";") %>% 
  group_by(p) %>% summarise_if(is.numeric, sum) %>% drop_na() %>% 
  filter(!p=="p__Fungi_phy_Incertae_sedis")%>% 
  filter(!p=="Unassigned") %>% nrow()}))
count_c<- suppressPackageStartupMessages(suppressWarnings(function(x,y){x %>% as.data.frame() %>% 
    rownames_to_column(var = "Feature.ID") %>% inner_join(
      y) %>% separate(
        "Taxon", c("k", "p", "c", "o", "f", "g", "s"), sep = ";") %>% 
    group_by(c) %>% summarise_if(is.numeric, sum) %>% drop_na() %>% 
    filter(!c=="Unassigned") %>% filter(!c=="c__Fungi_cls_Incertae_sedis")  %>% nrow() }))
count_o<- suppressPackageStartupMessages(suppressWarnings(function(x,y){x %>% as.data.frame() %>% 
    rownames_to_column(var = "Feature.ID") %>% inner_join(
      y) %>% separate(
        "Taxon", c("k", "p", "c", "o", "f", "g", "s"), sep = ";") %>% 
    group_by(o) %>% summarise_if(is.numeric, sum) %>% drop_na() %>% 
    filter(!o=="Unassigned") %>% filter(!o=="o__Fungi_ord_Incertae_sedis")  %>% nrow() }))
count_f<- suppressPackageStartupMessages(suppressWarnings(function(x,y){x %>% as.data.frame() %>% 
    rownames_to_column(var = "Feature.ID") %>% inner_join(
      y) %>% separate(
        "Taxon", c("k", "p", "c", "o", "f", "g", "s"), sep = ";") %>% 
    group_by(f) %>% summarise_if(is.numeric, sum) %>% drop_na() %>% 
    filter(!f=="Unassigned") %>% filter(!f=="f__Fungi_fam_Incertae_sedis")  %>% nrow() }))
count_g<- suppressPackageStartupMessages(suppressWarnings(function(x,y){x %>% as.data.frame() %>% 
    rownames_to_column(var = "Feature.ID") %>% inner_join(
      y) %>% separate(
        "Taxon", c("k", "p", "c", "o", "f", "g", "s"), sep = ";") %>% 
    group_by(g) %>% summarise_if(is.numeric, sum) %>% drop_na() %>% 
    filter(!g=="Unassigned") %>% filter(!g=="g__Fungi_gen_Incertae_sedis")  %>% nrow() }))
count_una<- suppressPackageStartupMessages(suppressWarnings(function(x,y){x %>% as.data.frame() %>% 
    rownames_to_column(var = "Feature.ID") %>% inner_join(
      y) %>% filter(Taxon=="k__Fungi") 
  
  }))



phyl<- mapply(count_p, list_table, list_taxa)
clas<- mapply(count_c, list_table, list_taxa)
ord<- mapply(count_o, list_table, list_taxa)
fam<- mapply(count_f, list_table, list_taxa)
gen<- mapply(count_g, list_table, list_taxa)
una<- mapply(count_una, list_table, list_taxa)

levels_df<- data.frame(phyla= phyl, classes=clas,
                       orders=ord, families=fam,
                       genera=gen, unassigned=lengths(una[1,]))

DT::datatable(
  levels_df %>% rownames_to_column(var = "method"),
  fillContainer = FALSE, options = list(pageLength = 10),
  filter=list(position = 'top', clear = TRUE, plain = FALSE)
)


