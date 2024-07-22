# Incidence-based Raup-Crick beta-diversity
library(qiime2R)
library(tidyverse)
library(vegan)
library(readxl)

abund_table=read_qza("Data/filt_table_oism.qza")$data # Change to "OTU16S_rar.csv" if you want to assess the bacterial dataset
#colnames(abund_table)<-stringr::str_extract(colnames(abund_table),  "^\\w..[^_]")
colnames(abund_table)<-stringr::str_extract(colnames(abund_table),  "[^_]+(?=_................\\w[^_]*$)")

dry<- read_excel("Data/Fisicoquimicos.xlsx", sheet = "metadata_secas_all") %>% mutate(Season="Dry")
wet<- read_excel("Data/Fisicoquimicos.xlsx", sheet = "metadata_lluvias_all") %>% mutate(Season="Wet")

drys<- dry %>% dplyr::select( Poligono,Sitio,Transecto,pH,MO,N,P,Season)
wets<- wet %>% dplyr::select(Poligono,Sitio,Transecto,pH,MO,N,P,Season)

metadata<- rbind(drys,wets) %>% mutate(sea=c(rep("D",36), rep("R",36))) %>% 
  unite("SampleID", c("Poligono", "Sitio", "Transecto", "sea"),remove = F, sep = "") %>% column_to_rownames(var = "SampleID")

meta_table=metadata
meta_table$Poligono<- as.factor(meta_table$Poligono)

tp1=meta_table[which(meta_table$Poligono=="1"), ]
tp2=meta_table[which(meta_table$Poligono=="2"), ]
tp3=meta_table[which(meta_table$Poligono=="3"), ]
tp4=meta_table[which(meta_table$Poligono=="4"), ]
tp5=meta_table[which(meta_table$Poligono=="5"), ]
tp6=meta_table[which(meta_table$Poligono=="6"), ]

tpd<-meta_table[which(meta_table$Season=="Dry"),]
tpr<-meta_table[which(meta_table$Season=="Wet"),]

tp1d<- meta_table[which(meta_table$Season=="Dry" & meta_table$Poligono=="1"),]
tp2d<- meta_table[which(meta_table$Season=="Dry" & meta_table$Poligono=="2"),]
tp3d<- meta_table[which(meta_table$Season=="Dry" & meta_table$Poligono=="3"),]
tp4d<- meta_table[which(meta_table$Season=="Dry" & meta_table$Poligono=="4"),]
tp5d<- meta_table[which(meta_table$Season=="Dry" & meta_table$Poligono=="5"),]
tp6d<- meta_table[which(meta_table$Season=="Dry" & meta_table$Poligono=="6"),]


tp1r<- meta_table[which(meta_table$Season=="Wet" & meta_table$Poligono=="1"),]
tp2r<- meta_table[which(meta_table$Season=="Wet" & meta_table$Poligono=="2"),]
tp3r<- meta_table[which(meta_table$Season=="Wet" & meta_table$Poligono=="3"),]
tp4r<- meta_table[which(meta_table$Season=="Wet" & meta_table$Poligono=="4"),]
tp5r<- meta_table[which(meta_table$Season=="Wet" & meta_table$Poligono=="5"),]
tp6r<- meta_table[which(meta_table$Season=="Wet" & meta_table$Poligono=="6"),]


tp1i=row.names(tp1)
tp2i=row.names(tp2)
tp3i=row.names(tp3)
tp4i=row.names(tp4)
tp5i=row.names(tp5)
tp6i=row.names(tp6)

tpdi<- rownames(tpd)
tpri<- rownames(tpr)

tp1di=row.names(tp1d)
tp2di=row.names(tp2d)
tp3di=row.names(tp3d)
tp4di=row.names(tp4d)
tp5di=row.names(tp5d)
tp6di=row.names(tp6d)

tp1ri=row.names(tp1r)
tp2ri=row.names(tp2r)
tp3ri=row.names(tp3r)
tp4ri=row.names(tp4r)
tp5ri=row.names(tp5r)
tp6ri=row.names(tp6r)

abund_table=t(abund_table)
row.names(abund_table)

abund_table1=subset(abund_table, row.names(abund_table) %in% tp1i)
abund_table2=subset(abund_table, row.names(abund_table) %in% tp2i)
abund_table3=subset(abund_table, row.names(abund_table) %in% tp3i)
abund_table4=subset(abund_table, row.names(abund_table) %in% tp4i)
abund_table5=subset(abund_table, row.names(abund_table) %in% tp5i)
abund_table6=subset(abund_table, row.names(abund_table) %in% tp6i)
abund_tabled= subset(abund_table, rownames(abund_table) %in% tpdi)
abund_tabler= subset(abund_table, rownames(abund_table) %in% tpri)


abund_table1d=subset(abund_table, row.names(abund_table) %in% tp1di)
abund_table2d=subset(abund_table, row.names(abund_table) %in% tp2di)
abund_table3d=subset(abund_table, row.names(abund_table) %in% tp3di)
abund_table4d=subset(abund_table, row.names(abund_table) %in% tp4di)
abund_table5d=subset(abund_table, row.names(abund_table) %in% tp5di)
abund_table6d=subset(abund_table, row.names(abund_table) %in% tp6di)

abund_table1r=subset(abund_table, row.names(abund_table) %in% tp1ri)
abund_table2r=subset(abund_table, row.names(abund_table) %in% tp2ri)
abund_table3r=subset(abund_table, row.names(abund_table) %in% tp3ri)
abund_table4r=subset(abund_table, row.names(abund_table) %in% tp4ri)
abund_table5r=subset(abund_table, row.names(abund_table) %in% tp5ri)
abund_table6r=subset(abund_table, row.names(abund_table) %in% tp6ri)

abund_table1 = abund_table1[,which(colSums(abund_table1) != 0)]
abund_table2 = abund_table2[,which(colSums(abund_table2) != 0)]
abund_table3 = abund_table3[,which(colSums(abund_table3) != 0)]
abund_table4 = abund_table4[,which(colSums(abund_table4) != 0)]
abund_table5 = abund_table5[,which(colSums(abund_table5) != 0)]
abund_table6 = abund_table6[,which(colSums(abund_table6) != 0)]

abund_tabled = abund_tabled[,which(colSums(abund_tabled) != 0)]
abund_tabler = abund_tabler[,which(colSums(abund_tabler) != 0)]

abund_table1d = abund_table1d[,which(colSums(abund_table1d) != 0)]
abund_table2d = abund_table2d[,which(colSums(abund_table2d) != 0)]
abund_table3d = abund_table3d[,which(colSums(abund_table3d) != 0)]
abund_table4d = abund_table4d[,which(colSums(abund_table4d) != 0)]
abund_table5d = abund_table5d[,which(colSums(abund_table5d) != 0)]
abund_table6d = abund_table6d[,which(colSums(abund_table6d) != 0)]

abund_table1r = abund_table1r[,which(colSums(abund_table1r) != 0)]
abund_table2r = abund_table2r[,which(colSums(abund_table2r) != 0)]
abund_table3r = abund_table3r[,which(colSums(abund_table3r) != 0)]
abund_table4r = abund_table4r[,which(colSums(abund_table4r) != 0)]
abund_table5r = abund_table5r[,which(colSums(abund_table5r) != 0)]
abund_table6r = abund_table6r[,which(colSums(abund_table6r) != 0)]


rc1=QsNullModels::raup_crick(abund_table1, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc2=QsNullModels::raup_crick(abund_table2, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc3=QsNullModels::raup_crick(abund_table3, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc4=QsNullModels::raup_crick(abund_table4, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc5=QsNullModels::raup_crick(abund_table5, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc6=QsNullModels::raup_crick(abund_table6, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)

rcd=QsNullModels::raup_crick(abund_tabled, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rcr=QsNullModels::raup_crick(abund_tabler, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
#rall=QsNullModels::raup_crick(abund_table, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)

rc1d=QsNullModels::raup_crick(abund_table1d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc2d=QsNullModels::raup_crick(abund_table2d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc3d=QsNullModels::raup_crick(abund_table3d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc4d=QsNullModels::raup_crick(abund_table4d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc5d=QsNullModels::raup_crick(abund_table5d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc6d=QsNullModels::raup_crick(abund_table6d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)

rc1r=QsNullModels::raup_crick(abund_table1r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc2r=QsNullModels::raup_crick(abund_table2r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc3r=QsNullModels::raup_crick(abund_table3r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc4r=QsNullModels::raup_crick(abund_table4r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc5r=QsNullModels::raup_crick(abund_table5r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc6r=QsNullModels::raup_crick(abund_table6r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)

library(reshape2)
#rc1m<- as.matrix(rc1)
#rc2m<- as.matrix(rc2)
#rc3m<- as.matrix(rc3)
#rc4m<- as.matrix(rc4)
#rc5m<- as.matrix(rc5)
#rc6m<- as.matrix(rc6)
#rc1md<- as.matrix(rc1d)
#rc2md<- as.matrix(rc2d)
#rc3md<- as.matrix(rc3d)
#rc4md<- as.matrix(rc4d)
#rc5md<- as.matrix(rc5d)
#rc6md<- as.matrix(rc6d)
#rc1mr<- as.matrix(rc1r)
#rc2mr<- as.matrix(rc2r)
#rc3mr<- as.matrix(rc3r)
#rc4mr<- as.matrix(rc4r)
#rc5mr<- as.matrix(rc5r)
#rc6mr<- as.matrix(rc6r)
#rc1m[upper.tri(rc1m)]<-NA

rcc1=subset(melt(as.matrix(rc1)),value!=0)
rcc2=subset(melt(as.matrix(rc2)),value!=0)
rcc3=subset(melt(as.matrix(rc3)),value!=0)
rcc4=subset(melt(as.matrix(rc4)),value!=0)
rcc5=subset(melt(as.matrix(rc5)),value!=0)
rcc6=subset(melt(as.matrix(rc6)),value!=0)
rccd=subset(melt(as.matrix(rcd)),value!=0)
rccr=subset(melt(as.matrix(rcr)),value!=0)

rcc1d=subset(melt(as.matrix(rc1d)),value!=0)
rcc2d=subset(melt(as.matrix(rc2d)),value!=0)
rcc3d=subset(melt(as.matrix(rc3d)),value!=0)
rcc4d=subset(melt(as.matrix(rc4d)),value!=0)
rcc5d=subset(melt(as.matrix(rc5d)),value!=0)
rcc6d=subset(melt(as.matrix(rc6d)),value!=0)

rcc1r=subset(melt(as.matrix(rc1r)),value!=0)
rcc2r=subset(melt(as.matrix(rc2r)),value!=0)
rcc3r=subset(melt(as.matrix(rc3r)),value!=0)
rcc4r=subset(melt(as.matrix(rc4r)),value!=0)
rcc5r=subset(melt(as.matrix(rc5r)),value!=0)
rcc6r=subset(melt(as.matrix(rc6r)),value!=0)

rcc1$pol=1
rcc2$pol=2
rcc3$pol=3
rcc4$pol=4
rcc5$pol=5
rcc6$pol=6
rccd$sea="dry"
rccr$sea="wet"

rcc1d$factor="P1-dry"
rcc2d$factor="P2-dry"
rcc3d$factor="P3-dry"
rcc4d$factor="P4-dry"
rcc5d$factor="P5-dry"
rcc6d$factor="P6-dry"

rcc1r$factor="P1-wet"
rcc2r$factor="P2-wet"
rcc3r$factor="P3-wet"
rcc4r$factor="P4-wet"
rcc5r$factor="P5-wet"
rcc6r$factor="P6-wet"


rcc_m=rbind(rcc1, rcc2,rcc3,rcc4,rcc5,rcc6)
rcc_dm=rbind(rccd, rccr)
rcc_dmp=rbind(rcc1r, rcc2r,rcc3r,rcc4r,rcc5r,rcc6r,
              rcc1d, rcc2d, rcc3d, rcc4d, rcc5d, rcc6d)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
library(ggpubr)
summary_se1=summarySE(rcc_m, measurevar = "value", groupvars = "pol")
summary_se2=summarySE(rcc_dm, measurevar = "value", groupvars = "sea")
summary_se3=summarySE(rcc_dmp, measurevar = "value", groupvars = "factor")

# Let's save the summary_se object
#write.csv(summary_se, "./RC_bact.csv")
betadiv=rcc_dmp %>% separate(factor, c("Pol", "Season"))


ppol<-ggerrorplot(rcc_m, x="pol", y="value", desc_stat = "mean_se", 
            error.plot = "errorbar",            # Change error plot type
            add = "mean")
psea<-ggerrorplot(rcc_dm, x="sea", y="value", desc_stat = "mean_se", 
            error.plot = "errorbar",            # Change error plot type
            add = "mean")
ggerrorplot(rcc_dmp, x="factor", y="value", desc_stat = "mean_se", 
            error.plot = "errorbar",            # Change error plot type
            add = "mean")

plot<-ggerrorplot(betadiv, x="Pol", y="value", desc_stat = "mean_se", 
            color="Season",facet.by = "Pol", size = 1,
            error.plot = "errorbar",            # Change error plot type
            add = "mean")+facet_grid(.~Pol, scales = "free")+
  geom_hline(yintercept = 0.95, linetype=2)+
  geom_hline(yintercept = -0.95, linetype=2)+theme_linedraw()+
  theme(#legend.position="none",
        #panel.border = element_blank(), 
        panel.spacing.x = unit(0.3,"line"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+ 
  stat_compare_means(aes(group = Season), label = "p.signif")+
  ylab("Incidence-based beta diversity - Raup-Crick metric (±SE)")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  scale_y_continuous(breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),
                     sec.axis = sec_axis(~.,
                                         breaks = c(-0.95,0.95),
                                         name=" ‘stochastic’ range (−0.95 < βRC < + 0.95)",
                                         labels = c("Deterministic \n more similar", 
                                                    "Deterministic \n dissimilar" )))+
  scale_color_manual(values = c("#FE922AFF","#4490FEFF"),labels=c("Secas", "Lluvias"))
plot
ggsave("Plot/comb-null-oism.png",width = 8, height = 5, dpi = 300, plot = plot, device = "png")

source("Scripts/raup_crick_abundance.R")
library(phyloseq)
library(picante)
library(ape)
rc1da=raup_crick_abundance(abund_table1d, plot_names_in_col1 = F, reps = 999,  set_all_species_equal = F)
rc2da=raup_crick_abundance(abund_table2d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc3da=raup_crick_abundance(abund_table3d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc4da=raup_crick_abundance(abund_table4d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc5da=raup_crick_abundance(abund_table5d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc6da=raup_crick_abundance(abund_table6d, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)

rc1ra=raup_crick_abundance(abund_table1r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc2ra=raup_crick_abundance(abund_table2r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc3ra=raup_crick_abundance(abund_table3r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc4ra=raup_crick_abundance(abund_table4r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc5ra=raup_crick_abundance(abund_table5r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
rc6ra=raup_crick_abundance(abund_table6r, plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)

library(reshape2)
rcc1da=subset(melt(as.matrix(rc1da)),value!=0)
rcc2da=subset(melt(as.matrix(rc2da)),value!=0)
rcc3da=subset(melt(as.matrix(rc3da)),value!=0)
rcc4da=subset(melt(as.matrix(rc4da)),value!=0)
rcc5da=subset(melt(as.matrix(rc5da)),value!=0)
rcc6da=subset(melt(as.matrix(rc6da)),value!=0)

rcc1ra=subset(melt(as.matrix(rc1ra)),value!=0)
rcc2ra=subset(melt(as.matrix(rc2ra)),value!=0)
rcc3ra=subset(melt(as.matrix(rc3ra)),value!=0)
rcc4ra=subset(melt(as.matrix(rc4ra)),value!=0)
rcc5ra=subset(melt(as.matrix(rc5ra)),value!=0)
rcc6ra=subset(melt(as.matrix(rc6ra)),value!=0)

rcc1da$factor="P1-dry"
rcc2da$factor="P2-dry"
rcc3da$factor="P3-dry"
rcc4da$factor="P4-dry"
rcc5da$factor="P5-dry"
rcc6da$factor="P6-dry"

rcc1ra$factor="P1-wet"
rcc2ra$factor="P2-wet"
rcc3ra$factor="P3-wet"
rcc4ra$factor="P4-wet"
rcc5ra$factor="P5-wet"
rcc6ra$factor="P6-wet"
rcc_dmpa=rbind(rcc1ra, rcc2ra,rcc3ra,rcc4ra,rcc5ra,rcc6ra,
              rcc1da, rcc2da, rcc3da, rcc4da, rcc5da, rcc6da)

summary_se4=summarySE(rcc_dmpa, measurevar = "value", groupvars = "factor")

library(tidyverse)
library(ggpubr)
# Let's save the summary_se object
#write.csv(summary_se, "./RC_bact.csv")
betadiva=rcc_dmpa %>% separate(factor, c("Pol", "Season"))
a<-ggerrorplot(betadiva, x="Pol", y="value", desc_stat = "mean_se", 
            color="Season",facet.by = "Pol", size = 1,
            error.plot = "errorbar",            # Change error plot type
            add = "mean")+facet_grid(.~Pol, scales = "free")+
  geom_hline(yintercept = 0.95, linetype=2)+
  geom_hline(yintercept = -0.95, linetype=2)+theme_linedraw()+
  theme(#legend.position="none",
    #panel.border = element_blank(), 
    panel.spacing.x = unit(0.3,"line"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"))+ 
  stat_compare_means(aes(group = Season), label = "p.signif")+
  ylab("Incidence-based beta diversity - Raup-Crick metric (±SE)")+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  scale_y_continuous(breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),
                     sec.axis = sec_axis(~.,
                                         breaks = c(-0.95,0.95),
                                         name=" ‘stochastic’ range (−0.95 < βRC < + 0.95)",
                                         labels = c("Deterministic \n more similar", 
                                                    "Deterministic \n dissimilar" )))
a
ggsave("Plot/comb-null-abund-oism.png",width = 8, height = 5, dpi = 300, plot = a, device = "png")

#nmds con todas
vare.mds=metaMDS(rall, distance = "jaccard", autotransform = T, k=3)

#dfs<-subset(df, row.names(df) %in% sitelist)


data.scores <- as.data.frame(scores(vare.mds,display = "sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
#data.scores$Season <- dfs$Season  #  add the grp variable created earlier
#data.scores$Poligono <- dfs$Poligono
#data.scores$date <- df$Date
#data.scores  #look at the data

data.scores<-merge(data.scores, metadata, by=0)

grp.a <- data.scores[data.scores$Season == "Dry", ][chull(data.scores[data.scores$Season == 
                                                                        "Dry", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$Season == "Wet", ][chull(data.scores[data.scores$Season == 
                                                                        "Wet", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data

# Figure S8 or S9 depending on the OTU data you loaded
nmds_otu<-ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,group=Season, fill=Season),alpha=0.2) +
  scale_fill_manual(values=c("Dry" = "#CC4800", "Wet" = "#3078CB")) +
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Season),size=2) + # add the point markers
  scale_shape_manual(values=c(19,21))+
  theme_bw() + 
  theme(axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=8), # remove x-axis labels
        axis.title.y = element_text(size=8), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "right")

nmds_otu
#ggsave('nmds_otu_season_j6.jpg', width =8, height = 7,  plot =nmds_otu)
tab<- read.csv("/home/yendi/Downloads/example_files/Consensus_table.csv", row.names = 1) %>% dplyr::select(-9:-15)

ver<- hillR::hill_taxa_parti_pairwise(comm =t(tab), q = 0, output = "matrix",  pairs = "full" )
