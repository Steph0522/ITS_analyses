
# Assessing environmental conditions and selecting structuring variables
library(reshape2)
library(picante)
library(readxl)
library(tidyverse)
dry<- read_excel("Data/Fisicoquimicos.xlsx", sheet = "metadata_secas_all") %>% mutate(Season="Dry")
wet<- read_excel("Data/Fisicoquimicos.xlsx", sheet = "metadata_lluvias_all") %>% mutate(Season="Wet")

drys<- dry %>% dplyr::select( Poligono,Sitio,Transecto,pH,MO,N,P,Season)
wets<- wet %>% dplyr::select(Poligono,Sitio,Transecto,pH,MO,N,P,Season)

env_season<- rbind(drys,wets) %>% mutate(sea=c(rep("D",36), rep("R",36))) %>% unite("SampleID", c("Poligono", "Sitio", "Transecto", "sea"),remove = F, sep = "")

df <-env_season  %>% column_to_rownames(var="SampleID") %>% dplyr::select(pH,MO,N,P)
dfs=data.frame(scale(df, center = T, scale = T)) #standardize environmental data

dfs$Poligono<-NA
dfs$Poligono<-env_season$Poligono
melted_df=melt(dfs, id="Poligono")
mean_table=Rmisc::summarySE(melted_df, measurevar = "value", groupvars=c("Poligono","variable"), na.rm = T)

# checking collinearity
cor.table(na.omit(dfs[,1:4]))

# RDA separately for 16S and 18S data
library(vegan)
library(qiime2R)
spp.ITS=read_qza("Data/filt_table_ojis.qza")$data
colnames(spp.ITS)<-stringr::str_extract(colnames(spp.ITS), '[^_]+(?=_............\\w[^_]*$)') 
env.ITS=env_season
spp.ITS=data.frame(t(spp.ITS), check.names = F)

mm=cbind(spp.ITS,env.ITS)
mm=na.omit(mm)
dim(mm)
head(mm)

env.ITS=mm[,14616:14619]
spp.ITS=mm[,1:14611]

dca=decorana(spp.ITS)
plot(dca) # RDA is okay to use

#select variables that did not show collinearity
spp.ITS_hell=decostand(spp.ITS, "hell") # Hellinger transformation
env.ITS_st=data.frame(scale(env.ITS, scale=T, center=F)) # standardize env. data

PC1=rda(spp.ITS_hell~., env.ITS_st)
plot(PC1, display = "sites" ,type = c("text"))
plot(PC1, type="n")
text(PC1, dis="cn", lwd=1, cex=1)
text(PC1, dis="sites", type=c("text"), cex=0.8)
points(PC1, dis="sites", lwd=1, cex=0.5)

#quartz.save("RDA0_16S.pdf", type="pdf")

# Forward selection procedure
PC1=rda(spp.ITS~., env.ITS)
plot(PC1, type="n")
text(PC1, dis="cn", lwd=2.5, cex=1.2)
text(PC1, dis="sites", type=c("text"), cex=0.8)
points(PC1, dis="sites", lwd=2.5)

cap.env=capscale(spp.ITS_hell~.,env.ITS_st,distance = "bray")
cap.env

mod0.env=capscale(spp.ITS_hell~1,env.ITS_st,distance = "bray")
step.env=ordistep(mod0.env,scope = formula(cap.env))
step.env$anova # These are the selected structuring variables that are included in the subsequent analyses for 16S


# Plot the temporal trends of the selected environmental variables over time
# Figure 2
#setwd("/Volumes/FREECOM HDD/Phd_project/R_copy/")
df <- env.ITS[,c("pH", "MO", "N", "P")]
df$Poligono<-NA
df$Poligono<-env_season$Poligono
melted_df=melt(df, id="Poligono")

df$Season<-NA
df$Season<-env_season$Season
melted_df=melt(df, id=c("Poligono", "Season"))

mean=Rmisc::summarySE(melted_df, measurevar = "value", groupvars=c("Poligono","variable", "Season"), na.rm = T)
mean$variable=factor(mean$variable, levels = c("pH", "MO", "N", "P"))

conservation_status <- c(
  pH = "pH",
  MO = "OM (%)",
  N = "N (mg/Kg)",
  P= "P (mg/Kg)")

library(plyr)
mean$conservation2 <- plyr::revalue(mean$variable, conservation_status)
mean$Poligon<- factor(mean$Poligono, levels = 1:6, labels=1:6)

melted_df$conservation2 <- plyr::revalue(melted_df$variable, conservation_status)
melted_df$Poligon<- factor(melted_df$Poligono, levels = 1:6, labels=1:6)

library(ggplot2)
library(ggpubr)
fq<-ggline(melted_df, x="Poligon", y="value", shape = "Season", 
       linetype="Season", facet.by = "conservation2", add = "mean_sd")+
  facet_grid(conservation2~., scales = "free_y")+labs(x="Poligon",y="Mean±SD")+
  theme_bw()+
    theme(legend.position="top")+
  theme(axis.text = element_text(size=8),
        plot.title = element_text(size=8, face="bold"))

ggsave('fisicoq-season.png', width =4, height = 7, dpi = 300, plot =fq)
#quartz.save("meanplot.pdf", type = "pdf") 
#quartz.save("meanplot.tiff", type = "tiff")

# It is clear from the Figure 2 plot that there are two distinct periods.
#But are they statistically different?
#### non-parametric test between dry and wet periods
# Table S1 in Supplementary material
library(car)
library(pgirmess)
env_season$Poligono<- as.factor(env_season$Poligono)

k.ph.s<-round(kruskal.test(env_season$pH~env_season$Season)$p.value, digits = 3)
k.ph.p<-round(kruskal.test(env_season$pH~env_season$Poligono)$p.value, digits = 3)
leveneTest(env_season$pH~env_season$Season)
leveneTest(env_season$pH~env_season$Poligono)

k.mo.s<-round(kruskal.test(env_season$MO~env_season$Season)$p.value, 3)
k.mo.p<-round(kruskal.test(env_season$MO~env_season$Poligono)$p.value,3)
leveneTest(env_season$MO~env_season$Season)
leveneTest(env_season$MO~env_season$Poligono)

k.n.s<-round(kruskal.test(env_season$N~env_season$Season)$p.value,3)
k.n.p<-round(kruskal.test(env_season$N~env_season$Poligono)$p.value,3)
leveneTest(env_season$N~env_season$Season)
leveneTest(env_season$N~env_season$Poligono)

k.p.s<-round(kruskal.test(env_season$P~env_season$Season)$p.value,3)
k.p.p<-round(kruskal.test(env_season$P~env_season$Poligono)$p.value,3)
leveneTest(env_season$P~env_season$Season)
leveneTest(env_season$P~env_season$Poligono)



lm.ph<-lm(pH~Season*Poligono ,data = env_season)
lm.ph.perm<-PermTest(lm.ph)
lm.ph.perm


lm.mo<-lm(MO~Season*Poligono ,data = env_season)
lm.mo.perm<-PermTest(lm.mo)
lm.mo.perm


lm.n<-lm(N~Season*Poligono ,data = env_season)
lm.n.perm<-PermTest(lm.n)
lm.n.perm


lm.p<-lm(P~Season*Poligono ,data = env_season)
lm.p.perm<-PermTest(lm.p)
lm.p.perm


all.lm<- cbind(lm.ph.perm$resultats, lm.mo.perm$resultats, lm.n.perm$resultats, lm.p.perm$resultats)
colnames(all.lm)<-c("pH", "OM", "N", "P")

all.kw<- data.frame(pH=c(k.ph.s,k.ph.p),
                    OM=c(k.mo.s, k.mo.p),
                    N=c(k.n.s, k.n.p),
                    P=c(k.p.s, k.p.p))
rownames(all.kw)<-c("Season", "Poligon")

write.csv(all.lm,"Data/all.lm.csv")
write.csv(all.kw,"Data/all.kw.csv")

#NMDS; assessing dry-wet periods 
library(ggfortify)
library(cluster)
library(vegan)

# NMDS plot showing the changes in environmental conditions in relation to the dry
# and wet period based on Euclidean distances




df <- env_season %>% column_to_rownames(var = "SampleID") %>% dplyr::select("pH","MO","N","P","Season","Poligono")
df=na.omit(df)

dff=df %>% dplyr::select("pH","MO","N","P")
dff=scale(dff, scale = T, center = F)
dff <- dff[,colSums(is.na(dff))<nrow(dff)]


dffd = (vegdist(dff, "euc"))
as<-anova(betadisper(dffd, df$Season)) 
ap<-anova(betadisper(dffd, df$Poligono)) 

write.csv(as, "Data/permisp_env_season.csv")
write.csv(ap, "Data/permisp_env_pol.csv")


adonis_location = adonis2(dff ~ df$Season*df$Poligono, method = "euc") #PERMANOVA
adonis_location

#write.csv(adonis_location, "Data/adonis_env.csv")
#dff1<- dff %>%as.data.frame() %>%  dplyr::select_at(vars(starts_with("1"))) %>% as.matrix()
dff<- dff %>%as.data.frame() 

dff = as.matrix((vegdist(dff, "euc")))
vare.mds=metaMDS(dff, distance = "euc", autotransform = F, na.rm=T, k=3)
data.scores <- as.data.frame(scores(vare.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$Season <- df$Season  #  add the grp variable created earlier
data.scores$Poligono<- df$Poligono
data.scores  #look at the data


grp.a <- data.scores[data.scores$Season == "Dry", ][chull(data.scores[data.scores$Season== 
                                                                        "Dry", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$Season == "Wet", ][chull(data.scores[data.scores$Season== 
                                                                        "Wet", c("NMDS1", "NMDS2")]), ]  # hull values for grp B


grp.1 <- data.scores[data.scores$Poligono == "1", ][chull(data.scores[data.scores$Poligono== 
                                                                        "1", c("NMDS1", "NMDS2")]), ]  # hull values for grp A

grp.2 <- data.scores[data.scores$Poligono == "2", ][chull(data.scores[data.scores$Poligono== 
                                                                        "2", c("NMDS1", "NMDS2")]), ] 


grp.3<- data.scores[data.scores$Poligono == "3", ][chull(data.scores[data.scores$Poligono== 
                                                                        "3", c("NMDS1", "NMDS2")]), ] 
grp.4 <- data.scores[data.scores$Poligono == "4", ][chull(data.scores[data.scores$Poligono== 
                                                                        "4", c("NMDS1", "NMDS2")]), ] 

grp.5 <- data.scores[data.scores$Poligono == "5", ][chull(data.scores[data.scores$Poligono== 
                                                                        "5", c("NMDS1", "NMDS2")]), ] 

grp.6 <- data.scores[data.scores$Poligono == "6", ][chull(data.scores[data.scores$Poligono== 
                                                                        "6", c("NMDS1", "NMDS2")]), ] 

hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
#hull.data <- rbind(grp.1, grp.2, grp.3, grp.4, grp.5, grp.6)  #combine grp.a and grp.b


hull.data

# Figure S7 in the Supplementary material
nmds_env<-ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,group=Season, fill=Season),alpha=0.2) +
  scale_fill_manual(values=c("Dry" = "#CC4800", "Wet" = "#3078CB")) +
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=Season),size=2) + # add the point markers
  scale_shape_manual(values=c(19,21))+
  theme_bw() + 
  theme(axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=10), # remove x-axis labels
        axis.title.y = element_text(size=10), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "right")# + annotate("text", x = 5, y = -4.5, label = "stress: 0.056")
#quartz.save("ENV_DIS.pdf", type="pdf") 
nmds_env

ggsave('nmds_env_season.jpg', width =8, height = 7,  plot =nmds_env)



# NMDS plot showing the changes in bacterial or microeukaryotic communities in relation to the dry and wet
# period based on Bray-Curtis dissimilarities 
#setwd("/Volumes/FREECOM HDD/Phd_project/R_copy")
otu.ITS=read_qza("Data/filt_table_ojis.qza")$data # Change to "OTU16S_rar.csv" if you want to assess the bacterial dataset
colnames(otu.ITS)<-stringr::str_extract(colnames(otu.ITS), '[^_]+(?=_............\\w[^_]*$)')
otu.ITS<- otu.ITS[,match(rownames(df), colnames(otu.ITS))]

otu.ITS=data.frame(t(otu.ITS))
otulist=row.names(otu.ITS)

otu.ITS_hell=decostand(otu.ITS, "hell")
sitelist=row.names(df)
sitelist<- sitelist[str_detect(sitelist, "^6")]

#otherwise run from this line:
#otu.ITS_hell=subset(otu.ITS_hell, row.names(otu.ITS_hell) %in% sitelist)
otu.ITS= subset(otu.ITS, row.names(otu.ITS) %in% sitelist)

vare.mds=metaMDS(otu.ITS, distance = "jaccard", autotransform = F, k=3)

dfs<-subset(df, row.names(df) %in% sitelist)

data.scores <- as.data.frame(scores(vare.mds,display = "sites"))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$Season <- dfs$Season  #  add the grp variable created earlier
data.scores$Poligono <- dfs$Poligono
#data.scores$date <- df$Date
data.scores  #look at the data

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
ggsave('nmds_otu_season_j6.jpg', width =8, height = 7,  plot =nmds_otu)

#P20+ annotate("text", x = -1.3, y = -1.5, label = "stress: 0.144")

#quartz.save("OTU18S_DIS.pdf", type="pdf")

# PERMANOVA and PERMDISP using Bray-Curtis distances for the bacterial or microeukaryotic dataset
# Table S2 and S3 in Supplementary material

adonis_location = adonis2(otu.ITS_hell ~ df$Season*df$Poligono, method = "bray")
adonis_location#PERMANOVA
dffd = (vegdist(otu.ITS_hell, "bray"))
bs<-anova(betadisper(dffd, df$Season)) #PERMDISP
bp<-anova(betadisper(dffd, df$Poligono)) #PERMDISP


ggsave('nmds_otu_season_bray.jpg', width =8, height = 7,  plot =nmds_otu)
write.csv(adonis_location, "Data/adonis_otu_bray.csv")
write.csv(bs, "Data/permdis_otu_season_bray.csv")
write.csv(bp, "Data/permdis_otu_poligon_bray.csv")
