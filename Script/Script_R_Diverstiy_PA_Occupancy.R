
# Script also pushed on Github: https://github.com/gonzaloprb/Diverstiy-with-Depth

# Install.packages
require(scatterplot3d); require(reshape); require (reshape2); require (data.table); require(RColorBrewer); require(dplyr); require (tidyr); require (plyr); require (ggplot2); require (vegan); require (goeveg);   require  (stats); require (dynRB); require (matrixStats); require(plyr); require (stringr);require (reshape2); require(gridExtra)
require (betapart); require (car); require (MASS); 
require (ggpubr); require (cowplot); require (patchwork)  


# Mid-domain effect
require(devtools)
library(reshape2)
install_github("cran/rangemodelR")
library(rangemodelR)

# Geospatial data
require (geodist)

# rm (list = ls()) 
# Set working directory etc.
setwd("~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/Diverstiy with Depth/Data")
PA_df <- read.csv(file = "photoquad_sppcount_DEEPHOPE_genus_form.csv", header = T, dec = ".", sep = ",", row.names = 1)

PA_df$Island <- gsub("Mangareva", "Gambier", PA_df$Island)
PA_df <- subset (PA_df,Coral_genus!="NA_Coral")

# All stats below come from PA_df. 
# I re-start every time for simplicity 

## Generic Richness profile 


# How many quadrats
nb_genera <- ddply(PA_df, ~ Island + Island_Site + Depth + Coral_genus ,function(x){
  c(nobserv=nrow(x)) }) 


# Transform values bigger than 30, because several "Porites branching" "Porites massive"....
nb_genera$nobserv[nb_genera$nobserv>30] <- 30

# Generic richness, unique number of genera
# length (unique (nb_genera$Coral_genus))
nb_genera <- ddply(nb_genera, ~ Island + Island_Site +  Depth  ,function(x){
  c(nobserv=nrow(unique (x))) })

# Complete rows for all Island, site, Depth, nbgenera = 0
nb_genera <- nb_genera %>% complete( Island,Island_Site,Depth,fill = list(nobserv = 0))

# Transform depth as a Qualitative variable  
nb_genera$Depth <- as.factor(as.character(nb_genera$Depth))
nb_genera$Depth = factor(nb_genera$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# ggplot with islands
ggplot(nb_genera, aes(x=Depth, y=nobserv)) +
  geom_boxplot() + geom_point(aes (colour = Island),size = 1) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Genus richness (n)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

## Mid-domain effect

# Ideally represent with Genus richness profile


all_short <- data.matrix(dcast(PA_df, Coral_genus ~ Depth))
all_short[which(all_short > 0)] = 1

all_short <- t(all_short[,2:7])

rm <- rangemod1d(all_short,cohesion = TRUE,var = NULL,rsize = "observed",reps = 10000)
ric <- rowSums(all_short) 
depth <- c(6,20,40,60,90,120)

plot(depth, ric, ylim = c(2,35), xlab = "Depth", ylab="Richness", type="n", main = "")
lines(depth, ric, lwd = 5, col="blue")
lines(depth, rm$mod.rich, lwd=1)
lines(depth, rm$q2.5, lty = 2)
lines(depth, rm$q97.5, lty = 2)

print (rm) 

## NMDS 
# Question 1: I cannot work with all Quadrats if we keep the PA-Occuppancy "Frequency = Nb Quadrats with genus / Total Nb Quadrats"



resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth + Coral_genus ,function(x){
  c(nobserv=nrow(x)) })

melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)

### Not the PA, the occupancy "Frequency" 
# Measure the Occupancy "Frecuency" f = (nº Quadrats gen (n) present / total n Quadrats)
melt_temp$value[melt_temp$value > 30] <- 30
melt_temp$value <- melt_temp$value / 30

cast_temp = dcast(melt_temp, Depth +  Island + Island_Site ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0

cast_temp$Island <- factor(cast_temp$Island, levels = c("Bora","Makatea","Gambier","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))
cast_database <- cast_temp

# Preparing them for format matrix
cast_temp$ID<- with(cast_temp, paste0(Depth, sep = "_",  Island,  sep = "_", Island_Site))
rownames (cast_temp) <- cast_temp$ID
# Re get database
final <- melt (cast_temp, id = c ("Island","Island_Site","Depth"), measure.vars=c(names(cast_temp[,-c(1,2,3,39)])
), na.rm = FALSE)

cast_temp <- cast_temp[,-c(1,2,3,38)]


# columns where sum is 0
cast_temp <- cast_temp %>% select_if(colSums(.) != 0)
# rows where sum is 0
cast_temp <- cast_temp[rowSums(cast_temp[,])>0, ]

# View (cast_temp)

# As I am working with the occupancy, "Frequency"
PA_NMDS <- metaMDS(cast_temp, k=2, trymax = 1000, distance = "bray") 

# Prepare to plot
treat=c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15)) # Alphabetic order for colours
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

par (mfrow =c(1,1))

ordiplot(PA_NMDS,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(PA_NMDS,display="sites",labels =T)
orditorp(PA_NMDS,display="sites",cex=0.4,air=0.01)
ordispider(PA_NMDS, groups=treat, col = colours,cex=0.6, lwd = 1.5)
ordihull(PA_NMDS,groups=treat,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
orditorp(PA_NMDS,display="species",col="red",air=0.01, cex =0.5)
legend(x="topright", y="top", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  pt.cex=0.8, cex = 0.70, horiz = F, text.width = 0.3, text.font = 20)
title ("Bray distance - Occupancy")
mysubtitle <- paste0("Stress = ", format(round(PA_NMDS$stress, 2)))
mtext(mysubtitle, side=1, line=-2, at=0, adj=0, cex=0.7)


## Beta.pair per depth
# Same question as 1: I cannot work with all Quadrats if we keep the PA-Occuppancy "Frequency = Nb Quadrats with genus / Total Nb Quadrats"


# Resume to keep only coral genus
resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth + Coral_genus ,function(x){
  c(nobserv=nrow(x)) })

# Complete rows for all genera with Coral cover = 0 for all genera
resume_df <- resume_df %>% complete( Island,Island_Site,Depth, Coral_genus,fill = list(nobserv = 0))

# Prepare columns and rows (by now keeping depth)
melt_depth = melt(resume_df, id=c("Depth","Island", "Island_Site","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)

### Not the PA, the occupancy "Frequency"
# Measure the Occupancy "Frecuency" f = (nº Quadrats gen (n) present / total n Quadrats)
melt_depth$value[melt_depth$value > 30] <- 30
melt_depth$value <- melt_depth$value / 30

melt_depth$ID<- with(melt_depth, paste0(Island, sep = "_", Island_Site))

cast_depth = dcast(melt_depth, Depth + Coral_genus ~ ID, mean, add.missing = T)
# cast_depth = dcast(melt_depth, Coral_genus ~ Depth, mean, add.missing = T)

### For 6m  #View(Depth_1) - Necessary rows as sites and columns as species (genera)
Depth_1 <- filter (cast_depth, Depth == 6)
rownames (Depth_1) <- Depth_1$Coral_genus
Depth_1 <- subset (Depth_1, select = - c(Depth, Coral_genus))

Depth_1_Beta <- as.data.frame (t(Depth_1))
#I think, I need to delete columns where all genus are 0 
Depth_1_Beta <- Depth_1_Beta %>% select_if(colSums(.) != 0) 

### Working with Occupancy - "Frequency"
# To do simply beta.pair, necessary to Really working with PA - Jaccard Sorensentransform to 0-1 . Otherwise, beta.pair.abund
# Depth_1_Beta[Depth_1_Beta > 0] <- 1

coral.matrices_Depth_1 <- beta.pair.abund(Depth_1_Beta, index.family = "bray")

mean (coral.matrices_Depth_1$beta.bray.bal) # This is species replacement 
mean (coral.matrices_Depth_1$beta.bray.gra) # This is species nestedness
mean (coral.matrices_Depth_1$beta.bray) # PArtially divided between species replacement and nestedness. Accounts for the two of them. 
# Turnover is higher than nestedness in the overall dissimilarity 


### For 20m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_2 <- filter (cast_depth, Depth == 20)
rownames (Depth_2) <- Depth_2$Coral_genus
Depth_2 <- subset (Depth_2, select = - c(Depth, Coral_genus))

Depth_2_Beta <- as.data.frame (t(Depth_2))
#I think, I need to delete columns where all genus are 0 
Depth_2_Beta <- Depth_2_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_2 <- beta.pair.abund(Depth_2_Beta, index.family = "bray")
mean (coral.matrices_Depth_2$beta.bray)

# So on with the remaining depths

## Betadisper 
# Here is where I am confused...
# I can make it for all depths together, considering groups as different depths. I obtain (1) the average distance to median ("b-dissimilarity" per depth ?), (2) anova and (3)permutest pair-wise differences between (depths)

# Working with all depths, either I work with the mother matrix ## 1 ##  or straight from beta.pair.abund matrix ## 2 ## 

# However, for making a betadisper per depth and check differences between islands is not possible. I need replicates. 
# - Valeriano: you said using all quadrats per depth. However, impossible if we work with the index occupancy-frequency. 

# betadisper using all depths


## 1 ## Betadisper from  the entire mother dissimilarity matrix...
resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth + Coral_genus ,function(x){
  c(nobserv=nrow(x)) })

melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)

### Not the PA, the occupancy "Frequency" 
# Measure the Occupancy "Frecuency" f = (nº Quadrats gen (n) present / total n Quadrats)
melt_temp$value[melt_temp$value > 30] <- 30
melt_temp$value <- melt_temp$value / 30

cast_temp = dcast(melt_temp, Depth +  Island + Island_Site ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0

cast_temp$Island <- factor(cast_temp$Island, levels = c("Bora","Makatea","Gambier","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))
cast_database <- cast_temp

# Preparing them for format matrix
cast_temp$ID<- with(cast_temp, paste0(Depth, sep = "_",  Island,  sep = "_", Island_Site))
rownames (cast_temp) <- cast_temp$ID

cast_temp <- cast_temp[,-c(1,2,3,38)]

# columns where sum is 0
cast_temp <- cast_temp %>% select_if(colSums(.) != 0)
# rows where sum is 0
cast_temp <- cast_temp[rowSums(cast_temp[,])>0, ]

# Distance of the entire dissimilarity matrix
dis <- vegdist(cast_temp)
# dis <- vegdist(cast_depth)
groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

# Calculate multivariate dispersions
mod <- betadisper(dis, groups)
mod

# Perform test
anova(mod)
# Permutation test for F
permutest(mod, pairwise = TRUE, permutations = 99)
# Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod))


## 2 ##  Betadisper from the beta.pair bray distance matrix 

resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth + Coral_genus ,function(x){
  c(nobserv=nrow(x)) })

melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)

### Not the PA, the occupancy "Frequency" 
# Measure the Occupancy "Frecuency" f = (nº Quadrats gen (n) present / total n Quadrats)
melt_temp$value[melt_temp$value > 30] <- 30
melt_temp$value <- melt_temp$value / 30

cast_temp = dcast(melt_temp, Depth +  Island + Island_Site ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0

cast_temp$Island <- factor(cast_temp$Island, levels = c("Bora","Makatea","Gambier","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))

cast_all_depth <- cast_temp
cast_all_depth$ID<- with(cast_all_depth, paste0(Depth, sep = "_",  Island,  sep = "_", Island_Site))
rownames (cast_all_depth) <- cast_all_depth$ID

cast_all_depth <- cast_all_depth[,-c(1,2,3,38)]
# columns where sum is 0
cast_all_depth <- cast_all_depth %>% select_if(colSums(.) != 0)
# rows where sum is 0
cast_all_depth <- cast_all_depth[rowSums(cast_all_depth[,])>0, ]

# Measure pair.abund
coral.matrices_Depth <- beta.pair.abund(cast_all_depth, index.family = "bray")
mean (coral.matrices_Depth$beta.bray)

Beta_Depth <- dist (coral.matrices_Depth$beta.bray)

mod_Depth <- betadisper (Beta_Depth,groups)
mod_Depth

# Perform test
anova(mod_Depth)
# Permutation test for F
permutest(mod_Depth, pairwise = TRUE, permutations = 99)
# Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod_Depth))


# betadisper separate per depths
# As I said above, this really doesn't work. I think I need replicates. Cannot use quadrats as replicates because we compute Occupancy - "Frequency"
# Either working from mother matrix (Depth_1_Beta) or from distance matrix from bray.pair.abund of Depth_1 (=6m)


### Depth_1
beta_Depth_1 <- dist (coral.matrices_Depth_1$beta.bray)
dis1 <- vegdist(Depth_1_Beta, method = "bray")

# Only with Depth 1
groups <- rownames (Depth_1_Beta)
# groups <- c(rep("6M",16))

mod1 <- betadisper (dis1,groups)
# mod1 <- betadisper (beta_Depth_1,groups)
mod1
# Perform test
anova(mod1)
# Permutation test for F
permutest(mod1, pairwise = TRUE, permutations = 99)

## Tukey's Honest Significant Differences
# (mod.HSD <- TukeyHSD(mod1))


### Depth_2
beta_Depth_2 <- dist (coral.matrices_Depth_2$beta.bray)
dis2 <- vegdist(Depth_2_Beta)

# Only with Depth 1
groups <- rownames (Depth_2_Beta)

mod2 <- betadisper (beta_Depth_2,groups, )
# mod2

str(beta_Depth_1)

library (usedist)

beta_Depth_1 <- dist_setNames(beta_Depth_1, paste0 ("6",sep = "_",rownames (Depth_1_Beta)))

beta_Depth_2 <- dist_setNames(beta_Depth_2, paste0 ("20",sep = "_",rownames (Depth_2_Beta)))


## Mantel tests with distance - per depth

# I can also plot bray-distance according to vertical depth distance measuring bray-distance per site and not per depth


# Mantel tests 
# It needs Beta-dissimilarity matrix and matrix distance

library (geodist)

Locations <- read.csv(file = "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/GIS_MAP/Deephope_sampling_locations_RAN*.csv", header = T, dec = ".", sep = ";", row.names = 1)
Locations$Island <- gsub("Mangareva", "Gambier", Locations$Island)
Locations$Island <- gsub("Bora Bora", "Bora", Locations$Island)

Locations$ID<- with(Locations, paste0(Island, sep = "_", Site))

# Keep only same ID as coral data 
unique (melt_depth$ID)
Locations <- Locations[Locations$ID %in% melt_depth$ID, ]
Locations <- Locations[order(Locations[,'ID']), ]

rownames (Locations) <- Locations$ID
Locations <- subset (Locations, select =  c(Latitude, Longitude))

dm <- geodist (Locations, measure = "geodesic", paired = T) /1000

dm <- as.data.frame(dm)
rownames (dm) <-  rownames (Locations)
names (dm) <- rownames (Locations) 

# Transform to object of class dist
dm <- dist (dm)

##### Make the mantel test ##### 
# For Depth 1 - 6m 
mantel(coral.matrices_Depth_1$beta.bray, dm)
mantel(coral.matrices_Depth_2$beta.bray, dm)

# and so on!
