
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
setwd("~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Data")
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

rm <- rangemod1d(all_short,cohesion = TRUE,var = NULL,rsize = "observed",reps = 10000, degen = T)
ric <- rowSums(all_short) 
depth <- c(6,20,40,60,90,120)

plot(depth, ric, ylim = c(2,35), xlab = "Depth", ylab="Richness", type="n", main = "")
lines(depth, ric, lwd = 5, col="blue")
lines(depth, rm$mod.rich, lwd=1)
lines(depth, rm$q2.5, lty = 2)
lines(depth, rm$q97.5, lty = 2)

print (rm) 


exp_mde <- rm$out.df

# Two ways of counting the generic richness diversity, Make it per island
nb_genera_all <- ddply(PA_df, ~ Depth + Coral_genus  ,function(x){
  c(nobserv=nrow(unique (x))) })
count <- ddply(nb_genera_all, ~ Depth ,function(x){
  c(nobserv=nrow(unique (x))) })

# Or straight
agg <- aggregate(Coral_genus ~ Depth, data=subset(PA_df), FUN=function(x) length(unique(x)))

exp_mde$obs <- count$nobserv
# VAlues are nearly the same

# Not sure if I have to make the model per island, archipalgo, island_site, etc. What is clear is that the expected results are 
# very close to the observed ones




# Mid-domain effect per island

# For each island
MDE_Island <- lapply(unique(PA_df$Island), function (x) {
  
  Island<- PA_df[PA_df$Island == x,]
  
  island_short <- data.matrix(dcast(Island, Coral_genus ~ Depth))
  island_short[which(island_short > 0)] = 1
  
  island_short <- t(island_short[,2:7])
  
  rm <- rangemod1d(island_short,cohesion = TRUE,var = NULL,rsize = "observed",reps = 10000)
  ric <- rowSums(island_short) 
  depth <- c(6,20,40,60,90,120)
  
  
  plot(depth, ric, ylim = c(2,32), xlab = "Depth", ylab="Richness", type="n", main = x)
  lines(depth, ric, lwd = 5, col="blue")
  lines(depth, rm$mod.rich, lwd=1)
  lines(depth, rm$q2.5, lty = 2)
  lines(depth, rm$q97.5, lty = 2)
  
  return(rm)
  return(ric)
  return(depth)
  
})

# MDE for different islands
MDE_Bora <- MDE_Island[[1]]
MDE_Makatea <- MDE_Island[[2]]
MDE_Gambier <- MDE_Island[[3]]
MDE_Moorea <- MDE_Island[[4]]
MDE_Rangiroa <- MDE_Island[[5]]
MDE_Raroia <- MDE_Island[[6]]
MDE_Tahiti <- MDE_Island[[7]]
MDE_Tikehau <- MDE_Island[[8]]


# Two ways of counting the generic richness diversity, - per island
nb_genera_all <- ddply(PA_df, ~ Island + Depth + Coral_genus  ,function(x){
  c(nobserv=nrow(unique (x))) })
count <- ddply(nb_genera_all, ~ Island + Depth ,function(x){
  c(nobserv=nrow(unique (x))) })

# Or straight
agg <- aggregate(Coral_genus ~ Island + Depth, data=subset(PA_df), FUN=function(x) length(unique(x)))

# Make a plot of expected results vs observations
exp_mde$obs <- count$nobserv




############ PA (0 or 1) ###################

# Working with PA only 

## NMDS 
# Question 1: I cannot work with all Quadrats if we keep the PA-Occuppancy "Frequency = Nb Quadrats with genus / Total Nb Quadrats"

resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth + Coral_genus ,function(x){
  c(nobserv=nrow(x)) })

melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)

### PA
melt_temp$value[melt_temp$value > 1] <- 1

cast_temp = dcast(melt_temp, Depth +  Island + Island_Site ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0

cast_temp$Island <- factor(cast_temp$Island, levels = c("Bora","Gambier","Makatea","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))
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
PA_NMDS <- metaMDS(cast_temp, k=2, trymax = 1000, distance = "jaccard") 

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
title ("Jaccard distance - PA")
mysubtitle <- paste0("Stress = ", format(round(PA_NMDS$stress, 2)))
mtext(mysubtitle, side=1, line=-2, at=0, adj=0, cex=0.7)


## Beta.pair per depth


# Resume to keep only coral genus
resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth + Coral_genus ,function(x){
  c(nobserv=nrow(x)) })

# Complete rows for all genera with Coral cover = 0 for all genera
resume_df <- resume_df %>% complete( Island,Island_Site,Depth, Coral_genus,fill = list(nobserv = 0))

# Prepare columns and rows (by now keeping depth)
melt_depth = melt(resume_df, id=c("Depth","Island", "Island_Site","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)

### PA (0 or 1)
melt_depth$value[melt_depth$value > 1] <- 1


cast_temp = dcast(melt_temp, Depth +  Island + Island_Site ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0


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

### Working with PA
# To do simply beta.pair, necessary to Really working with PA - Jaccard Sorensentransform to 0-1 . Otherwise, beta.pair.abund
# Depth_1_Beta[Depth_1_Beta > 0] <- 1

coral.matrices_Depth_1 <- beta.pair(Depth_1_Beta, index.family = "jaccard")
# coral.matrices_Depth_1 <- beta.pair(Depth_1_Beta, index.family = "sorensen")

mean (coral.matrices_Depth_1$beta.jtu) # This is species replacement 
mean (coral.matrices_Depth_1$beta.jne) # This is species nestedness
mean (coral.matrices_Depth_1$beta.jac)
# Turnover is higher than nestedness in the overall dissimilarity 


### For 20m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_2 <- filter (cast_depth, Depth == 20)
rownames (Depth_2) <- Depth_2$Coral_genus
Depth_2 <- subset (Depth_2, select = - c(Depth, Coral_genus))

Depth_2_Beta <- as.data.frame (t(Depth_2))
#I think, I need to delete columns where all genus are 0 
Depth_2_Beta <- Depth_2_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_2 <- beta.pair(Depth_2_Beta, index.family = "jaccard")
mean (coral.matrices_Depth_2$beta.jac)

### For 40m  #View(Depth_3) - Necessary rows as sites and columns as species (genera)
Depth_3 <- filter (cast_depth, Depth == 40)
rownames (Depth_3) <- Depth_3$Coral_genus
Depth_3 <- subset (Depth_3, select = - c(Depth, Coral_genus))

Depth_3_Beta <- as.data.frame (t(Depth_3))
#I think, I need to delete columns where all genus are 0 
Depth_3_Beta <- Depth_3_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_3 <- beta.pair(Depth_3_Beta, index.family = "jaccard")
mean (coral.matrices_Depth_3$beta.jac)

### For 60m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_4 <- filter (cast_depth, Depth == 60)
rownames (Depth_4) <- Depth_4$Coral_genus
Depth_4 <- subset (Depth_4, select = - c(Depth, Coral_genus))

Depth_4_Beta <- as.data.frame (t(Depth_4))
#I think, I need to delete columns where all genus are 0 
Depth_4_Beta <- Depth_4_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_4 <- beta.pair(Depth_4_Beta, index.family = "jaccard")
mean (coral.matrices_Depth_4$beta.jac)

### For 90m  #View(Depth_5) - Necessary rows as sites and columns as species (genera)
Depth_5 <- filter (cast_depth, Depth == 90)
rownames (Depth_5) <- Depth_5$Coral_genus
Depth_5 <- subset (Depth_5, select = - c(Depth, Coral_genus))

Depth_5_Beta <- as.data.frame (t(Depth_5))
#I think, I need to delete columns where all genus are 0 
Depth_5_Beta <- Depth_5_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_5 <- beta.pair(Depth_5_Beta, index.family = "jaccard")
mean (coral.matrices_Depth_5$beta.jac)

### For 120m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_6 <- filter (cast_depth, Depth == 120)
rownames (Depth_6) <- Depth_6$Coral_genus
Depth_6 <- subset (Depth_6, select = - c(Depth, Coral_genus))

Depth_6_Beta <- as.data.frame (t(Depth_6))
#I think, I need to delete columns where all genus are 0 
Depth_6_Beta <- Depth_6_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_6 <- beta.pair(Depth_6_Beta, index.family = "jaccard")
mean (coral.matrices_Depth_6$beta.jac)


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
mantel(coral.matrices_Depth_1$beta.jac, dm) # 6m
mantel(coral.matrices_Depth_2$beta.jac, dm) # 20m
mantel(coral.matrices_Depth_3$beta.jac, dm) # 40m 
mantel(coral.matrices_Depth_4$beta.jac, dm) # 60m 
mantel(coral.matrices_Depth_5$beta.jac, dm) #90m #
mantel(coral.matrices_Depth_6$beta.jac, dm) #120m # No significance
# This means that the beta.jac distance is not correlated to distance


### Betadisper 

# I can make it for all depths together, considering groups as different depths. I obtain (1) the average distance to median ("b-dissimilarity" per depth ?), (2) anova and (3)permutest pair-wise differences between (depths)

# 1st betadisper using all depths
##   Betadisper from the beta.pair bray distance matrix 

resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth + Coral_genus ,function(x){
  c(nobserv=nrow(x)) })

melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)

### PA (0 or 1)
melt_temp$value[melt_temp$value > 1] <- 1


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
coral.matrices_Depth <- beta.pair(cast_all_depth, index.family = "jaccard")
mean (coral.matrices_Depth$beta.jac)

# Beta_Depth <- dist (coral.matrices_Depth$beta.jac)
Beta_Depth <- coral.matrices_Depth$beta.jac

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

mod_Depth <- betadisper (Beta_Depth,groups)
mod_Depth

plot (mod_Depth)
boxplot (mod_Depth)

# Perform test
anova(mod_Depth)
# Permutation test for F
permutest(mod_Depth, pairwise = TRUE, permutations = 99)
# Tukey's Honest Significant Differences
plot (mod.HSD <- TukeyHSD(mod_Depth))


distances_centroid <- as.data.frame(mod_Depth$distances)
colnames (distances_centroid) <- "Distances"
distances_centroid$Depth <- sapply(strsplit(rownames(distances_centroid), "_"), "[", 1)
mean <- aggregate (Distances ~ Depth, distances_centroid,"mean")
sd <- aggregate (Distances ~ Depth, distances_centroid,"sd")
mean$sd <- sd [2]
mean$Distances<- round(mean$Distances, digits = 3)
mean$sd<- round(mean$sd, digits = 3)

# Or 
with(mod_Depth, tapply(distances, groups, "mean"))
with(mod_Depth, tapply(distances, groups, "median"))

distances_centroid$Depth = factor(distances_centroid$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median") + xlab ("Depth (m)") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 


####################################################################################


# ############### "Frequency = Nb Quadrats with genus / Total Nb Quadrats" ###################
## NMDS 

resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth + Coral_genus ,function(x){
  c(nobserv=nrow(x)) })

melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)

### Not the PA, the occupancy "Frequency" 
# Measure the Occupancy "Frecuency" f = (nº Quadrats gen (n) present / total n Quadrats)
melt_temp$value[melt_temp$value > 30] <- 30
melt_temp$value <- melt_temp$value / 30

cast_temp = dcast(melt_temp, Depth +  Island + Island_Site ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0

cast_temp$Island <- factor(cast_temp$Island, levels = c("Bora","Gambier","Makatea","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))
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

### For 40m  #View(Depth_3) - Necessary rows as sites and columns as species (genera)
Depth_3 <- filter (cast_depth, Depth == 40)
rownames (Depth_3) <- Depth_3$Coral_genus
Depth_3 <- subset (Depth_3, select = - c(Depth, Coral_genus))

Depth_3_Beta <- as.data.frame (t(Depth_3))
#I think, I need to delete columns where all genus are 0 
Depth_3_Beta <- Depth_3_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_3 <- beta.pair.abund(Depth_3_Beta, index.family = "bray")
mean (coral.matrices_Depth_3$beta.bray)

### For 60m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_4 <- filter (cast_depth, Depth == 60)
rownames (Depth_4) <- Depth_4$Coral_genus
Depth_4 <- subset (Depth_4, select = - c(Depth, Coral_genus))

Depth_4_Beta <- as.data.frame (t(Depth_4))
#I think, I need to delete columns where all genus are 0 
Depth_4_Beta <- Depth_4_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_4 <- beta.pair.abund(Depth_4_Beta, index.family = "bray")
mean (coral.matrices_Depth_4$beta.bray)

### For 90m  #View(Depth_5) - Necessary rows as sites and columns as species (genera)
Depth_5 <- filter (cast_depth, Depth == 90)
rownames (Depth_5) <- Depth_5$Coral_genus
Depth_5 <- subset (Depth_5, select = - c(Depth, Coral_genus))

Depth_5_Beta <- as.data.frame (t(Depth_5))
#I think, I need to delete columns where all genus are 0 
Depth_5_Beta <- Depth_5_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_5 <- beta.pair.abund(Depth_5_Beta, index.family = "bray")
mean (coral.matrices_Depth_5$beta.bray)

### For 120m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_6 <- filter (cast_depth, Depth == 120)
rownames (Depth_6) <- Depth_6$Coral_genus
Depth_6 <- subset (Depth_6, select = - c(Depth, Coral_genus))

Depth_6_Beta <- as.data.frame (t(Depth_6))
#I think, I need to delete columns where all genus are 0 
Depth_6_Beta <- Depth_6_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_6 <- beta.pair.abund(Depth_6_Beta, index.family = "bray")
mean (coral.matrices_Depth_6$beta.bray)


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
mantel(coral.matrices_Depth_1$beta.bray, dm) # 6m
mantel(coral.matrices_Depth_2$beta.bray, dm) # 20m
mantel(coral.matrices_Depth_3$beta.bray, dm) # 40m # Significance close to 0, it means no correlation
mantel(coral.matrices_Depth_4$beta.bray, dm) # 60m # Significance distante to 0, it means correlation
mantel(coral.matrices_Depth_5$beta.bray, dm) #90m # Significance close to 0, it means no correlation
mantel(coral.matrices_Depth_6$beta.bray, dm) #120m # No significance




### Betadisper 

# I can make it for all depths together, considering groups as different depths. I obtain (1) the average distance to median ("b-dissimilarity" per depth ?), (2) anova and (3)permutest pair-wise differences between (depths)

# 1st betadisper using all depths
##   Betadisper from the beta.pair bray distance matrix 

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

# Beta_Depth <- dist (coral.matrices_Depth$beta.bray)
Beta_Depth <- coral.matrices_Depth$beta.bray

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

mod_Depth <- betadisper (Beta_Depth,groups)
mod_Depth

plot (mod_Depth)
boxplot (mod_Depth)

# Perform test
anova(mod_Depth)
# Permutation test for F
permutest(mod_Depth, pairwise = TRUE, permutations = 99)
# Tukey's Honest Significant Differences
plot (mod.HSD <- TukeyHSD(mod_Depth))


distances_centroid <- as.data.frame(mod_Depth$distances)
colnames (distances_centroid) <- "Distances"
distances_centroid$Depth <- sapply(strsplit(rownames(distances_centroid), "_"), "[", 1)
mean <- aggregate (Distances ~ Depth, distances_centroid,"mean")
sd <- aggregate (Distances ~ Depth, distances_centroid,"sd")
mean$sd <- sd [2]
mean$Distances<- round(mean$Distances, digits = 3)
mean$sd<- round(mean$sd, digits = 3)

# Or 
with(mod_Depth, tapply(distances, groups, "mean"))
with(mod_Depth, tapply(distances, groups, "median"))

distances_centroid$Depth = factor(distances_centroid$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median") + xlab ("Depth (m)") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 





