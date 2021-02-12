
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

par(mfrow=c(2,4))

# For each island
MDE_Island <- lapply(unique(PA_df$Island), function (x) {
  
  Island<- PA_df[PA_df$Island == x,]
  
  island_short <- data.matrix(dcast(Island, Coral_genus ~ Depth))
  island_short[which(island_short > 0)] = 1
  
  island_short <- t(island_short[,2:7])
  
  rm <- rangemod1d(island_short,cohesion = TRUE,var = NULL,rsize = "observed",reps = 10000)
  ric <- rowSums(island_short) 
  depth <- c(6,20,40,60,90,120)
  
  
  plot(depth, ric, ylim = c(2,32), xlab = "Depth", ylab="Generic richness", type="n", main = x)
  lines(depth, ric, lwd = 5, col="blue")
  lines(depth, rm$mod.rich, lwd=1)
  lines(depth, rm$q2.5, lty = 2)
  lines(depth, rm$q97.5, lty = 2)

  
  
  return(rm)
  return(ric)
  return(depth)
  
})

# Two ways of counting the generic richness diversity, - per island
nb_genera_all <- ddply(PA_df, ~ Island + Depth + Coral_genus  ,function(x){
  c(nobserv=nrow(unique (x))) })
count <- ddply(nb_genera_all, ~ Island + Depth ,function(x){
  c(nobserv=nrow(unique (x))) })

# Or straight
agg <- aggregate(Coral_genus ~ Island + Depth, data=subset(PA_df), FUN=function(x) length(unique(x)))
colnames (agg) [3]<- "Richness"

# MDE for different islands from the generated list
MDE_Bora <- cbind (MDE_Island[[1]], subset(agg, Island == 'Bora'))
MDE_Makatea <- cbind (MDE_Island[[2]], subset(agg, Island == 'Makatea'))
MDE_Gambier <- cbind (MDE_Island[[3]], subset(agg, Island == 'Gambier'))
MDE_Moorea <- cbind (MDE_Island[[4]], subset(agg, Island == 'Moorea'))
MDE_Rangiroa <- cbind (MDE_Island[[5]], subset(agg, Island == 'Rangiroa'))
MDE_Raroia <- cbind (MDE_Island[[6]], subset(agg, Island == 'Raroia'))
MDE_Tahiti <- cbind (MDE_Island[[7]], subset(agg, Island == 'Tahiti'))
MDE_Tikehau <- cbind (MDE_Island[[8]], subset(agg, Island == 'Tikehau'))

MDE_all <- rbind (MDE_Bora,MDE_Makatea,MDE_Gambier,MDE_Moorea,MDE_Rangiroa,MDE_Raroia,MDE_Tahiti,MDE_Tikehau)

# Maybe use some tests here? 
cor.test(MDE_all$mod.rich, MDE_all$Richness,  method = "pearson", use = "complete.obs")

RMSE = function(obs, pred){
  sqrt(mean((obs - pred)^2))
}

pred <- MDE_all$mod.rich
obs <- MDE_all$Richness

RMSE (obs, pred)
bias <- mean (pred)-mean (obs)


## ggplot 
ggplot(MDE_all, aes(x=Depth, y=Richness)) + 
  geom_point(aes(y=Richness,colour = Island),size = 0.8, shape = 19, alpha = 0.8) +
 # geom_line(aes(y=Richness,colour = Island),size = 0.5,linetype="solid", alpha = 0.8)  + 
  geom_line(aes(y=mod.rich,colour = Island),size = 0.8,linetype="dashed", alpha = 0.8)  +  
  geom_line(aes(y=q2.5,colour = Island), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_line(aes(y=q97.5,colour = Island), size=0.3, linetype="dotted",alpha = 0.8) +
  stat_summary(fun=mean, geom="point", shape=18, color="black", size=4) + 
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", color="black", size=0.5, width=2) +
  scale_x_continuous(name ="Depth (m)", limits=c(3,122), breaks = c(6,20,40,60,90,120)) +
  scale_y_continuous(name ="Generic richness (%)", limits=c(0,30), breaks = c(0,10,20,30)) +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 


# Separate per islands
ggplot(MDE_all, aes(x=Depth, y=Richness)) + 
  geom_point(aes(y=Richness),size = 1, shape = 19) +
  geom_line(aes(y=Richness),size = 1,linetype="solid", alpha = 0.8, colour = "blue")  + 
  geom_line(aes(y=mod.rich),size = 1,linetype="dotted", alpha = 1)  +  
  geom_line(aes(y=q2.5), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_line(aes(y=q97.5), size=0.3, linetype="dotted",alpha = 0.8) +
  facet_wrap(~Island, ncol = 2)  +
  scale_x_continuous(name ="Depth (m)", limits=c(3,122), breaks = c(6,20,40,60,90,120)) +
  scale_y_continuous(name ="Generic richness (%)", limits=c(0,30), breaks = c(0,10,20,30)) +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 



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

# Working with PA only 
PA_NMDS <- metaMDS(cast_temp, k=2, trymax = 1000, distance = "jaccard") 

# Prepare to plot
treat=c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15)) # Alphabetic order for colours
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

par (mfrow =c(1,1))

pdf("~/Desktop/NMDS_Diversity_PA_Jaccard.pdf", bg = "white") # starts writing a PDF to file
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
dev.off()


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
beta_1 <- mean (coral.matrices_Depth_1$beta.jac)
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
beta_2 <- mean (coral.matrices_Depth_2$beta.jac)
mean (coral.matrices_Depth_2$beta.jac)

### For 40m  #View(Depth_3) - Necessary rows as sites and columns as species (genera)
Depth_3 <- filter (cast_depth, Depth == 40)
rownames (Depth_3) <- Depth_3$Coral_genus
Depth_3 <- subset (Depth_3, select = - c(Depth, Coral_genus))

Depth_3_Beta <- as.data.frame (t(Depth_3))
#I think, I need to delete columns where all genus are 0 
Depth_3_Beta <- Depth_3_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_3 <- beta.pair(Depth_3_Beta, index.family = "jaccard")
beta_3 <- mean (coral.matrices_Depth_3$beta.jac)
mean (coral.matrices_Depth_3$beta.jac)

### For 60m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_4 <- filter (cast_depth, Depth == 60)
rownames (Depth_4) <- Depth_4$Coral_genus
Depth_4 <- subset (Depth_4, select = - c(Depth, Coral_genus))

Depth_4_Beta <- as.data.frame (t(Depth_4))
#I think, I need to delete columns where all genus are 0 
Depth_4_Beta <- Depth_4_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_4 <- beta.pair(Depth_4_Beta, index.family = "jaccard")
beta_4 <- mean (coral.matrices_Depth_4$beta.jac)
mean (coral.matrices_Depth_4$beta.jac)

### For 90m  #View(Depth_5) - Necessary rows as sites and columns as species (genera)
Depth_5 <- filter (cast_depth, Depth == 90)
rownames (Depth_5) <- Depth_5$Coral_genus
Depth_5 <- subset (Depth_5, select = - c(Depth, Coral_genus))

Depth_5_Beta <- as.data.frame (t(Depth_5))
#I think, I need to delete columns where all genus are 0 
Depth_5_Beta <- Depth_5_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_5 <- beta.pair(Depth_5_Beta, index.family = "jaccard")
beta_5 <- mean (coral.matrices_Depth_5$beta.jac)
mean (coral.matrices_Depth_5$beta.jac)

### For 120m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_6 <- filter (cast_depth, Depth == 120)
rownames (Depth_6) <- Depth_6$Coral_genus
Depth_6 <- subset (Depth_6, select = - c(Depth, Coral_genus))

Depth_6_Beta <- as.data.frame (t(Depth_6))
#I think, I need to delete columns where all genus are 0 
Depth_6_Beta <- Depth_6_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_6 <- beta.pair(Depth_6_Beta, index.family = "jaccard")
beta_6 <- mean (coral.matrices_Depth_6$beta.jac)
mean (coral.matrices_Depth_6$beta.jac)


# Create data.frame
boxplot (coral.matrices_Depth_1$beta.jac,coral.matrices_Depth_2$beta.jac, coral.matrices_Depth_3$beta.jac,coral.matrices_Depth_4$beta.jac,coral.matrices_Depth_5$beta.jac, coral.matrices_Depth_6$beta.jac, 
         xlab = "Depth (m)",
         ylab = "Beta.bray - Jaccard",
         names = c("6","20","40","60","90", "120"), 
         main = "Jaccard distance - PA")

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                        beta_jac=c(beta_1, beta_2, beta_3, beta_4, beta_5, beta_6))

beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_jac~ Depth,beta_div_depth ))



############ I think I can delete all this ############
# 
# beta_div_depth$Depth = factor(beta_div_depth$Depth,levels = c ("6", "20", "40", "60", "90", "120"))
# 
# ggplot (beta_div_depth, aes (x = Depth, y = beta_jac)) + geom_col()
# 
# 
# 
# 
# Beta_Depth_1_matrix <- melt(as.matrix(coral.matrices_Depth_1$beta.jac), varnames = c("row", "col"))
# Beta_Depth_1_matrix$Depth <- "6"
# Beta_Depth_3_matrix <- melt(as.matrix(coral.matrices_Depth_3$beta.jac), varnames = c("row", "col"))
# Beta_Depth_3_matrix$Depth <- "40"
# Beta_Depth_6_matrix <- melt(as.matrix(coral.matrices_Depth_6$beta.jac), varnames = c("row", "col"))
# Beta_Depth_6_matrix$Depth <- "120"
# 
# Beta_Depth_All_matrix <- rbind (Beta_Depth_1_matrix,Beta_Depth_6_matrix)
# Beta_Depth_All_matrix$Depth <- as.numeric (Beta_Depth_All_matrix$Depth)
# 
# summary (lm (value ~ Depth, Beta_Depth_All_matrix))
# 
# MRM(dist(value) ~ dist(sitelocation) + dist(forestpct), data=graze, nperm=10)
# 
# # I can also plot bray-distance according to vertical depth distance measuring bray-distance per site and not per depth
############ I think I can delete all this ############

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
# Necessary for distance afterwards
cast_temp$ID <- with(cast_temp, paste0(Island,  sep = "_", Island_Site))

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

# Beta jaccard
# Beta_Depth <- dist (coral.matrices_Depth$beta.jac)
Beta_Depth <- coral.matrices_Depth$beta.jac

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

mod_Depth <- betadisper (Beta_Depth,groups)
mod_Depth

# Here make the NMDS for Jaccard separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
pdf("~/Desktop/NMDS_Diversity_PA_Jaccard.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  horiz = F, text.width = 0.2, text.font = 4)
title ("Jaccard distance - PA")
dev.off()


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
with(mod_Depth, tapply(distances, groups, "mean")) ## This is the beta-diversity value for the table
beta_jacc_depth <- with(mod_Depth, tapply(distances, groups, "mean"))
with(mod_Depth, tapply(distances, groups, "median"))


distances_centroid$Depth = factor(distances_centroid$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# Display same colours as ordination plots
colours <- "black"
# colours <- c("wheat", "aquamarine2","deepskyblue","blue","navyblue","black")

# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot(colour = colours) + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median") + xlab ("Depth (m)") + ggtitle ("Jaccard PA") + 
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 



### Beta jaccard-turnover
# Beta_Depth <- dist (coral.matrices_Depth$beta.jac)
Beta_jtu_Depth <- coral.matrices_Depth$beta.jtu

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

mod_Depth <- betadisper (Beta_jtu_Depth,groups)
mod_Depth

# Here make the NMDS for Jaccard separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
pdf("~/Desktop/NMDS_Diversity_PA_Jaccard_turn.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  horiz = F, text.width = 0.2, text.font = 4)
title ("Jaccard distance (turnover) - PA")
dev.off()

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
with(mod_Depth, tapply(distances, groups, "mean")) ## This is the beta-div turnover value
with(mod_Depth, tapply(distances, groups, "median"))

distances_centroid$Depth = factor(distances_centroid$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median - jturnover") + xlab ("Depth (m)") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 


### Beta jaccard-nestedness
# Beta_Depth <- dist (coral.matrices_Depth$beta.jac)
Beta_nte_Depth <- coral.matrices_Depth$beta.jne

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

mod_Depth <- betadisper (Beta_nte_Depth,groups)
mod_Depth

# Here make the NMDS for Jaccard separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
pdf("~/Desktop/NMDS_Diversity_PA_Jaccard_nest.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  horiz = F, text.width = 0.2, text.font = 4)
title ("Jaccard distance (nestedness) - PA")
dev.off()

boxplot (mod_Depth)

# Perform test
anova(mod_Depth)
# Permutation test for F
permutest(mod_Depth, pairwise = TRUE, permutations = 99)
# Tukey's Honest Significant Differences
plot (mod.HSD <- TukeyHSD(mod_Depth))

print (mod_Depth)


distances_centroid <- as.data.frame(mod_Depth$distances)
colnames (distances_centroid) <- "Distances"
distances_centroid$Depth <- sapply(strsplit(rownames(distances_centroid), "_"), "[", 1)
mean <- aggregate (Distances ~ Depth, distances_centroid,"mean")
sd <- aggregate (Distances ~ Depth, distances_centroid,"sd")
mean$sd <- sd [2]
mean$Distances<- round(mean$Distances, digits = 3)
mean$sd<- round(mean$sd, digits = 3)

# Or 
with(mod_Depth, tapply(distances, groups, "mean"))  ## This is the beta-div nestedness value
with(mod_Depth, tapply(distances, groups, "median"))

distances_centroid$Depth = factor(distances_centroid$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median - jnestedness") + xlab ("Depth (m)") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 


groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

# Mantel tests just for the Jaccard distance from the beta.pair
#####################################
###### First measure distances ######
library (geodist)

Locations <- read.csv(file = "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/GIS_MAP/Deephope_sampling_locations_RAN*.csv", header = T, dec = ".", sep = ";", row.names = 1)
Locations$Island <- gsub("Mangareva", "Gambier", Locations$Island)
Locations$Island <- gsub("Bora Bora", "Bora", Locations$Island)

Locations$ID<- with(Locations, paste0(Island, sep = "_", Site))

# Keep only same ID as coral data 
unique (cast_temp$ID)
Locations <- Locations[Locations$ID %in% cast_temp$ID, ]
Locations <- Locations[order(Locations[,'ID']), ]

rownames (Locations) <- Locations$ID
Locations <- subset (Locations, select =  c(Latitude, Longitude))

dm <- geodist (Locations, measure = "geodesic", paired = T) /1000

dm <- as.data.frame(dm)
rownames (dm) <-  rownames (Locations)
names (dm) <- rownames (Locations) 
###### First measure distances ######
#####################################

Beta_Depth <- coral.matrices_Depth$beta.jac
# Beta_Depth <- coral.matrices_Depth$beta.jtu

Beta_Depth_Matrix <- as.matrix (Beta_Depth)
# Separate Depth 1 = 6m
Beta_Depth_Matrix_Depth_1 <- Beta_Depth_Matrix[c(1:16),c(1:16)]
Beta_Depth_Matrix_Depth_1 <- as.dist(Beta_Depth_Matrix_Depth_1)
mantel(Beta_Depth_Matrix_Depth_1, dm) # 6m
# Separate Depth 2 = 20m 
Beta_Depth_Matrix_Depth_2 <- Beta_Depth_Matrix[c(17:32),c(17:32)]
Beta_Depth_Matrix_Depth_2 <- as.dist(Beta_Depth_Matrix_Depth_2)
mantel(Beta_Depth_Matrix_Depth_2, dm) # 20m
# Separate Depth 3 = 40m
Beta_Depth_Matrix_Depth_3 <- Beta_Depth_Matrix[c(33:48),c(33:48)]
Beta_Depth_Matrix_Depth_3 <- as.dist(Beta_Depth_Matrix_Depth_3)
mantel(Beta_Depth_Matrix_Depth_3, dm) # 40m
# Separate Depth 4 = 60m 
Beta_Depth_Matrix_Depth_4 <- Beta_Depth_Matrix[c(49:64),c(49:64)]
Beta_Depth_Matrix_Depth_4 <- as.dist(Beta_Depth_Matrix_Depth_4)
mantel(Beta_Depth_Matrix_Depth_4, dm) # 60m
# Separate Depth 5 = 90m
Beta_Depth_Matrix_Depth_5 <- Beta_Depth_Matrix[c(65:80),c(65:80)]
Beta_Depth_Matrix_Depth_5 <- as.dist(Beta_Depth_Matrix_Depth_5)
mantel(Beta_Depth_Matrix_Depth_5, dm) # 90m
# Separate Depth 6 = 120m 
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:95),c(81:95)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)
# For 120m, Necessary to change the dimesnsion as I had to remove one site (no corals at all)
# Necessary to create a new for 120 m because in Tahiti 2 there aren's corals
dm_120 <- dm
# Delete columns
dm_120 <- subset(dm_120, select=-c(Tahiti_2))
# Delete rows
row.names.remove <- c("Tahiti_2")
dm_120 <- dm_120[!(row.names(dm_120) %in% row.names.remove), ]
# Create distance object again
dm_120 <- dist (dm_120)
# Finally the test
mantel(Beta_Depth_Matrix_Depth_6, dm_120) #



# Plot of the lm for measuring increase of beta.jac diversity with depth
boxplot (Beta_Depth_Matrix_Depth_1,Beta_Depth_Matrix_Depth_2, Beta_Depth_Matrix_Depth_3,Beta_Depth_Matrix_Depth_4,Beta_Depth_Matrix_Depth_5, Beta_Depth_Matrix_Depth_6, 
         xlab = "Depth (m)",
         ylab = "Beta.bray - Jaccard",
         names = c("6","20","40","60","90", "120"), 
         main = "Jaccard distance - PA")


# Introduce values manually: 
beta_jacc_depth # Very important to verify they come from the model considering jaccard (not turnover or nestedness)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_jac=c(0.2698211,0.2491936, 0.2538985, 0.2751023, 0.4701546, 0.3677171))


beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_jac~ Depth,beta_div_depth ))





####################################################################################
####################################################################################
####################################################################################
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
beta_1 <- mean (coral.matrices_Depth_1$beta.bray) # PArtially divided between species replacement and nestedness. Accounts for the two of them. 
mean (coral.matrices_Depth_1$beta.bray)
# Turnover is higher than nestedness in the overall dissimilarity 


### For 20m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_2 <- filter (cast_depth, Depth == 20)
rownames (Depth_2) <- Depth_2$Coral_genus
Depth_2 <- subset (Depth_2, select = - c(Depth, Coral_genus))

Depth_2_Beta <- as.data.frame (t(Depth_2))
#I think, I need to delete columns where all genus are 0 
Depth_2_Beta <- Depth_2_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_2 <- beta.pair.abund(Depth_2_Beta, index.family = "bray")
beta_2 <- mean (coral.matrices_Depth_2$beta.bray)
mean (coral.matrices_Depth_2$beta.bray)

### For 40m  #View(Depth_3) - Necessary rows as sites and columns as species (genera)
Depth_3 <- filter (cast_depth, Depth == 40)
rownames (Depth_3) <- Depth_3$Coral_genus
Depth_3 <- subset (Depth_3, select = - c(Depth, Coral_genus))

Depth_3_Beta <- as.data.frame (t(Depth_3))
#I think, I need to delete columns where all genus are 0 
Depth_3_Beta <- Depth_3_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_3 <- beta.pair.abund(Depth_3_Beta, index.family = "bray")
beta_3 <- mean (coral.matrices_Depth_3$beta.bray)
mean (coral.matrices_Depth_3$beta.bray)

### For 60m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_4 <- filter (cast_depth, Depth == 60)
rownames (Depth_4) <- Depth_4$Coral_genus
Depth_4 <- subset (Depth_4, select = - c(Depth, Coral_genus))

Depth_4_Beta <- as.data.frame (t(Depth_4))
#I think, I need to delete columns where all genus are 0 
Depth_4_Beta <- Depth_4_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_4 <- beta.pair.abund(Depth_4_Beta, index.family = "bray")
beta_4 <- mean (coral.matrices_Depth_4$beta.bray)
mean (coral.matrices_Depth_4$beta.bray)

### For 90m  #View(Depth_5) - Necessary rows as sites and columns as species (genera)
Depth_5 <- filter (cast_depth, Depth == 90)
rownames (Depth_5) <- Depth_5$Coral_genus
Depth_5 <- subset (Depth_5, select = - c(Depth, Coral_genus))

Depth_5_Beta <- as.data.frame (t(Depth_5))
#I think, I need to delete columns where all genus are 0 
Depth_5_Beta <- Depth_5_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_5 <- beta.pair.abund(Depth_5_Beta, index.family = "bray")
beta_5 <- mean (coral.matrices_Depth_5$beta.bray)
mean (coral.matrices_Depth_5$beta.bray)

### For 120m  #View(Depth_2) - Necessary rows as sites and columns as species (genera)
Depth_6 <- filter (cast_depth, Depth == 120)
rownames (Depth_6) <- Depth_6$Coral_genus
Depth_6 <- subset (Depth_6, select = - c(Depth, Coral_genus))

Depth_6_Beta <- as.data.frame (t(Depth_6))
#I think, I need to delete columns where all genus are 0 
Depth_6_Beta <- Depth_6_Beta %>% select_if(colSums(.) != 0) 
coral.matrices_Depth_6 <- beta.pair.abund(Depth_6_Beta, index.family = "bray")
beta_6 <- mean (coral.matrices_Depth_6$beta.bray)
mean (coral.matrices_Depth_6$beta.bray)

#### I think I have it: 

# Create data.frame
boxplot (coral.matrices_Depth_1$beta.bray,coral.matrices_Depth_2$beta.bray, coral.matrices_Depth_3$beta.bray,coral.matrices_Depth_4$beta.bray,coral.matrices_Depth_5$beta.bray, coral.matrices_Depth_6$beta.bray,
         xlab = "Depth (m)",
         ylab = "Beta.bray Diversity",
         names = c("6","20","40","60","90", "120"), 
         main = "Bray distance - Occupancy")

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray=c(beta_1, beta_2, beta_3, beta_4, beta_5, beta_6))

beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_bray~ Depth,beta_div_depth ))


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

cast_all_depth <- cast_temp
# Necessary for distance afterwards
cast_temp$ID <- with(cast_temp, paste0(Island,  sep = "_", Island_Site))

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


### Beta bray-
# Beta_Depth <- dist (coral.matrices_Depth$beta.bray)
Beta_Depth <- coral.matrices_Depth$beta.bray

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

mod_Depth <- betadisper (Beta_Depth,groups)
mod_Depth

# Here make the NMDS for Bray separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
pdf("~/Desktop/NMDS_Diversity_Occupancy_Bray.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  horiz = F, text.width = 0.2, text.font = 4)
title ("Bray distance - Occupancy")
dev.off()

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
beta_bray_depth <- with(mod_Depth, tapply(distances, groups, "mean"))
with(mod_Depth, tapply(distances, groups, "median"))

distances_centroid$Depth = factor(distances_centroid$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# Display same colours as ordination plots
colours <- "black"
# colours <- c("wheat", "aquamarine2","deepskyblue","blue","navyblue","black")


# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median - Bray") + xlab ("Depth (m)") + ggtitle ("Occupancy") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"))  


### Beta balance turnover
Beta_Depth <- coral.matrices_Depth$beta.bray.bal

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

mod_Depth <- betadisper (Beta_Depth,groups)
mod_Depth

# Here make the NMDS for Bray separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
pdf("~/Desktop/NMDS_Diversity_Occupancy_Bray_balance.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  horiz = F, text.width = 0.2, text.font = 4)
title ("Bray balance distance - Occupancy")
dev.off()

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

# Display same colours as ordination plots
colours <- "black"
# colours <- c("wheat", "aquamarine2","deepskyblue","blue","navyblue","black")


# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median - Balance") + xlab ("Depth (m)") + ggtitle ("Bray Balance Occupancy") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"))  


### Beta gradient nestedness
Beta_Depth <- coral.matrices_Depth$beta.bray.gra

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

mod_Depth <- betadisper (Beta_Depth,groups)
mod_Depth

# Here make the NMDS for Bray separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
pdf("~/Desktop/NMDS_Diversity_Occupancy_Bray_gradient.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  horiz = F, text.width = 0.2, text.font = 4)
title ("Bray gradient distance - Occupancy")
dev.off()

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

# Display same colours as ordination plots
colours <- "black"
# colours <- c("wheat", "aquamarine2","deepskyblue","blue","navyblue","black")


# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median - Gradient") + xlab ("Depth (m)") + ggtitle ("Bray Gradient Occupancy") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"))


groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

# Mantel tests just for the Bray distance from the beta.pair
#####################################
###### First measure distances ######
library (geodist)

Locations <- read.csv(file = "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/GIS_MAP/Deephope_sampling_locations_RAN*.csv", header = T, dec = ".", sep = ";", row.names = 1)
Locations$Island <- gsub("Mangareva", "Gambier", Locations$Island)
Locations$Island <- gsub("Bora Bora", "Bora", Locations$Island)

Locations$ID<- with(Locations, paste0(Island, sep = "_", Site))

# Keep only same ID as coral data 
unique (cast_temp$ID)
Locations <- Locations[Locations$ID %in% cast_temp$ID, ]
Locations <- Locations[order(Locations[,'ID']), ]

rownames (Locations) <- Locations$ID
Locations <- subset (Locations, select =  c(Latitude, Longitude))

dm <- geodist (Locations, measure = "geodesic", paired = T) /1000

dm <- as.data.frame(dm)
rownames (dm) <-  rownames (Locations)
names (dm) <- rownames (Locations) 
###### First measure distances ######
#####################################


Beta_Depth <- coral.matrices_Depth$beta.bray

Beta_Depth_Matrix <- as.matrix (Beta_Depth)
# Separate Depth 1 = 6m
Beta_Depth_Matrix_Depth_1 <- Beta_Depth_Matrix[c(1:16),c(1:16)]
Beta_Depth_Matrix_Depth_1 <- as.dist(Beta_Depth_Matrix_Depth_1)
mantel(Beta_Depth_Matrix_Depth_1, dm) # 6m
# Separate Depth 2 = 20m 
Beta_Depth_Matrix_Depth_2 <- Beta_Depth_Matrix[c(17:32),c(17:32)]
Beta_Depth_Matrix_Depth_2 <- as.dist(Beta_Depth_Matrix_Depth_2)
mantel(Beta_Depth_Matrix_Depth_2, dm) # 20m
# Separate Depth 3 = 40m
Beta_Depth_Matrix_Depth_3 <- Beta_Depth_Matrix[c(33:48),c(33:48)]
Beta_Depth_Matrix_Depth_3 <- as.dist(Beta_Depth_Matrix_Depth_3)
mantel(Beta_Depth_Matrix_Depth_3, dm) # 40m
# Separate Depth 4 = 60m 
Beta_Depth_Matrix_Depth_4 <- Beta_Depth_Matrix[c(49:64),c(49:64)]
Beta_Depth_Matrix_Depth_4 <- as.dist(Beta_Depth_Matrix_Depth_4)
mantel(Beta_Depth_Matrix_Depth_4, dm) # 60m
# Separate Depth 5 = 90m
Beta_Depth_Matrix_Depth_5 <- Beta_Depth_Matrix[c(65:80),c(65:80)]
Beta_Depth_Matrix_Depth_5 <- as.dist(Beta_Depth_Matrix_Depth_5)
mantel(Beta_Depth_Matrix_Depth_5, dm) # 90m
# Separate Depth 6 = 120m 
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:95),c(81:95)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)
# For 120m, Necessary to change the dimesnsion as I had to remove one site (no corals at all)
# Necessary to create a new for 120 m because in Tahiti 2 there aren's corals
dm_120 <- dm
# Delete columns
dm_120 <- subset(dm_120, select=-c(Tahiti_2))
# Delete rows
row.names.remove <- c("Tahiti_2")
dm_120 <- dm_120[!(row.names(dm_120) %in% row.names.remove), ]
# Create distance object again
dm_120 <- dist (dm_120)
# Finally the test
mantel(Beta_Depth_Matrix_Depth_6, dm_120) #



# Plot of the lm for measuring increase of beta.bray diversity with depth
boxplot (Beta_Depth_Matrix_Depth_1,Beta_Depth_Matrix_Depth_2, Beta_Depth_Matrix_Depth_3,Beta_Depth_Matrix_Depth_4,Beta_Depth_Matrix_Depth_5, Beta_Depth_Matrix_Depth_6, 
         xlab = "Depth (m)",
         ylab = "Beta.bray",
         names = c("6","20","40","60","90", "120"), 
         main = "Bray distance - Occupancy")


# Introduce values manually: 
beta_bray_depth # Very important to verify they come from the model considering jaccard (not turnover or nestedness)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray=c(0.1973027,0.2242932, 0.3139426, 0.2987933, 0.3681515, 0.4483688))


beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_bray~ Depth,beta_div_depth ))



