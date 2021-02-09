# rm (list = ls())

#   title: "Script_Diversity_Stats"


# Script also pushed on Github: https://github.com/gonzaloprb/Diverstiy-with-Depth

# Install.packages
library(scatterplot3d); library(reshape); library (reshape2); library (data.table); library(RColorBrewer); library(dplyr); library (tidyr); library (plyr); library (ggplot2); require (vegan); require (goeveg);   require  (stats); require (dynRB); require (matrixStats); require(plyr); require (stringr);require (reshape2); require(gridExtra)
require (betapart); library (car); library (MASS); library (glmm); library (glmnet); library (glmmTMB); library (lme4)
require (ggpubr); library (cowplot); library (patchwork)


# Mid-domain effect (there's no mid-domain effect)
require(devtools)
library(reshape2) 
install_github("cran/rangemodelR")
library(rangemodelR)

# Geospatial data
require (geodist)

# rm (list = ls()) 
# Set working directory etc.
setwd("~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Data")
coral_data <- read.csv(file = "coral_data.csv", header = T, dec = ".", sep = ",", row.names = 1)
rndpts_df <- read.csv(file = "photoquad_rndpts_DEEPHOPE.csv", header = T, dec = ".", sep = ",", row.names = 1)

# Correct the island name, Michel's suggestion 
rndpts_df$Island <- gsub("Mangareva", "Gambier", rndpts_df$Island)
coral_data$Island <- gsub("Mangareva", "Gambier", coral_data$Island)

# DElete the unnecessary columns for coral_data
coral_data <- subset (coral_data, select = - c(Archipelago, Cattegory))

# Keep only corals in rndpts_df / fix the genus/form issue, aggregate
rndpts_df <- filter (rndpts_df, Cattegory == "Scleractinian")
rndpts_df <- rndpts_df %>% separate(Coral_Genus_Form,  c("Coral_genus", "Coral_form"), " ", extra = "merge")
rndpts_df <- subset (rndpts_df, select = - c(Cattegory, Coral_form))
rndpts_df <- aggregate (Cover ~ Unique_Image_ID + Date + Site + Archipelago + Island + Island_Site + Depth + Quadrat + Coral_genus, rndpts_df , sum)



# All stats below come from rndpts_df or pooled together quadrats coral_data
# I re-start every time for simplicity 


## Coral cover profile 

# Plot considering the effect of Island and Island_Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
coral_cover <- aggregate (Cover ~ Island + Island_Site + Depth, coral_data , sum)


# Transform depth as a Qualitative variable  
coral_cover$Depth <- as.factor(as.character(coral_cover$Depth))
coral_cover$Depth = factor(coral_cover$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# ggplot with island sites
ggplot(coral_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (colour = Island),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 




# no-mid-domain effect. 

#### Mixed-Model for showing cover decreases with depth in all sites
# Just very quickly linear model. I don't think I need to complicate myself with glm or lmer

coral_cover2 <- coral_cover
coral_cover2$Depth <- as.numeric (as.character(coral_cover2$Depth))
coral_cover2$Island_Site <- as.character(coral_cover2$Island_Site)
coral_cover2$Island_Site_2 <- paste(coral_cover2$Island, "_", coral_cover2$Island_Site)

# lm_1 <- lm(Cover ~  Depth, data = coral_cover2) 
# summary (lm_1)

lm_2 <- lm(Cover ~  Depth:Island_Site_2, data = coral_cover2) 
summary (lm_2)

# lm_3 <- lm(Cover ~ Depth + Depth*Island_Site_2, data = coral_cover2) 
# summary (lm_3)
### ? 

# Mixed model with depth as fixed factor and Island_Site as Random
# Convergence code: 0 means no error

mixed.lmer <- lmer(Cover ~ Depth + (1|Island_Site_2), data = coral_cover2)
summary(mixed.lmer)
# No variance as Site is only considered as random intercept
plot(mixed.lmer) 
qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer))
# Check singularity
# diag.vals <- getME(mixed.lmer,"theta")[getME(mixed.lmer,"lower") == 0]
# any(diag.vals < 1e-6) # FALSE

# P-values og glmer
library (nlme)
nlme_1 <- lme(Cover ~ Depth, random = ~1|Island_Site, coral_cover2)
anova(nlme_1)

mixed.lmer_2 <- lmer(Cover ~ Depth + (Depth | Island_Site_2), coral_cover2)
summary(mixed.lmer_2)
# 1.766e+01/(1.766e+01 +1.804e+02) --> Island_Site_2 explains 8% of variance
plot(mixed.lmer_2) 
qqnorm(resid(mixed.lmer_2))
qqline(resid(mixed.lmer_2))

mixed.lmer_3 <- lmer(Cover ~ Depth + (Depth || Island_Site_2), coral_cover2)
summary(mixed.lmer_3)
plot(mixed.lmer_3) 
qqnorm(resid(mixed.lmer_3))
qqline(resid(mixed.lmer_3))

mixed.lmer_4 <- lmer(Cover ~ 1 + Depth + (1 + Depth | Island_Site_2), coral_cover2)
summary(mixed.lmer_4)
# 1.766e+01/(1.766e+01 +1.804e+02) --> Intercept of Island_Site_2 explains 8% of variance
# 1.726e-03 /(1.726e-03  +1.804e+02) --> Intercept of Island_Site_2 explains 8% of variance
plot(mixed.lmer_4) 
qqnorm(resid(mixed.lmer_4))
qqline(resid(mixed.lmer_4))

mixed.lmer_5 <- lmer(Cover ~ 1 + Depth + (1 + Depth | Island_Site_2) +  (1|Island_Site_2) , coral_cover2)
summary(mixed.lmer_5)
plot(mixed.lmer_5) 
qqnorm(resid(mixed.lmer_5))
qqline(resid(mixed.lmer_5))

anova (mixed.lmer,mixed.lmer_2,mixed.lmer_3,mixed.lmer_4,mixed.lmer_5)

# Smaller AIC and BIC is mixed.lmer     Cover = Depth * -0.31038 + 38.82288 

# Fitted plot from model 

coral_cover2$fit <- predict(mixed.lmer) # Add model to the dataframe

ggplot(coral_cover2,aes(Depth, Cover,  col=Island )) + geom_smooth(span = 0.9, method = "loess", size = 0.2, colour = "transparent") +
  geom_line(aes(y=fit), size=0.8, linetype="dashed", col = "red") + 
  geom_point(alpha = 0.3) + 
  theme_bw()

coral_cover$fit <- predict(mixed.lmer)
# Make the plot again
coral_cover$Depth <- as.numeric (as.character(coral_cover$Depth))

ggplot(coral_cover, aes(x=Depth, y=Cover)) + geom_point(aes (colour = Island),size = 0.5, alpha = 0.8)  + geom_smooth(span = 0.8, method = "loess", se = F) + 
  geom_line(aes(y=fit), size=0.7, linetype="dashed", col = "red") +   stat_summary(fun=mean, geom="point", shape=18, color="red", size=3) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 




ggplot(coral_cover, aes(x=Depth, y=Cover)) + geom_point(aes(colour = Island),size = 0.5, alpha = 0.8)  + 
  geom_ribbon(stat = "smooth",method = "loess", se = T, alpha = 0, colour = "black", linetype = "dotted")+
  stat_summary(fun=mean, geom="point", shape=18, color="black", size=4) + 
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", color="black", size=0.5, width=2) +
  geom_line(aes(y=fit), size=1, linetype="dashed", col = "blue") +
  scale_x_continuous(name ="Depth (m)", limits=c(-2,130), breaks = c(6,20,40,60,90,120)) +
  scale_y_continuous(name ="Coral cover (%)", limits=c(-5,80), breaks = c(0,20,40,60,80)) +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                  axis.text = element_text(size=10, colour="black"),
                  axis.title = element_text(size=11, face="bold", colour="black")) 




# NMDS for cover,
# 1- working with all quadrats pooled together first (coral_data)

# Complete rows for all genera with Coral cover = for all equal to 0
# coral_data_all <- coral_data %>% complete(Depth, Island,Island_Site, Coral_genus,fill = list(Cover = 0))

melt_island = melt(coral_data, id=c("Depth","Island", "Island_Site","Coral_genus"), measure.vars="Cover", na.rm=FALSE)
cast_island = dcast(melt_island, Depth + Island + Island_Site ~ Coral_genus, mean, add.missing = T)

# Missing values to 0
cast_island[is.na(cast_island)] <- 0

# Islands as factors to choose the order
cast_island$Island <- factor(cast_island$Island, levels = c("Bora","Gambier","Makatea","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))

# Generation of new database for adonis stats
cast_database <- cast_island
# Preparing unique rownames for format matrix
cast_island$ID<- with(cast_island, paste0(Depth, sep = "_",Island, sep = "_", Island_Site))
rownames (cast_island) <- cast_island$ID

# Delete some genera
cast_island <- subset (cast_island, select = - c(Non_Id_Coral))

# Delete unnecessary columns
cast_island <- cast_island[,-c(1,2,3,36)]

# As I am working with the Abundance, "Coral cover"
# Compute distances - Bray-curtis
dis <- vegdist (cast_island)
total_NMDS <- metaMDS(cast_island, k=2, trymax = 1000, distance = "bray")     # With two dimensions we are already below 0.2  

# Prepare for plot 
treat=c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",14)) # Alphabetic order for colours
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")


par (mfrow =c(1,1))
# 120, 20, 40, 60, 6, 90)
pdf("~/Desktop/NMDS_Diversity_Coral_Cover_Bray.pdf", bg = "white") # starts writing a PDF to file
ordiplot(total_NMDS,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(total_NMDS,display="sites",labels =T, cex = 0.1)
orditorp(total_NMDS,display="sites",cex=0.3,air=0.01)
ordispider(total_NMDS, groups=treat, col = colours,cex=0.6, lwd = 1.5)
ordihull(total_NMDS,groups=treat,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
orditorp(total_NMDS,display="species",col="red",air=0.01, cex =0.5)
legend(x="topright", y="top", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  pt.cex=0.8, cex = 0.70, horiz = F, text.width = 0.3, text.font = 20)
title ("Bray distance - Coral cover")
mysubtitle <- paste0("Stress = ", format(round(total_NMDS$stress, 2)))
mtext(mysubtitle, side=1, line=-2, at=0, adj=0, cex=0.7) 
dev.off()


########## Maybe here trying to make the NMDS with all quadrats,(rndpts_df) and not quadrats pooled together #########

## Beta.pair per depth
# Prepare database for all quadrats together (coral_data)

# Complete rows for all genera with Coral cover = 0 for all genera
coral_data_all <- coral_data %>% complete(Depth, Island,Island_Site, Coral_genus,fill = list(Cover = 0))

melt_depth = melt(coral_data_all, id=c("Depth","Island", "Island_Site","Coral_genus"), measure.vars="Cover", na.rm=FALSE)
melt_depth$ID<- with(melt_depth, paste0(Island, sep = "_", Island_Site))
cast_depth = dcast(melt_depth, Depth + Coral_genus ~ ID, mean, add.missing = T)

cast_depth <- cast_depth[!grepl("Non_Id_Coral", cast_depth$Coral_genus),]



# Doing same as above in a different way and from the entire quadrats database, gives the same 
# Resume to keep only coral genus
resume_df <- ddply(rndpts_df, ~ Island + Island_Site + Depth + Coral_genus ,
                   function(x){c(Cover= sum(x$Cover)/30) })

# Complete rows for all genera with Coral cover = 0 for all genera
resume_df <- resume_df %>% complete( Island,Island_Site,Depth, Coral_genus,fill = list(Cover = 0))
# Prepare columns and rows (by now keeping depth)
melt_depth = melt(resume_df, id=c("Depth","Island", "Island_Site","Coral_genus"), measure.vars="Cover", na.rm=FALSE)
melt_depth$ID<- with(melt_depth, paste0(Island, sep = "_", Island_Site))
cast_depth = dcast(melt_depth, Depth + Coral_genus ~ ID, mean, add.missing = T)


# Same as in Diversity_PA 

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
# Turnover is slightly smaller than nestedness in the overall dissimilarity 


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

### For 60m  #View(Depth_4) - Necessary rows as sites and columns as species (genera)
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

### For 120m  #View(Depth_6) - Necessary rows as sites and columns as species (genera)
Depth_6 <- filter (cast_depth, Depth == 120)
rownames (Depth_6) <- Depth_6$Coral_genus
Depth_6 <- subset (Depth_6, select = - c(Depth, Coral_genus))

Depth_6_Beta <- as.data.frame (t(Depth_6))
#I think, I need to delete columns where all genus are 0 
Depth_6_Beta <- Depth_6_Beta %>% select_if(colSums(.) != 0) 
Depth_6_Beta <- Depth_6_Beta[rowSums(Depth_6_Beta[,])>0, ]
coral.matrices_Depth_6 <- beta.pair.abund(Depth_6_Beta, index.family = "bray")
beta_6 <- mean (coral.matrices_Depth_6$beta.bray) 
mean (coral.matrices_Depth_6$beta.bray)


#### I think I have it: 

# Create data.frame
boxplot (coral.matrices_Depth_1$beta.bray,coral.matrices_Depth_2$beta.bray, coral.matrices_Depth_3$beta.bray,coral.matrices_Depth_4$beta.bray,coral.matrices_Depth_5$beta.bray, coral.matrices_Depth_6$beta.bray,
         xlab = "Depth (m)",
         ylab = "Beta.bray Diversity",
         names = c("6","20","40","60","90", "120"), 
         main = "Bray distance - Coral cover")

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray=c(beta_1, beta_2, beta_3, beta_4, beta_5, beta_6))

beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_bray~ Depth,beta_div_depth ))


######## I think I can delete ##############

# Beta_Depth_1_matrix <- melt(as.matrix(coral.matrices_Depth_1$beta.bray), varnames = c("row", "col"))
# Beta_Depth_1_matrix$Depth <- "6"
# Beta_Depth_2_matrix <- melt(as.matrix(coral.matrices_Depth_2$beta.bray), varnames = c("row", "col"))
# Beta_Depth_2_matrix$Depth <- "20"
# Beta_Depth_3_matrix <- melt(as.matrix(coral.matrices_Depth_3$beta.bray), varnames = c("row", "col"))
# Beta_Depth_3_matrix$Depth <- "40"
# Beta_Depth_5_matrix <- melt(as.matrix(coral.matrices_Depth_5$beta.bray), varnames = c("row", "col"))
# Beta_Depth_5_matrix$Depth <- "90"
# Beta_Depth_6_matrix <- melt(as.matrix(coral.matrices_Depth_6$beta.bray), varnames = c("row", "col"))
# Beta_Depth_6_matrix$Depth <- "120"
# 
# Beta_Depth_All_matrix <- rbind (Beta_Depth_1_matrix,Beta_Depth_2_matrix,Beta_Depth_3_matrix,Beta_Depth_5_matrix,Beta_Depth_6_matrix)
# Beta_Depth_All_matrix$Depth <- as.numeric (Beta_Depth_All_matrix$Depth)
# 
# summary (lm (value ~ Depth, Beta_Depth_All_matrix))
# 
# library (usedist)
# 
# beta_Depth_1 <- dist_setNames(beta_Depth_1, paste0 ("6",sep = "_",rownames (Depth_1_Beta)))
# 
# beta_Depth_2 <- dist_setNames(beta_Depth_2, paste0 ("20",sep = "_",rownames (Depth_2_Beta)))

######## I think I can delete ##############


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
# Necessary to create a new for 120 m where in Moorea 1 and Tahiti 2 there aren's corals
dm_120 <- dm

# Transform to object of class dist
dm <- dist (dm)

##### Make the mantel test ##### (the more permutations, more significant)
# For Depth 1 - 6m 
mantel(coral.matrices_Depth_1$beta.bray, dm, method = "pearson", permutations = 10000) # 6m , strong + correlation (significant)
# Bray distance 0.6224 correlated and significant
mantel(coral.matrices_Depth_2$beta.bray, dm) # 20m # Bray distance not correlated and not significan with distance
mantel(coral.matrices_Depth_3$beta.bray, dm) # 40m
mantel(coral.matrices_Depth_4$beta.bray, dm) # 60m
mantel(coral.matrices_Depth_5$beta.bray, dm,method = "pearson",permutations = 10000 ) #90m

# For 120m, Necessary to change the dimesnsion as I had to remove two sites (no corals at all)
# Delete columns
dm_120 <- subset(dm_120, select=-c(Moorea_1,Tahiti_2))
# Delete rows
row.names.remove <- c("Moorea_1", "Tahiti_2")
dm_120 <- dm_120[!(row.names(dm_120) %in% row.names.remove), ]
# Create distance object
dm_120 <- dist (dm_120)
mantel(coral.matrices_Depth_6$beta.bray, dm_120) #120m weak + correlation (significant), close to 0 indicates no correlation






### Betadisper 
# I can make it for all depths together, considering groups as different depths. I obtain (1) the average distance to median ("b-dissimilarity" per depth ?), (2) anova and (3)permutest pair-wise differences between (depths)


# 1st betadisper using all depths

##   Betadisper from the beta.pair bray distance matrix, this is the good one!

resume_df <- ddply(rndpts_df, ~ Island + Island_Site + Depth + Coral_genus ,
                   function(x){c(Cover= sum(x$Cover)/30) })

melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Coral_genus"), measure.vars="Cover", na.rm=FALSE)
cast_temp = dcast(melt_temp, Depth +  Island + Island_Site ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0

cast_temp$Island <- factor(cast_temp$Island, levels = c("Bora","Gambier","Makatea","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))

cast_all_depth <- cast_temp
# Necessary for distance afterwards
cast_temp$ID <- with(cast_temp, paste0(Island,  sep = "_", Island_Site))

cast_all_depth$ID<- with(cast_all_depth, paste0(Depth, sep = "_",  Island,  sep = "_", Island_Site))
rownames (cast_all_depth) <- cast_all_depth$ID

cast_all_depth <- cast_all_depth[,-c(1,2,3,37)]
# columns where sum is 0
cast_all_depth <- cast_all_depth %>% select_if(colSums(.) != 0)
# rows where sum is 0
cast_all_depth <- cast_all_depth[rowSums(cast_all_depth[,])>0, ]

# Measure pair.abund
coral.matrices_Depth <- beta.pair.abund(cast_all_depth, index.family = "bray")
mean (coral.matrices_Depth$beta.bray)



### Beta.bray 
# Beta_Depth <- dist (coral.matrices_Depth$beta.bray)
Beta_Depth <- coral.matrices_Depth$beta.bray

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",14))

mod_Depth <- betadisper (Beta_Depth,groups, type = "median", bias.adjust=TRUE)
mod_Depth

# Here make the NMDS for Bray separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
pdf("~/Desktop/NMDS_Diversity_Coralcover_Bray.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  horiz = F, text.width = 0.2, text.font = 4)
title ("Bray distance - Coral cover")
dev.off()


boxplot (mod_Depth)

# Perform test
anova(mod_Depth)
# Permutation test for F
permutest(mod_Depth, pairwise = TRUE, permutations = 999)
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
  ylab ("Distance to median - Bray") + xlab ("Depth (m)") + ggtitle ("Bray Coral cover") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


### Beta. balance turnover
# Beta_Depth <- dist (coral.matrices_Depth$beta.bray)
Beta_Depth <- coral.matrices_Depth$beta.bray.bal

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",14))

mod_Depth <- betadisper (Beta_Depth,groups, type = "median", bias.adjust=TRUE)
mod_Depth

# Here make the NMDS for Bray separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
pdf("~/Desktop/NMDS_Diversity_Coral_cover_Bray.Balance.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  horiz = F, text.width = 0.2, text.font = 4)
title ("Bray distance balance - Coral cover")
dev.off()


boxplot (mod_Depth)

# Perform test
anova(mod_Depth)
# Permutation test for F
permutest(mod_Depth, pairwise = TRUE, permutations = 999)
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
  ylab ("Distance to median - Balance") + xlab ("Depth (m)") + ggtitle ("Bray Balance Coral cover") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 



### Beta. gradient nestedness
# Beta_Depth <- dist (coral.matrices_Depth$beta.bray)
Beta_Depth <- coral.matrices_Depth$beta.bray.gra

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",14))

mod_Depth <- betadisper (Beta_Depth,groups, type = "median", bias.adjust=TRUE)
mod_Depth

# Here make the NMDS for Bray separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
pdf("~/Desktop/NMDS_Diversity_Coral_cover_Bray.Gradient.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  horiz = F, text.width = 0.2, text.font = 4)
title ("Bray distance gradient - Coral cover")
dev.off()


boxplot (mod_Depth)

# Perform test
anova(mod_Depth)
# Permutation test for F
permutest(mod_Depth, pairwise = TRUE, permutations = 999)
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
  ylab ("Distance to median - Gradient") + xlab ("Depth (m)") + ggtitle ("Bray Gradient - Coral cover") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 



groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",14))

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
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:94),c(81:94)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)
# For 120m, Necessary to change the dimesnsion as I had to remove one site (no corals at all)
# Necessary to create a new for 120 m because in Moorea 1 and Tahiti 2 there aren't corals
dm_120 <- dm
# Delete columns
dm_120 <- subset(dm_120, select=-c(Moorea_1,Tahiti_2))
# Delete rows
row.names.remove <- c("Moorea_1", "Tahiti_2")
dm_120 <- dm_120[!(row.names(dm_120) %in% row.names.remove), ]
# Create distance object
dm_120 <- dist (dm_120)
# Finally the test
mantel(Beta_Depth_Matrix_Depth_6, dm_120) #



#### Plot of the lm for measuring increase of beta.bray diversity with depth
boxplot (Beta_Depth_Matrix_Depth_1,Beta_Depth_Matrix_Depth_2, Beta_Depth_Matrix_Depth_3,Beta_Depth_Matrix_Depth_4,Beta_Depth_Matrix_Depth_5, Beta_Depth_Matrix_Depth_6, 
         xlab = "Depth (m)",
         ylab = "Beta.bray",
         names = c("6","20","40","60","90", "120"), 
         main = "Bray distance - Coral cover")


# Introduce values manually: 
beta_bray_depth # Very important to verify they come from the model considering jaccard (not turnover or nestedness)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray=c(0.3148875,0.4232184, 0.5143043, 0.5256691, 0.5515951, 0.5754556))


beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_bray~ Depth,beta_div_depth ))





