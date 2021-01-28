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
  theme_bw()  + ylab ("Cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# no-mid-domain effect. No test necessary. 



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
mean (coral.matrices_Depth_1$beta.bray) # PArtially divided between species replacement and nestedness. Accounts for the two of them. 
# Turnover is slightly smaller than nestedness in the overall dissimilarity 


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

### For 60m  #View(Depth_4) - Necessary rows as sites and columns as species (genera)
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

### For 120m  #View(Depth_6) - Necessary rows as sites and columns as species (genera)
Depth_6 <- filter (cast_depth, Depth == 120)
rownames (Depth_6) <- Depth_6$Coral_genus
Depth_6 <- subset (Depth_6, select = - c(Depth, Coral_genus))

Depth_6_Beta <- as.data.frame (t(Depth_6))
#I think, I need to delete columns where all genus are 0 
Depth_6_Beta <- Depth_6_Beta %>% select_if(colSums(.) != 0) 
Depth_6_Beta <- Depth_6_Beta[rowSums(Depth_6_Beta[,])>0, ]
coral.matrices_Depth_6 <- beta.pair.abund(Depth_6_Beta, index.family = "bray")
mean (coral.matrices_Depth_6$beta.bray)












