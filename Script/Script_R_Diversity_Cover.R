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

# Test of different models 
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
plot(mixed.lmer) # Data a bit structured. I am not sure if this is good. 
qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer))
# Check singularity
# diag.vals <- getME(mixed.lmer,"theta")[getME(mixed.lmer,"lower") == 0]
# any(diag.vals < 1e-6) # FALSE

# P-values og glmer
library (nlme)
nlme_1 <- lme(Cover ~ Depth, random = ~1|Island_Site, coral_cover2)
anova(nlme_1)

# Extra mixed models
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

### This is not the good one ####
# ggplot(coral_cover2,aes(Depth, Cover,  col=Island )) + geom_smooth(span = 0.9, method = "loess", size = 0.2, colour = "transparent") +
#   geom_line(aes(y=fit), size=0.8, linetype="dashed", col = "red") + 
#   geom_point(alpha = 0.3) + 
#   theme_bw()
### This is not the good one ####

coral_cover$fit <- predict(mixed.lmer)
# Make the plot again
coral_cover$Depth <- as.numeric (as.character(coral_cover$Depth))

### This is not the good one ####
# ggplot(coral_cover, aes(x=Depth, y=Cover)) + geom_point(aes (colour = Island),size = 0.5, alpha = 0.8)  + geom_smooth(span = 0.8, method = "loess", se = F) + 
#   geom_line(aes(y=fit), size=0.7, linetype="dashed", col = "red") +   stat_summary(fun=mean, geom="point", shape=18, color="red", size=3) + 
#   theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
#   theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
#         axis.text = element_text(size=10, colour="black"),
#         axis.title = element_text(size=11, face="bold", colour="black")) 
### This is not the good one ####

# This is the good one!
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




# NMDS from community matrix of coral cover
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
# pdf("~/Desktop/NMDS_Diversity_Coral_Cover_Bray.pdf", bg = "white") # starts writing a PDF to file
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
# dev.off()


################## Vertical beta diversity ##########################3

#### Now I want to check the vertical dissimilarity within the Island, site (or even Island changing ddply)
#### After meeting 2 keeping all quadrats

##   Betadisper from the beta.pair bray distance matrix, this is the good one!
resume_df <- ddply(rndpts_df, ~ Island + Island_Site +  Depth + Quadrat   +  Coral_genus ,
                   function(x){c(Cover= sum(x$Cover)) })
melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Quadrat","Coral_genus"), measure.vars="Cover", na.rm=FALSE)
cast_temp = dcast(melt_temp,Island + Island_Site +  Depth + Quadrat ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0
cast_temp$Island <- factor(cast_temp$Island, levels = c("Bora","Gambier","Makatea","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))
cast_all_depth <- cast_temp
# Necessary for distance afterwards

# complete to have all quadrats although empty 
cast_all_depth <- cast_all_depth %>% complete( Island,Island_Site,Depth, Quadrat,fill = list(Cover = 0))
cast_all_depth[is.na(cast_all_depth)] <- 0

cast_all_depth$ID<- with(cast_all_depth, paste0(Island, sep = "_",  Island_Site, sep = "_Depth_",Depth, sep = "_", Quadrat ))
cast_all_depth <- as.matrix (cast_all_depth)
row.names(cast_all_depth) <- cast_all_depth[,"ID"]

cast_all_depth <- cast_all_depth[,-c(1,2,3,4,38)]

class(cast_all_depth) <- "numeric"
### Unnecessary now #####
# columns where sum is 0
cast_all_depth <- cast_all_depth[,colSums(cast_all_depth[,])>0]
# rows where sum is 0
cast_all_depth <- cast_all_depth[rowSums(cast_all_depth[,])>0, ]
### Unnecessary now #####
# Measure pair.abund
coral.matrices_Depth <- beta.pair.abund(cast_all_depth, index.family = "bray")
mean (coral.matrices_Depth$beta.bray)

# betadisper per island site, betadisper(island, Prof) and permanova
## Beta-pair in depth per island using 6 m as reference

# For doing it use the function below. It works for all of them, except Moorea_1 and Tahiti_2 because at 120m no corals 

# The function needs that you define: 
# Island_Site: (Moorea_2)
# beta_type : (e.g., beta.bray, beta.bray.bal, or beta.bray.gra)

Div_profile <-function(island_site,beta_type){
  
  island <-  island_site
  beta_div <- beta_type
  
  Beta_Depth <- coral.matrices_Depth [[beta_div]]
  Beta_Depth_Matrix <- as.matrix (Beta_Depth)
  
  Beta_Depth_island <- Beta_Depth_Matrix[grepl(island, rownames(Beta_Depth_Matrix)),grepl(island, colnames(Beta_Depth_Matrix))]
  Beta_Depth_island_dist <- as.dist(Beta_Depth_island)
  
  Depth_6 <- count(rownames(Beta_Depth_island) %like% "Depth_6_")
  Depth_20 <- count(rownames(Beta_Depth_island) %like% "Depth_20_")
  Depth_40 <- count(rownames(Beta_Depth_island) %like% "Depth_40_")
  Depth_60 <- count(rownames(Beta_Depth_island) %like% "Depth_60_")
  Depth_90 <- count(rownames(Beta_Depth_island) %like% "Depth_90_")
  Depth_120 <- count(rownames(Beta_Depth_island) %like% "Depth_120_")
  
  # betadisper for Island at the different depths
  groups <- c(rep("6M",Depth_6),rep("20M",Depth_20), rep("40M",Depth_40), rep("60M",Depth_60), rep("90M",Depth_90), rep("120M",Depth_120)) # In case I keep all quadrats
  
  mod_Depth <- betadisper (Beta_Depth_island_dist,groups)
  
  # plot (mod_Depth, main = island)
  
  # Prepare to plot
  colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
  # colours <- c( "aquamarine2","deepskyblue","blue","wheat","navyblue") # For Moorea_1 and Tahiti_2
  
  # 120, 20, 40, 60, 6, 90)
  
  # par (mfrow =c(1,1))
  # pdf(paste0("~/Desktop/NMDS_",island,sep = "_", beta_div,".pdf"), bg = "white") # starts writing a PDF to file
  ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
  orditorp(mod_Depth,display="sites",labels =T)
  ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
  ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
  # orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
  legend(x="topleft", y="topleft", bg = "transparent",bty="n", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), cex= 0.5, horiz = F)
  title (paste0(island,sep = "_", beta_div))
  #  dev.off()
  
  # boxplot (mod_Depth)
  # Perform test
  anova(mod_Depth)
  # Permutation test for F
  permutest(mod_Depth, pairwise = TRUE, permutations = 99)
  
  distances_centroid <- as.data.frame(mod_Depth$distances)
  colnames (distances_centroid) <- "Distances"
  distances_centroid$Depth <- sapply(strsplit(rownames(distances_centroid), "_"), "[", 4)
  mean <- aggregate (Distances ~ Depth, distances_centroid,"mean")
  sd <- aggregate (Distances ~ Depth, distances_centroid,"sd")
  mean$sd <- sd [2]
  mean$Distances<- round(mean$Distances, digits = 3)
  mean$sd<- round(mean$sd, digits = 3)
  
  distances_centroid$Depth = factor(distances_centroid$Depth,levels = c ("6", "20", "40", "60", "90", "120"))
  
  # Red points are average distance to median of betadisper
  boxplot_betadisper <- ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
    geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
    ylab ("Distance to median") + xlab ("Depth (m)") + ggtitle (paste0(island,sep = "_", beta_div)) + 
    theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                        axis.text = element_text(size=10, colour="black"),
                        axis.title = element_text(size=11, face="bold", colour="black")) 
  
  # Compare beta in reference to 6m (6m will be 0, 20m will be 6 and 20, 40m will be 6 and 40)
  
  # reference 6m 
  Beta_Depth_island_ref <- Beta_Depth_island[c(1:Depth_6),1]
  # Beta_Depth_island_ref <- as.dist(Beta_Depth_island_ref)
  mean (Beta_Depth_island_ref)
  
  # 6m vs 20 m
  Beta_Depth_island_ref_20 <- Beta_Depth_island[c(1:Depth_6),c((Depth_6+1):(Depth_6+Depth_20))]
  # Beta_Depth_island_ref_20 <- as.dist(Beta_Depth_island_ref_20)
  mean (Beta_Depth_island_ref_20)
  
  # 6m vs 40 m
  Beta_Depth_island_ref_40 <- Beta_Depth_island[c(1:Depth_6),c((Depth_6+Depth_20+1):(Depth_6+Depth_20+Depth_40))]
  # Beta_Depth_island_ref_40 <- as.dist(Beta_Depth_island_ref_40)
  mean (Beta_Depth_island_ref_40)
  
  # 6m vs 60 m
  Beta_Depth_island_ref_60 <- Beta_Depth_island[c(1:Depth_6),c((Depth_6+Depth_20+Depth_40+1):(Depth_6+Depth_20+Depth_40+Depth_60))]
  # Beta_Depth_island_ref_60 <- as.dist(Beta_Depth_island_ref_60)
  mean (Beta_Depth_island_ref_60)
  
  # 6m vs 90 m
  Beta_Depth_island_ref_90 <- Beta_Depth_island[c(1:Depth_6),c((Depth_6+Depth_20+Depth_40+Depth_60+1):(Depth_6+Depth_20+Depth_40+Depth_60+Depth_90))]
  # Beta_Depth_island_ref_90 <- as.dist(Beta_Depth_island_ref_90)
  mean (Beta_Depth_island_ref_90)
  
  # 6m vs 120 m
  Beta_Depth_island_ref_120 <- Beta_Depth_island[c(1:Depth_6),c((Depth_6+Depth_20+Depth_40+Depth_60+Depth_90+1):(Depth_6+Depth_20+Depth_40+Depth_60+Depth_90+Depth_120))]
  # Beta_Depth_island_ref_120 <- as.dist(Beta_Depth_island_ref_120)
  mean (Beta_Depth_island_ref_120)
  
  beta_vert_prof <<- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), Island_Site = c(island,island,island,island,island,island),
                                beta=c(0,mean (Beta_Depth_island_ref_20), mean (Beta_Depth_island_ref_40), mean (Beta_Depth_island_ref_60), mean (Beta_Depth_island_ref_90), mean (Beta_Depth_island_ref_120)))
  
  # beta_vert_prof <- data.frame(Depth=c("6", "20", "40", "60", "90"),beta=c(0,mean (Beta_Depth_island_ref_20), mean (Beta_Depth_island_ref_40), mean (Beta_Depth_island_ref_60), mean (Beta_Depth_island_ref_90)))
  # For Moorea_1 and Tahiti_2 
  
  print (island)
  print (beta_div)
  print ("- BETADISPER MOD -")
  print (mod_Depth)
  print ("- BETADISPER ANOVA -")
  print (anova(mod_Depth))
  print ("BETA RELATIVE TO 6m")
  
  print(beta_vert_prof)
  print(boxplot_betadisper)
  
}

# Let's start using the function:

# Bora 1
Bora_1_bray <- Div_profile("Bora_1", "beta.bray")
Bora_1_bray_bal <- Div_profile("Bora_1", "beta.bray.bal")
Bora_1_bray_gra <- Div_profile("Bora_1", "beta.bray.gra")

# Bora 2
Bora_2_bray <- Div_profile("Bora_2", "beta.bray")
Bora_2_bray_bal <- Div_profile("Bora_2", "beta.bray.bal")
Bora_2_bray_gra <- Div_profile("Bora_2", "beta.bray.gra")
# Gambier 1
Gambier_1_bray <- Div_profile("Gambier_1", "beta.bray")
Gambier_1_bray_bal <- Div_profile("Gambier_1", "beta.bray.bal")
Gambier_1_bray_gra <- Div_profile("Gambier_1", "beta.bray.gra")# It does not work
island_site <- "Gambier_1"
beta_type <- "beta.bray.gra"
# Gambier 2
Gambier_2_bray <- Div_profile("Gambier_2", "beta.bray")
Gambier_2_bray_bal <- Div_profile("Gambier_2", "beta.bray.bal")
Gambier_2_bray_gra <- Div_profile("Gambier_2", "beta.bray.gra")
# Makatea 1
Makatea_1_bray <- Div_profile("Makatea_1", "beta.bray")
Makatea_1_bray_bal <- Div_profile("Makatea_1", "beta.bray.bal")
Makatea_1_bray_gra <- Div_profile("Makatea_1", "beta.bray.gra")
# Makatea 2
Makatea_2_bray <- Div_profile("Makatea_2", "beta.bray")
Makatea_2_bray_bal <- Div_profile("Makatea_2", "beta.bray.bal")
Makatea_2_bray_gra <- Div_profile("Makatea_2", "beta.bray.gra")
# Moorea 1 it doesn't work because at 120m no corals at all, need to do it manually
island_site <- "Moorea_1"
beta_type <- "beta.bray"
beta_type <- "beta.bray.bal"
beta_type <- "beta.bray.gra"
# Moorea 2
Moorea_2_bray <- Div_profile("Moorea_2", "beta.bray")
Moorea_2_bray_bal <- Div_profile("Moorea_2", "beta.bray.bal")
Moorea_2_bray_gra <- Div_profile("Moorea_2", "beta.bray.gra")
# Rangiroa 1
Rangiroa_1_bray <- Div_profile("Rangiroa_1", "beta.bray")
Rangiroa_1_bray_bal <- Div_profile("Rangiroa_1", "beta.bray.bal")
Rangiroa_1_bray_gra <- Div_profile("Rangiroa_1", "beta.bray.gra")
# Rangiroa 2
Rangiroa_2_bray <- Div_profile("Rangiroa_2", "beta.bray")
Rangiroa_2_bray_bal <- Div_profile("Rangiroa_2", "beta.bray.bal")
Rangiroa_2_bray_gra <- Div_profile("Rangiroa_2", "beta.bray.gra")
# Raroia 1
Raroia_1_bray <- Div_profile("Raroia_1", "beta.bray")
Raroia_1_bray_bal <- Div_profile("Raroia_1", "beta.bray.bal")
Raroia_1_bray_gra <- Div_profile("Raroia_1", "beta.bray.gra")
# Raroia 2
Raroia_2_bray <- Div_profile("Raroia_2", "beta.bray")
Raroia_2_bray_bal <- Div_profile("Raroia_2", "beta.bray.bal")
Raroia_2_bray_gra <- Div_profile("Raroia_2", "beta.bray.gra")
# Tahiti 1
Tahiti_1_bray <- Div_profile("Tahiti_1", "beta.bray")
Tahiti_1_bray_bal <- Div_profile("Tahiti_1", "beta.bray.bal")
Tahiti_1_bray_gra <- Div_profile("Tahiti_1", "beta.bray.gra")
# Tahiti 2, it doesn't work because at 120m no corals at all, need to do it manually
island_site <- "Tahiti_2"
beta_type <- "beta.bray"
beta_type <- "beta.bray.bal"
beta_type <- "beta.bray.gra"
# Tikehau 1
Tikehau_1_bray <- Div_profile("Tikehau_1", "beta.bray")
Tikehau_1_bray_bal <- Div_profile("Tikehau_1", "beta.bray.bal")
Tikehau_1_bray_gra <- Div_profile("Tikehau_1", "beta.bray.gra")# It does not work
island_site <- "Tikehau_1"
beta_type <- "beta.bray.gra"
# Tikehau 2
Tikehau_2_bray <- Div_profile("Tikehau_2", "beta.bray")
Tikehau_2_bray_bal <- Div_profile("Tikehau_2", "beta.bray.bal")
Tikehau_2_bray_gra <- Div_profile("Tikehau_2", "beta.bray.gra")

################## Vertical beta diversity ##########################3


################## Spatial beta diversity ##########################3
# betadisper
# I do ti for all sites and depths together, considering groups as different depths. I obtain (1) the average distance to median (2) anova and (3)permutest pair-wise differences between (depths)

### Betadisper from the beta.pair bray distance matrix, this is the good one!

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

# plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
# pdf("~/Desktop/NMDS_Diversity_Coralcover_Bray.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  cex= 0.5, horiz = F)
title ("Bray distance - Coral cover")
# dev.off()

# boxplot (mod_Depth)

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


# Measure of mean per depth of beta.bray
Beta_Depth_Matrix <- as.matrix (Beta_Depth)
# Separate Depth 1 = 6m
Beta_Depth_Matrix_Depth_1 <- Beta_Depth_Matrix[c(1:16),c(1:16)]
Beta_Depth_Matrix_Depth_1 <- as.dist(Beta_Depth_Matrix_Depth_1)

# Separate Depth 2 = 20m 
Beta_Depth_Matrix_Depth_2 <- Beta_Depth_Matrix[c(17:32),c(17:32)]
Beta_Depth_Matrix_Depth_2 <- as.dist(Beta_Depth_Matrix_Depth_2)

# Separate Depth 3 = 40m
Beta_Depth_Matrix_Depth_3 <- Beta_Depth_Matrix[c(33:48),c(33:48)]
Beta_Depth_Matrix_Depth_3 <- as.dist(Beta_Depth_Matrix_Depth_3)

# Separate Depth 4 = 60m 
Beta_Depth_Matrix_Depth_4 <- Beta_Depth_Matrix[c(49:64),c(49:64)]
Beta_Depth_Matrix_Depth_4 <- as.dist(Beta_Depth_Matrix_Depth_4)

# Separate Depth 5 = 90m
Beta_Depth_Matrix_Depth_5 <- Beta_Depth_Matrix[c(65:80),c(65:80)]
Beta_Depth_Matrix_Depth_5 <- as.dist(Beta_Depth_Matrix_Depth_5)

# Separate Depth 6 = 120m 
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:94),c(81:94)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth # These are the values for the table


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

# plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
# pdf("~/Desktop/NMDS_Diversity_Coral_cover_Bray.Balance.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  cex= 0.5, horiz = F)
title ("Bray distance balance - Coral cover")
# dev.off()

# boxplot (mod_Depth)

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

# Measure of mean per depth of beta.bray.bal
Beta_Depth_Matrix <- as.matrix (Beta_Depth)
# Separate Depth 1 = 6m
Beta_Depth_Matrix_Depth_1 <- Beta_Depth_Matrix[c(1:16),c(1:16)]
Beta_Depth_Matrix_Depth_1 <- as.dist(Beta_Depth_Matrix_Depth_1)

# Separate Depth 2 = 20m 
Beta_Depth_Matrix_Depth_2 <- Beta_Depth_Matrix[c(17:32),c(17:32)]
Beta_Depth_Matrix_Depth_2 <- as.dist(Beta_Depth_Matrix_Depth_2)

# Separate Depth 3 = 40m
Beta_Depth_Matrix_Depth_3 <- Beta_Depth_Matrix[c(33:48),c(33:48)]
Beta_Depth_Matrix_Depth_3 <- as.dist(Beta_Depth_Matrix_Depth_3)

# Separate Depth 4 = 60m 
Beta_Depth_Matrix_Depth_4 <- Beta_Depth_Matrix[c(49:64),c(49:64)]
Beta_Depth_Matrix_Depth_4 <- as.dist(Beta_Depth_Matrix_Depth_4)

# Separate Depth 5 = 90m
Beta_Depth_Matrix_Depth_5 <- Beta_Depth_Matrix[c(65:80),c(65:80)]
Beta_Depth_Matrix_Depth_5 <- as.dist(Beta_Depth_Matrix_Depth_5)

# Separate Depth 6 = 120m 
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:94),c(81:94)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray_bal=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth

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

# plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
# pdf("~/Desktop/NMDS_Diversity_Coral_cover_Bray.Gradient.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  cex= 0.5, horiz = F)
title ("Bray distance gradient - Coral cover")
# dev.off()

# boxplot (mod_Depth)

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

# Measure of mean per depth of beta.bray.gra
Beta_Depth_Matrix <- as.matrix (Beta_Depth)
# Separate Depth 1 = 6m
Beta_Depth_Matrix_Depth_1 <- Beta_Depth_Matrix[c(1:16),c(1:16)]
Beta_Depth_Matrix_Depth_1 <- as.dist(Beta_Depth_Matrix_Depth_1)

# Separate Depth 2 = 20m 
Beta_Depth_Matrix_Depth_2 <- Beta_Depth_Matrix[c(17:32),c(17:32)]
Beta_Depth_Matrix_Depth_2 <- as.dist(Beta_Depth_Matrix_Depth_2)

# Separate Depth 3 = 40m
Beta_Depth_Matrix_Depth_3 <- Beta_Depth_Matrix[c(33:48),c(33:48)]
Beta_Depth_Matrix_Depth_3 <- as.dist(Beta_Depth_Matrix_Depth_3)

# Separate Depth 4 = 60m 
Beta_Depth_Matrix_Depth_4 <- Beta_Depth_Matrix[c(49:64),c(49:64)]
Beta_Depth_Matrix_Depth_4 <- as.dist(Beta_Depth_Matrix_Depth_4)

# Separate Depth 5 = 90m
Beta_Depth_Matrix_Depth_5 <- Beta_Depth_Matrix[c(65:80),c(65:80)]
Beta_Depth_Matrix_Depth_5 <- as.dist(Beta_Depth_Matrix_Depth_5)

# Separate Depth 6 = 120m 
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:94),c(81:94)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray_gra=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth # These are the values of the table



# Mantel tests just for the Bray distance from the beta.pair
#####################################
###### First measure distances ######
library (geodist)
setwd("~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Data")
Locations <- read.csv(file = "Deephope_sampling_locations_RAN*.csv", header = T, dec = ".", sep = ";", row.names = 1)
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

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",14))

# I re-open the beta.bray
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

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth

beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_bray~ Depth,beta_div_depth ))

################## Spatial beta diversity ##########################3






