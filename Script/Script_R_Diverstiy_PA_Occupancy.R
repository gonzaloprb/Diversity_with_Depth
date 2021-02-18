
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

rm (list = ls()) 
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


## ggplot of generic richness vs predicted richness of the model for all islands.
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


# Separate per islands if you prefer
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



## NMDS from community matrix (slide 33 skipped in the presentation)

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

# pdf("~/Desktop/NMDS_Diversity_PA_Jaccard.pdf", bg = "white") # starts writing a PDF to file
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
# dev.off()



################## Vertical beta diversity ##########################3

#### Now I want to check the vertical dissimilarity within the Island site 
# After meeting 2 trying to keep all quadrats

##   Betadisper from the beta.pair jac distance matrix 
resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth + Quadrat +  Coral_genus ,function(x){
  c(nobserv=nrow(x)) })
melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Quadrat","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)
### PA (0 or 1)
melt_temp$value[melt_temp$value > 1] <- 1
cast_temp = dcast(melt_temp, Island + Island_Site + Depth + Quadrat ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0
cast_temp$Island <- factor(cast_temp$Island, levels = c("Bora","Makatea","Gambier","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))
cast_all_depth <- cast_temp
# Necessary for distance afterwards
cast_temp$ID <- with(cast_temp, paste0(Island,  sep = "_", Island_Site))

# complete to have all quadrats although empty 
cast_all_depth <- cast_all_depth %>% complete( Island,Island_Site,Depth, Quadrat,fill = list(nobserv = 0))
cast_all_depth[is.na(cast_all_depth)] <- 0

cast_all_depth$ID<- with(cast_all_depth,paste0(Island, sep = "_",  Island_Site, sep = "_Depth_",Depth, sep = "_", Quadrat ))
cast_all_depth <- as.matrix (cast_all_depth)
row.names(cast_all_depth) <- cast_all_depth[,"ID"]

cast_all_depth <- cast_all_depth[,-c(1,2,3,4,39)]

class(cast_all_depth) <- "numeric"

# columns where sum is 0
cast_all_depth <- cast_all_depth[,colSums(cast_all_depth[,])>0]
# rows where sum is 0
cast_all_depth <- cast_all_depth[rowSums(cast_all_depth[,])>0, ]

# Measure pair.abund
coral.matrices_Depth <- beta.pair(cast_all_depth, index.family = "jaccard")
mean (coral.matrices_Depth$beta.jac)


# betadisper per island, betadisper(island, Prof) and permanova
# Beta-pair in depth per island using 6 m as reference

# For doing it use the function below. It works for all of them, except Tahiti_2 because at 120m no corals 

# The function needs that you define: 
# Island_Site: (Moorea_2)
# beta_type : (e.g., beta.jac, beta.jtu, or beta.jne)

Div_profile_jac <-function(island_site,beta_type){
  
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
Bora_1_bray <- Div_profile_jac("Bora_1", "beta.jac")
Bora_1_bray_bal <- Div_profile_jac("Bora_1", "beta.jtu")
Bora_1_bray_gra <- Div_profile_jac("Bora_1", "beta.jne")

# Bora 2
Bora_2_bray <- Div_profile_jac("Bora_2", "beta.jac")
Bora_2_bray_bal <- Div_profile_jac("Bora_2", "beta.jtu")
Bora_2_bray_gra <- Div_profile_jac("Bora_2", "beta.jne")
# Gambier 1
Gambier_1_bray <- Div_profile_jac("Gambier_1", "beta.jac")
Gambier_1_bray_bal <- Div_profile_jac("Gambier_1", "beta.jtu")
Gambier_1_bray_gra <- Div_profile_jac("Gambier_1", "beta.jne")
# Gambier 2
Gambier_2_bray <- Div_profile_jac("Gambier_2", "beta.jac")
Gambier_2_bray_bal <- Div_profile_jac("Gambier_2", "beta.jtu")
Gambier_2_bray_gra <- Div_profile_jac("Gambier_2", "beta.jne")
# Makatea 1
Makatea_1_bray <- Div_profile_jac("Makatea_1", "beta.jac")
Makatea_1_bray_bal <- Div_profile_jac("Makatea_1", "beta.jtu")
Makatea_1_bray_gra <- Div_profile_jac("Makatea_1", "beta.jne")
# Makatea 2
Makatea_2_bray <- Div_profile_jac("Makatea_2", "beta.jac")
Makatea_2_bray_bal <- Div_profile_jac("Makatea_2", "beta.jtu")
Makatea_2_bray_gra <- Div_profile_jac("Makatea_2", "beta.jne")
# Moorea 1 
Moorea_1_jac <- Div_profile_jac("Moorea_1", "beta.jac")
Moorea_1_bray_bal <- Div_profile_jac("Moorea_1", "beta.jtu")
Moorea_1_bray_gra <- Div_profile_jac("Moorea_1", "beta.jne")
# Moorea 2
Moorea_2_bray <- Div_profile_jac("Moorea_2", "beta.jac")
Moorea_2_bray_bal <- Div_profile_jac("Moorea_2", "beta.jtu")
Moorea_2_bray_gra <- Div_profile_jac("Moorea_2", "beta.jne")
# Rangiroa 1
Rangiroa_1_bray <- Div_profile_jac("Rangiroa_1", "beta.jac")
Rangiroa_1_bray_bal <- Div_profile_jac("Rangiroa_1", "beta.jtu")
Rangiroa_1_bray_gra <- Div_profile_jac("Rangiroa_1", "beta.jne") # It does not work
island_site <- "Rangiroa_1"
beta_type <- "beta.jne"
# Rangiroa 2
Rangiroa_2_bray <- Div_profile_jac("Rangiroa_2", "beta.jac")
Rangiroa_2_bray_bal <- Div_profile_jac("Rangiroa_2", "beta.jtu")
Rangiroa_2_bray_gra <- Div_profile_jac("Rangiroa_2", "beta.jne")
# Raroia 1
Raroia_1_bray <- Div_profile_jac("Raroia_1", "beta.jac")
Raroia_1_bray_bal <- Div_profile_jac("Raroia_1", "beta.jtu")
Raroia_1_bray_gra <- Div_profile_jac("Raroia_1", "beta.jne")
# Raroia 2
Raroia_2_bray <- Div_profile_jac("Raroia_2", "beta.jac")
Raroia_2_bray_bal <- Div_profile_jac("Raroia_2", "beta.jtu")
Raroia_2_bray_gra <- Div_profile_jac("Raroia_2", "beta.jne")
# Tahiti 1
Tahiti_1_bray <- Div_profile_jac("Tahiti_1", "beta.jac")
Tahiti_1_bray_bal <- Div_profile_jac("Tahiti_1", "beta.jtu")
Tahiti_1_bray_gra <- Div_profile_jac("Tahiti_1", "beta.jne")
# Tahiti 2, it doesn't work because at 120m no corals at all, need to do it manually
island_site <- "Tahiti_2"
beta_type <- "beta.jac"
beta_type <- "beta.jtu"
beta_type <- "beta.jne"
# Tikehau 1
Tikehau_1_bray <- Div_profile_jac("Tikehau_1", "beta.jac")
Tikehau_1_bray_bal <- Div_profile_jac("Tikehau_1", "beta.jtu")
Tikehau_1_bray_gra <- Div_profile_jac("Tikehau_1", "beta.jne")
# Tikehau 2
Tikehau_2_bray <- Div_profile_jac("Tikehau_2", "beta.jac")
Tikehau_2_bray_bal <- Div_profile_jac("Tikehau_2", "beta.jtu")
Tikehau_2_bray_gra <- Div_profile_jac("Tikehau_2", "beta.jne")

################## Vertical beta diversity ##########################3


################## Spatial beta diversity ##########################3
### Betadisper 

# I do it for all sites and depths together, considering groups as different depths. I obtain (1) the average distance to median, (2) anova and (3)permutest pair-wise differences between (depths)
# Each row is a Depth_Island_Site

##   Betadisper from the beta.pair jac distance matrix of PA

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

### Beta jaccard
# Beta_Depth <- dist (coral.matrices_Depth$beta.jac)
Beta_Depth <- coral.matrices_Depth$beta.jac

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

mod_Depth <- betadisper (Beta_Depth,groups)
mod_Depth # Average distance to median is not the beta.jac and turn, nest

# Here make the NMDS for Jaccard separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
# pdf("~/Desktop/NMDS_Diversity_PA_Jaccard.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  cex= 0.5, horiz = F)
title ("Jaccard distance - PA")
# dev.off()

# boxplot (mod_Depth)

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
with(mod_Depth, tapply(distances, groups, "mean")) ## This is NOT the beta-diversity value for the table. We want the mean of the beta pair for the depth
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


# Measure of mean per depth of beta.jac
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
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:95),c(81:95)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_jac=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth # This is the value for the table


### Beta jaccard-turnover

# Beta_Depth <- dist (coral.matrices_Depth$beta.jac)
Beta_Depth <- coral.matrices_Depth$beta.jtu

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
# pdf("~/Desktop/NMDS_Diversity_PA_Jaccard_turn.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), cex= 0.5, horiz = F)
title ("Jaccard distance (turnover) - PA")
# dev.off()

# boxplot (mod_Depth)

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
with(mod_Depth, tapply(distances, groups, "mean")) ## This is NOT the beta-div turnover value
with(mod_Depth, tapply(distances, groups, "median"))

distances_centroid$Depth = factor(distances_centroid$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median - jturnover") + xlab ("Depth (m)") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 

# Measure of mean per depth of beta.jac turn
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
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:95),c(81:95)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_jac_turn=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth # This is the beta.turn value

### Beta jaccard-nestedness
# Beta_Depth <- dist (coral.matrices_Depth$beta.jac)
Beta_Depth <- coral.matrices_Depth$beta.jne

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

mod_Depth <- betadisper (Beta_Depth,groups)
mod_Depth

# Here make the NMDS for Jaccard separating by depths
plot (mod_Depth)
# Prepare to plot
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

# plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
# pdf("~/Desktop/NMDS_Diversity_PA_Jaccard_nest.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  cex= 0.5, horiz = F)
title ("Jaccard distance (nestedness) - PA")
# dev.off()

# boxplot (mod_Depth)

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
with(mod_Depth, tapply(distances, groups, "mean"))  ## This is not the beta-div nestedness value
with(mod_Depth, tapply(distances, groups, "median"))

distances_centroid$Depth = factor(distances_centroid$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median - jnestedness") + xlab ("Depth (m)") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 

# Measure of mean per depth of beta.jac_nes
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
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:95),c(81:95)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_jac_nes=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth   # This is the value for the table



# Mantel tests just for the Jaccard distance from the beta.pair
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

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

### Re-open the beta.jac dis matrix
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
                             beta_jac=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))


beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_jac~ Depth,beta_div_depth ))
################## Spatial beta diversity ##########################3




####################################################################################
####################################################################################
####################################################################################
####################################################################################

# ############### "Frequency = Nb Quadrats with genus / Total Nb Quadrats" ###################
## NMDS 

# rm (list = ls()) 
# Set working directory etc.
setwd("~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Data")
PA_df <- read.csv(file = "photoquad_sppcount_DEEPHOPE_genus_form.csv", header = T, dec = ".", sep = ",", row.names = 1)

PA_df$Island <- gsub("Mangareva", "Gambier", PA_df$Island)
PA_df <- subset (PA_df,Coral_genus!="NA_Coral")


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

# NMDS from the community matrix of occupancy, "Frequency"

PA_NMDS <- metaMDS(cast_temp, k=2, trymax = 1000, distance = "bray") 

# Prepare to plot
treat=c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15)) # Alphabetic order for colours
colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 120, 20, 40, 60, 6, 90)

par (mfrow =c(1,1))
# pdf("~/Desktop/NMDS_Diversity_Frequency/Occupancy_Bray.pdf", bg = "white") # starts writing a PDF to file
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
# dev.off()

################## Vertical beta diversity cannot be done keepign all quadrats because of the nature of the database ##########################3
# Check at the end of the script, but certainly not good #



################## Spatial beta diversity ##########################3
### Betadisper 

# I do it for all sites and depths together, considering groups as different depths. I obtain (1) the average distance to median, (2) anova and (3)permutest pair-wise differences between (depths)
# Each row is a Depth_Island_Site

##   Betadisper from the beta.pair bray distance matrix of Frequency

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

# plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
# pdf("~/Desktop/NMDS_Diversity_Occupancy_Bray.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  cex= 0.5, horiz = F)
title ("Bray distance - Occupancy")
# dev.off()

# boxplot (mod_Depth)

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
  ylab ("Distance to median") + xlab ("Depth (m)") + ggtitle ("Occupancy - Bray") +
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
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:95),c(81:95)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth


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

# plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
# pdf("~/Desktop/NMDS_Diversity_Occupancy_Bray_balance.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), cex= 0.5, horiz = F)
# dev.off()

# boxplot (mod_Depth)

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
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:95),c(81:95)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray_bal=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth



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

# plot(mod_Depth, hull = FALSE, ellipse = F, labels =F)

par (mfrow =c(1,1))
# pdf("~/Desktop/NMDS_Diversity_Occupancy_Bray_gradient.pdf", bg = "white") # starts writing a PDF to file
ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  cex= 0.5, horiz = F)
title ("Bray gradient distance - Occupancy")
# dev.off()

# boxplot (mod_Depth)

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

# Measure of mean per depth of beta.bray_nes
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
Beta_Depth_Matrix_Depth_6 <- Beta_Depth_Matrix[c(81:95),c(81:95)]
Beta_Depth_Matrix_Depth_6 <- as.dist(Beta_Depth_Matrix_Depth_6)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray_nes=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth



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

groups <- c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15))

# Re-open the beta.bray distance matrix
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
                             beta_bray=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))
beta_div_depth

beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_bray~ Depth,beta_div_depth ))


################## Spatial beta diversity ##########################3





######################## From here below, unnecessary ################################
# This is the try of vertical beta-diversity
# No replicates because no quadrats. Only one value from each row which is depth, Island, site
## Beta-pair in depth per island site using 6 m as reference
# betadisper per island, betadisper(island, Prof) and permanova

# Choose the island, beta or component. Need to do it manually!!!
Beta_Depth <- coral.matrices_Depth$beta.bray
beta_div <- "beta.bray"
# Beta_Depth <- coral.matrices_Depth$beta.bray.bal
# beta_div <- "beta.bray.bal"
# Beta_Depth <- coral.matrices_Depth$beta.bray.gra
# beta_div <- "beta.bray.gra"

 

Beta_Depth_Matrix <- as.matrix (Beta_Depth)


# Choose one of the islands ("Bora","Makatea","Gambier","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))
island <-  "Moorea"
# island <-  "Tahiti"

Beta_Depth_island <- Beta_Depth_Matrix[grepl(island, rownames(Beta_Depth_Matrix)),grepl(island, colnames(Beta_Depth_Matrix))]
Beta_Depth_island_dist <- as.dist(Beta_Depth_island)

# betadisper for Island at the different depths - Moorea
groups <- c(rep("6M",2),rep("20M",2), rep("40M",2), rep("60M",2), rep("90M",2), rep("120M",2))
# groups <- c(rep("6M",2),rep("20M",2), rep("40M",2), rep("60M",2), rep("90M",2), rep("120M",1)) ### For Tahiti missing depths

mod_Depth <- betadisper (Beta_Depth_island_dist,groups)
mod_Depth
plot (mod_Depth, main = island)
# boxplot (mod_Depth)
# Perform test
anova(mod_Depth)
# Permutation test for F
permutest(mod_Depth, pairwise = TRUE, permutations = 99)
print (mod_Depth)

distances_centroid <- as.data.frame(mod_Depth$distances)
colnames (distances_centroid) <- "Distances"
distances_centroid$Depth <- sapply(strsplit(rownames(distances_centroid), "_"), "[", 1)
mean <- aggregate (Distances ~ Depth, distances_centroid,"mean")
sd <- aggregate (Distances ~ Depth, distances_centroid,"sd")
mean$sd <- sd [2]
mean$Distances<- round(mean$Distances, digits = 3)
mean$sd<- round(mean$sd, digits = 3)

distances_centroid$Depth = factor(distances_centroid$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# Red points are average distance to median of betadisper
ggplot(data = distances_centroid, aes(y = Distances, x = Depth)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Depth), colour = "red", size = 3) + 
  ylab ("Distance to median") + xlab ("Depth (m)") + ggtitle (paste0(island,sep = "_", beta_div)) + 
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 

# Compare beta in reference to 6m 

# reference 6m 
Beta_Depth_island_ref <- Beta_Depth_island[c(1:2),c(1:2)]
Beta_Depth_island_ref <- as.dist(Beta_Depth_island_ref)
mean (Beta_Depth_island_ref)

# 6m vs 20 m
Beta_Depth_island_ref_20 <- Beta_Depth_island[c(1:4),c(1:4)]
Beta_Depth_island_ref_20 <- as.dist(Beta_Depth_island_ref_20)
mean (Beta_Depth_island_ref_20)

# 6m vs 40 m
Beta_Depth_island_ref_40 <- Beta_Depth_island[c(1,2,5,6),c(1,2,5,6)]
Beta_Depth_island_ref_40 <- as.dist(Beta_Depth_island_ref_40)
mean (Beta_Depth_island_ref_40)

# 6m vs 60 m
Beta_Depth_island_ref_60 <- Beta_Depth_island[c(1,2,7,8),c(1,2,7,8)]
Beta_Depth_island_ref_60 <- as.dist(Beta_Depth_island_ref_60)
mean (Beta_Depth_island_ref_60)

# 6m vs 90 m
Beta_Depth_island_ref_90 <- Beta_Depth_island[c(1,2,9,10),c(1,2,9,10)]
Beta_Depth_island_ref_90 <- as.dist(Beta_Depth_island_ref_90)
mean (Beta_Depth_island_ref_90)

# 6m vs 120 m
Beta_Depth_island_ref_120 <- Beta_Depth_island[c(1,2,11,12),c(1,2,11,12)]
# Beta_Depth_island_ref_120 <- Beta_Depth_island[c(1,2,11),c(1,2,11)] ### For Tahiti, missing corals in one depth
Beta_Depth_island_ref_120 <- as.dist(Beta_Depth_island_ref_120)
mean (Beta_Depth_island_ref_120)


