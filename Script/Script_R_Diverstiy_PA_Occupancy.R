
# Script also pushed on Github: https://github.com/gonzaloprb/PhD_Diversity_Depth

# Install.packages
require(scatterplot3d); require(reshape); require (reshape2); require (data.table); require(RColorBrewer); require(dplyr); require (tidyr); require (plyr); require (ggplot2); require (vegan); require (goeveg);   require  (stats); require (dynRB); require (matrixStats); require(plyr); require (stringr);require (reshape2); require(gridExtra)
require (betapart); require (car); require (MASS); 
require (ggpubr); require (cowplot); require (patchwork)  ;library(viridisLite); library(fishualize);library(tidyverse) 



# Mid-domain effect
require(devtools)    
library(reshape2)
install_github("cran/rangemodelR")
library(rangemodelR) 

# Geospatial data
require (geodist)

rm (list = ls()) 
# Set working directory etc.
PA_df <- read.csv(file = "Data/photoquad_sppcount_DEEPHOPE_genus_form.csv", header = T, dec = ".", sep = ",", row.names = 1)

PA_df$Island <- gsub("Mangareva", "Gambier", PA_df$Island)
PA_df <- subset (PA_df,Coral_genus!="NA_Coral")

# Necessary to correct several genus typing mistakes
PA_df$Coral_Genus_Form <- gsub("Dipsastrea massive", "Dipsastraea massive", PA_df$Coral_Genus_Form)
PA_df$Coral_genus <- gsub("Dipsastrea", "Dipsastraea", PA_df$Coral_genus)
PA_df$Coral_Genus_Form <- gsub("Echynophillia encrusting", "Echinophyllia encrusting", PA_df$Coral_Genus_Form)
PA_df$Coral_Genus_Form <- gsub("Echynophillia laminar", "Echinophyllia laminar", PA_df$Coral_Genus_Form)
PA_df$Coral_genus <- gsub("Echynophillia", "Echinophyllia", PA_df$Coral_genus)
PA_df$Coral_Genus_Form <- gsub("Coscinarea encrusting", "Coscinaraea encrusting", PA_df$Coral_Genus_Form)
PA_df$Coral_Genus_Form <- gsub("Coscinarea laminar", "Coscinaraea laminar", PA_df$Coral_Genus_Form)
PA_df$Coral_genus <- gsub("Coscinarea", "Coscinaraea", PA_df$Coral_genus)

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


# Measure mean and standard error
summary <- ddply(MDE_all, .(Depth), summarize, Mean=mean(Richness), Richness_se=sd(Richness) / sqrt(length(Richness)),
                 Mean_model=mean(mod.rich), Q2.5_model=mean(q2.5), Q97.5_model=mean(q97.5))
MDE_all <- merge (MDE_all,summary)

cols <- brewer.pal(8, "Dark2")
## ggplot of generic richness vs predicted richness of the model for all islands.
Fig_1B <- ggplot(MDE_all, aes(x=Depth, y=Richness)) + 
 # geom_point(aes(y=Richness,colour = Island),size = 0.8, shape = 19, alpha = 0.8) +
 # geom_point(aes(fill = Island, colour = Island), shape = 21, size = 0.8) +
 # geom_line(aes(y=Richness,colour = Island),size = 0.5,linetype="solid", alpha = 0.8)  + 
  geom_line(aes(y=Mean_model),size = 1,linetype="dashed", alpha = 1)  +  
  geom_line(aes(y=Q2.5_model), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_line(aes(y=Q97.5_model), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean - Richness_se, ymax = Mean + Richness_se), color="black", size=1, width=2) +
  scale_x_continuous(name ="Depth (m)", limits=c(5,122), breaks = c(6,20,40,60,90,120)) +
  scale_y_continuous(name ="Generic richness (Nb. of genera)", limits=c(0,30), breaks = c(0,10,20,30)) +
  scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
  geom_point(aes(y=Mean), shape=21, size=4, fill = "white")+
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "bottom")

Fig_1B # Mean values of the model and IQ 2.5 and 97.5 among islands - to see differences between islands check Supp. 
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_1B.pdf", Fig_1B, width = 4, height = 3.5)

# Plot both databases together. Coral Cover and PA. It is necessary to add Fig. 1C still
(Figure_1 = Fig_1A + Fig_1B + plot_layout(guides = 'collect', ncol = 1)  & theme(legend.position='bottom'))

ggsave( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Figure_1.pdf", Figure_1,  width = 4, height = 8)

# MEasure generic richness per site and not by island
agg2 <- aggregate(Coral_genus ~ Island + Island_Site +  Depth, data=subset(PA_df), FUN=function(x) length(unique(x)))
colnames (agg2) [4]<- "Richness_Site"

MDE_all <- merge (MDE_all, agg2)

MDE_all$Island_Site <- as.factor (MDE_all$Island_Site)
# Separate per islands if you prefer supplementary figure
Fig_S2 <-ggplot(MDE_all, aes(x=Depth, y=Richness)) + 
  facet_wrap(~Island, ncol = 4, scales = "free")  +
  #geom_line(aes(y=Richness),size = 1,linetype="solid", alpha = 0.8, colour = "blue")  + 
  geom_line(aes(y=Mean_model),size = 0.5,linetype="dashed", alpha = 1)  +  
  geom_line(aes(y=Q2.5_model), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_line(aes(y=Q97.5_model), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_point(aes(y=Richness), shape=21, size=1.5, fill = "grey")+
  geom_point(aes(y=Richness_Site, shape = Island_Site), size=0.5, fill = "grey")+ # Variability of sites
  scale_x_continuous(name ="Depth (m)", limits=c(3,122), breaks = c(6,20,40,60,90,120)) +
  scale_y_continuous(name ="Generic richness (%)", limits=c(0,30), breaks = c(0,10,20,30)) +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"), 
                      strip.text = element_text(size = 11,face="bold", colour="black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(), legend.position = "none")
Fig_S2  

ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_S2.pdf", Fig_S2, width = 6, height = 5)




# Depth distribution of all genera!

nb_genera <- ddply(PA_df, ~ Coral_genus + Depth + Site ,function(x){
  c(nobserv=nrow(x)) }) 

nb_genera <- ddply(PA_df, ~ Coral_genus + Depth  ,function(x){
  c(nobserv=nrow(x)) }) 

ggplot(nb_genera, aes(Coral_genus, Depth, fill = factor(Coral_genus))) +  geom_violin(scale = "width",position = "dodge",draw_quantiles = c(0.5)) + 
  scale_y_reverse(name ="Depth (m)",lim=c(120,0), breaks = c(6,20,40,60,90,120)) +
  scale_x_discrete(position = "top", name ="Coral genus") +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text.x = element_text(angle = 90, size = 8),
                          axis.text.y = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none") 





nb_genera$Coral_genus = factor(nb_genera$Coral_genus,levels = c ("Plesiastrea","Hydnophora", "Platygyra",
                                                                 "Pocillopora","Acanthastrea","Lobophyllia","Gardinoseris","Fungia",  "Lobactis", "Napopora",
                                                                 "Herpolitha","Cyphastrea" ,"Paragoniastrea","Cantharellus",
                                                                 "Astrea","Leptoria", "Acropora","Dipsastraea","Leptastrea","Psammocora", "Pavona","Stylocoeniella",
                                                                 "Porites", "Montipora",  "Astreopora","Sandalolitha","Pleuractis","Pachyseris","Leptoseris", 
                                                                  "Echinophyllia","Cycloseris","Coscinaraea","Alveopora","Turbinaria"))
                                                                 

# nb_genera <- subset(nb_genera, nb_genera$nobserv > 50)   # Apply subset function

 
Fig_Depth_Distribution <- ggplot(nb_genera, aes(Coral_genus, fill = "black", Depth, width = nobserv)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 30, ymax = 120),fill = "grey", alpha = 0.03) +
  geom_path(aes (colour = "black"),size = 1, lineend = "round") +
  geom_violin(aes (fill = "black"),scale = "width" ,position = "dodge",draw_quantiles = c(0.5)) +
  scale_y_reverse(name ="Depth (m)",lim=c(120,5), breaks = c(6,20,40,60,90,120)) +
  scale_x_discrete(position = "top", name ="Coral genus") +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text.x = element_text(angle = 90, size = 10),
                          axis.text.y = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none") 
Fig_Depth_Distribution 
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_Depth_Distribution_PA.pdf", Fig_Depth_Distribution, width = 8, height = 6)



# nb_genera <- subset(nb_genera, nb_genera$nobserv > 50)   # Apply subset function



Fig_Depth_Distribution_Main <- ggplot(nb_genera, aes(Coral_genus, Depth, width = nobserv)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 30, ymax = 120),fill = "grey", alpha = 0.03) +
  geom_path( aes (colour = "black"),size = 1, lineend = "round") +
  geom_violin(aes (fill = "black"),scale = "width",position = "dodge",draw_quantiles = c(0.5), trim=T) +
  scale_y_reverse(name ="Depth (m)",lim=c(120,5), breaks = c(6,20,40,60,90,120)) +
  scale_x_discrete(position = "top", name ="Coral genus") +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text.x = element_text(angle = 90, size = 10),
                          axis.text.y = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none") 
Fig_Depth_Distribution_Main 
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_Depth_Distribution_50.pdf", Fig_Depth_Distribution_Main, width = 8, height = 6)



nb_genera <- ddply(PA_df, ~ Coral_genus + Depth + Site ,function(x){
  c(nobserv=nrow(x)) }) 

# Put names of sites more short and beautiful and in order
nb_genera$Site <- gsub('SOCMOO1', 'MOO1', nb_genera$Site)
nb_genera$Site <- gsub('SOCMOO2', 'MOO2', nb_genera$Site)
nb_genera$Site <- gsub('SOCTAH1', 'TAH1', nb_genera$Site)
nb_genera$Site <- gsub('SOCTAH2', 'TAH2', nb_genera$Site)
nb_genera$Site <- gsub('SOCBOR1', 'BOR1', nb_genera$Site)
nb_genera$Site <- gsub('SOCBOR2', 'BOR2', nb_genera$Site)
nb_genera$Site <- gsub('TUATIK1', 'TIK1', nb_genera$Site)
nb_genera$Site <- gsub('TUATIK2', 'TIK2', nb_genera$Site)
nb_genera$Site <- gsub('TUARAN1', 'RAN1', nb_genera$Site)
nb_genera$Site <- gsub('TUARAN2', 'RAN2', nb_genera$Site)
nb_genera$Site <- gsub('TUARAR1', 'RAR1', nb_genera$Site)
nb_genera$Site <- gsub('TUARAR2', 'RAR2', nb_genera$Site)
nb_genera$Site <- gsub('TUAMAK1', 'MAK1', nb_genera$Site)
nb_genera$Site <- gsub('TUAMAK2', 'MAK2', nb_genera$Site)
nb_genera$Site <- gsub('GAMMAN1', 'GAM1', nb_genera$Site)
nb_genera$Site <- gsub('GAMMAN2', 'GAM2', nb_genera$Site)

nb_genera$Site <- factor(nb_genera$Site, levels = c("MOO1","MOO2","TAH1","TAH2","BOR1","BOR2","TIK1","TIK2","RAN1","RAN2","RAR1","RAR2","MAK1","MAK2","GAM1","GAM2"))


nb_genera$Coral_genus = factor(nb_genera$Coral_genus,levels = c ("Plesiastrea","Hydnophora", "Platygyra",
                                                                 "Pocillopora","Acanthastrea","Lobophyllia","Gardinoseris","Fungia",  "Lobactis", "Napopora",
                                                                 "Herpolitha","Cyphastrea" ,"Paragoniastrea","Cantharellus",
                                                                 "Astrea","Leptoria", "Acropora","Dipsastraea","Leptastrea","Psammocora", "Pavona","Stylocoeniella",
                                                                 "Porites", "Montipora",  "Astreopora","Sandalolitha","Pleuractis","Pachyseris","Leptoseris", 
                                                                 "Echinophyllia","Cycloseris","Coscinaraea","Alveopora","Turbinaria"))




Fig_Site_Distribution <- ggplot(nb_genera, aes(Coral_genus, Depth, width = nobserv)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 30, ymax = 120),fill = "grey", alpha = 0.03) +
  geom_path( aes (colour = "black"),size = 1, lineend = "round") +
  #geom_violin(scale = "width",position = "dodge",draw_quantiles = c(0.5), trim=T) +
  facet_wrap(~Site, ncol = 2)+
  scale_y_reverse(name ="Depth (m)",lim=c(120,5), breaks = c(6,20,40,60,90,120)) +
  scale_x_discrete(position = "top", name ="Coral genus") +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text.x = element_text(angle = 90, size = 12),
                          axis.text.y = element_text(size=14, colour="black"),
                          axis.title = element_text(size=14, face="bold", colour="black"), legend.position = "none") 
Fig_Site_Distribution
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_Site_Distribution.pdf", Fig_Site_Distribution, width = 16, height = 18)




## NMDS from community matrix (slide 33 skipped in the presentation)

resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth + Coral_genus ,function(x){
  c(nobserv=nrow(x)) })

melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)

### PA        

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

# Prepare for ggplot 
PA_NMDS_Species <- PA_NMDS[['species']]
PA_NMDS_Species <- as.data.frame(PA_NMDS_Species)
PA_NMDS_Species$Species <- rownames (PA_NMDS_Species)

PA_NMDS_coordinates <- PA_NMDS[['points']]
PA_NMDS_coordinates <- as.data.frame(PA_NMDS_coordinates)
PA_NMDS_coordinates$ID <- rownames (PA_NMDS_coordinates)
# Extract depths
PA_NMDS_coordinates$Depth <- sub("\\_.*", "", PA_NMDS_coordinates$ID)
# Extract Island Sites
PA_NMDS_coordinates$Island <- sub("^[^_]*_", "", PA_NMDS_coordinates$ID)
# Remove ID column 
PA_NMDS_coordinates <- subset(PA_NMDS_coordinates, select=-c(ID))


#### Necessary theme for NMDS ####
theme_blank = function(base_size = 12, base_family = "") { 
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(), 
      axis.text.x = element_text(size = base_size*0.7, color = "black", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.7, color = "black", lineheight = 0.9),  
      axis.ticks = NULL,  
      axis.title.x = element_text(size = base_size*0.9, color = "black", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size*0.9, color = "black", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "black"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "black"),  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "white", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "black"),  
      panel.grid.major = element_line(color = "white"),  
      panel.grid.minor = element_line(color = "white"),  
      # Specify facetting options
      strip.background = NULL,  
      strip.text.x = element_text(size = base_size*0.7, color = "black"),  
      strip.text.y = element_text(size = base_size*0.7, color = "black",angle = -90),  
      # Specify plot options
      plot.title = element_text(size = base_size*1.2, color = "black"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
}
#### Necessary theme for NMDS ####
 
# Plot NMDS
PA_NMDS_coordinates = PA_NMDS_coordinates %>% mutate(Depth = factor(Depth, levels = c(6, 20, 40, 60, 90, 120)))
hull_PA_NMDS <- PA_NMDS_coordinates %>% group_by(Depth) %>% slice(chull(MDS1,MDS2))

keep_species <-  c ("Acropora","Astreopora","Montipora","Leptastrea","Leptoseris","Pachyseris","Pavona","Echynophillia","Fungia","Pocillopora","Porites")
PA_NMDS_Species <- PA_NMDS_Species[PA_NMDS_Species$Species %in% keep_species, ]

# colours of Fishualize 
# scales::show_col(fish(n = 6, option = "Ostracion_whitleyi", direction = -1))
colours <- c("#CBEAF1FF", "#73C2E1FF", "#3F73B1FF", "#28407AFF", "#222E47FF", "#191819FF")


Fig_3A = PA_NMDS_coordinates %>% ggplot() + theme_blank() + 
  geom_point(aes(x = MDS1, y = MDS2, fill = Depth), shape = 21, show.legend = F) +
  geom_polygon(data = hull_PA_NMDS, aes(x = MDS1, y = MDS2, fill = Depth), alpha = 0.75, color = "black", show.legend = T) +
 # geom_text(data = PA_NMDS_Species, aes(x = MDS1, y = MDS2, label = Species), size = 2) +
  geom_label(data = PA_NMDS_Species, aes(x = MDS1, y = MDS2, label = Species), size = 2, alpha = .75) +
  scale_x_continuous(name ="NMDS 1") + scale_y_continuous(name ="NMDS 2") + ggtitle("Jaccard - Presence/Absence") +
  fishualize::scale_fill_fish_d(option = "Ostracion_whitleyi", direction = -1) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=11, face="bold")) +
  annotate(geom = 'text', label = 'Stress = 0.11', 
           x = range(PA_NMDS_coordinates$MDS1)[1], y = range(PA_NMDS_Species$MDS2)[1], hjust = 0, vjust = 7)



# Old plot without ggplot
# treat=c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",15)) # Alphabetic order for colours
# colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# # 120, 20, 40, 60, 6, 90)
# 
# par (mfrow =c(1,1))
# 
# pdf("~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/NMDS_Diversity_PA_Jaccard.pdf", bg = "white") # starts writing a PDF to file
# ordiplot(PA_NMDS,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
# orditorp(PA_NMDS,display="sites",labels =T)
# orditorp(PA_NMDS,display="sites",cex=0.4,air=0.01)
# ordispider(PA_NMDS, groups=treat, col = colours,cex=0.6, lwd = 1.5)
# ordihull(PA_NMDS,groups=treat,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(PA_NMDS,display="species",col="red",air=0.01, cex =0.5)
# legend(x="topright", y="top", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  pt.cex=0.8, cex = 0.70, horiz = F, text.width = 0.3, text.font = 20)
# title ("Jaccard distance - PA")
# mysubtitle <- paste0("Stress = ", format(round(PA_NMDS$stress, 2)))
# mtext(mysubtitle, side=1, line=-2, at=0, adj=0, cex=0.7)
# dev.off()


############# Heatmap of beta-diversity at the site level - PA Jaccard - Figure 3 A and B #############


resume_df <- ddply(PA_df, ~ Island + Island_Site + Depth +  Coral_genus ,function(x){
  c(nobserv=nrow(x)) })
melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Coral_genus"), measure.vars="nobserv", na.rm=FALSE)
### PA (0 or 1)
melt_temp$value[melt_temp$value > 1] <- 1
cast_temp = dcast(melt_temp, Island + Island_Site + Depth  ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0
cast_temp$Island <- factor(cast_temp$Island, levels = c("Bora","Makatea","Gambier","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))
cast_all_depth <- cast_temp
# Necessary for distance afterwards
cast_temp$ID <- with(cast_temp, paste0(Island,  sep = "_", Island_Site))

# complete to have all quadrats although empty 
cast_all_depth <- cast_all_depth %>% complete( Island,Island_Site,Depth,fill = list(nobserv = 0))
cast_all_depth[is.na(cast_all_depth)] <- 0

cast_all_depth$ID<- with(cast_all_depth,paste0(Island, sep = "_",  Island_Site, sep = "_Depth_",Depth ))
cast_all_depth <- as.matrix (cast_all_depth)
row.names(cast_all_depth) <- cast_all_depth[,"ID"]

cast_all_depth <- cast_all_depth[,-c(1,2,3,4,38)] 

class(cast_all_depth) <- "numeric"

# columns where sum is 0
# cast_all_depth <- cast_all_depth[,colSums(cast_all_depth[,])>0]
# rows where sum is 0
# cast_all_depth <- cast_all_depth[rowSums(cast_all_depth[,])>0, ]

# Measure pair.abund
coral.matrices_Depth <- beta.pair(cast_all_depth, index.family = "jaccard")


island_site <-  c("Bora_1", "Bora_2", "Gambier_1", "Gambier_2", "Makatea_1", "Makatea_2", "Moorea_1", "Moorea_2", "Rangiroa_1",
                  "Rangiroa_2", "Raroia_1", "Raroia_2", "Tahiti_1", "Tahiti_2", "Tikehau_1", "Tikehau_2")

depth <- c(6,20,40,60,90,120)

summary_final <- data.frame()
Jaccard_final  <- data.frame()
for (i in unique (island_site)){
  Beta_Depth <-  coral.matrices_Depth$beta.jac
  mat_jac <- as.matrix (Beta_Depth)
  mat_jac_Island <- mat_jac[grepl(i, rownames(mat_jac)),grepl(i, colnames(mat_jac))]
  
  mat_jac = melt(mat_jac_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames (mat_jac) <- c("Var1", "Var2", "nobserv")
  
  mat_jac$Depth_Var1 <- as.numeric (sub('.*\\_', '', mat_jac$Var1))
  mat_jac$Depth_Var2 <- as.numeric (sub('.*\\_', '', mat_jac$Var2))
  
  # Comparing the overlap of 40 or the other depths, with the rest of depths
  # AS 40m is the reference
  mat_jac<-mat_jac[(mat_jac$Depth_Var1==20 ),]
  
  mat_jac$SDist <- (mat_jac$Depth_Var1 - mat_jac$Depth_Var2)
  
  # # All pairs between themselves out
  mat_jac<-mat_jac[!(mat_jac$nobserv==0),]
  # # All pairs between the same depth factors
  mat_jac<-mat_jac[!(mat_jac$SDist==0),]
  
  colnames(mat_jac)[3] <- "Jac_dissimilarity"
  
  
  # Add the other two components
  # Turnover
  Turn_Depth <- coral.matrices_Depth$beta.jtu
  mat_jac_turn <- as.matrix (Turn_Depth)
  mat_jac_turn_Island <- mat_jac_turn[grepl(i, rownames(mat_jac_turn)),grepl(i, colnames(mat_jac_turn))]
  mat_jac_turn = melt(mat_jac_turn_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames (mat_jac_turn) <- c("Var1", "Var2", "Jac_turnover")
  # The merge keeping only all mat_bray already cleans and uses 6m as reference 
  mat_jac <- merge(mat_jac, mat_jac_turn, by=c("Var1","Var2"), all.x = T, all.y = F)
  
  # Nestedness
  Nest_Depth <- coral.matrices_Depth$beta.jne
  mat_jac_nest <- as.matrix (Nest_Depth)
  mat_jac_nest_Island <- mat_jac_nest[grepl(i, rownames(mat_jac_nest)),grepl(i, colnames(mat_jac_nest))]
  mat_jac_nest = melt(mat_jac_nest_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames(mat_jac_nest) <- c("Var1", "Var2", "Jac_nestedness")
  # The merge keeping only all mat_jac already cleans and uses 6m as reference 
  mat_jac <- merge(mat_jac, mat_jac_nest, by=c("Var1","Var2"), all.x = T, all.y = F)
  
  
  # Delete all quadrats without corals 
  mat_jac <- na.omit (mat_jac)
  
  # Compute the index suggested by Seb Villeger - indicator of importance of turnover "Toussaint et al 2014"
  #  mat_jac$pturn = mat_jac$Turn_Similarity / mat_jac$Jac_Similarity
  mat_jac$pturn_Diss = mat_jac$Jac_turnover / mat_jac$Jac_dissimilarity
  
  # Means and Sd in funciton of SDist
  summary_mat_jac <- ddply(mat_jac, .(Depth_Var1,Depth_Var2,SDist), summarize, Jaccard=mean(Jac_dissimilarity), Jac_se=sd(Jac_dissimilarity) / sqrt(length(Jac_dissimilarity)),
                           Turnover=mean(Jac_turnover), Turn_se=sd(Jac_turnover) / sqrt(length(Jac_turnover)),
                           Nestedness=mean(Jac_nestedness), Nest_se=sd(Jac_nestedness) / sqrt(length(Jac_nestedness)), 
                           # P_Turn=mean(pturn, na.omit=T), P_Turn_se=sd(pturn, na.omit=T) / sqrt(length(pturn, na.omit=T)),
                           P_Turn_Diss=mean(pturn_Diss), P_Turn_Diss_se=sd(pturn_Diss) / sqrt(length(pturn_Diss)))
  
  
  summary_2 <- melt(summary_mat_jac, id.vars=c("Depth_Var1","Depth_Var2","SDist"), measure.vars = c("Jaccard", "Turnover", "Nestedness",  "P_Turn_Diss"))
  colnames (summary_2) <- c("Depth_Var1","Depth_Var2","SDist","Component", "Dissimilarity")
  
  summary_se <- melt(summary_mat_jac, id.vars=c("Depth_Var1","Depth_Var2","SDist"), measure.vars = c("Jac_se", "Turn_se", "Nest_se",  "P_Turn_Diss_se"))
  colnames (summary_se) <- c("Depth_Var1","Depth_Var2","SDist","Component", "Se")
  
  summary_se$Component <- gsub('Jac_se', 'Jaccard',  summary_se$Component)
  summary_se$Component <- gsub('Turn_se', 'Turnover',  summary_se$Component)
  summary_se$Component<- gsub('Nest_se', 'Nestedness',  summary_se$Component)
  # summary_se$Component<- gsub('P_Turn_se', 'P_Turn',  summary_se$Component)
  summary_se$Component<- gsub('P_Turn_Diss_se', 'P_Turn_Diss',  summary_se$Component)
  
  summary <- merge (summary_2,  summary_se)
  
  keep <- c("Jaccard", "Turnover","P_Turn_Diss")
  summary <- filter(summary, Component %in% keep)
  
  summary$Island_Site <- i
  
  summary_final = rbind (summary_final, summary)
  
}


# Measure mean and standard error
summary <- ddply(summary_final, ~ SDist + Component, summarize, Sim_mean=mean(Dissimilarity), Sim_se=sd(Dissimilarity) / sqrt(length(Dissimilarity)))
summary_final <- merge (summary_final,summary)


# Change the names on axis
summary_final$SDist<- paste(summary_final$Depth_Var1,"vs",summary_final$Depth_Var2)

# Run the respective row for each depth only
# summary_final$SDist <- factor(summary_final$SDist, levels = c("6 vs 20","6 vs 40","6 vs 60","6 vs 90","6 vs 120"))
# summary_final$SDist <- factor(summary_final$SDist, levels = c("20 vs 6","20 vs 40","20 vs 60","20 vs 90","20 vs 120"))
summary_final$SDist <- factor(summary_final$SDist, levels = c("40 vs 6","40 vs 20","40 vs 60","40 vs 90","40 vs 120"))
# summary_final$SDist <- factor(summary_final$SDist, levels = c("60 vs 6","60 vs 20","60 vs 40","60 vs 90","60 vs 120"))
# summary_final$SDist <- factor(summary_final$SDist, levels = c("90 vs 6","90 vs 20","90 vs 40","90 vs 60","90 vs 120"))
# summary_final$SDist <- factor(summary_final$SDist, levels = c("120 vs 6","120 vs 20","120 vs 40","120 vs 60","120 vs 90"))


# Measure the mean fo all locations
pturn <- aggregate(Dissimilarity ~ SDist + Component, summary_final, mean )

pturn$Dissimilarity[pturn$Component=="P_Turn_Diss"] 

pturn$Dissimilarity[pturn$Component=="Jaccard"] 

pturn$Dissimilarity[pturn$Component=="Turnover"] 


# Necessary to fill a table manually!!!! 
# Do the same with Script_R_Diversity for abundance and Bray Curtis
# And make the calculation for each depth separately! later open the New_database

pturn <- read.csv(file = "Data/Pturn.csv", header = T, dec = ".", sep = ";", row.names = 1)

pturn$Depth1 <- as.factor (pturn$Depth1)
pturn$Depth2 <- as.factor (pturn$Depth2)

pturn$Depth2 <- factor(pturn$Depth2, levels = c("6","20","40","60","90", "120"))


pturn$Quart_Pturn_Jac <- ntile(pturn$Pturn_Jac_Site, 3)

pturn$Quart_Pturn_Jac <- factor(pturn$Pturn_Jac_Site, levels=c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1") )

# Jaccard
jac <- ggplot(pturn, aes(Depth2, Depth1, fill=Jaccard_Site))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.5, limit=c(0,1), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), space="Lab", name="Dissimilarity")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="Depth (m)",y="Depth (m)", title = "Jaccard") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))


# Pturnover Jaccard
pturn_jac <- ggplot(pturn, aes(Depth2, Depth1, fill=Pturn_Jac_Site))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.5, limit=c(0,1), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), space="Lab", name="P Turnover")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="",y="Depth (m)", title = "Replacement") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))


# Bray
bray <- ggplot(pturn, aes(Depth2, Depth1, fill=Bray_Site))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.5, limit=c(0,1), breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), space="Lab", name="Dissimilarity")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="Depth (m)",y="Depth (m)", title = "Bray-Curtis") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))

# Pturnover Bray
pturn_bray <- ggplot(pturn, aes(Depth2, Depth1, fill=Pturn_Bray_Site))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.5, limit=c(0,1),breaks = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), space="Lab", name="P Turnover")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="",y="Depth (m)", title = "Substitution") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))


# Figure heat map
# Plot both graphs together. Jaccard and Bray
(Figure_heatmap_Site = jac +  bray  + pturn_jac + pturn_bray +  plot_layout(guides = 'collect' )) 


ggsave( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Figure_heatmap_Site.pdf", Figure_heatmap_Site,  width = 8, height = 8)

############# Heatmap of beta-diversity at the site level - PA Jaccard - Figure 3 A and B #############

########################################## End of heat map ############################################

# !!!!!! FROM HERE REDUCE EVERYTHING (with the three points) UNTIL Spatial beta diversity !!!!!! #

################## Vertical beta diversity ##########################

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
  print ("- PERMUTEST -")
  print (permutest(mod_Depth, pairwise = TRUE, permutations = 99))
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

################## Distance decay dissimilarity ###################### 
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
# cast_all_depth <- cast_all_depth[,colSums(cast_all_depth[,])>0]
# rows where sum is 0
# cast_all_depth <- cast_all_depth[rowSums(cast_all_depth[,])>0, ]

# Measure pair.abund
coral.matrices_Depth <- beta.pair(cast_all_depth, index.family = "jaccard")



### Trying to automatize the code below ######

island_site <-  c("Bora_1", "Bora_2", "Gambier_1", "Gambier_2", "Makatea_1", "Makatea_2", "Moorea_1", "Moorea_2", "Rangiroa_1",
                  "Rangiroa_2", "Raroia_1", "Raroia_2", "Tahiti_1", "Tahiti_2", "Tikehau_1", "Tikehau_2")

summary_final <- data.frame()
Jaccard_final  <- data.frame()
for (i in unique (island_site)){
  Beta_Depth <-  coral.matrices_Depth$beta.jac
  mat_jac <- as.matrix (Beta_Depth)
  mat_jac_Island <- mat_jac[grepl(i, rownames(mat_jac)),grepl(i, colnames(mat_jac))]
  
  mat_jac = melt(mat_jac_Island, measure.vars="nobserv", na.rm=FALSE)
  
  colnames (mat_jac) <- c("Var1", "Var2", "nobserv")
  
  mat_jac$Depth_Var1 <- as.numeric (gsub("(.*_){2}(\\d+)_.+", "\\2", mat_jac$Var1))
  mat_jac$Depth_Var2 <- as.numeric (gsub("(.*_){2}(\\d+)_.+", "\\2", mat_jac$Var2))

  
  ####### adding this new# Make the comparison against shallow waters! 6 and 20m together grouped as 15 m. Not 6 m anymore
  # mat_jac$Depth_Var1[mat_jac$Depth_Var1 == 6] <- 15
  # mat_jac$Depth_Var1[mat_jac$Depth_Var1 == 20] <- 15

  # mat_jac$Depth_Var2[mat_jac$Depth_Var2 == 6] <- 15
  # mat_jac$Depth_Var2[mat_jac$Depth_Var2 == 20] <- 15
  #
  #
  mat_jac$SDist <- abs(mat_jac$Depth_Var1 - mat_jac$Depth_Var2)
  #
  ####### adding this new#  Against 15 m because 6 and 20m (as surface reference) instead of 6 m
  # mat_jac<-mat_jac[(mat_jac$Depth_Var1==15),]

  # All pairs between themselves out
  # mat_jac<-mat_jac[!(mat_jac$value==0),]
  # All pairs between the same depth factors
  mat_jac<-mat_jac[!(mat_jac$SDist==0),]

 
  
  # as 6m is reference keep only the comparisons against 6m (we do not want comparisons between 90 and 120m)
  # mat_jac<-mat_jac[(mat_jac$Depth_Var1==6 ),]
  
  # as 20m reference - necessary to delete the comparisons against 6 m also 
  # mat_jac<-mat_jac[(mat_jac$Depth_Var1==20 ),]
  # mat_jac<-mat_jac[!(mat_jac$Depth_Var2==6),]
  
  ###### LAetitia request Only compare dissimilariry between 40 and 60, and 90 and 120m #######
  # mat_jac$Depth_Var1[mat_jac$Depth_Var1 == 40] <- 50
  # mat_jac$Depth_Var1[mat_jac$Depth_Var1 == 60] <- 50
  # mat_jac$Depth_Var2[mat_jac$Depth_Var2 == 40] <- 50
  # mat_jac$Depth_Var2[mat_jac$Depth_Var2 == 60] <- 50
  # 
  # # # Depth 90 and 120
  # mat_jac$Depth_Var1[mat_jac$Depth_Var1 == 90] <- 100
  # mat_jac$Depth_Var1[mat_jac$Depth_Var1 == 120] <- 100
  # mat_jac$Depth_Var2[mat_jac$Depth_Var2 == 90] <- 100
  # mat_jac$Depth_Var2[mat_jac$Depth_Var2 == 120] <- 100
  # # 
  # mat_jac<-mat_jac[(mat_jac$Depth_Var1==50),]
  # mat_jac<-mat_jac[(mat_jac$Depth_Var2==100),]
  # # 
  # mat_jac$SDist <- abs(mat_jac$Depth_Var1 - mat_jac$Depth_Var2)
  # # 
  # # # All pairs between themselves out
  # mat_jac<-mat_jac[!(mat_jac$value==0),]
  # # # All pairs between the same depth factors
  # mat_jac<-mat_jac[!(mat_jac$SDist==0),]
  # 
  ###### LAetitia request Only compare dissimilariry between 40 and 60, and 90 and 120m #######
  
  
  # Comparing the overlap of 40 or the other depths, with the rest of depths
  # AS 40m is the reference
  # mat_jac<-mat_jac[(mat_jac$Depth_Var1==6 ),]
  # 
  # mat_jac$SDist <- abs(mat_jac$Depth_Var1 - mat_jac$Depth_Var2)
  # 
  # # # All pairs between themselves out
  # mat_jac<-mat_jac[!(mat_jac$value==0),]
  # # # All pairs between the same depth factors
  # mat_jac<-mat_jac[!(mat_jac$SDist==0),]
  
  colnames(mat_jac)[3] <- "Jac_dissimilarity"
  
  
  linm <- lm(Jac_dissimilarity ~ SDist, data = mat_jac)
  linm_resume <- summary(linm)
  plot (Jac_dissimilarity ~ SDist, data = mat_jac, ylab = "Dissimilarity", xlab = "Vertical distance", cex = 0.5, main = paste0(i,sep = "_", "Jaccard"))
  mtext(paste0("R squared: ",round(linm_resume$r.squared,5)),adj = 0)
  mtext(paste0("P-value: ", format.pval(pf(linm_resume$fstatistic[1], # F-statistic
                                           linm_resume$fstatistic[2], # df
                                           linm_resume$fstatistic[3], # df
                                           lower.tail = FALSE))))
  mtext(paste0("y = ",round(linm_resume$coefficients[1],2)," + ",
               round(linm_resume$coefficients[2],5),"x"),adj = 1)
  # 
  # Add the other two components
  # Turnover
  Turn_Depth <- coral.matrices_Depth$beta.jtu
  mat_jac_turn <- as.matrix (Turn_Depth)
  mat_jac_turn_Island <- mat_jac_turn[grepl(i, rownames(mat_jac_turn)),grepl(i, colnames(mat_jac_turn))]
  mat_jac_turn = melt(mat_jac_turn_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames(mat_jac_turn) <- c("Var1","Var2","Jac_turnover")
  # The merge keeping only all mat_bray already cleans and uses 6m as reference 
  mat_jac <- merge(mat_jac, mat_jac_turn, by=c("Var1","Var2"), all.x = T, all.y = F)
  
  # Nestedness
  Nest_Depth <- coral.matrices_Depth$beta.jne
  mat_jac_nest <- as.matrix (Nest_Depth)
  mat_jac_nest_Island <- mat_jac_nest[grepl(i, rownames(mat_jac_nest)),grepl(i, colnames(mat_jac_nest))]
  mat_jac_nest = melt(mat_jac_nest_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames(mat_jac_nest) <- c("Var1","Var2","Jac_nestedness")
  
  # The merge keeping only all mat_jac already cleans and uses 6m as reference 
  mat_jac <- merge(mat_jac, mat_jac_nest, by=c("Var1","Var2"), all.x = T, all.y = F)
  
  
  # Delete all quadrats without corals 
  mat_jac <- na.omit (mat_jac)
  
  # Transform from Dissimilarity to Similarity
  mat_jac$Jac_Similarity <- 1 - mat_jac$Jac_dissimilarity
  # Necessary to measure the respective turnover and nestedness to the new Similarity (inverse dissimilarity)
  mat_jac$Turn_Similarity <- (mat_jac$Jac_turnover * mat_jac$Jac_Similarity) / mat_jac$Jac_dissimilarity
  mat_jac$Nest_Similarity <- (mat_jac$Jac_nestedness * mat_jac$Jac_Similarity) / mat_jac$Jac_dissimilarity
  
  # Compute the index suggested by Seb Villeger - indicator of importance of turnover "Toussaint et al 2014"
  #  mat_jac$pturn = mat_jac$Turn_Similarity / mat_jac$Jac_Similarity
  mat_jac$pturn_Diss = mat_jac$Jac_turnover / mat_jac$Jac_dissimilarity
  
  # Means and Sd in funciton of SDist
  summary_mat_jac <- ddply(mat_jac, .(SDist), summarize, Jaccard=mean(Jac_Similarity), Jac_se=sd(Jac_Similarity) / sqrt(length(Jac_Similarity)),
                           Turnover=mean(Turn_Similarity), Turn_se=sd(Turn_Similarity) / sqrt(length(Turn_Similarity)),
                           Nestedness=mean(Nest_Similarity), Nest_se=sd(Nest_Similarity) / sqrt(length(Nest_Similarity)), 
                           # P_Turn=mean(pturn, na.omit=T), P_Turn_se=sd(pturn, na.omit=T) / sqrt(length(pturn, na.omit=T)),
                           P_Turn_Diss=mean(pturn_Diss), P_Turn_Diss_se=sd(pturn_Diss) / sqrt(length(pturn_Diss)))
  
  
  summary_2 <- melt(summary_mat_jac, id.vars=c("SDist"), measure.vars = c("Jaccard", "Turnover", "Nestedness",  "P_Turn_Diss"))
  colnames (summary_2) <- c("SDist","Component", "Similarity")
  
  summary_se <- melt(summary_mat_jac, id.vars=c("SDist"), measure.vars = c("Jac_se", "Turn_se", "Nest_se",  "P_Turn_Diss_se"))
  colnames (summary_se) <- c("SDist","Component", "Se")
  
  summary_se$Component <- gsub('Jac_se', 'Jaccard',  summary_se$Component)
  summary_se$Component <- gsub('Turn_se', 'Turnover',  summary_se$Component)
  summary_se$Component<- gsub('Nest_se', 'Nestedness',  summary_se$Component)
 # summary_se$Component<- gsub('P_Turn_se', 'P_Turn',  summary_se$Component)
  summary_se$Component<- gsub('P_Turn_Diss_se', 'P_Turn_Diss',  summary_se$Component)
  
  summary <- merge (summary_2,  summary_se)
  
  summary$Island_Site <- i
  
  summary_final = rbind (summary_final, summary)
  
  # For the model 
  mat_jac_all <- melt(mat_jac, id.vars=c("SDist"), measure.vars = c("Jac_Similarity", "Turn_Similarity", "Nest_Similarity"))
  colnames (mat_jac_all) <- c("SDist","Component", "Similarity")
  
  mat_jac_all$Island_Site <- i
  
  Jaccard_final = rbind (Jaccard_final, mat_jac_all)

}


# Measure mean and standard error
summary <- ddply(summary_final, ~ SDist + Component, summarize, Sim_mean=mean(Similarity), Sim_se=sd(Similarity) / sqrt(length(Similarity)))
summary_final <- merge (summary_final,summary)

m1 <- lm(Similarity ~ SDist*Component, data=summary_final)

summary(m1)
a1 <- m1$coef[2] #slope JAccard
a2 <- m1$coef[2] + m1$coef[6] # turnover slope
a3 <- m1$coef[2] + m1$coef[7] # Nestedness slope

b1 <- m1$coef[1]
b2 <- b1  + m1$coef[3] 

# Geom point with colour for island and shape for component.

Fig_2A = summary_final %>% filter(., Component %in% c("Jaccard", "Turnover")) %>% 
  ggplot(., aes(x=SDist, y=Similarity, colour = Component)) +
  scale_color_manual(values = c("#DC143C","#6495ED")) + 
  geom_point (size = 1, shape = 21, alpha = .75) +
  geom_errorbar(aes(ymin = Sim_mean - Sim_se, ymax = Sim_mean + Sim_se), size=1, width=2) +
 # geom_abline(intercept = b1, slope = a1,  colour ="#DC143C") +
 # geom_abline(intercept = b2, slope = a2,  colour ="#6495ED") +
  scale_x_continuous(name ="Vertical distance (m)", limits=c(5,120), breaks = c(25,50,75,100)) +
  scale_y_continuous(name ="Similarity", limits=c(0,0.75), breaks = c(0,0.25,0.5,0.75)) +
  ggtitle("Jaccard - PA") +
  geom_point(aes(y=Sim_mean, shape = Component, colour = Component), size=3, shape = 21, fill = "white") + 
  theme_classic() + theme(axis.text = element_text(size=10), axis.title = element_text(size=11, face="bold"), plot.title = element_text(hjust = 0.5, size=12, face="bold"))
Fig_2A 
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_2A_bis20.pdf", Fig_2A, width = 4, height = 5)

# Test of piecewise/segmented regressions

m1 <- lm(
  Similarity ~ ifelse(SDist<60,SDist-60,0) 
  + ifelse(SDist>=60,SDist-60,0) * Component, data=summary_final)
summary(m1)

qplot(SDist, Similarity, group = SDist > 60, geom = c('point', 'smooth'), 
      method = 'lm', se = F, data = summary_final)

#### Model for predicting a bit the three components ####
# 
# # Tried lmer and brm
# vertical.lmer <- lmer(Dissimilarity ~ SDist*Component + (1|Island_Site), data = Jaccard_final)
# summary(vertical.lmer)
# 
# # Jac = 3.688e-03 * (Vert_dist) + 6.212e-01
# # Turn = 2.026e-03 * (Vert_dist) + (6.212e-01 + (-2.227e-01))
# # Nest = (-5.546e-03) * (Vert_dist) + (6.212e-01 +(- 4.039e-01))
# 
# # No variance as Site is only considered as random intercept
# plot(vertical.lmer) # Data a bit structured. I am not sure if this is good. 
# qqnorm(resid(vertical.lmer))
# qqline(resid(vertical.lmer))
# # Check singularity
# diag.vals <- getME(vertical.lmer,"theta")[getME(vertical.lmer,"lower") == 0]
# any(diag.vals < 1e-6) # FALSE
# 
# # P-values og glmer
# library (nlme)
# nlme_vertical <- lme(Dissimilarity ~ SDist*Component, random = ~1|Island_Site, Jaccard_final)
# anova(nlme_vertical)
# 
# Jaccard_final$fit <- predict(vertical.lmer)
# 
# ggplot(data= Jaccard_final, aes(x=SDist, y=fit,colour = Component))  +
#   geom_point(shape=8)+ facet_wrap(~Island_Site, nrow = 2, strip.position = "top") +
#   scale_x_continuous(name ="Vertical distance (m)", limits=c(5,120), breaks = c(0,20,40,60,80,100,120)) +
#   scale_y_continuous(name ="Jaccard dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
#   theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
#                       axis.text = element_text(size=10, colour="black"),
#                       axis.title = element_text(size=11, face="bold", colour="black")) 
# 
# resume_Jaccard_final <- ddply(Jaccard_final, ~ Island_Site + Component + SDist ,function(x){
#   c(fit=mean(x$fit)) }) 
# 
# ggplot(data= resume_Jaccard_final, aes(x=SDist, y=fit,colour = Component))  +
#   geom_point(shape=8)+ facet_wrap(~Island_Site, nrow = 2, strip.position = "top") +
#   scale_x_continuous(name ="Vertical distance (m)", limits=c(5,120), breaks = c(0,20,40,60,80,100,120)) +
#   scale_y_continuous(name ="Jaccard dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
#   theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
#                       axis.text = element_text(size=10, colour="black"),
#                       axis.title = element_text(size=11, face="bold", colour="black")) 

#### Model for predicting a bit the three components ####

# Change the names on axis
summary_final$SDist<- gsub('114', '6vs120m', summary_final$SDist)
summary_final$SDist<- gsub('84', '6vs90m', summary_final$SDist)
summary_final$SDist<- gsub('54', '6vs60m', summary_final$SDist)
summary_final$SDist<- gsub('34', '6.vs40m', summary_final$SDist)
summary_final$SDist<- gsub('14', '6vs20m', summary_final$SDist)


summary_final$SDist <- factor(summary_final$SDist, levels = c("6vs20m","6.vs40m","6vs60m","6vs90m","6vs120m"))


SFig_3A <- summary_final %>% filter(., Component %in% c("Jaccard", "Turnover")) %>% 
  ggplot(., aes(x=SDist, y=Similarity, fill=Component)) + 
  scale_color_manual(values = c("#DC143C","#6495ED")) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.9) + facet_wrap(~Island_Site, nrow = 2, strip.position = "top") +
  geom_errorbar(aes(ymin = Similarity - Se, ymax = Similarity + Se),size =0.3, width = 0.3, position = position_dodge(0.9))+
  scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.1))+
  ylab ("Similarity") + xlab ("Vertical distance") + 
  theme_bw()  + theme(
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"), 
                      strip.text = element_text(size = 10, colour="black"),
                      panel.spacing = unit(0.5, "lines"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank()) 
SFig_3A
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/SFig_3A.pdf", SFig_3A,width = 14, height = 4)


# Measure of the mean of the index of Villeger 
pturn <- aggregate( Similarity ~ SDist + Component, summary_final, mean )

pturn$Similarity[pturn$Component=="P_Turn_Diss"] 

#Jaccard dissimilarity
1 - pturn$Similarity[pturn$Component=="Jaccard"] 


##################

#### Getting the heatmap at the Quadrats level ####
island_site <-  c("Bora_1", "Bora_2", "Gambier_1", "Gambier_2", "Makatea_1", "Makatea_2", "Moorea_1", "Moorea_2", "Rangiroa_1",
                  "Rangiroa_2", "Raroia_1", "Raroia_2", "Tahiti_1", "Tahiti_2", "Tikehau_1", "Tikehau_2")

depth <- c(6,20,40,60,90,120)

summary_final <- data.frame()
Jaccard_final  <- data.frame()
for (i in unique (island_site)){
  Beta_Depth <-  coral.matrices_Depth$beta.jac
  mat_jac <- as.matrix (Beta_Depth)
  mat_jac_Island <- mat_jac[grepl(i, rownames(mat_jac)),grepl(i, colnames(mat_jac))]
  
  mat_jac = melt(mat_jac_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames (mat_jac) <- c("Var1", "Var2", "nobserv")
  
  mat_jac$Depth_Var1 <- as.numeric (gsub("(.*_){2}(\\d+)_.+", "\\2", mat_jac$Var1))
  mat_jac$Depth_Var2 <- as.numeric (gsub("(.*_){2}(\\d+)_.+", "\\2", mat_jac$Var2))
  
  
  # Comparing the overlap of 40 or the other depths, with the rest of depths
  # AS 40m is the reference
  mat_jac<-mat_jac[(mat_jac$Depth_Var1==40 ),]
  
  mat_jac$SDist <- (mat_jac$Depth_Var1 - mat_jac$Depth_Var2)
  
  # # All pairs between themselves out
  mat_jac<-mat_jac[!(mat_jac$nobserv==0),]
  # # All pairs between the same depth factors
  mat_jac<-mat_jac[!(mat_jac$SDist==0),]
  
  colnames(mat_jac)[3] <- "Jac_dissimilarity"
  
  
  # Add the other two components
  # Turnover
  Turn_Depth <- coral.matrices_Depth$beta.jtu
  mat_jac_turn <- as.matrix (Turn_Depth)
  mat_jac_turn_Island <- mat_jac_turn[grepl(i, rownames(mat_jac_turn)),grepl(i, colnames(mat_jac_turn))]
  mat_jac_turn = melt(mat_jac_turn_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames (mat_jac_turn) <- c("Var1", "Var2", "Jac_turnover")
  # The merge keeping only all mat_bray already cleans and uses 6m as reference 
  mat_jac <- merge(mat_jac, mat_jac_turn, by=c("Var1","Var2"), all.x = T, all.y = F)
  
  # Nestedness
  Nest_Depth <- coral.matrices_Depth$beta.jne
  mat_jac_nest <- as.matrix (Nest_Depth)
  mat_jac_nest_Island <- mat_jac_nest[grepl(i, rownames(mat_jac_nest)),grepl(i, colnames(mat_jac_nest))]
  mat_jac_nest = melt(mat_jac_nest_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames(mat_jac_nest) <- c("Var1", "Var2", "Jac_nestedness")
  # The merge keeping only all mat_jac already cleans and uses 6m as reference 
  mat_jac <- merge(mat_jac, mat_jac_nest, by=c("Var1","Var2"), all.x = T, all.y = F)
  
  
  # Delete all quadrats without corals 
  mat_jac <- na.omit (mat_jac)
  
  # Compute the index suggested by Seb Villeger - indicator of importance of turnover "Toussaint et al 2014"
  #  mat_jac$pturn = mat_jac$Turn_Similarity / mat_jac$Jac_Similarity
  mat_jac$pturn_Diss = mat_jac$Jac_turnover / mat_jac$Jac_dissimilarity
  
  # Means and Sd in funciton of SDist
  summary_mat_jac <- ddply(mat_jac, .(Depth_Var1,Depth_Var2,SDist), summarize, Jaccard=mean(Jac_dissimilarity), Jac_se=sd(Jac_dissimilarity) / sqrt(length(Jac_dissimilarity)),
                           Turnover=mean(Jac_turnover), Turn_se=sd(Jac_turnover) / sqrt(length(Jac_turnover)),
                           Nestedness=mean(Jac_nestedness), Nest_se=sd(Jac_nestedness) / sqrt(length(Jac_nestedness)), 
                           # P_Turn=mean(pturn, na.omit=T), P_Turn_se=sd(pturn, na.omit=T) / sqrt(length(pturn, na.omit=T)),
                           P_Turn_Diss=mean(pturn_Diss), P_Turn_Diss_se=sd(pturn_Diss) / sqrt(length(pturn_Diss)))
  
  
  summary_2 <- melt(summary_mat_jac, id.vars=c("Depth_Var1","Depth_Var2","SDist"), measure.vars = c("Jaccard", "Turnover", "Nestedness",  "P_Turn_Diss"))
  colnames (summary_2) <- c("Depth_Var1","Depth_Var2","SDist","Component", "Dissimilarity")
  
  summary_se <- melt(summary_mat_jac, id.vars=c("Depth_Var1","Depth_Var2","SDist"), measure.vars = c("Jac_se", "Turn_se", "Nest_se",  "P_Turn_Diss_se"))
  colnames (summary_se) <- c("Depth_Var1","Depth_Var2","SDist","Component", "Se")
  
  summary_se$Component <- gsub('Jac_se', 'Jaccard',  summary_se$Component)
  summary_se$Component <- gsub('Turn_se', 'Turnover',  summary_se$Component)
  summary_se$Component<- gsub('Nest_se', 'Nestedness',  summary_se$Component)
  # summary_se$Component<- gsub('P_Turn_se', 'P_Turn',  summary_se$Component)
  summary_se$Component<- gsub('P_Turn_Diss_se', 'P_Turn_Diss',  summary_se$Component)
  
  summary <- merge (summary_2,  summary_se)
  
  keep <- c("Jaccard", "Turnover","P_Turn_Diss")
  summary <- filter(summary, Component %in% keep)
  
  summary$Island_Site <- i
  
  summary_final = rbind (summary_final, summary)
  
}


# Measure mean and standard error
summary <- ddply(summary_final, ~ SDist + Component, summarize, Sim_mean=mean(Dissimilarity), Sim_se=sd(Dissimilarity) / sqrt(length(Dissimilarity)))
summary_final <- merge (summary_final,summary)


# Change the names on axis
summary_final$SDist<- paste(summary_final$Depth_Var1,"vs",summary_final$Depth_Var2)

# Run the respective row for each depth only
# summary_final$SDist <- factor(summary_final$SDist, levels = c("6 vs 20","6 vs 40","6 vs 60","6 vs 90","6 vs 120"))
# summary_final$SDist <- factor(summary_final$SDist, levels = c("20 vs 6","20 vs 40","20 vs 60","20 vs 90","20 vs 120"))
summary_final$SDist <- factor(summary_final$SDist, levels = c("40 vs 6","40 vs 20","40 vs 60","40 vs 90","40 vs 120"))
# summary_final$SDist <- factor(summary_final$SDist, levels = c("60 vs 6","60 vs 20","60 vs 40","60 vs 90","60 vs 120"))
# summary_final$SDist <- factor(summary_final$SDist, levels = c("90 vs 6","90 vs 20","90 vs 40","90 vs 60","90 vs 120"))
# summary_final$SDist <- factor(summary_final$SDist, levels = c("120 vs 6","120 vs 20","120 vs 40","120 vs 60","120 vs 90"))


# Measure the mean fo all locations
pturn <- aggregate(Dissimilarity ~ SDist + Component, summary_final, mean )

pturn$Dissimilarity[pturn$Component=="P_Turn_Diss"] 

pturn$Dissimilarity[pturn$Component=="Jaccard"] 

pturn$Dissimilarity[pturn$Component=="Turnover"] 


### Necessary to fill a table manually!!!! And make the calculation for each depth separately! later open the New_database

pturn <- read.csv(file = "Data/Pturn.csv", header = T, dec = ".", sep = ";", row.names = 1)

pturn$Depth1 <- as.factor (pturn$Depth1)
pturn$Depth2 <- as.factor (pturn$Depth2)

pturn$Depth2 <- factor(pturn$Depth2, levels = c("6","20","40","60","90", "120"))


pturn$Quart_Pturn_Jac <- ntile(pturn$Pturn_Jac, 3)

pturn$Quart_Pturn_Jac <- factor(pturn$Pturn_Jac, levels=c("0.5","0.6","0.7","0.8","0.9","1") )

# Jaccard
jac <- ggplot(pturn, aes(Depth2, Depth1, fill=Jaccard))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.79, limit=c(0.6,1), space="Lab", name="Dissimilarity")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="Depth (m)",y="Depth (m)", title = "Jaccard") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))


# Pturnover Jaccard
pturn_jac <- ggplot(pturn, aes(Depth2, Depth1, fill=Turnover_Jac))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.69, limit=c(0.35,1), breaks = c(0.4,0.5,0.6,0.7,0.8,0.9,1), space="Lab", name="P Turnover")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="",y="Depth (m)", title = "Replacement") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))


# Bray
bray <- ggplot(pturn, aes(Depth2, Depth1, fill=Bray))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.79, limit=c(0.6,1), space="Lab", name="Dissimilarity")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="Depth (m)",y="Depth (m)", title = "Bray-Curtis") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))

# Pturnover Bray
pturn_bray <- ggplot(pturn, aes(Depth2, Depth1, fill=Pturn_Bray))+
  geom_tile(color="white")+
  scale_fill_gradient2(low="lemonchiffon", mid = "orange", high="red",  midpoint = 0.69, limit=c(0.35,1),breaks = c(0.4,0.5,0.6,0.7,0.8,0.9,1), space="Lab", name="P Turnover")+
  scale_x_discrete(position = "top")+scale_y_discrete(limits=rev)+
  theme_bw()+
  coord_fixed()+
  labs(x="",y="Depth (m)", title = "Substitution") +
  guides(fill = guide_legend(reverse = FALSE))+
  theme(panel.grid = element_blank(), axis.title = element_text(size=11, face="bold", colour="black"), 
        strip.text = element_text(size = 10, colour="black"))


# Figure heat map
# Plot both graphs together. Jaccard and Bray
(Figure_heatmap = jac +  bray  + pturn_jac + pturn_bray +  plot_layout(guides = 'collect' )) 


ggsave( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Figure_heatmap.pdf", Figure_heatmap,  width = 8, height = 8)

#### Getting the heatmap at the Quadrats level ####










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


ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  cex= 0.5, horiz = F)

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
                             beta_jac=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)),
                             Se = c(sd(Beta_Depth_Matrix_Depth_1) / sqrt(length(Beta_Depth_Matrix_Depth_1)), sd(Beta_Depth_Matrix_Depth_2) / sqrt(length(Beta_Depth_Matrix_Depth_2)),sd(Beta_Depth_Matrix_Depth_3) / sqrt(length(Beta_Depth_Matrix_Depth_3)),
                                    sd(Beta_Depth_Matrix_Depth_4) / sqrt(length(Beta_Depth_Matrix_Depth_4)),sd(Beta_Depth_Matrix_Depth_5) / sqrt(length(Beta_Depth_Matrix_Depth_5)),sd(Beta_Depth_Matrix_Depth_6) / sqrt(length(Beta_Depth_Matrix_Depth_6))))
beta_div_depth # These are the values for the table


# check linear model for Jaccard
beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_jac ~ Depth,beta_div_depth ))

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

ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), cex= 0.5, horiz = F)
title ("Jaccard distance (turnover) - PA")


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

ordiplot(mod_Depth,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(mod_Depth,display="species",col="red",air=0.01, cex =1)
legend(x="topright", y="topleft", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  cex= 0.5, horiz = F)
title ("Jaccard distance (nestedness) - PA")


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

Locations <- read.csv(file = "Data/Deephope_sampling_locations_RAN*.csv", header = T, dec = ".", sep = ";", row.names = 1)
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


# Create dataframe of three columns: Site1, Site2, SDist
Spa_Dist <- as.matrix (dm)
Spa_Dist = melt(Spa_Dist, measure.vars="Dist", na.rm=FALSE)
colnames(Spa_Dist) <- c("Site1", "Site2", "Spa_Distance")
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
pdf("~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Boxplot_Jac.pdf", bg = "white", width = 6, height = 4) 
boxplot (Beta_Depth_Matrix_Depth_1,Beta_Depth_Matrix_Depth_2, Beta_Depth_Matrix_Depth_3,Beta_Depth_Matrix_Depth_4,Beta_Depth_Matrix_Depth_5, Beta_Depth_Matrix_Depth_6, 
         xlab = "Depth (m)",
         ylab = "Beta.jaccard",
         names = c("6","20","40","60","90", "120"), 
         main = "Jaccard - Composition")
dev.off()



# Introduce values manually: 
beta_jacc_depth # Very important to verify they come from the model considering jaccard (not turnover or nestedness)

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_jac=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))


beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_jac~ Depth,beta_div_depth ))


############################################
# Get the entire beta.pair for plot
Beta_Depth <- coral.matrices_Depth$beta.jac
Beta_Depth_Matrix <- as.matrix (Beta_Depth)

Jac_spa = melt(Beta_Depth_Matrix, measure.vars="Jac", na.rm=FALSE)

Jac_spa$Depth_Var1 <- as.numeric (str_extract(Jac_spa$Var1, "[^_]+"))
Jac_spa$Depth_Var2 <- as.numeric (str_extract(Jac_spa$Var2, "[^_]+"))

# Keep only pairs between same depths
Jac_spa$VDist <- abs(Jac_spa$Depth_Var1 - Jac_spa$Depth_Var2)
Jac_spa<-Jac_spa[(Jac_spa$VDist==0),]

# Delete pairs that are the same "6_Bora_1" vs "6_Bora_1"
same <- which(!(Jac_spa$Var1==Jac_spa$Var2))
Jac_spa<-Jac_spa[same, ]

# Change names and only keep the necessary columns 
colnames(Jac_spa)[3] <- "Jac_dissimilarity"
colnames(Jac_spa)[4] <- "Depth"
Jac_spa <- Jac_spa[c("Var1","Var2","Jac_dissimilarity","Depth")]

ggplot(data = Jac_spa, aes(y = Jac_dissimilarity, x = as.factor(Depth))) + 
  geom_boxplot()  + 
  ylab ("Jac_dissimilarity") + xlab ("Depth (m)") + 
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 

# Add the other two components
# Turnover
Turn_Depth <- coral.matrices_Depth$beta.jtu
Turn_Depth_Matrix <- as.matrix (Turn_Depth)

Turn_spa = melt(Turn_Depth_Matrix, measure.vars="Turn", na.rm=FALSE)
colnames(Turn_spa)[3] <- "Jac_turnover"
# The merge keeps already only the rows we want based on "Var1" and "Var2"
Jac_spa <- merge(Jac_spa, Turn_spa, by=c("Var1","Var2"), all.x = T, all.y = F)

# Nestedness
Nest_Depth <- coral.matrices_Depth$beta.jne
Nest_Depth_Matrix <- as.matrix (Nest_Depth)

Nest_spa = melt(Nest_Depth_Matrix, measure.vars="Nest", na.rm=FALSE)
colnames(Nest_spa)[3] <- "Jac_nestedness"
# The merge keeps already only the rows we want based on "Var1" and "Var2"
Jac_spa <- merge(Jac_spa, Nest_spa, by=c("Var1","Var2"), all.x = T, all.y = F)



# Means and Se in funciton of Depth
summary_Jac_spa <- ddply(Jac_spa, .(Depth), summarize, Jac=mean(Jac_dissimilarity), Jac_se=sd(Jac_dissimilarity) / sqrt(length(Jac_dissimilarity)),
                          Turn=mean(Jac_turnover), Turn_se=sd(Jac_turnover) / sqrt(length(Jac_turnover)),
                          Nest=mean(Jac_nestedness), Nest_se=sd(Jac_nestedness) / sqrt(length(Jac_nestedness)))

summary_3 <- melt(summary_Jac_spa, id.vars=c("Depth"), measure.vars = c("Jac", "Turn", "Nest"))
colnames (summary_3) <- c("Depth","Component", "Dissimilarity")

summary_se_3 <- melt(summary_Jac_spa, id.vars=c("Depth"), measure.vars = c("Jac_se", "Turn_se", "Nest_se"))
colnames (summary_se_3) <- c("Depth","Component", "Se")

summary_se_3$Component <- gsub('Jac_se', 'Jac', summary_se_3$Component)
summary_se_3$Component <- gsub('Turn_se', 'Turn', summary_se_3$Component)
summary_se_3$Component<- gsub('Nest_se', 'Nest', summary_se_3$Component)

summary_spa <- merge (summary_3, summary_se_3)


# Plot of the three components from summary 
dis_vs_spa <- ggplot(data= summary_spa, aes(x=Depth, y=Dissimilarity, shape = Component, colour = Component))  +
  geom_point(aes( shape = Component, colour = Component), size=4) + 
  geom_errorbar(aes(ymin = Dissimilarity - Se, ymax = Dissimilarity + Se), size=0.5, width=2, colour = "black") +
  scale_x_continuous(name ="Depth (m)", limits=c(5,120), breaks = c(0,20,40,60,80,100,120)) +
  scale_y_continuous(name ="Jac dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"))
dis_vs_spa 


ggplot(data= summary_spa, aes(x=Depth, y=Dissimilarity, shape = Component, colour = Component))  +
  geom_point(aes( shape = Component, colour = Component), size=4) + facet_wrap(~Component, nrow = 3)+
  geom_errorbar(aes(ymin = Dissimilarity - Se, ymax = Dissimilarity + Se), size=1, width=2, colour = "black") +
  scale_x_continuous(name ="Depth (m)", limits=c(5,120), breaks = c(0,20,40,60,80,100,120)) +
  scale_y_continuous(name ="Jac dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"))


# Working with the entire database pairs for boxplot and mean and se
Jac_spa_all <- melt(Jac_spa, id.vars=c("Depth"), measure.vars = c("Jac_dissimilarity", "Jac_turnover"))
colnames (Jac_spa_all) <- c("Depth","Component", "Dissimilarity")

# Measure mean and standard error
summary <- ddply(Jac_spa_all, ~Depth + Component, summarize, Dissimilarity_mean=mean(Dissimilarity), Se=sd(Dissimilarity) / sqrt(length(Dissimilarity)))
Jac_spa_all <- merge (Jac_spa_all,summary)

SFig_6A <- ggplot(data = Jac_spa_all, aes(y = Dissimilarity, x = as.factor(Depth), colour = Component)) + 
  geom_boxplot()  + 
  scale_color_manual(values = c("#DC143C","#6495ED", "#6B8E23")) + 
  geom_errorbar(aes(ymin = Dissimilarity_mean - Se, ymax = Dissimilarity_mean + Se), size=1, width=0.1, colour = "black") +
  geom_point(aes(y = Dissimilarity_mean), size=1.5, fill = "white", shape = 21) +
  facet_wrap(~Component, nrow = 3, scales = "free")+
  scale_y_continuous(name ="Dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  ylab ("Dissimilarity") + xlab ("Depth (m)") + ggtitle ("Jaccard - Presence/Absence") +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), 
                          strip.background = element_blank(),
                          strip.text.x = element_blank(), 
                          legend.position = "right") 

SFig_6A
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/SFig_6A.pdf", SFig_6A ,width = 6, height = 6)


summary.lm(aov(Dissimilarity_mean ~ Depth*Component, Jac_spa_all))

spatial.lmer <- lmer(Dissimilarity ~ Depth*Component  + (1|Component) , data = Jac_spa_all)
summary (spatial.lmer)



# Make a final figure: beta vs geographical distance (x=geo_dist, y=beta, points= pairs for all different given depths to see if it's correlated to 
# geographical distance. We know they are not because 

# From Jac_Spa, all pairs

Jac_spa$Site1 <- sub("^[^_]*_", "", Jac_spa$Var1)
Jac_spa$Site2 <- sub("^[^_]*_", "", Jac_spa$Var2)


# Merge the two dataframes of pairs and geographical distance
Spatial_Jaccard <- merge (Jac_spa, Spa_Dist )

# Same values as Beta_Div_Depth that are the ones I used for Mantel tests
summary_Spatial_Jaccard <- ddply(Spatial_Jaccard, ~Depth, summarize, Jaccard_mean=mean(Jac_dissimilarity), Se=sd(Jac_dissimilarity) / sqrt(length(Jac_dissimilarity)))
Spatial_Jaccard <- merge (Spatial_Jaccard,summary_Spatial_Jaccard)

# Colours from fishualize as NMDS
colours <- c("#CBEAF1FF", "#73C2E1FF", "#3F73B1FF", "#28407AFF", "#222E47FF", "#191819FF")

SFig_5A <- ggplot(data=Spatial_Jaccard,aes(x=Spa_Distance, y=Jac_dissimilarity, colour = as.factor(Depth)))  +
  geom_point (size = 0.5)+
  facet_wrap(~Depth, nrow = 1, scales = "free")+
  # geom_line(aes(y=Jaccard_mean,colour = as.factor(Depth)), size=1) + 
  # geom_line(aes(y=Jaccard_mean-Se,colour = as.factor(Depth)), size=0.3, linetype="dotted",alpha = 0.8) +
  # geom_line(aes(y=Jaccard_mean+Se,colour = as.factor(Depth)), size=0.3, linetype="dotted",alpha = 0.8) +
  scale_colour_manual(values = c("#CBEAF1FF", "#73C2E1FF", "#3F73B1FF", "#28407AFF", "#222E47FF", "#191819FF")) +
  scale_x_continuous(name ="Geographical distance (km)", limits=c(0,2000), breaks = c(0,500,1000,1500,2000)) +
  scale_y_continuous(name ="Jaccard dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  theme_classic()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      strip.text.x = element_text(size = 12, color = "black", face = "bold"),
                      strip.background = element_rect(fill="grey"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none") 

SFig_5A
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/SFig_5A.pdf", SFig_5A ,width = 12, height = 3)


################## Spatial beta diversity ##########################3
 









