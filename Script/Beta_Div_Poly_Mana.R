

###########################  Start of script   ################################

rm (list = ls()) 

par(mfrow=c(1,1))

library (ggplot2); library (data.table); library  (stats); library (dynRB); library(dplyr); library (matrixStats); library(plyr); library (stringr)
library (vegan); library (goeveg); library (reshape2); library (scatterplot3d); library(betapart)

# From the beggining once again 
setwd("~/Documents/AAASea_Science/Publications/Coral Dynamics of French Polynesia/Final_Data_Scripts_Figures_Gonzalo")
poly_mana <- read.csv(file = "Benthos_Polynesia_Mana_new_June.csv", header = T, dec = ".", sep = ";")
poly_mana <- data.table(poly_mana, keep.rownames=FALSE, check.names=FALSE, key=NULL)

# Keep only islands we want
colnames (poly_mana) = c("Date","Year","Island","Site","Island_Site","Observer","Quadrat","Genus","Coral_cover")

keep_sites <-  c ("Tahiti Faaa", "Moorea Entre 2 Baies", "Nengo Nengo", "Raiatea", "Tetiaroa", "Tikehau", "Tubuai")
poly_mana <- poly_mana[poly_mana$Island_Site %in% keep_sites, ]


# Change some names because of typing mistakes. Here we do not transform in others
summary (as.factor (poly_mana$Genus))
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Synarea ", "Porites")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Synarea", "Porites")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Porites nigrescens", "Porites")     # Remember to mention mostly Porites nigrescens
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Porites rus", "Porites")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Poritesrus", "Porites")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Dipsastraea rotumana", "Dipsastrea")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Dipsastraea", "Dipsastrea")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Goniastrea stelligera (ex Favia stelligera)", "Goniastrea")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Porites ", "Porites")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Porites ", "Porites")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Leptatstrea", "Leptastrea")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Astrea curta", "Astrea")
poly_mana$Genus <- str_replace_all(poly_mana$Genus, "Autres", "Others")




# We still have the problem with some locations (now it doesn't change anything)
poly_mana <- setDT(poly_mana)[Coral_cover < 0, Coral_cover := 0]






# Make aggregate for each of the 20 quadrats
poly_mana_genus_sites <-  aggregate(Coral_cover ~ Island_Site + Genus + Year , poly_mana, mean)

# Reshape data with melt data so that each row is a unique id-variable combination f
poly_mana_genus_sites_Final <- melt (poly_mana_genus_sites, id.vars = c ("Genus", "Island_Site", "Year"), na.rm = TRUE, measure.vars = c("Coral_cover"), value.name = c("Coral_cover"))

colnames (poly_mana_genus_sites_Final) <- c("Genus","Island_Site","Year","variable","Coral_cover")
# Cast the data for having details by genus 
poly_mana_genus_sites_Final <- dcast (poly_mana_genus_sites_Final, Island_Site + Year ~ Genus, fun=sum, value.var = c("Coral_cover"))

poly_mana_genus_sites_Final[is.na(poly_mana_genus_sites_Final)] <- 0



# Different time periods now
# Select years of the study for the before-after 
df <- poly_mana_genus_sites_Final
df$Period = df$Year
for(i in 1:length(df$Period)){
  # cat(df$Period[i])
  if (df$Period [i] > 1993 & df$Period [i] < 1996){
    df$Period[i] <- "1st_Period"}
  else if (df$Period [i] > 2003 & df$Period [i] < 2006) {
    df$Period[i] <- "2nd_Period"}
  else if (df$Period [i] > 2015 & df$Period [i] < 2018) {
    df$Period[i] <- "3rd_Period"}
}
poly_mana_final_2 <- df



# Keep only rows of periods
keep_rows <-  c ("1st_Period", "2nd_Period", "3rd_Period")
poly_mana_final_2 <- poly_mana_final_2[poly_mana_final_2$Period %in% keep_rows, ]

# Remove Tubuai because we do not have the three periods
poly_mana_final_2 <- poly_mana_final_2[- grep("Tubuai", poly_mana_final_2$Island_Site),]
poly_mana_final_2 <- poly_mana_final_2[- grep("Nengo Nengo", poly_mana_final_2$Island_Site),]

poly_mana_genus_sites_Final <- poly_mana_final_2 

poly_mana_genus_sites_Final$Island_Site <- str_replace_all (poly_mana_genus_sites_Final$Island_Site, "Moorea Entre 2 Baies", "MOOREA")
poly_mana_genus_sites_Final$Island_Site <- str_replace_all (poly_mana_genus_sites_Final$Island_Site, "Tahiti Faaa", "TAHITI")
poly_mana_genus_sites_Final$Island_Site <- str_replace_all (poly_mana_genus_sites_Final$Island_Site, "Tetiaroa", "TETIAROA")
poly_mana_genus_sites_Final$Island_Site <- str_replace_all (poly_mana_genus_sites_Final$Island_Site, "Raiatea", "RAIATEA")
poly_mana_genus_sites_Final$Island_Site <- str_replace_all (poly_mana_genus_sites_Final$Island_Site, "Tikehau", "TIKEHAU")



poly_mana_genus_sites_Final$ID<- with(poly_mana_genus_sites_Final, paste0(Island_Site, sep = "_",Year))
rownames (poly_mana_genus_sites_Final) <- poly_mana_genus_sites_Final$ID


### Plot by periods

poly_mana_genus_sites_Final$ID<- with(poly_mana_genus_sites_Final, paste0(Island_Site, sep = "_",Period))
rownames (poly_mana_genus_sites_Final) <- poly_mana_genus_sites_Final$ID


poly_mana_genus_sites_Final <- subset(poly_mana_genus_sites_Final, select= -c(Island_Site,Year,Period,ID))


# Start calculation of nmds or beta-diversity
dis <- vegdist (poly_mana_genus_sites_Final)
total_NMDS <- metaMDS(poly_mana_genus_sites_Final, k=2, trymax = 1000, distance = "bray")     # With two dimensions we are already below 0.2  

# Prepare for ggplot 
CC_NMDS_Species <- total_NMDS[['species']]
CC_NMDS_Species <- as.data.frame(CC_NMDS_Species)
CC_NMDS_Species$Species <- rownames (CC_NMDS_Species)

CC_NMDS_coordinates <- total_NMDS[['points']]
CC_NMDS_coordinates <- as.data.frame(CC_NMDS_coordinates)
CC_NMDS_coordinates$ID <- rownames (CC_NMDS_coordinates)
# Extract Island
CC_NMDS_coordinates$Island <- sub("\\_.*", "", CC_NMDS_coordinates$ID)
# Extract Periods
CC_NMDS_coordinates$Period <- sub("^[^_]*_", "", CC_NMDS_coordinates$ID)
# Remove ID column 
CC_NMDS_coordinates <- subset(CC_NMDS_coordinates, select=-c(ID))





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


# Plot of beta-diversity - NMDS by periods
CC_NMDS_coordinates = CC_NMDS_coordinates %>% mutate(Period = factor(Period, levels = c("1st_Period","2nd_Period","3rd_Period")))
hull_CC_NMDS <- CC_NMDS_coordinates %>% group_by(Period) %>% slice(chull(MDS1,MDS2))


Bray_Curtis_Poly_Mana <- CC_NMDS_coordinates %>% ggplot() + theme_blank() + 
  geom_point(aes(x = MDS1, y = MDS2, fill = Period), shape = 21, show.legend = F) +
  geom_polygon(data = hull_CC_NMDS, aes(x = MDS1, y = MDS2, fill = Period), alpha = 0.75, color = "black", show.legend = T) +
  #geom_text(data = CC_NMDS_Species, aes(x = MDS1, y = MDS2, label = Species), size = 2) +
  geom_label(data = CC_NMDS_Species, aes(x = MDS1, y = MDS2, label = Species), size = 2, alpha = .75) +
  scale_x_continuous(name ="NMDS 1") + scale_y_continuous(name ="NMDS 2") + ggtitle("Bray-Curtis - Coral cover") +
  fishualize::scale_fill_fish_d(option = "Ostracion_whitleyi", direction = -1) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=11, face="bold")) +
  annotate(geom = 'text', label = 'Stress = 0.15', 
           x = range(CC_NMDS_coordinates$MDS1)[1], y = range(CC_NMDS_Species$MDS2)[1], hjust = 0, vjust = 6)

 


##############

# Now with Jaccard Presence/Absence
poly_mana_genus_sites_Final_PA <- poly_mana_genus_sites_Final
poly_mana_genus_sites_Final_PA[poly_mana_genus_sites_Final_PA > 0] <- 1 



# Working with PA only 
PA_NMDS <- metaMDS(poly_mana_genus_sites_Final_PA, k=2, trymax = 1000, distance = "jaccard") 

# Prepare for ggplot 
PA_NMDS_Species <- PA_NMDS[['species']]
PA_NMDS_Species <- as.data.frame(PA_NMDS_Species)
PA_NMDS_Species$Species <- rownames (PA_NMDS_Species)

PA_NMDS_coordinates <- PA_NMDS[['points']]
PA_NMDS_coordinates <- as.data.frame(PA_NMDS_coordinates)
PA_NMDS_coordinates$ID <- rownames (PA_NMDS_coordinates)
# Extract Island
PA_NMDS_coordinates$Island <- sub("\\_.*", "", PA_NMDS_coordinates$ID)
# Extract Periods
PA_NMDS_coordinates$Period <- sub("^[^_]*_", "", PA_NMDS_coordinates$ID)
# Remove ID column 
PA_NMDS_coordinates <- subset(PA_NMDS_coordinates, select=-c(ID))



PA_NMDS_coordinates = PA_NMDS_coordinates %>% mutate(Period = factor(Period, levels = c("1st_Period","2nd_Period","3rd_Period")))
hull_PA_NMDS <- PA_NMDS_coordinates %>% group_by(Period) %>% slice(chull(MDS1,MDS2))


Jaccard_Poly_Mana <- PA_NMDS_coordinates %>% ggplot() + theme_blank() + 
  geom_point(aes(x = MDS1, y = MDS2, fill = Period), shape = 21, show.legend = F) +
  geom_polygon(data = hull_PA_NMDS, aes(x = MDS1, y = MDS2, fill = Period), alpha = 0.75, color = "black", show.legend = T) +
  # geom_text(data = PA_NMDS_Species, aes(x = MDS1, y = MDS2, label = Species), size = 2) +
  geom_label(data = PA_NMDS_Species, aes(x = MDS1, y = MDS2, label = Species), size = 2, alpha = .75) +
  scale_x_continuous(name ="NMDS 1") + scale_y_continuous(name ="NMDS 2") + ggtitle("Jaccard - Presence/Absence") +
  fishualize::scale_fill_fish_d(option = "Ostracion_whitleyi", direction = -1) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=11, face="bold")) +
  annotate(geom = 'text', label = 'Stress = 0.11', 
           x = range(PA_NMDS_coordinates$MDS1)[1], y = range(PA_NMDS_Species$MDS2)[1], hjust = 0, vjust = 7)


(Sup_Fig_7 = Jaccard_Poly_Mana + Bray_Curtis_Poly_Mana + plot_layout(guides = 'collect'))

ggsave( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Sup_Fig_7.pdf", Sup_Fig_7,  width = 10, height = 5)


# Alpha richness

poly_mana_genus_sites_Final_PA$Richness <- rowSums(poly_mana_genus_sites_Final_PA)

poly_mana_genus_sites_Final_PA <- subset(poly_mana_genus_sites_Final_PA, select=c(Richness))



# Decompose Bray_Curtis
Bray_Decomposition <- beta.pair.abund(poly_mana_genus_sites_Final, index.family = "bray")

mean (Bray_Decomposition$beta.bray)
mean (Bray_Decomposition$beta.bray.bal)
mean (Bray_Decomposition$beta.bray.gra)



# Bray
Beta <- Bray_Decomposition$beta.bray
groups <- c(rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1))

mod_Depth <- betadisper (Beta,groups, type = "median", bias.adjust=TRUE)
mod_Depth
plot (mod_Depth)

colours <- c("blue", "red","green")

ordiplot(mod_Depth,type="n", choices = c(1,2)) 
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
legend(x="topright", y="topleft", legend=c ("1st", "2nd", "3rd"), col=c("blue", "red","green"), fill = c("blue", "red","green"),  cex= 0.5, horiz = F)


distances_centroid <- as.data.frame(mod_Depth$distances)
colnames (distances_centroid) <- "Distances"
distances_centroid$Period <- sapply(strsplit(rownames(distances_centroid), "_"), "[", 2)
mean <- aggregate (Distances ~ Period, distances_centroid,"mean")

ggplot(data = distances_centroid, aes(y = Distances, x = Period)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Period), colour = "red", size = 3) + 
  ylab ("Distance to median - Bray") + xlab ("Period") + ggtitle ("Bray Coral cover") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 



# Balance 
Beta_bal <- Bray_Decomposition$beta.bray.bal
groups <- c(rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1))

mod_Depth <- betadisper (Beta_bal,groups, type = "median", bias.adjust=TRUE)
mod_Depth
plot (mod_Depth)

ordiplot(mod_Depth,type="n", choices = c(1,2)) 
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
legend(x="topright", y="topleft", legend=c ("1st", "2nd", "3rd"), col=c("blue", "red","green"), fill = c("blue", "red","green"),  cex= 0.5, horiz = F)


distances_centroid <- as.data.frame(mod_Depth$distances)
colnames (distances_centroid) <- "Distances"
distances_centroid$Period <- sapply(strsplit(rownames(distances_centroid), "_"), "[", 2)
mean <- aggregate (Distances ~ Period, distances_centroid,"mean")

ggplot(data = distances_centroid, aes(y = Distances, x = Period)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Period), colour = "red", size = 3) + 
  ylab ("Distance to median - Bray Balance") + xlab ("Period") + ggtitle ("Bray Coral cover") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 


# Gradient
Beta_gra <- Bray_Decomposition$beta.bray.gra
groups <- c(rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1),rep("1st",1),rep("2nd",1), rep("3rd",1))

mod_Depth <- betadisper (Beta_gra,groups, type = "median", bias.adjust=TRUE)
mod_Depth
plot (mod_Depth)

ordiplot(mod_Depth,type="n", choices = c(1,2)) 
orditorp(mod_Depth,display="sites",labels =T)
ordispider(mod_Depth, groups=groups, col = colours,cex=0.6, lwd = 1.5)
ordihull(mod_Depth,groups=groups,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
legend(x="topright", y="topleft", legend=c ("1st", "2nd", "3rd"), col=c("blue", "red","green"), fill = c("blue", "red","green"),  cex= 0.5, horiz = F)

distances_centroid <- as.data.frame(mod_Depth$distances)
colnames (distances_centroid) <- "Distances"
distances_centroid$Period <- sapply(strsplit(rownames(distances_centroid), "_"), "[", 2)
mean <- aggregate (Distances ~ Period, distances_centroid,"mean")

ggplot(data = distances_centroid, aes(y = Distances, x = Period)) + 
  geom_boxplot() + geom_point (data = mean, aes(y = Distances, x = Period), colour = "red", size = 3) + 
  ylab ("Distance to median - Bray Gradient") + xlab ("Period") + ggtitle ("Bray Coral cover") +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 




