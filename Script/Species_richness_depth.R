
### Demonstrating correlation between generic and species richness"
# "Gonzalo"


  

require(reshape); require (reshape2); require (data.table); require(dplyr); require (tidyr); require (plyr); require (ggplot2);  require  (stats); require (dynRB); require (matrixStats); require (stringr); require(gridExtra); require (patchwork)



## R Markdown to show that genus is correlated to species level

### From Rocha database (Coral trait - database)
 

 
rm(list = ls ())

Max_Min <- read.csv(file = "Data/External_data_genus_species/Rocha_Max_Min_Coral_Trait_Database.csv", header = T, dec = ".", sep = ";")

colnames (Max_Min) <- c("Species", "Min_Depth", "Max_Depth")
# to numeric
cols.num <- c("Min_Depth", "Max_Depth")
Max_Min[cols.num] <- sapply(Max_Min[cols.num],as.numeric)

str(Max_Min)


depth_diversity <- data.frame(
  Depth=seq(from = 5, to = 170, by = 5),
  Diversity=NA
)

# Measure generic richness for each depth, based on Max and Min Depths
for(i in unique(depth_diversity$Depth)){
  new <- filter(Max_Min, Min_Depth <= i)
  new2 <- filter(new, Max_Depth >= i)
  depth_diversity$Diversity[depth_diversity$Depth == i] <- nrow (new2)}

depth_diversity$Diversity <- as.numeric (depth_diversity$Diversity)

depth_diversity$Depth <- as.numeric (depth_diversity$Depth)

# Plot of species diversity
species_plot <- ggplot(depth_diversity, aes(x=Depth, y=Diversity)) +
  geom_point(shape = 19, color = "black", fill = "red", size = 2) + 
  theme_bw()  + ylab ("Species richness (n)") + xlab ("Depth (m)") + theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=8, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# Separate specie in genus also 
Max_Min_Depth_2 <-  Max_Min %>% separate(Species,  c("Coral_genus", "Coral_species"), " ", extra = "merge")

# Measure generic richness for each depth, based on Max and Min Depths
Max_Min_Depth_Genus <- data_frame(Coral_genus = NA, Max_Depth = NA, Min_Depth= NA)
for(i in unique(Max_Min_Depth_2$Coral_genus)){
  new <- filter(Max_Min_Depth_2, Coral_genus == i)
  new$Max_Depth <- max(new$Max_Depth)
  new$Min_Depth <- min(new$Min_Depth)
  new <- new[c("Coral_genus","Max_Depth","Min_Depth")]
  new <- new[!duplicated(new), ]
  Max_Min_Depth_Genus <- rbind (Max_Min_Depth_Genus,new)}


Genus_Diversity <- data.frame(
  Depth=seq(from = 5, to = 170, by = 5),
  Diversity=NA)

# Measure generic richness for each depth, based on Max and Min Depths
for(i in unique(Genus_Diversity$Depth)){
  new <- filter(Max_Min_Depth_Genus, Min_Depth <= i)
  new2 <- filter(new, Max_Depth >= i)
  Genus_Diversity$Diversity[Genus_Diversity$Depth == i] <- nrow (new2)}


# Plot of genus diversity
genus_plot <- ggplot(Genus_Diversity, aes(x=Depth, y=Diversity)) +
  geom_point(shape = 19, color = "black", fill = "red", size = 2) +
  theme_bw()  + ylab ("Genus richness (n)") + xlab ("Depth (m)") + theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=8, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# Correlate species vs genus 
cor.test (depth_diversity$Diversity, Genus_Diversity$Diversity)
cor.test(depth_diversity$Diversity, Genus_Diversity$Diversity, method ="spearman")

# Create a database together for plots 

depth_diversity$Genus <- Genus_Diversity$Diversity

linm <- lm(Diversity ~ Genus, data = depth_diversity)
linm_resume <- summary(linm)
a1 <- linm_resume$coef[2]
b1 <- linm_resume$coef[1]
plot (Diversity ~ Genus, data = depth_diversity, ylab = "Species richness", xlab = "Generic richness")
mtext(paste0("R squared: ",round(linm_resume$r.squared,5)),adj = 1)
# mtext(paste0("P-value: ", format.pval(pf(linm_resume$fstatistic[1], # F-statistic
#                                          linm_resume$fstatistic[2], # df
#                                          linm_resume$fstatistic[3], # df
#                                          lower.tail = FALSE))))
# mtext(paste0("y = ",round(linm_resume$coefficients[1],2)," + ", 
#              round(linm_resume$coefficients[2],2),"x"),adj = 1)
# abline(lm(depth_diversity$Diversity ~ depth_diversity$Genus), col = "black")


Correlation_Rocha <- ggplot(depth_diversity, aes(x=Genus, y=Diversity)) + 
  geom_point(shape = 21, size = 2) +
  # geom_abline(intercept = b1, slope = a1, colour ="black") +
  scale_x_continuous(name ="Genus richness (n)") +
  scale_y_continuous(name ="Species richness (n)") +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour= "black"), legend.position = "none")


(Figure_Rocha = species_plot + genus_plot + Correlation_Rocha + plot_layout(guides = 'collect'))

ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Figure_Rocha.pdf", Figure_Rocha, width =12, height= 4)



### From Roberts et al 2019 database Testing Biodiversity Theory
rm(list = ls ())

Max_Min_Depth <- read.csv(file = "Data/External_data_genus_species/Roberts et al 2019_Data_Testing Biodiversity Theory.csv", header = T, dec = ".", sep = ";")

#Separate max and min depth for each genera
Max_Min_Depth <- ddply(Max_Min_Depth, .(ID), summarize, Max_Depth=max(Depth), Min_Depth=min(Depth))

Max_Min_Depth <- Max_Min_Depth [c("ID", "Min_Depth","Max_Depth")]

colnames (Max_Min_Depth) <- c("Coral_species","Min_Depth", "Max_Depth")

# to numeric
cols.num <- c("Min_Depth", "Max_Depth")
Max_Min_Depth[cols.num] <- sapply(Max_Min_Depth[cols.num],as.numeric)

str(Max_Min_Depth)

# There a lot of duplicated species from the different sites, etc
# Remove the duplicated species and keep only the wider upper and lower boundaries
Max_Min_Depth_Species <- data_frame(Coral_species = NA, Max_Depth = NA, Min_Depth= NA)
for(i in unique(Max_Min_Depth$Coral_species)){
  new <- filter(Max_Min_Depth, Coral_species == i)
  new$Max_Depth <- max(new$Max_Depth)
  new$Min_Depth <- min(new$Min_Depth)
  new <- new[c("Coral_species","Max_Depth","Min_Depth")]
  new <- new[!duplicated(new), ]
  Max_Min_Depth_Species <- rbind (Max_Min_Depth_Species,new)}


depth_diversity <- data.frame(
  Depth=seq(from = 1, to = 45, by = 1),
  Diversity=NA)


# Measure generic richness for each depth, based on Max and Min Depths
for(i in unique(depth_diversity$Depth)){
  new <- filter(Max_Min_Depth_Species, Min_Depth <= i)
  new2 <- filter(new, Max_Depth >= i)
  depth_diversity$Diversity[depth_diversity$Depth == i] <- nrow (new2)}

depth_diversity$Diversity <- as.numeric (depth_diversity$Diversity)
depth_diversity$Depth <- as.numeric (depth_diversity$Depth)

# Plot of species diversity
species_plot <- ggplot(depth_diversity, aes(x=Depth, y=Diversity)) +
  geom_point(shape = 19, color = "black", fill = "red", size = 2) +
  theme_bw()  + ylab ("Species richness (n)") + xlab ("Depth (m)") + theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=8, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# Separate species in genus also 
Max_Min_Depth_2 <-  Max_Min_Depth %>% separate(Coral_species,  c("Coral_genus", "Coral_species"), " ", extra = "merge")

# Measure generic richness for each depth, based on Max and Min Depths
Max_Min_Depth_Genus <- data_frame(Coral_genus = NA, Max_Depth = NA, Min_Depth= NA)
for(i in unique(Max_Min_Depth_2$Coral_genus)){
  new <- filter(Max_Min_Depth_2, Coral_genus == i)
  new$Max_Depth <- max(new$Max_Depth)
  new$Min_Depth <- min(new$Min_Depth)
  new <- new[c("Coral_genus","Max_Depth","Min_Depth")]
  new <- new[!duplicated(new), ]
  Max_Min_Depth_Genus <- rbind (Max_Min_Depth_Genus,new)}


Genus_Diversity <- data.frame(
  Depth=seq(from = 1, to = 45, by = 1),
  Diversity=NA)

# Measure generic richness for each depth, based on Max and Min Depths
for(i in unique(Genus_Diversity$Depth)){
  new <- filter(Max_Min_Depth_Genus, Min_Depth <= i)
  new2 <- filter(new, Max_Depth >= i)
  Genus_Diversity$Diversity[Genus_Diversity$Depth == i] <- nrow (new2)}

# Plot of genus diversity
genus_plot <- ggplot(Genus_Diversity, aes(x=Depth, y=Diversity)) +
  geom_point(shape = 19, color = "black", fill = "red", size = 2) +
  theme_bw()  + ylab ("Genus richness (n)") + xlab ("Depth (m)") + theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=8, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 



# Correlate species vs genus 
cor.test (depth_diversity$Diversity, Genus_Diversity$Diversity)
cor.test(depth_diversity$Diversity, Genus_Diversity$Diversity, method ="spearman")

# Create a database together for plots 

depth_diversity$Genus <- Genus_Diversity$Diversity

linm <- lm(Diversity ~ Genus, data = depth_diversity)
linm_resume <- summary(linm)
a1 <- linm_resume$coef[2]
b1 <- linm_resume$coef[1]
plot (Diversity ~ Genus, data = depth_diversity, ylab = "Species richness", xlab = "Generic richness")
mtext(paste0("R squared: ",round(linm_resume$r.squared,5)),adj = 1)
# mtext(paste0("P-value: ", format.pval(pf(linm_resume$fstatistic[1], # F-statistic
#                                          linm_resume$fstatistic[2], # df
#                                          linm_resume$fstatistic[3], # df
#                                          lower.tail = FALSE))))
# mtext(paste0("y = ",round(linm_resume$coefficients[1],2)," + ", 
#              round(linm_resume$coefficients[2],2),"x"),adj = 1)
# abline(lm(depth_diversity$Diversity ~ depth_diversity$Genus), col = "black")


Correlation_Roberts_TestingTheory <- ggplot(depth_diversity, aes(x=Genus,y=Diversity)) +
  geom_point(shape = 21, size = 2) +
  # geom_abline(intercept = b1, slope = a1, colour ="black") +
  scale_x_continuous(name ="Genus richness (n)") +
  scale_y_continuous(name ="Species richness (n)") +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none")

(Figure_Roberts_TestingTheory = species_plot + genus_plot + Correlation_Roberts_TestingTheory + plot_layout(guides = 'collect'))

ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Figure_Roberts_TestingTheory.pdf",Figure_Roberts_TestingTheory, width = 12, height = 4)


### From Roberts et al 2019 Paradox diversity - From DRRHData 


rm(list = ls ())

Max_Min_Depth <- read.csv(file = "Data/External_data_genus_species/Roberts et al 2019_ParadoxDRRHData.csv", header = T, dec = ".", sep = ";")
# length (unique (Max_Min_Depth$Species))

Max_Min_Depth <- Max_Min_Depth [c("Species", "Min.Depth", "Max.Depth")]

colnames (Max_Min_Depth) <- c("Coral_species", "Min_Depth", "Max_Depth")
# to numeric
cols.num <- c("Min_Depth", "Max_Depth")
Max_Min_Depth[cols.num] <- sapply(Max_Min_Depth[cols.num],as.numeric)

str(Max_Min_Depth)

# In this case, there are not duplicated species from the different sites, etc

depth_diversity <- data.frame(
  Depth=seq(from = 1, to = 45, by = 1),
  Diversity=NA)


# Measure generic richness for each depth, based on Max and Min Depths
for(i in unique(depth_diversity$Depth)){
  new <- filter(Max_Min_Depth, Min_Depth <= i)
  new2 <- filter(new, Max_Depth >= i)
  depth_diversity$Diversity[depth_diversity$Depth == i] <- nrow (new2)}

depth_diversity$Diversity <- as.numeric (depth_diversity$Diversity)
depth_diversity$Depth <- as.numeric (depth_diversity$Depth)


# Plot of species diversity
species_plot <- ggplot(depth_diversity, aes(x=Depth, y=Diversity)) +
  geom_point(shape = 19, color = "black", fill = "red", size = 2) +
  theme_bw()  + ylab ("Species richness (n)") + xlab ("Depth (m)") + theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=8, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# Separate specie in genus also 
Max_Min_Depth_2 <-  Max_Min_Depth %>% separate(Coral_species,  c("Coral_genus", "Coral_species"), " ", extra = "merge")

# Measure generic richness for each depth, based on Max and Min Depths
Max_Min_Depth_Genus <- data_frame(Coral_genus = NA, Max_Depth = NA, Min_Depth= NA)
for(i in unique(Max_Min_Depth_2$Coral_genus)){
  new <- filter(Max_Min_Depth_2, Coral_genus == i)
  new$Max_Depth <- max(new$Max_Depth)
  new$Min_Depth <- min(new$Min_Depth)
  new <- new[c("Coral_genus","Max_Depth","Min_Depth")]
  new <- new[!duplicated(new), ]
  Max_Min_Depth_Genus <- rbind (Max_Min_Depth_Genus,new)}

# Measure generic richness for each depth with another methodology. To make sure our computation was correct
Genus_Diversity <- data.frame(
  Depth=seq(from = 1, to = 45, by = 1),
  Diversity=NA)

for(i in unique(Genus_Diversity$Depth)){
  new <- filter(Max_Min_Depth_Genus, Min_Depth <= i)
  new2 <- filter(new, Max_Depth >= i)
  Genus_Diversity$Diversity[Genus_Diversity$Depth == i] <- nrow (new2)}


# Plot of genus diversity
genus_plot <- ggplot(Genus_Diversity, aes(x=Depth, y=Diversity)) +
  geom_point(shape = 19, color = "black", fill = "red", size = 2) +
  theme_bw()  + ylab ("Genus richness (n)") + xlab ("Depth (m)") + theme_classic() +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=8, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

#
# Correlate species vs genus 
cor.test (depth_diversity$Diversity, Genus_Diversity$Diversity)
cor.test(depth_diversity$Diversity, Genus_Diversity$Diversity, method ="spearman")

# Create a database together for plots 

depth_diversity$Genus <- Genus_Diversity$Diversity


linm <- lm(Diversity ~ Genus, data = depth_diversity)
linm_resume <- summary(linm)
a1 <- linm_resume$coef[2]
b1 <- linm_resume$coef[1]

plot (Diversity ~ Genus, data = depth_diversity, ylab = "Species richness", xlab = "Generic richness")
mtext(paste0("R squared: ",round(linm_resume$r.squared,5)),adj = 0)
# mtext(paste0("P-value: ", format.pval(pf(linm_resume$fstatistic[1], # F-statistic
#                                          linm_resume$fstatistic[2], # df
#                                          linm_resume$fstatistic[3], # df
#                                          lower.tail = FALSE))))
# mtext(paste0("y = ",round(linm_resume$coefficients[1],2)," + ", 
#              round(linm_resume$coefficients[2],2),"x"),adj = 1)
# abline(lm(depth_diversity$Diversity ~ depth_diversity$Genus), col = "black")

Correlation_Roberts_Paradox <- ggplot(depth_diversity, aes(x=Genus, y=Diversity)) +
  geom_point(shape = 21, size = 2) +
  # geom_abline(intercept = b1, slope = a1, colour ="black") +
  scale_x_continuous(name ="Genus richness (n)") +
  scale_y_continuous(name ="Species richness (n)") +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none")

(Figure_Roberts_Paradox = species_plot + genus_plot + Correlation_Roberts_Paradox + plot_layout(guides = 'collect'))


ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Figure_Roberts_Paradox.pdf",Figure_Roberts_Paradox,width = 12, height = 4)



