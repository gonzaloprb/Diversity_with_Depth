# rm (list = ls())

#   title: "Script_Diversity_Cover"   

# Script also pushed on Github: https://github.com/gonzaloprb/PhD_Diversity_Depth

 
# Install.packages
library(scatterplot3d); library(reshape); library (reshape2); library (data.table); library(RColorBrewer); library(dplyr); library (tidyr); library (plyr); library (ggplot2); require (vegan); require (goeveg);   require  (stats); require (dynRB); require (matrixStats); require(plyr); require (stringr);require (reshape2); require(gridExtra)
require (betapart); library (car); library (MASS); library (glmm); library (glmnet); library (glmmTMB); library (lme4); library(randomcoloR)
require (ggpubr); library (cowplot); library (patchwork); library (scales);library(viridisLite) ; library(patchwork); library(fishualize);library(tidyverse) 
require (DHARMa); require (MuMIn); require (merTools)

# Mid-domain effect (there's no mid-domain effect)
require(devtools)
library(reshape2)

# Geospatial data
require (geodist)

# rm (list = ls()) 
# Set working directory etc.
rndpts_df <- read.csv(file = "Data/photoquad_rndpts_DEEPHOPE.csv", header = T, dec = ".", sep = ",", row.names = 1)

rndpts_df$Island <- gsub("Mangareva", "Gambier", rndpts_df$Island)

# Keep only corals in rndpts_df / fix the genus/form issue, aggregate
rndpts_df <- filter (rndpts_df, Cattegory == "Scleractinian")
rndpts_df <- rndpts_df %>% separate(Coral_Genus_Form,  c("Coral_genus", "Coral_form"), " ", extra = "merge")
rndpts_df <- subset (rndpts_df, select = - c(Cattegory, Coral_form))
rndpts_df <- aggregate (Cover ~ Unique_Image_ID + Date + Site + Archipelago + Island + Island_Site + Depth + Quadrat + Coral_genus, rndpts_df , sum)


data <- read.csv(file = "Data/photoquad_rndpts_DEEPHOPE.csv", header = T, dec = ".", sep = ",", row.names = 1)

data$Island <- gsub("Mangareva", "Gambier", data$Island)
# We want to keep all quadrats
# Keep only corals in rndpts_df / fix the genus/form issue, aggregate
data <- subset (data, Cattegory == "Scleractinian")
data <- data %>% separate(Coral_Genus_Form,  c("Coral_genus", "Coral_form"), " ", extra = "merge")

data <- aggregate (Cover ~ Unique_Image_ID + Date + Site + Archipelago + Island + Island_Site + Depth + Quadrat + Cattegory +  Coral_genus, data , sum)

data_df <- data

data <- ddply(data, ~ Archipelago + Island +Island_Site + Depth + Cattegory + Coral_genus, function(x){c(Cover = sum(x$Cover)/30) })

coral_data <- data


# All stats below come from rndpts_df or pooled together quadrats coral_data
# I re-start every time for simplicity 

## Coral cover profile 
# Plot considering the effect of Island and Island_Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
coral_cover <- aggregate (Cover ~ Island + Island_Site + Depth, coral_data , sum)
# For keeping quadrats - perhaps unnecessary
data_df <- aggregate (Cover ~ Island + Island_Site + Quadrat +  Depth, data_df , sum)

# Transform depth as a Qualitative variable  
coral_cover$Depth <- as.factor(as.character(coral_cover$Depth))
coral_cover$Depth = factor(coral_cover$Depth,levels = c ("6", "20", "40", "60", "90", "120"))

# ggplot with island sites
ggplot(coral_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# no-mid-domain effect. 

#### Mixed-Model for showing cover decreases with depth in all sites
# Just very quickly linear model. I don't think I need to complicate myself with glm or lmer

coral_cover2 <- coral_cover
coral_cover2 <- coral_cover2 %>% complete(Island,Island_Site, Depth,fill = list(Cover = 0))
coral_cover2$Depth <- as.numeric (as.character(coral_cover2$Depth))
coral_cover2$Island_Site <- as.character(coral_cover2$Island_Site)
coral_cover2$Island_Site_2 <- paste(coral_cover2$Island, "_", coral_cover2$Island_Site)

#### Keeping quadrats - perhaps unnecessary ####
data_df2 <- data_df
data_df2$Depth <- as.numeric (as.character(data_df2$Depth))
data_df2$Island_Site <- as.character(data_df2$Island_Site)
data_df2$Island_Site_2 <- paste(data_df2$Island, "_", data_df2$Island_Site)
#### Keeping quadrats - perhaps unnecessary ####

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

mixed.lmer.full <- lmer(Cover ~ Depth +   (1|Island_Site_2), data = coral_cover2, REML = TRUE)
mixed.lmer <- lmer(Cover ~  (1|Island_Site_2),data = coral_cover2, REML = TRUE)


#### Keeping quadrats - it gives the same ####
mixed.lmer.quadrats.full <- lmer(Cover ~ Depth +   (1|Island_Site_2), data = data_df2, REML = TRUE)
mixed.lmer.quadrats <- lmer(Cover ~  (1|Island_Site_2),data = data_df2, REML = TRUE)
#### Keeping quadrats - it gives the same ####


# Testing nested structure! According to AIC better without the nested structure. 
AIC(lmer(Cover ~ (1|Island_Site_2/Island), data = coral_cover2, REML = TRUE)) # Nested effect. Higher AIC and give a warning! Better without nested structure
AIC(lmer(Cover ~ (1|Island_Site_2), data = coral_cover2, REML = TRUE))

# Getting the p-value for depth!
anova(mixed.lmer.full,mixed.lmer)

# Check normality and assumptions
shapiro.test(residuals(mixed.lmer.full))

leveneTest(residuals(mixed.lmer.full) ~ as.factor(coral_cover2$Depth)) # Levenetest does not accept quantitive data so I made Depth as factor.
# p-value is less than 0.05, which means we accept the null hypothesis that there is no variance among the groups
# Since the result is not significant, the assumption of equal variances (homoscedasticity) is met.

boxplot(residuals(mixed.lmer.full) ~ coral_cover2$Depth)


plot(mixed.lmer.full)


summary(mixed.lmer.full)
confint(mixed.lmer.full)
anova (mixed.lmer.full)
coef(mixed.lmer.full)





# No variance as Site is only considered as random intercept
plot(mixed.lmer.full) # Data is okay! The deviaition of the residuals stays kind of constant
qqnorm(resid(mixed.lmer.full))
qqline(resid(mixed.lmer.full))
# Check singularity
# diag.vals <- getME(mixed.lmer.full,"theta")[getME(mixed.lmer.full,"lower") == 0]
# any(diag.vals < 1e-6) # FALSE

# With all quadrats
plot(mixed.lmer.quadrats.full) # Data is not okay! The deviaition of the residuals increases a lot!
qqnorm(resid(mixed.lmer.quadrats.full))

# Plots of model with DHARMa package
testDispersion(mixed.lmer.full)
simulationOutput <- simulateResiduals(fittedModel = mixed.lmer.full, plot = T)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plotQQunif(simulationOutput) # left plot in plot.DHARMa()
plotResiduals(simulationOutput) # right plot in plot.DHARMa()


# Plots of model with DHARMa package - but with all quadrats
testDispersion(mixed.lmer.quadrats.full)
simulationOutput <- simulateResiduals(fittedModel = mixed.lmer.quadrats.full, plot = T)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plotQQunif(simulationOutput) # left plot in plot.DHARMa()
plotResiduals(simulationOutput) # right plot in plot.DHARMa()


# Get the marginal and conditional R2
r.squaredGLMM(mixed.lmer.full)
r.squaredGLMM(mixed.lmer)
# Marginal is explained by the fixed effects. 
# Conditional is explained by the full model



# Extra mixed models
mixed.lmer_2 <- lmer(Cover ~ Depth + (Depth | Island_Site_2), coral_cover2)
summary(mixed.lmer_2)
r.squaredGLMM(mixed.lmer_2)
# 1.766e+01/(1.766e+01 +1.804e+02) --> Island_Site_2 explains 8% of variance
plot(mixed.lmer_2) 
qqnorm(resid(mixed.lmer_2))
qqline(resid(mixed.lmer_2))

mixed.lmer_3 <- lmer(Cover ~ Depth + (Depth || Island_Site_2), coral_cover2)
summary(mixed.lmer_3)
r.squaredGLMM(mixed.lmer_3)
plot(mixed.lmer_3) 
qqnorm(resid(mixed.lmer_3))
qqline(resid(mixed.lmer_3))

mixed.lmer_4 <- lmer(Cover ~ 1 + Depth + (1 + Depth | Island_Site_2), coral_cover2)
summary(mixed.lmer_4)
r.squaredGLMM(mixed.lmer_4)
# 1.766e+01/(1.766e+01 +1.804e+02) --> Intercept of Island_Site_2 explains 8% of variance
# 1.726e-03 /(1.726e-03  +1.804e+02) --> Intercept of Island_Site_2 explains 8% of variance
plot(mixed.lmer_4) 
qqnorm(resid(mixed.lmer_4))
qqline(resid(mixed.lmer_4))

mixed.lmer_5 <- lmer(Cover ~ 1 + Depth + (1 + Depth | Island_Site_2) +  (1|Island_Site_2) , coral_cover2)
summary(mixed.lmer_5)
r.squaredGLMM(mixed.lmer_5)
plot(mixed.lmer_5) 
qqnorm(resid(mixed.lmer_5))
qqline(resid(mixed.lmer_5))

mixed.lmer_nested <- lmer(Cover ~  (1|Island_Site_2/Island), data = coral_cover2, REML = TRUE)
mixed.lmer_nested_full <- lmer(Cover ~ Depth +  (1|Island_Site_2/Island), data = coral_cover2, REML = TRUE)

r.squaredGLMM(mixed.lmer_nested)
r.squaredGLMM(mixed.lmer_nested_full)

anova(mixed.lmer_nested_full,mixed.lmer_nested)

# mixed.log model
mixed.lmer.log <- lmer(Cover ~ log(Depth) + (1|Island_Site_2), data = coral_cover2) 



anova (mixed.lmer, mixed.lmer.full,mixed.lmer_2,mixed.lmer_3,mixed.lmer_4,mixed.lmer_5, mixed.lmer_nested, mixed.lmer_nested_full, mixed.lmer.log)








#################################### new #############################



### Use a beta distribution to avoid problem of normality 

coral_cover2$Cover_Beta <- coral_cover2$Cover/100


coral_cover2$Tot_Points <- 75
coral_cover2$Coral_points <- (coral_cover2$Cover * coral_cover2$Tot_Points) / 100
coral_cover2$Coral_points  <- round(coral_cover2$Coral_points,0)
coral_cover2$NonCoral_points <- coral_cover2$Tot_Points - coral_cover2$Coral_points 

coral_cover2$Proportion = coral_cover2$Coral_points / (coral_cover2$Coral_points + coral_cover2$NonCoral_points)
coral_cover2$Proportion [coral_cover2$Proportion == 0] <-  .001 # Otherwise beta distribution with 0 it does not work

mixed.lmer.full_beta <- lmer(Proportion ~ Depth +   (1|Island_Site_2), data = coral_cover2, REML = TRUE)
r.squaredGLMM(mixed.lmer.full_beta)

# Check normality and assumptions
shapiro.test(residuals(mixed.lmer.full_beta)) # Smaller than 0.05, accept null hypothesis, reject normality

leveneTest(residuals(mixed.lmer.full_beta) ~ as.factor(coral_cover2$Depth)) # Levenetest does not accept quantitive data so I made Depth as factor.
# p-value is smaller than 0.05, which means we accept the null hypothesis that there is  variance among the groups
# Since the result is not significant, the assumption of equal variances (homoscedasticity) is met.

boxplot(residuals(mixed.lmer.full_beta) ~ coral_cover2$Depth)

plot(mixed.lmer.full_beta)


# With beta distribution and glmmTMB package
library (glmmTMB)

mixed_beta_full <- glmmTMB (Proportion ~ Depth + (1|Island_Site_2), data = coral_cover2, REML = TRUE, family = beta_family())
summary (mixed_beta_full)

mixed_beta_null <- glmmTMB (Proportion ~  (1|Island_Site_2), data = coral_cover2, REML = TRUE, family = beta_family())
summary (mixed_beta_null)

mixed_beta_2 <- glmmTMB (Proportion ~ Depth + (Depth | Island_Site_2) , data = coral_cover2, REML = TRUE, family = beta_family())
# No convergence      summary (mixed_beta_2)

mixed_beta_3 <- glmmTMB (Proportion ~ 1 + Depth + (1 + Depth | Island_Site_2), data = coral_cover2, REML = TRUE, family = beta_family())
# No convergence      summary (mixed_beta_3)

mixed_beta_4 <- glmmTMB (Proportion ~ 1 + Depth + (1 + Depth | Island_Site_2) +  (1|Island_Site_2), data = coral_cover2, REML = TRUE, family = beta_family())
# No convergence      summary (mixed_beta_4)

mixed.lmer_beta_nested <- glmmTMB(Proportion ~  (1|Island_Site_2/Island), data = coral_cover2, REML = TRUE, family = beta_family())
summary (mixed.lmer_beta_nested)
mixed.lmer_beta_nested_full <- glmmTMB(Proportion ~ Depth +  (1|Island_Site_2/Island), data = coral_cover2, REML = TRUE, family = beta_family())
summary (mixed.lmer_beta_nested_full)


anova(mixed_beta_null,mixed.lmer_beta_nested) # Better without nested effect

anova(mixed_beta_full,mixed_beta_2,mixed_beta_3,mixed_beta_4)
# Most did not converge beta 2, 3, 4


shapiro.test(residuals(mixed_beta_full)) # Normality rejected again

leveneTest(residuals(mixed_beta_full) ~ as.factor(coral_cover2$Depth)) # Levenetest does not accept quantitive data so I made Depth as factor.
# normality rejected again 

# Plots of model with DHARMa package 
testDispersion(mixed_beta_full)
simulationOutput <- simulateResiduals(fittedModel = mixed_beta_full, plot = T) # It looks like okay but we still have problem of normality

# The only good one, mixed_beta_full but we have the problem of normality
summary (mixed_beta_full)




# With another package
library (GLMMadaptive)

mixed_model_beta_full <- mixed_model(fixed = Proportion ~  Depth, random = ~ 1 | Island_Site_2, data = coral_cover2,
            family = beta.fam())

mixed_model_full_2 <- mixed_model(fixed = Proportion ~ 1 + Depth , random = ~ 1 + Depth | Island_Site_2, data = coral_cover2,
                               family = beta.fam())

shapiro.test(residuals(mixed_model_beta_full)) # does not accept normality

leveneTest(residuals(mixed_model_beta_full) ~ as.factor(coral_cover2$Depth)) # Levenetest does not accept quantitive data so I made Depth as factor.

testDispersion(mixed_model_beta_full)
simulationOutput <- simulateResiduals(fittedModel = mixed_model_beta_full, plot = T) # Significant deviations within groups


summary (mixed_model_beta_full)

anova (mixed_model_beta_full,mixed_model_full_2)


# Try with Binomial distribution

glmer_binomial_model <- glmer(cbind(Coral_points, Tot_Points - Coral_points) ~ Depth + (1 | Island_Site_2),
      data = coral_cover2, family = binomial)
summary (glmer_binomial_model)

glmer_binomial_model_null <- glmer(cbind(Coral_points, Tot_Points - Coral_points) ~  (1 | Island_Site_2),
                              data = coral_cover2, family = binomial)
summary (glmer_binomial_model_null)

r.squaredGLMM(glmer_binomial_model) # The binomial has the highest R2
anova (glmer_binomial_model,glmer_binomial_model_null, mixed.lmer.full) # but the AIC is much lower 


# The best solution is to use Bayesian modelling.

library(brms);library('rstan'); library("stam");library("parallel"); library ("performance")
Binomial_Cover_model <- brm(Coral_points | trials(Tot_Points) ~  Depth + (1 | Island_Site_2),
                                  data = coral_cover2, family = binomial(),  # prior = my_priors,
                                  control = list(adapt_delta = 0.9, max_treedepth = 11),
                                  iter = 4000, warmup = 1000, chains = 2, cores = 2) 
save(Binomial_Cover_model, file="Data/Binomial_Cover_model_Diversity.RData")
load("Data/Binomial_Cover_model_Diversity.RData") 


Beta_Cover_model <- brm(Cover_Beta ~  Depth + (1 | Island_Site_2),
                            data = coral_cover2, family = Beta(),  # prior = my_priors,
                            control = list(adapt_delta = 0.9, max_treedepth = 11),
                            iter = 4000, warmup = 1000, chains = 2, cores = 2) 
save(Beta_Cover_model, file="Data/Beta_Cover_model_Diversity.RData")
load("Data/Beta_Cover_model_Diversity.RData") 

Beta_Cover_model_zero <- brm(Cover_Beta ~  Depth + (1 | Island_Site_2),
                        data = coral_cover2, family = zero_inflated_beta(),  # prior = my_priors,
                        control = list(adapt_delta = 0.9, max_treedepth = 11),
                        iter = 4000, warmup = 1000, chains = 2, cores = 2) 
save(Beta_Cover_model_zero, file="Data/Beta_Cover_model_zero_Diversity.RData")
load("Data/Beta_Cover_model_zero_Diversity.RData") 


Student_model <- brm(Cover ~  Depth + (1 | Island_Site_2),
                             data = coral_cover2, family = student(),  # prior = my_priors,
                             control = list(adapt_delta = 0.9, max_treedepth = 11),
                             iter = 2000, warmup = 1000, chains = 2, cores = 2) 
save(Student_model, file="Data/Student_model.RData")
load("Data/Student_model.RData") 

# The binomial is better than the beta model 
summary (Beta_Cover_model)
summary (Binomial_Cover_model)
summary (Beta_Cover_model_zero)
summary (Student_model)

bayes_R2(Beta_Cover_model)
bayes_R2(Binomial_Cover_model)
bayes_R2(Beta_Cover_model_zero)
bayes_R2(Student_model)

# model_weights(Binomial_Cover_model, Beta_Cover_model, Beta_Cover_model_zero, weights = "waic")

# The best is to use the Binomial_Cover_Model
plot(Binomial_Cover_model)
pp_check(Binomial_Cover_model, type = "scatter_avg") # Not structured data
bayes_R2(Binomial_Cover_model) # 
r2_bayes(Binomial_Cover_model)

me_null <- conditional_effects(Binomial_Cover_model, nsamples = 1000, probs = c(0.05, 0.95), spaghetti = F) # Default is 0.95
plot(me_null, ask = FALSE, points = F) # Probability scale!


# Work with posterior predict 
# Create the reference dataframe - newdata (Here for all depths)
Depth <- unique (Binomial_Cover_model$data$Depth)
Island_Site_2 <- unique (Binomial_Cover_model$data$Island_Site_2)
Coral_points <- unique (Binomial_Cover_model$data$Coral_points)
Tot_Points <- unique (Binomial_Cover_model$data$Tot_Points)
#
ref_data <- crossing(Depth, Island_Site_2, Tot_Points)

# Forest plots with the posteriors
fitted_values <- posterior_epred(Binomial_Cover_model, newdata = ref_data, re_formula = 'Coral_points | trials(Tot_Points) ~  Depth + (1 | Island_Site_2)')

# Number of rows equals to ref_data and number of dimensions is equal to (Island_Island_Site*Depths)
# str (fitted_values)
# dim (fitted_values)

# Necessary to traspose
fitted_values <- t(fitted_values)

# Create combination of ref_data
ref_data_fitted <- cbind (ref_data [c(1,2)],fitted_values)

# Multiple columns of predictions into a single one. 
ref_data_fitted <- melt (ref_data_fitted, id.vars = c ("Depth", "Island_Site_2"), na.rm = F, measure.vars = c(3:6002), value.name = c("Posterior_Prob"))
ref_data_fitted$Depth = factor (ref_data_fitted$Depth, levels = c ("120", "90","60", "40","20", "6"))


# Summary depth and island among all iterations 

summary <- ddply(ref_data_fitted, .(Depth,Island_Site_2), summarize, Post_Mean=mean(Posterior_Prob), Post_Sd = sd(Posterior_Prob), Post_se=sd(Posterior_Prob) / sqrt(length(Posterior_Prob)), 
                 Post_Margin.error = qt(p=0.05/2, df=length (Posterior_Prob)-1,lower.tail=F) * Post_se)

coral_cover2 <- merge (coral_cover2, summary)



# Extract the number of points from cover / Same as done before
coral_cover2$Tot_Points <- 75
coral_cover2$Coral_points <- (coral_cover2$Cover * coral_cover2$Tot_Points) / 100
# Round the points
coral_cover2$Coral_points <-round(coral_cover2$Coral_points,0)
coral_cover2$NonCoral_points <-abs (coral_cover2$Coral_points - coral_cover2$Tot_Points)


# Transform the posterior predict of binomial from points proportions (out of 75) to coral cover (%)
coral_cover2$Post_Mean_Cover <- (coral_cover2$Post_Mean * 100) / coral_cover2$Tot_Points
coral_cover2$Post_Sd_Cover <- (coral_cover2$Post_Sd * 100) / coral_cover2$Tot_Points
coral_cover2$Post_se_Cover <- (coral_cover2$Post_se * 100) / coral_cover2$Tot_Points
coral_cover2$Post_Margin.error_Cover <- (coral_cover2$Post_Margin.error * 100) / coral_cover2$Tot_Points


#################################### new #############################


# Measure mean and standard error and confidence intervals from raw cover values

coral_cover <- coral_cover2
summary <- ddply(coral_cover, .(Depth), summarize, Mean=mean(Cover), Cover_se=sd(Cover) / sqrt(length(Cover)), 
                 Margin.error = qt(p=0.05/2, df=length (Cover)-1,lower.tail=F) * Cover_se)

# Confidence Interval manually would be
# alpha = 0.05
# degrees.freedom = 16 - 1
# t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
# print(t.score)
# margin.error <- t.score * summary$Cover_se
# # Confidence Interval

# Combine dataframes
coral_cover <- merge (coral_cover,summary)

# Make a single value for each depth of Post Bayesian
summary2 <- ddply(coral_cover, .(Depth), summarize, Post_Cover_Depth=mean(Post_Mean_Cover), Post_Cover_sd_Depth = sd (Post_Mean_Cover))
coral_cover <- merge (coral_cover,summary2)


# Make the plot again
coral_cover$Depth <- as.numeric (as.character(coral_cover$Depth))




# Plot with Bayesian Binomial model (with standard deviation of the model and with coral cover mean values with CI

Fig_1A <- ggplot(coral_cover, aes(x=Depth, y=Cover)) +
  # geom_point(aes(fill = Island, colour = Island),shape = 21, size = 0.8)  + 
  geom_errorbar(aes(ymin = Mean - Margin.error, ymax = Mean + Margin.error), color="black", size=1, width=2) +
  geom_line(aes(y=Post_Cover_Depth), size=1, linetype="dashed", alpha = 1) +
  geom_line(aes(y=Post_Cover_Depth + Post_Cover_sd_Depth), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_line(aes(y=Post_Cover_Depth - Post_Cover_sd_Depth), size=0.3, linetype="dotted",alpha = 0.8) +
  scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
  geom_point(aes(y=Mean), shape=21, fill="white", size=4) +  
  scale_x_continuous(name ="Depth (m)", limits=c(5,122), breaks = c(6,20,40,60,90,120)) +
  scale_y_continuous(name ="Coral cover (%)", limits=c(-5,85), breaks = c(0,20,40,60,80)) +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none") 
Fig_1A
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_1A.pdf", Fig_1A,width = 4, height = 3.5)




# Measure the mean of cover for each island (mean of site 1 and 2)

Island_cover <- aggregate (Cover ~ Island  + Depth, coral_cover , mean)

colnames (Island_cover) <- c("Island", "Depth","Island_cover")

coral_cover <- merge (coral_cover,Island_cover)

Island_sd <- aggregate (Cover ~ Island  + Depth, coral_cover , sd)
colnames (Island_sd) <- c("Island", "Depth","Island_sd")

coral_cover <- merge (coral_cover,Island_sd)

# Measure mean and standard error
summary_island <- ddply(coral_cover, ~ Island + Depth, summarize, Island_cover=mean(Cover), Island_se=sd(Cover) / sqrt(length(Cover)))

coral_cover <- merge (coral_cover,summary_island)

# Separate per islands if finally keep a supplementary figure 

# Byesian has not interval confidences, it has sd of the model isntead

Fig_SXX <- ggplot(coral_cover, aes(x=Depth, y=Cover)) + 
  facet_wrap(~Island, ncol = 4, scales = "free")  +
  geom_line(aes(y=Post_Cover_Depth), size=1, linetype="dashed", alpha = 1) +
  geom_line(aes(y=Post_Cover_Depth + Post_Cover_sd_Depth), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_line(aes(y=Post_Cover_Depth - Post_Cover_sd_Depth), size=0.3, linetype="dotted",alpha = 0.8) +
  # geom_errorbar(aes(ymin = Island_cover - Island_se, ymax = Island_cover + Island_se), color="black", size=0.5, width=0.5) +
  geom_point(aes(y=Island_cover), shape=21, size=1.5, fill = "grey")+
  geom_point(aes(y=Cover, shape = Island_Site), size=0.5, fill = "grey")+ # Variability of sites
  scale_x_continuous(name ="Depth (m)", limits=c(3,122), breaks = c(6,20,40,60,90,120)) +
  scale_y_continuous(name ="Coral cover (%)", limits=c(-5,85), breaks = c(0,20,40,60,80)) +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"), 
                      strip.text = element_text(size = 11,face="bold", colour="black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(), legend.position = "none") 
Fig_SXX

ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_SCover_Island.pdf", Fig_SXX, width = 6, height = 5)





# Stack Area plot
# Relative contribution of each genera to visually see if it changes or not

coral_data_main_Genus <- aggregate(Cover ~  Island + Island_Site + Depth + Coral_genus,  coral_data, sum)
coral_data_main_Genus <- coral_data_main_Genus %>% complete(Island,Island_Site, Depth,Coral_genus,fill = list(Cover = 0))

coral_data_main_Total <- aggregate(Cover ~  Depth +  Island + Island_Site,  coral_data, sum)
coral_data_main <- merge (coral_data_main_Genus,coral_data_main_Total, by=c("Island", "Island_Site","Depth"))
names (coral_data_main) <- c("Island", "Island_Site","Depth", "Coral_genus","Genus_cover", "Total_cover")
# Calculation of Relative contribution to Total coral cover
coral_data_main$Relative <- (coral_data_main$Genus_cover/ coral_data_main$Total_cover)*100
rm (coral_data_main_Genus, coral_data_main_Total)

coral_data_main <- coral_data_main %>% unite(Island_Island_Site, Island, Island_Site)


# Stack area plot
# colours <- rainbow(33)
colours <- distinctColorPalette(33)
# brewer.pal(33, "BrBG")
Sup_Fig_X <- ggplot(coral_data_main, aes(x=Depth, y=Relative, fill=Coral_genus)) +
  scale_fill_manual(values= colours) +
  geom_area(alpha=1 , size=0.5, colour="black") + facet_wrap(~Island_Island_Site,ncol = 8) + 
  scale_x_reverse(lim=c(120,0), breaks = c(6,20,40,60,90,120)) +    coord_flip() + 
  ggtitle ("Relative percentatge of generic composition") + xlab ("Depth (m)") + ylab ("") + labs(fill='Coral genera', size = 17) +
  theme_classic() + theme(axis.title = element_text(size=17),
                          plot.title = element_text(size=20),
                          axis.text.y =  element_text(size=15),
                          axis.text.x =  element_blank(),
                          axis.ticks.x=element_blank(),
                          strip.text.x =  element_text(size = 15,colour = "black", angle = 0),
                          strip.text.y = element_text(size = 6),strip.text.y.left = element_text(angle=0, size = 15),
                          strip.background = element_rect(fill="grey"),
                          legend.position = "bottom") +
                          guides(fill = guide_legend(ncol = 8, title.hjust = 0.4))
Sup_Fig_X



# Possible Fig 1C 
# Stack Area plot - Composition!
# Relative contribution of each genera to visually see if it changes or not

coral_data_main_Genus <- aggregate(Cover ~  Island + Island_Site + Depth + Coral_genus,  coral_data, sum)
coral_data_main_Genus <- coral_data_main_Genus %>% complete(Island,Island_Site, Depth,Coral_genus,fill = list(Cover = 0))

# Group the very rare corals together! Otherwise, impossible to see the plot. At least display:
unique (coral_data_main_Genus$Coral_genus)

keep_species <-  c ("Acropora","Astreopora","Montipora","Leptoseris","Leptastrea","Pachyseris","Pavona","Echynophillia","Pocillopora","Porites")

# Transform rest of genus to others
coral_data_main_Genus$Coral_genus <- gsub("Dipsastrea|Astrea|Leptoria|Plesiastrea|Cyphastrea|Paragoniastrea|Cycloseris|Cantharellus|Pleuractis|Lobactis|Herpolitha|Acanthastrea|Alveopora|Coscinarea|Gardinoseris|Lobophyllia|Napopora|Non_Id_Coral|Psammocora|Stylocoeniella|Turbinaria|Sandalolitha|Fungia", "Others", coral_data_main_Genus$Coral_genus)

######### This is measuring relative composition for each site and later making the mean of the relatives. This is good!
coral_data_compo <- aggregate(Cover ~   Island + Island_Site + Depth + Coral_genus,  coral_data_main_Genus, mean)
coral_data_compo_total <- aggregate(Cover ~  Island + Island_Site + Depth,  coral_data_compo, sum)

coral_data_compo <- merge (coral_data_compo,coral_data_compo_total, by=c("Island","Island_Site","Depth"))
names (coral_data_compo) <- c("Island","Island_Site","Depth", "Coral_genus","Genus_cover", "Total_cover")
# Calculation of Relative contribution to Total coral cover
coral_data_compo$Relative <- (coral_data_compo$Genus_cover/ coral_data_compo$Total_cover)*100

coral_data_compo <- aggregate(Relative ~  Depth + Coral_genus,  coral_data_compo, mean)

rm (coral_data_main_Genus, coral_data_compo_total)

#########

### This is pooling all sites together and later measuring relative. This is not good!
# coral_data_compo <- aggregate(Cover ~  Depth + Coral_genus,  coral_data_main_Genus, mean)
# coral_data_compo_total <- aggregate(Cover ~  Depth,  coral_data_compo, sum)
# 
# coral_data_compo <- merge (coral_data_compo,coral_data_compo_total, by=c("Depth"))
# names (coral_data_compo) <- c("Depth", "Coral_genus","Genus_cover", "Total_cover")
# # Calculation of Relative contribution to Total coral cover
# coral_data_compo$Relative <- (coral_data_compo$Genus_cover/ coral_data_compo$Total_cover)*100
# rm (coral_data_main_Genus, coral_data_compo_total)
###


# As factor to plot in the order we want
coral_data_compo$Coral_genus = factor(coral_data_compo$Coral_genus,levels = c ("Pocillopora", "Pachyseris", "Leptoseris","Montipora","Porites", "Acropora", "Pavona",  "Echynophillia","Leptastrea","Astreopora","Others"))

# Set good colours for final figure 1C

# colours_corals <- distinctColorPalette(13)
# save(colours_corals, file = "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/colours_corals.Rdata")

# Load colours_corals
 load(file = "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/colours_corals.Rdata")


# Vertical bars plot 
Fig_1C <- ggplot(coral_data_compo, aes(x=Depth, y=Relative)) +
  geom_bar(aes(fill=factor(Coral_genus)), stat="identity") + 
  scale_fill_manual(values= colours_corals) +
  # geom_line(aes(y=fit_se), size=1, linetype="dotted",alpha = 0.8, col = "blue") +
  # geom_line(aes(y=fit_seneg), size=1, linetype="dotted",alpha = 0.8, col = "blue") +
  scale_x_continuous(name ="Depth (m)", breaks = c(6,20,40,60,90,120)) +
  scale_y_continuous(name ="Dominance relative composition (%)", limits=c(0,101), breaks = c(0,20,40,60,80,100)) +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "bottom") +
                          guides(fill = guide_legend(ncol = 4, title = NULL))
Fig_1C
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_1C.pdf", Fig_1C,width = 4, height = 4)




# Stack area plot

Sup_Fig_4_1 <- ggplot(coral_data_compo, aes(x=Depth, y=Relative, fill=Coral_genus)) +
  scale_fill_manual(values= colours_corals) +
  geom_area(alpha=1 , size=0.5, colour="black") + 
  scale_x_reverse(lim=c(120,0), breaks = c(6,20,40,60,90,120)) +    coord_flip() + 
  ggtitle ("Domminance community structure") + xlab ("Depth (m)") + ylab ("Coral genera") + 
  theme_classic() + theme(axis.title = element_text(size=17),
                          plot.title = element_text(size=20),
                          axis.text.y =  element_text(size=15),
                          axis.text.x =  element_blank(),
                          axis.ticks.x=element_blank(),
                          strip.background = element_rect(fill="grey"),
                          legend.text=element_text(size=4),
                          legend.title=element_blank(), 
                          legend.position = "bottom") 
Sup_Fig_4_1



# Plot for supplementary - Composition by the different sites
coral_data_main_Genus <- aggregate(Cover ~  Island + Island_Site + Depth + Coral_genus,  coral_data, sum)
coral_data_main_Genus <- coral_data_main_Genus %>% complete(Island,Island_Site, Depth,Coral_genus,fill = list(Cover = 0))

# Group the very rare corals together! Otherwise, impossible to see the plot. At least display:
unique (coral_data_main_Genus$Coral_genus)

keep_species <-  c ("Acropora","Astreopora","Montipora","Leptoseris","Leptastrea","Pachyseris","Pavona","Echynophillia","Pocillopora","Porites")

# Transform rest of genus to others
coral_data_main_Genus$Coral_genus <- gsub("Dipsastrea|Astrea|Leptoria|Plesiastrea|Cyphastrea|Paragoniastrea|Cycloseris|Cantharellus|Pleuractis|Lobactis|Herpolitha|Acanthastrea|Alveopora|Coscinarea|Gardinoseris|Lobophyllia|Napopora|Non_Id_Coral|Psammocora|Stylocoeniella|Turbinaria|Sandalolitha|Fungia", "Others", coral_data_main_Genus$Coral_genus)


coral_data_compo <- aggregate(Cover ~  Island + Island_Site + Depth + Coral_genus ,  coral_data_main_Genus, mean)
coral_data_compo_total <- aggregate(Cover ~ Island + Island_Site + Depth,  coral_data_compo, sum)

coral_data_compo <- merge (coral_data_compo,coral_data_compo_total, by=c("Island", "Island_Site","Depth"))
names (coral_data_compo) <- c("Island","Island_Site","Depth", "Coral_genus","Genus_cover", "Total_cover")
# Calculation of Relative contribution to Total coral cover
coral_data_compo$Relative <- (coral_data_compo$Genus_cover/ coral_data_compo$Total_cover)*100
rm (coral_data_main_Genus, coral_data_compo_total)

# As factor to plot in the order we want
coral_data_compo$Coral_genus = factor(coral_data_compo$Coral_genus,levels = c ("Pocillopora", "Pachyseris", "Leptoseris","Montipora","Porites", "Acropora", "Pavona",  "Echynophillia","Leptastrea","Astreopora","Others"))

coral_data_compo <- coral_data_compo %>% unite(Island_Island_Site, Island, Island_Site)


coral_data_compo$Island_Island_Site <- gsub('Moorea_1', 'MOO1', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Moorea_2', 'MOO2', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Tahiti_1', 'TAH1', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Tahiti_2', 'TAH2', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Bora_1', 'BOR1', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Bora_2', 'BOR2', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Tikehau_1', 'TIK1', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Tikehau_2', 'TIK2', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Rangiroa_1', 'RAN1', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Rangiroa_2', 'RAN2', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Raroia_1', 'RAR1', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Raroia_2', 'RAR2', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Makatea_1', 'MAK1', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Makatea_2', 'MAK2', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Gambier_1', 'GAM1', coral_data_compo$Island_Island_Site)
coral_data_compo$Island_Island_Site <- gsub('Gambier_2', 'GAM2', coral_data_compo$Island_Island_Site)

# coral_data_main_Island_Site$Island_Island_Site <- factor(coral_data_main_Island_Site$Island_Island_Site, levels = c("Moorea_1","Moorea_2","Tahiti_1","Tahiti_2","Bora_1","Bora_2","Tikehau_1","Tikehau_2","Rangiroa_2","Rangiroa_3","Raroia_1","Raroia_2","Makatea_1","Makatea_2","Gambier_1","Gambier_2"))

coral_data_compo$Island_Island_Site <- factor(coral_data_compo$Island_Island_Site, levels = c("MOO1","MOO2","TAH1","TAH2","BOR1","BOR2","TIK1","TIK2","RAN1","RAN2","RAR1","RAR2","MAK1","MAK2","GAM1","GAM2"))


Fig_Sup_4 <- ggplot(coral_data_compo, aes(x=Island_Island_Site, y=Relative)) +
  geom_bar(aes(fill=factor(Coral_genus)), stat="identity", width = 0.3) +  facet_grid (rows = vars(Depth), switch = "both") +
  scale_fill_manual(values= colours_corals) +
  # geom_line(aes(y=fit_se), size=1, linetype="dotted",alpha = 0.8, col = "blue") +
  # geom_line(aes(y=fit_seneg), size=1, linetype="dotted",alpha = 0.8, col = "blue") +
  scale_y_continuous(position = "right",breaks = c(0,50,100))+ scale_x_discrete ()+ ylab ("Percentatge (%)") + xlab ("Islands") + 
  ggtitle("Relative community composition")+
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=16, face="bold"),
                          axis.text.y = element_text(size=14, colour="black"),
                          axis.text.x = element_text(size=12, colour="black"),
                          strip.text.y = element_text(size=16, colour="black"),
                          axis.title = element_text(size=14, face="bold", colour="black"), legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 14, title = NULL))
Fig_Sup_4
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_Rel_Comm_Domin.pdf", Fig_Sup_4, width = 20, height = 10)









### NMDS from community matrix of coral cover
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


# Prepare for ggplot 
CC_NMDS_Species <- total_NMDS[['species']]
CC_NMDS_Species <- as.data.frame(CC_NMDS_Species)
CC_NMDS_Species$Species <- rownames (CC_NMDS_Species)

CC_NMDS_coordinates <- total_NMDS[['points']]
CC_NMDS_coordinates <- as.data.frame(CC_NMDS_coordinates)
CC_NMDS_coordinates$ID <- rownames (CC_NMDS_coordinates)
# Extract depths
CC_NMDS_coordinates$Depth <- sub("\\_.*", "", CC_NMDS_coordinates$ID)
# Extract Island Sites
CC_NMDS_coordinates$Island <- sub("^[^_]*_", "", CC_NMDS_coordinates$ID)
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

# Plot NMDS
CC_NMDS_coordinates = CC_NMDS_coordinates %>% mutate(Depth = factor(Depth, levels = c(6, 20, 40, 60, 90, 120)))
hull_CC_NMDS <- CC_NMDS_coordinates %>% group_by(Depth) %>% slice(chull(MDS1,MDS2))

keep_species <-  c ("Acropora","Astreopora","Montipora","Leptastrea","Leptoseris","Pachyseris","Pavona","Echynophillia","Fungia","Pocillopora","Porites")
CC_NMDS_Species <- CC_NMDS_Species[CC_NMDS_Species$Species %in% keep_species, ]

Fig_3B = CC_NMDS_coordinates %>% ggplot() + theme_blank() + 
  geom_point(aes(x = MDS1, y = MDS2, fill = Depth), shape = 21, show.legend = F) +
  geom_polygon(data = hull_CC_NMDS, aes(x = MDS1, y = MDS2, fill = Depth), alpha = 0.75, color = "black", show.legend = T) +
  #geom_text(data = CC_NMDS_Species, aes(x = MDS1, y = MDS2, label = Species), size = 2) +
  geom_label(data = CC_NMDS_Species, aes(x = MDS1, y = MDS2, label = Species), size = 2, alpha = .75) +
  scale_x_continuous(name ="NMDS 1") + scale_y_continuous(name ="NMDS 2") + ggtitle("Bray-Curtis - Coral cover") +
  fishualize::scale_fill_fish_d(option = "Ostracion_whitleyi", direction = -1) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=11, face="bold")) +
  annotate(geom = 'text', label = 'Stress = 0.15', 
           x = range(CC_NMDS_coordinates$MDS1)[1], y = range(CC_NMDS_Species$MDS2)[1], hjust = 0, vjust = 6)


(Figure_3 = Fig_3A + Fig_3B + plot_layout(guides = 'collect'))

ggsave( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Figure_3.pdf", Figure_3,  width = 10, height = 5)
  
# Old plot without ggplot
# treat=c(rep("6M",16),rep("20M",16), rep("40M",16), rep("60M",16), rep("90M",16), rep("120M",14)) # Alphabetic order for colours
# colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# 
# 
# par (mfrow =c(1,1))
# # 120, 20, 40, 60, 6, 90)
# pdf("~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/NMDS_Diversity_Coral_Cover_Bray.pdf", bg = "white") # starts writing a PDF to file
# ordiplot(total_NMDS,type="n", choices = c(1,2)) # Check with numbers and see if I can get three axes
# orditorp(total_NMDS,display="sites",labels =T, cex = 0.1)
# orditorp(total_NMDS,display="sites",cex=0.3,air=0.01)
# ordispider(total_NMDS, groups=treat, col = colours,cex=0.6, lwd = 1.5)
# ordihull(total_NMDS,groups=treat,border = colours, col = colours, draw="polygon",label=F, lwd = 1.5)
# orditorp(total_NMDS,display="species",col="red",air=0.01, cex =0.5)
# legend(x="topright", y="top", legend=c ("120m","90m","60m", "40m","20m","6m"), col=c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"), fill = c("black", "navyblue","blue","deepskyblue","aquamarine2","wheat"),  pt.cex=0.8, cex = 0.70, horiz = F, text.width = 0.3, text.font = 20)
# title ("Bray distance - Coral cover")
# mysubtitle <- paste0("Stress = ", format(round(total_NMDS$stress, 2)))
# mtext(mysubtitle, side=1, line=-2, at=0, adj=0, cex=0.7) 
# # dev.off()




############# Heatmap of beta-diversity at the site level - Abundance Bray Curtis - Figure 3 A and B####


resume_df <- ddply(rndpts_df, ~ Island + Island_Site +  Depth  +  Coral_genus ,
                   function(x){c(Cover= sum(x$Cover)) })
melt_temp = melt(resume_df, id=c("Island","Island_Site","Depth","Coral_genus"), measure.vars="Cover", na.rm=FALSE)
cast_temp = dcast(melt_temp,Island + Island_Site +  Depth  ~ Coral_genus, mean, add.missing = T)
cast_temp[is.na(cast_temp)] <- 0
cast_temp$Island <- factor(cast_temp$Island, levels = c("Bora","Gambier","Makatea","Moorea","Rangiroa","Raroia","Tahiti","Tikehau"))
cast_all_depth <- cast_temp
# Necessary for distance afterwards

# complete to have all quadrats although empty 
cast_all_depth <- cast_all_depth %>% complete( Island,Island_Site,Depth,fill = list(Cover = 0))
cast_all_depth[is.na(cast_all_depth)] <- 0

cast_all_depth$ID<- with(cast_all_depth, paste0(Island, sep = "_",  Island_Site, sep = "_Depth_",Depth))
cast_all_depth <- as.matrix (cast_all_depth)
row.names(cast_all_depth) <- cast_all_depth[,"ID"]

cast_all_depth <- cast_all_depth[,-c(1,2,3,4,37)]

class(cast_all_depth) <- "numeric"
### Unnecessary now #####
# columns where sum is 0
# cast_all_depth <- cast_all_depth[,colSums(cast_all_depth[,])>0]
# rows where sum is 0
# cast_all_depth <- cast_all_depth[rowSums(cast_all_depth[,])>0, ]
### Unnecessary now #####
# Measure pair.abund
coral.matrices_Depth <- beta.pair.abund(cast_all_depth, index.family = "bray")


island_site <-  c("Bora_1", "Bora_2", "Gambier_1", "Gambier_2", "Makatea_1", "Makatea_2", "Moorea_1", "Moorea_2", "Rangiroa_1",
                  "Rangiroa_2", "Raroia_1", "Raroia_2", "Tahiti_1", "Tahiti_2", "Tikehau_1", "Tikehau_2")

depth <- c(6,20,40,60,90,120)

summary_final <- data.frame()
Bray_final  <- data.frame()
for (i in unique (island_site)){
  Beta_Depth <-  coral.matrices_Depth$beta.bray
  mat_bray <- as.matrix (Beta_Depth)
  mat_bray_Island <- mat_bray[grepl(i, rownames(mat_bray)),grepl(i, colnames(mat_bray))]
  
  mat_bray = melt(mat_bray_Island, measure.vars="nobserv", na.rm=FALSE)
  
  colnames (mat_bray) <- c("Var1", "Var2", "nobserv")
  
  mat_bray$Depth_Var1 <- as.numeric (sub('.*\\_', '', mat_bray$Var1))
  mat_bray$Depth_Var2 <- as.numeric (sub('.*\\_', '', mat_bray$Var2))
  
  
  # Comparing the overlap of 40 or the other depths, with the rest of depths
  # AS 40m is the reference
  # Necessary to change manually
  mat_bray<-mat_bray[(mat_bray$Depth_Var1==40 ),]
  
  mat_bray$SDist <- (mat_bray$Depth_Var1 - mat_bray$Depth_Var2)
  
  # # All pairs between themselves out
  mat_bray<-mat_bray[!(mat_bray$nobserv==0),]
  # # All pairs between the same depth factors
  mat_bray<-mat_bray[!(mat_bray$SDist==0),]
  
  colnames(mat_bray)[3] <- "Bray_dissimilarity"
  
  
  # Add the other two components
  # Turnover
  Turn_Depth <- coral.matrices_Depth$beta.bray.bal
  mat_bray_turn <- as.matrix (Turn_Depth)
  mat_bray_turn_Island <- mat_bray_turn[grepl(i, rownames(mat_bray_turn)),grepl(i, colnames(mat_bray_turn))]
  mat_bray_turn = melt(mat_bray_turn_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames(mat_bray_turn) <- c("Var1", "Var2","Bray_turnover")
  # The merge keeping only all mat_bray already cleans and uses 6m as reference 
  mat_bray <- merge(mat_bray, mat_bray_turn, by=c("Var1","Var2"), all.x = T, all.y = F)
  
  # Nestedness
  Nest_Depth <- coral.matrices_Depth$beta.bray.gra
  mat_bray_nest <- as.matrix (Nest_Depth)
  mat_bray_nest_Island <- mat_bray_nest[grepl(i, rownames(mat_bray_nest)),grepl(i, colnames(mat_bray_nest))]
  mat_bray_nest = melt(mat_bray_nest_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames(mat_bray_nest) <- c("Var1", "Var2","Bray_nestedness")
  # The merge keeping only all mat_bray already cleans and uses 6m as reference 
  mat_bray <- merge(mat_bray, mat_bray_nest, by=c("Var1","Var2"), all.x = T, all.y = F)
  
  
  # Delete all quadrats without corals 
  mat_bray <- na.omit (mat_bray)
  
  # Compute the index suggested by Seb Villeger - indicator of importance of turnover "Toussaint et al 2014"
  #  mat_bray$pturn = mat_bray$Turn_Similarity / mat_bray$Jac_Similarity
  mat_bray$pturn_Diss = mat_bray$Bray_turnover / mat_bray$Bray_dissimilarity
  
  # Means and Sd in funciton of SDist
  summary_mat_bray <- ddply(mat_bray, .(Depth_Var1,Depth_Var2,SDist), summarize, Bray=mean(Bray_dissimilarity), Bray_se=sd(Bray_dissimilarity) / sqrt(length(Bray_dissimilarity)),
                            Turnover=mean(Bray_turnover), Turn_se=sd(Bray_turnover) / sqrt(length(Bray_turnover)),
                            Nestedness=mean(Bray_nestedness), Nest_se=sd(Bray_nestedness) / sqrt(length(Bray_nestedness)), 
                            # P_Turn=mean(pturn, na.omit=T), P_Turn_se=sd(pturn, na.omit=T) / sqrt(length(pturn, na.omit=T)),
                            P_Turn_Diss=mean(pturn_Diss), P_Turn_Diss_se=sd(pturn_Diss) / sqrt(length(pturn_Diss)))
  
  
  summary_2 <- melt(summary_mat_bray, id.vars=c("Depth_Var1","Depth_Var2","SDist"), measure.vars = c("Bray", "Turnover", "Nestedness",  "P_Turn_Diss"))
  colnames (summary_2) <- c("Depth_Var1","Depth_Var2","SDist","Component", "Dissimilarity")
  
  summary_se <- melt(summary_mat_bray, id.vars=c("Depth_Var1","Depth_Var2","SDist"), measure.vars = c("Bray_se", "Turn_se", "Nest_se",  "P_Turn_Diss_se"))
  colnames (summary_se) <- c("Depth_Var1","Depth_Var2","SDist","Component", "Se")
  
  summary_se$Component <- gsub('Bray_se', 'Bray',  summary_se$Component)
  summary_se$Component <- gsub('Turn_se', 'Turnover',  summary_se$Component)
  summary_se$Component<- gsub('Nest_se', 'Nestedness',  summary_se$Component)
  # summary_se$Component<- gsub('P_Turn_se', 'P_Turn',  summary_se$Component)
  summary_se$Component<- gsub('P_Turn_Diss_se', 'P_Turn_Diss',  summary_se$Component)
  
  summary <- merge (summary_2,  summary_se)
  
  keep <- c("Bray", "Turnover","P_Turn_Diss")
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

pturn$Dissimilarity[pturn$Component=="Bray"] 

pturn$Dissimilarity[pturn$Component=="Turnover"] 

############# Heatmap of beta-diversity at the site level - Abundance Bray Curtis - Figure 3 A and B####

# Necessary to fill a table manually!!!! And make the calculation for each depth separately! later open the New_database
# Go to the Script of PA!!! Figure of heatmap is made there!

########################################## End of heat map ############################################

# !!!!!! FROM HERE REDUCE EVERYTHING (with the three points) UNTIL Spatial beta diversity !!!!!! #


################## Vertical beta diversity ##########################

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






################## Distance decay dissimilarity ######################
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
# cast_all_depth <- cast_all_depth[,colSums(cast_all_depth[,])>0]
# rows where sum is 0
# cast_all_depth <- cast_all_depth[rowSums(cast_all_depth[,])>0, ]
### Unnecessary now #####
# Measure pair.abund
coral.matrices_Depth <- beta.pair.abund(cast_all_depth, index.family = "bray")



### Trying to automatize the code below ######

island_site <-  c("Bora_1", "Bora_2", "Gambier_1", "Gambier_2", "Makatea_1", "Makatea_2", "Moorea_1", "Moorea_2", "Rangiroa_1",
                  "Rangiroa_2", "Raroia_1", "Raroia_2", "Tahiti_1", "Tahiti_2", "Tikehau_1", "Tikehau_2")

summary_final <- data.frame()
Bray_final <- data.frame ()
for (i in unique (island_site)){
  Beta_Depth <-  coral.matrices_Depth$beta.bray
  mat_bray <- as.matrix (Beta_Depth)
    mat_bray_Island <- mat_bray[grepl(i, rownames(mat_bray)),grepl(i, colnames(mat_bray))]
    
    mat_bray = melt(mat_bray_Island, measure.vars="Cover", na.rm=FALSE)
    
    mat_bray$Depth_Var1 <- as.numeric (gsub("(.*_){2}(\\d+)_.+", "\\2", mat_bray$Var1))
    mat_bray$Depth_Var2 <- as.numeric (gsub("(.*_){2}(\\d+)_.+", "\\2", mat_bray$Var2))
    
    ####### adding this new# Make the comparison against shallow waters! 6 and 20m together grouped as 15 m. Not 6 m anymore
    mat_bray$Depth_Var1[mat_bray$Depth_Var1 == 6] <- 15
    mat_bray$Depth_Var1[mat_bray$Depth_Var1 == 20] <- 15
    #
    mat_bray$Depth_Var2[mat_bray$Depth_Var2 == 6] <- 15
    mat_bray$Depth_Var2[mat_bray$Depth_Var2 == 20] <- 15
    #
    mat_bray$SDist <- abs(mat_bray$Depth_Var1 - mat_bray$Depth_Var2)

    ####### adding this new#  Against 15 m because 6 and 20m (as surface reference) instead of 6 m
    mat_bray<-mat_bray[(mat_bray$Depth_Var1==15),]

    # All pairs between themselves out
    mat_bray<-mat_bray[!(mat_bray$value==0),]
    # All pairs between the same depth factors
    mat_bray<-mat_bray[!(mat_bray$SDist==0),]

    
    # 
    # as 6m is reference keep only the compairsons against 6m (we do not want comparisons between 90 and 120m)
    # mat_bray<-mat_bray[(mat_bray$Depth_Var1==6 ),]
    
    # as 20m reference - necessary to delete the comparisons against 6 m also 
    # mat_bray<-mat_bray[(mat_bray$Depth_Var1==20 ),]
    # mat_bray<-mat_bray[!(mat_bray$Depth_Var2==6),]
    
    ###### LAetitia request Only compare dissimilariry between 40 and 60, and 90 and 120m #######
    # mat_bray$Depth_Var1[mat_bray$Depth_Var1 == 40] <- 50
    # mat_bray$Depth_Var1[mat_bray$Depth_Var1 == 60] <- 50
    # mat_bray$Depth_Var2[mat_bray$Depth_Var2 == 40] <- 50
    # mat_bray$Depth_Var2[mat_bray$Depth_Var2 == 60] <- 50
    # 
    # # Depth 90 and 120
    # mat_bray$Depth_Var1[mat_bray$Depth_Var1 == 90] <- 100
    # mat_bray$Depth_Var1[mat_bray$Depth_Var1 == 120] <- 100
    # mat_bray$Depth_Var2[mat_bray$Depth_Var2 == 90] <- 100
    # mat_bray$Depth_Var2[mat_bray$Depth_Var2 == 120] <- 100
    # 
    # mat_bray<-mat_bray[(mat_bray$Depth_Var1==50),]
    # mat_bray<-mat_bray[(mat_bray$Depth_Var2==100),]
    # 
    # mat_bray$SDist <- abs(mat_bray$Depth_Var1 - mat_bray$Depth_Var2)
    # 
    # # All pairs between themselves out
    # mat_bray<-mat_bray[!(mat_bray$value==0),]
    # # # All pairs between the same depth factors
    # mat_bray<-mat_bray[!(mat_bray$SDist==0),]
    # 
    ###### LAetitia request Only compare dissimilariry between 40 and 60, and 90 and 120m #######
    
    
    # Comparing the overlap of 40 or the other depths, m with the rest of depths
    # # AS 40m is the reference
    # mat_bray<-mat_bray[(mat_bray$Depth_Var1==120 ),]
    # 
    # mat_bray$SDist <- (mat_bray$Depth_Var1 - mat_bray$Depth_Var2)
    # 
    # # # All pairs between themselves out
    # mat_bray<-mat_bray[!(mat_bray$value==0),]
    # # # All pairs between the same depth factors
    # mat_bray<-mat_bray[!(mat_bray$SDist==0),]
    
    
    colnames(mat_bray)[3] <- "Bray_dissimilarity"
    
    
    linm <- lm(Bray_dissimilarity ~ SDist, data = mat_bray)
    linm_resume <- summary(linm)
    plot (Bray_dissimilarity ~ SDist, data = mat_bray, ylab = "Dissimilarity", xlab = "Vertical distance", cex = 0.5, main = paste0(i,sep = "_", "Bray"))
    mtext(paste0("R squared: ",round(linm_resume$r.squared,5)),adj = 0)
    mtext(paste0("P-value: ", format.pval(pf(linm_resume$fstatistic[1], # F-statistic
                                             linm_resume$fstatistic[2], # df
                                             linm_resume$fstatistic[3], # df
                                             lower.tail = FALSE))))
    mtext(paste0("y = ",round(linm_resume$coefficients[1],2)," + ",
                 round(linm_resume$coefficients[2],5),"x"),adj = 1)
    
    # Add the other two components
    # Turnover
    Turn_Depth <- coral.matrices_Depth$beta.bray.bal
    mat_bray_turn <- as.matrix (Turn_Depth)
    mat_bray_turn_Island <- mat_bray_turn[grepl(i, rownames(mat_bray_turn)),grepl(i, colnames(mat_bray_turn))]
    mat_bray_turn = melt(mat_bray_turn_Island, measure.vars="Cover", na.rm=FALSE)
    colnames(mat_bray_turn)[3] <- "Bray_turnover"
    # The merge keeping only all mat_bray already cleans and uses 6m as reference 
    mat_bray <- merge(mat_bray, mat_bray_turn, by=c("Var1","Var2"), all.x = T, all.y = F)
    
    # Nestedness
    Nest_Depth <- coral.matrices_Depth$beta.bray.gra
    mat_bray_nest <- as.matrix (Nest_Depth)
    mat_bray_nest_Island <- mat_bray_nest[grepl(i, rownames(mat_bray_nest)),grepl(i, colnames(mat_bray_nest))]
    mat_bray_nest = melt(mat_bray_nest_Island, measure.vars="Cover", na.rm=FALSE)
    colnames(mat_bray_nest)[3] <- "Bray_nestedness"
    # The merge keeping only all mat_bray already cleans and uses 6m as reference 
    mat_bray <- merge(mat_bray, mat_bray_nest, by=c("Var1","Var2"), all.x = T, all.y = F)
    
    # Delete all quadrats without corals 
    mat_bray <- na.omit (mat_bray)

      
    # Transform from Dissimilarity to Similarity
    mat_bray$Bray_Similarity <- 1 - mat_bray$Bray_dissimilarity
    # Necessary to measure the respective turnover and nestedness to the new Similarity (inverse dissimilarity)
    mat_bray$Turn_Similarity <- (mat_bray$Bray_turnover * mat_bray$Bray_Similarity) / mat_bray$Bray_dissimilarity
    mat_bray$Nest_Similarity <- (mat_bray$Bray_nestedness * mat_bray$Bray_Similarity) / mat_bray$Bray_dissimilarity
    
    # Compute the index suggested by Seb Villeger - indicator of importance of turnover "Toussaint et al 2014"
   # mat_bray$pturn = mat_bray$Turn_Similarity / mat_bray$Bray_Similarity
    mat_bray$pturn_Diss = mat_bray$Bray_turnover / mat_bray$Bray_dissimilarity
    
    # Means and Sd in funciton of SDist
    summary_mat_bray <- ddply(mat_bray, .(SDist), summarize, Bray=mean(Bray_Similarity), Bray_se=sd(Bray_Similarity) / sqrt(length(Bray_Similarity)),
                              Turnover=mean(Turn_Similarity), Turn_se=sd(Turn_Similarity) / sqrt(length(Turn_Similarity)),
                              Nestedness=mean(Nest_Similarity), Nest_se=sd(Nest_Similarity) / sqrt(length(Nest_Similarity)), 
                              P_Turn_Diss=mean(pturn_Diss), P_Turn_Diss_se=sd(pturn_Diss) / sqrt(length(pturn_Diss)))
    
    
    summary_2 <- melt(summary_mat_bray, id.vars=c("SDist"), measure.vars = c("Bray", "Turnover", "Nestedness", "P_Turn_Diss"))
    colnames (summary_2) <- c("SDist","Component", "Similarity")
    
    summary_se <- melt(summary_mat_bray, id.vars=c("SDist"), measure.vars = c("Bray_se", "Turn_se", "Nest_se", "P_Turn_Diss_se"))
    colnames (summary_se) <- c("SDist","Component", "Se")
    
    summary_se$Component <- gsub('Bray_se', 'Bray', summary_se$Component)
    summary_se$Component <- gsub('Turn_se', 'Turnover', summary_se$Component)
    summary_se$Component<- gsub('Nest_se', 'Nestedness', summary_se$Component)
    summary_se$Component<- gsub('P_Turn_Diss_se', 'P_Turn_Diss',  summary_se$Component)    
    
    summary <- merge (summary_2, summary_se)
    
    summary$Island_Site <- i
    
    summary_final = rbind (summary_final, summary)
    
    # For the model 
    mat_bray_all <- melt(mat_bray, id.vars=c("SDist"), measure.vars = c("Bray_Similarity", "Turn_Similarity", "Nest_Similarity"))
    colnames (mat_bray_all) <- c("SDist","Component", "Similarity")
    
    mat_bray_all$Island_Site <- i
    
    Bray_final = rbind (Bray_final, mat_bray_all)

}



# Measure mean and standard error
summary <- ddply(summary_final, ~ SDist + Component, summarize, Sim_mean=mean(Similarity), Sim_se=sd(Similarity) / sqrt(length(Similarity)))
summary_final <- merge (summary_final,summary)

m1 <- lm(Similarity ~ SDist*Component, data=summary_final)
summary(m1)
a1 <- m1$coef[2] # Bray slope
a2 <- m1$coef[2] + m1$coef[6] # Balanced variation (substitution) slope
a3 <- m1$coef[2] + m1$coef[7] # Gradient slope

b1 <- m1$coef[1]
b2 <- b1  + m1$coef[3] 

# Change to real names!
summary_final$Component <- gsub('Bray', 'Bray-Curtis', summary_final$Component)
summary_final$Component <- gsub('Turnover', 'Balanced variation', summary_final$Component)
summary_final$Component <- gsub('Nestedness', 'Abundance gradient', summary_final$Component)
summary_final$Component <- gsub('P_Turn_Diss', 'P_Turn_Diss', summary_final$Component)

summary_final$Component <- factor(summary_final$Component, levels = c("Bray-Curtis","Balanced variation","Abundance gradient", "P_Turn_Diss"))



# Geom point with colour for island and shape for component.
  
Fig_2B = summary_final %>% filter(., Component %in% c("Bray-Curtis", "Balanced variation")) %>% 
  ggplot(., aes(x=SDist, y=Similarity, colour = Component)) +
  scale_color_manual(values = c("#DC143C","#6495ED")) + 
  geom_point (size = 1, shape = 21, alpha = .75) +
  geom_errorbar(aes(ymin = Sim_mean - Sim_se, ymax = Sim_mean + Sim_se), size=1, width=2) +
#  geom_abline(intercept = b1, slope = a1,  colour ="#DC143C") +
#  geom_abline(intercept = b2, slope = a2,  colour ="#6495ED") +
  scale_x_continuous(name ="Vertical distance (m)", limits=c(15,110), breaks = c(25,50,75,100)) +
  scale_y_continuous(name ="Similarity", limits=c(0,0.75), breaks = c(0,0.25,0.5,0.75)) +
  ggtitle("Bray-Curtis - Coral cover") +
  geom_point(aes(y=Sim_mean, shape = Component, colour = Component), size=3, shape = 21, fill = "white") + 
  theme_classic() + theme(axis.text = element_text(size=10), axis.title = element_text(size=11, face="bold"), plot.title = element_text(hjust = 0.5, size=12, face="bold"))
Fig_2B
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Fig_2B_bis2.pdf", Fig_2B, width = 4, height = 5)


# Plot both graphs together. Jaccard and Bray
(Figure_2 = Fig_2A + Fig_2B + plot_layout(guides = 'collect'))

ggsave( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Figure_2.pdf", Figure_2,  width = 9, height = 4)
  

#### Model for predicting a bit the three components ####
# # Tried lmer and brm
# vertical.lmer <- lmer(Dissimilarity ~ SDist*Component  + (1|Island_Site) , data = Bray_final)
# # vertical.lmer <- glmer(Dissimilarity ~ SDist*Component  + (1|Island_Site) , data = Bray_final, family = poisson()) # Log model
# summary(vertical.lmer)
# 
# summary.lm(aov(Dissimilarity ~ SDist*Component , Bray_final))
# 
# 
# # These are values for jaccard
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
# # library (nlme)
# # nlme_vertical <- lme(Dissimilarity ~ SDist*Component, random = ~1|Island_Site, Bray_final)
# # anova(nlme_vertical)
# 
# Bray_final$fit <- predict(vertical.lmer)
# 
# ggplot(data= Bray_final, aes(x=SDist, y=fit,colour = Component))  +
#   geom_point(shape=8)+ facet_wrap(~Island_Site, nrow = 2, strip.position = "top") +
#   scale_x_continuous(name ="Vertical distance (m)", limits=c(5,120), breaks = c(0,20,40,60,80,100,120)) +
#   scale_y_continuous(name ="Bray dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
#   theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
#                       axis.text = element_text(size=10, colour="black"),
#                       axis.title = element_text(size=11, face="bold", colour="black")) 
# 
# resume_Bray_final <- ddply(Bray_final, ~ Island_Site + Component + SDist ,function(x){
#   c(fit=mean(x$fit)) }) 
# 
# ggplot(data= resume_Bray_final, aes(x=SDist, y=fit,colour = Component))  +
#   geom_point(shape=8)+ facet_wrap(~Island_Site, nrow = 2, strip.position = "top") +
#   scale_x_continuous(name ="Vertical distance (m)", limits=c(5,120), breaks = c(0,20,40,60,80,100,120)) +
#   scale_y_continuous(name ="Bray dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
#   theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
#                       axis.text = element_text(size=10, colour="black"),
#                       axis.title = element_text(size=11, face="bold", colour="black")) 


#### Model for predicting a bit the three components ####



# Change the names on axis
summary_final$SDist<- gsub('105', 'Surf.vs120m', summary_final$SDist)
summary_final$SDist<- gsub('75', 'Surf.vs90m', summary_final$SDist)
summary_final$SDist<- gsub('45', 'Surf.vs60m', summary_final$SDist)
summary_final$SDist<- gsub('25', 'Surf.vs40m', summary_final$SDist)

summary_final$SDist <- factor(summary_final$SDist, levels = c("Surf.vs40m","Surf.vs60m","Surf.vs90m","Surf.vs120m"))

SFig_3B <- summary_final %>% filter(., Component %in% c("Bray-Curtis", "Balanced variation")) %>% 
  ggplot(., aes(x=SDist, y=Similarity, fill=Component)) + 
  scale_color_manual(values = c("#DC143C","#6495ED")) + 
  geom_bar(stat="identity", position=position_dodge(), width = 0.9) + facet_wrap(~Island_Site, nrow = 2, strip.position = "top") +
  geom_errorbar(aes(ymin = Similarity - Se, ymax = Similarity + Se),size =0.3, width = 0.3, position = position_dodge(0.9))+
  scale_y_continuous(limits = c(0, 0.67), breaks = seq(0, 0.7, by = 0.1))+
  ylab ("Similarity") + xlab ("Vertical distance") + 
  theme_bw()  + theme(
    axis.text = element_text(size=10, colour="black"),
    axis.title = element_text(size=11, face="bold", colour="black"), 
    strip.text = element_text(size = 10, colour="black"),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), legend.position = "none") 
SFig_3B
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/SFig_3B.pdf", SFig_3B,width = 14, height = 4)



# Measure of the mean of the index of Villeger 
pturn <- aggregate( Similarity ~ SDist + Component, summary_final, mean )

pturn$Similarity[pturn$Component=="P_Turn_Diss"] 





###############

# Trying to get the database for heatmap
island_site <-  c("Bora_1", "Bora_2", "Gambier_1", "Gambier_2", "Makatea_1", "Makatea_2", "Moorea_1", "Moorea_2", "Rangiroa_1",
                  "Rangiroa_2", "Raroia_1", "Raroia_2", "Tahiti_1", "Tahiti_2", "Tikehau_1", "Tikehau_2")

depth <- c(6,20,40,60,90,120)

summary_final <- data.frame()
Bray_final  <- data.frame()
for (i in unique (island_site)){
  Beta_Depth <-  coral.matrices_Depth$beta.bray
  mat_bray <- as.matrix (Beta_Depth)
  mat_bray_Island <- mat_bray[grepl(i, rownames(mat_bray)),grepl(i, colnames(mat_bray))]
  
  mat_bray = melt(mat_bray_Island, measure.vars="nobserv", na.rm=FALSE)
  
  mat_bray$Depth_Var1 <- as.numeric (gsub("(.*_){2}(\\d+)_.+", "\\2", mat_bray$Var1))
  mat_bray$Depth_Var2 <- as.numeric (gsub("(.*_){2}(\\d+)_.+", "\\2", mat_bray$Var2))
  
  
  # Comparing the overlap of 40 or the other depths, with the rest of depths
  # AS 40m is the reference
  mat_bray<-mat_bray[(mat_bray$Depth_Var1==40 ),]
  
  mat_bray$SDist <- (mat_bray$Depth_Var1 - mat_bray$Depth_Var2)
  
  # # All pairs between themselves out
  mat_bray<-mat_bray[!(mat_bray$value==0),]
  # # All pairs between the same depth factors
  mat_bray<-mat_bray[!(mat_bray$SDist==0),]
  
  colnames(mat_bray)[3] <- "Bray_dissimilarity"
  
  
  # Add the other two components
  # Turnover
  Turn_Depth <- coral.matrices_Depth$beta.bray.bal
  mat_bray_turn <- as.matrix (Turn_Depth)
  mat_bray_turn_Island <- mat_bray_turn[grepl(i, rownames(mat_bray_turn)),grepl(i, colnames(mat_bray_turn))]
  mat_bray_turn = melt(mat_bray_turn_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames(mat_bray_turn)[3] <- "Bray_turnover"
  # The merge keeping only all mat_bray already cleans and uses 6m as reference 
  mat_bray <- merge(mat_bray, mat_bray_turn, by=c("Var1","Var2"), all.x = T, all.y = F)
  
  # Nestedness
  Nest_Depth <- coral.matrices_Depth$beta.bray.gra
  mat_bray_nest <- as.matrix (Nest_Depth)
  mat_bray_nest_Island <- mat_bray_nest[grepl(i, rownames(mat_bray_nest)),grepl(i, colnames(mat_bray_nest))]
  mat_bray_nest = melt(mat_bray_nest_Island, measure.vars="nobserv", na.rm=FALSE)
  colnames(mat_bray_nest)[3] <- "Bray_nestedness"
  # The merge keeping only all mat_bray already cleans and uses 6m as reference 
  mat_bray <- merge(mat_bray, mat_bray_nest, by=c("Var1","Var2"), all.x = T, all.y = F)
  
  
  # Delete all quadrats without corals 
  mat_bray <- na.omit (mat_bray)
  
  # Compute the index suggested by Seb Villeger - indicator of importance of turnover "Toussaint et al 2014"
  #  mat_bray$pturn = mat_bray$Turn_Similarity / mat_bray$Jac_Similarity
  mat_bray$pturn_Diss = mat_bray$Bray_turnover / mat_bray$Bray_dissimilarity
  
  # Means and Sd in funciton of SDist
  summary_mat_bray <- ddply(mat_bray, .(Depth_Var1,Depth_Var2,SDist), summarize, Bray=mean(Bray_dissimilarity), Bray_se=sd(Bray_dissimilarity) / sqrt(length(Bray_dissimilarity)),
                            Turnover=mean(Bray_turnover), Turn_se=sd(Bray_turnover) / sqrt(length(Bray_turnover)),
                            Nestedness=mean(Bray_nestedness), Nest_se=sd(Bray_nestedness) / sqrt(length(Bray_nestedness)), 
                            # P_Turn=mean(pturn, na.omit=T), P_Turn_se=sd(pturn, na.omit=T) / sqrt(length(pturn, na.omit=T)),
                            P_Turn_Diss=mean(pturn_Diss), P_Turn_Diss_se=sd(pturn_Diss) / sqrt(length(pturn_Diss)))
  
  
  summary_2 <- melt(summary_mat_bray, id.vars=c("Depth_Var1","Depth_Var2","SDist"), measure.vars = c("Bray", "Turnover", "Nestedness",  "P_Turn_Diss"))
  colnames (summary_2) <- c("Depth_Var1","Depth_Var2","SDist","Component", "Dissimilarity")
  
  summary_se <- melt(summary_mat_bray, id.vars=c("Depth_Var1","Depth_Var2","SDist"), measure.vars = c("Bray_se", "Turn_se", "Nest_se",  "P_Turn_Diss_se"))
  colnames (summary_se) <- c("Depth_Var1","Depth_Var2","SDist","Component", "Se")
  
  summary_se$Component <- gsub('Bray_se', 'Bray',  summary_se$Component)
  summary_se$Component <- gsub('Turn_se', 'Turnover',  summary_se$Component)
  summary_se$Component<- gsub('Nest_se', 'Nestedness',  summary_se$Component)
  # summary_se$Component<- gsub('P_Turn_se', 'P_Turn',  summary_se$Component)
  summary_se$Component<- gsub('P_Turn_Diss_se', 'P_Turn_Diss',  summary_se$Component)
  
  summary <- merge (summary_2,  summary_se)
  
  keep <- c("Bray", "Turnover","P_Turn_Diss")
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
summary_final$SDist <- factor(summary_final$SDist, levels = c("6 vs 20","6 vs 40","6 vs 60","6 vs 90","6 vs 120"))
summary_final$SDist <- factor(summary_final$SDist, levels = c("20 vs 6","20 vs 40","20 vs 60","20 vs 90","20 vs 120"))
summary_final$SDist <- factor(summary_final$SDist, levels = c("40 vs 6","40 vs 20","40 vs 60","40 vs 90","40 vs 120"))
summary_final$SDist <- factor(summary_final$SDist, levels = c("60 vs 6","60 vs 20","60 vs 40","60 vs 90","60 vs 120"))
summary_final$SDist <- factor(summary_final$SDist, levels = c("90 vs 6","90 vs 20","90 vs 40","90 vs 60","90 vs 120"))
summary_final$SDist <- factor(summary_final$SDist, levels = c("120 vs 6","120 vs 20","120 vs 40","120 vs 60","120 vs 90"))


# Measure the mean fo all locations
pturn <- aggregate(Dissimilarity ~ SDist + Component, summary_final, mean )

pturn$Dissimilarity[pturn$Component=="P_Turn_Diss"] 

pturn$Dissimilarity[pturn$Component=="Bray"] 

pturn$Dissimilarity[pturn$Component=="Turnover"] 

# Necessary to fill a table manually!!!! And make the calculation for each depth separately! later open the New_database

########################################################















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
                             beta_bray=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)),
                             Se = c(sd(Beta_Depth_Matrix_Depth_1) / sqrt(length(Beta_Depth_Matrix_Depth_1)), sd(Beta_Depth_Matrix_Depth_2) / sqrt(length(Beta_Depth_Matrix_Depth_2)),sd(Beta_Depth_Matrix_Depth_3) / sqrt(length(Beta_Depth_Matrix_Depth_3)),
                                    sd(Beta_Depth_Matrix_Depth_4) / sqrt(length(Beta_Depth_Matrix_Depth_4)),sd(Beta_Depth_Matrix_Depth_5) / sqrt(length(Beta_Depth_Matrix_Depth_5)),sd(Beta_Depth_Matrix_Depth_6) / sqrt(length(Beta_Depth_Matrix_Depth_6))))
beta_div_depth # These are the values for the table

# check linear model for Bray
beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_bray ~ Depth,beta_div_depth ))

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
                             beta_bray_bal=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)),
                             Se = c(sd(Beta_Depth_Matrix_Depth_1) / sqrt(length(Beta_Depth_Matrix_Depth_1)), sd(Beta_Depth_Matrix_Depth_2) / sqrt(length(Beta_Depth_Matrix_Depth_2)),sd(Beta_Depth_Matrix_Depth_3) / sqrt(length(Beta_Depth_Matrix_Depth_3)),
                                    sd(Beta_Depth_Matrix_Depth_4) / sqrt(length(Beta_Depth_Matrix_Depth_4)),sd(Beta_Depth_Matrix_Depth_5) / sqrt(length(Beta_Depth_Matrix_Depth_5)),sd(Beta_Depth_Matrix_Depth_6) / sqrt(length(Beta_Depth_Matrix_Depth_6))))
beta_div_depth

# check linear model for turnover
beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_bray_bal ~ Depth,beta_div_depth ))

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
                             beta_bray_gra=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)),
                             Se = c(sd(Beta_Depth_Matrix_Depth_1) / sqrt(length(Beta_Depth_Matrix_Depth_1)), sd(Beta_Depth_Matrix_Depth_2) / sqrt(length(Beta_Depth_Matrix_Depth_2)),sd(Beta_Depth_Matrix_Depth_3) / sqrt(length(Beta_Depth_Matrix_Depth_3)),
                                    sd(Beta_Depth_Matrix_Depth_4) / sqrt(length(Beta_Depth_Matrix_Depth_4)),sd(Beta_Depth_Matrix_Depth_5) / sqrt(length(Beta_Depth_Matrix_Depth_5)),sd(Beta_Depth_Matrix_Depth_6) / sqrt(length(Beta_Depth_Matrix_Depth_6))))
beta_div_depth # These are the values of the table

# check linear model for gradient
beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_bray_gra ~ Depth,beta_div_depth ))

# Mantel tests just for the Bray distance from the beta.pair
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
pdf("~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/Boxplot_Bray.pdf", bg = "white", width = 6, height = 4) 
boxplot (Beta_Depth_Matrix_Depth_1,Beta_Depth_Matrix_Depth_2, Beta_Depth_Matrix_Depth_3,Beta_Depth_Matrix_Depth_4,Beta_Depth_Matrix_Depth_5, Beta_Depth_Matrix_Depth_6, 
         xlab = "Depth (m)",
         ylab = "Beta.bray",
         names = c("6","20","40","60","90", "120"), 
         main = "Bray Curtis - Dominance")
dev.off()



# Introduce values manually: 

beta_div_depth <- data.frame(Depth=c("6", "20", "40", "60", "90", "120"), 
                             beta_bray=c(mean(Beta_Depth_Matrix_Depth_1),mean(Beta_Depth_Matrix_Depth_2), mean(Beta_Depth_Matrix_Depth_3), mean(Beta_Depth_Matrix_Depth_4), mean(Beta_Depth_Matrix_Depth_5), mean(Beta_Depth_Matrix_Depth_6)))


beta_div_depth$Depth <- as.numeric (beta_div_depth$Depth)
summary (lm(beta_bray~ Depth,beta_div_depth ))


############################################
# Get the entire beta.pair for plot
Beta_Depth <- coral.matrices_Depth$beta.bray
Beta_Depth_Matrix <- as.matrix (Beta_Depth)

Bray_spa = melt(Beta_Depth_Matrix, measure.vars="Bray", na.rm=FALSE)

Bray_spa$Depth_Var1 <- as.numeric (str_extract(Bray_spa$Var1, "[^_]+"))
Bray_spa$Depth_Var2 <- as.numeric (str_extract(Bray_spa$Var2, "[^_]+"))

# Keep only pairs between same depths
Bray_spa$VDist <- abs(Bray_spa$Depth_Var1 - Bray_spa$Depth_Var2)
Bray_spa<-Bray_spa[(Bray_spa$VDist==0),]

# Delete pairs that are the same "6_Bora_1" vs "6_Bora_1"
same <- which(!(Bray_spa$Var1==Bray_spa$Var2))
Bray_spa<-Bray_spa[same, ]

# Change names and only keep the necessary columns 
colnames(Bray_spa)[3] <- "Bray_dissimilarity"
colnames(Bray_spa)[4] <- "Depth"
Bray_spa <- Bray_spa[c("Var1","Var2","Bray_dissimilarity","Depth")]

ggplot(data = Bray_spa, aes(y = Bray_dissimilarity, x = as.factor(Depth))) + 
  geom_boxplot()  + 
  ylab ("Bray_dissimilarity") + xlab ("Depth (m)") + 
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black")) 

# Add the other two components
# Turnover
Turn_Depth <- coral.matrices_Depth$beta.bray.bal
Turn_Depth_Matrix <- as.matrix (Turn_Depth)

Turn_spa = melt(Turn_Depth_Matrix, measure.vars="Turn", na.rm=FALSE)
colnames(Turn_spa)[3] <- "Bray_turnover"
# The merge keeps already only the rows we want based on "Var1" and "Var2"
Bray_spa <- merge(Bray_spa, Turn_spa, by=c("Var1","Var2"), all.x = T, all.y = F)

# Nestedness
Nest_Depth <- coral.matrices_Depth$beta.bray.gra
Nest_Depth_Matrix <- as.matrix (Nest_Depth)

Nest_spa = melt(Nest_Depth_Matrix, measure.vars="Nest", na.rm=FALSE)
colnames(Nest_spa)[3] <- "Bray_nestedness"
# The merge keeps already only the rows we want based on "Var1" and "Var2"
Bray_spa <- merge(Bray_spa, Nest_spa, by=c("Var1","Var2"), all.x = T, all.y = F)



# Means and Se in funciton of Depth
summary_Bray_spa <- ddply(Bray_spa, .(Depth), summarize, Bray=mean(Bray_dissimilarity), Bray_se=sd(Bray_dissimilarity) / sqrt(length(Bray_dissimilarity)),
                          Turn=mean(Bray_turnover), Turn_se=sd(Bray_turnover) / sqrt(length(Bray_turnover)),
                          Nest=mean(Bray_nestedness), Nest_se=sd(Bray_nestedness) / sqrt(length(Bray_nestedness)))

summary_3 <- melt(summary_Bray_spa, id.vars=c("Depth"), measure.vars = c("Bray", "Turn", "Nest"))
colnames (summary_3) <- c("Depth","Component", "Dissimilarity")

summary_se_3 <- melt(summary_Bray_spa, id.vars=c("Depth"), measure.vars = c("Bray_se", "Turn_se", "Nest_se"))
colnames (summary_se_3) <- c("Depth","Component", "Se")

summary_se_3$Component <- gsub('Bray_se', 'Bray', summary_se_3$Component)
summary_se_3$Component <- gsub('Turn_se', 'Turn', summary_se_3$Component)
summary_se_3$Component<- gsub('Nest_se', 'Nest', summary_se_3$Component)

summary_spa <- merge (summary_3, summary_se_3)


# Plot of the three components from summary 
dis_vs_spa <- ggplot(data= summary_spa, aes(x=Depth, y=Dissimilarity, shape = Component, colour = Component))  +
  geom_point(aes( shape = Component, colour = Component), size=4) + 
  geom_errorbar(aes(ymin = Dissimilarity - Se, ymax = Dissimilarity + Se), size=0.5, width=2, colour = "black") +
  scale_x_continuous(name ="Depth (m)", limits=c(5,120), breaks = c(0,20,40,60,80,100,120)) +
  scale_y_continuous(name ="Bray dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"))
dis_vs_spa 


ggplot(data= summary_spa, aes(x=Depth, y=Dissimilarity, shape = Component, colour = Component))  +
  geom_point(aes( shape = Component, colour = Component), size=4) + facet_wrap(~Component, nrow = 3)+
  geom_errorbar(aes(ymin = Dissimilarity - Se, ymax = Dissimilarity + Se), size=1, width=2, colour = "black") +
  scale_x_continuous(name ="Depth (m)", limits=c(5,120), breaks = c(0,20,40,60,80,100,120)) +
  scale_y_continuous(name ="Bray dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  theme_bw()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                      axis.text = element_text(size=10, colour="black"),
                      axis.title = element_text(size=11, face="bold", colour="black"))


# Working with the entire database pairs for boxplot and mean and se
Bray_spa_all <- melt(Bray_spa, id.vars=c("Depth"), measure.vars = c("Bray_dissimilarity", "Bray_turnover"))
colnames (Bray_spa_all) <- c("Depth","Component", "Dissimilarity")

# Measure mean and standard error
summary <- ddply(Bray_spa_all, ~Depth + Component, summarize, Dissimilarity_mean=mean(Dissimilarity), Se=sd(Dissimilarity) / sqrt(length(Dissimilarity)))
Bray_spa_all <- merge (Bray_spa_all,summary)

SFig_6B <- ggplot(data = Bray_spa_all, aes(y = Dissimilarity, x = as.factor(Depth), colour = Component)) + 
  geom_boxplot()  + 
  scale_color_manual(values = c("#DC143C","#6495ED")) + 
  geom_errorbar(aes(ymin = Dissimilarity_mean - Se, ymax = Dissimilarity_mean + Se), size=1, width=0.1, colour = "black") +
  geom_point(aes(y = Dissimilarity_mean), size=1.5, fill = "white", shape = 21) +
  facet_wrap(~Component, nrow = 3, scales = "free")+
  ylab ("Dissimilarity") + xlab ("Depth (m)") + ggtitle ("Bray-Curtis - Coral cover") +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), 
                          strip.background = element_blank(),
                          strip.text.x = element_blank(), 
                          legend.position = "right") 
SFig_6B
ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/SFig_6B.pdf", SFig_6B ,width = 6, height = 6)


summary.lm(aov(Dissimilarity_mean ~ Depth*Component, Bray_spa_all))

spatial.lmer <- lmer(Dissimilarity ~ Depth*Component  + (1|Component) , data = Bray_spa_all)
summary (spatial.lmer)



# Make a final figure: beta vs geographical distance (x=geo_dist, y=beta, points= pairs for all different given depths to see if it's correlated to 
# geographical distance. We know they are not because 

# From Bray_Spa, all pairs

Bray_spa$Site1 <- sub("^[^_]*_", "", Bray_spa$Var1)
Bray_spa$Site2 <- sub("^[^_]*_", "", Bray_spa$Var2)


# Merge the two dataframes of pairs and geographical distance
Spatial_Bray <- merge (Bray_spa, Spa_Dist)

# Same values as Beta_Div_Depth that are the ones I used for Mantel tests
summary_Spatial_Bray <- ddply(Spatial_Bray, ~Depth, summarize, Bray_mean=mean(Bray_dissimilarity), Se=sd(Bray_dissimilarity) / sqrt(length(Bray_dissimilarity)))
Spatial_Bray <- merge (Spatial_Bray,summary_Spatial_Bray)

# colours <- c("black", "aquamarine2","deepskyblue","blue","wheat","navyblue")
# Colours from fishualize as NMDS
colours <- c("#CBEAF1FF", "#73C2E1FF", "#3F73B1FF", "#28407AFF", "#222E47FF", "#191819FF")

SFig_5B <- ggplot(data=Spatial_Bray,aes(x=Spa_Distance, y=Bray_dissimilarity, colour = as.factor(Depth)))  +
  geom_point (size = 0.5)+
  facet_wrap(~Depth, nrow = 1, scales = "free")+
 # geom_line(aes(y=Bray_mean,colour = as.factor(Depth)), size=1) + 
 # geom_line(aes(y=Bray_mean-Se,colour = as.factor(Depth)), size=0.3, linetype="dotted",alpha = 0.8) +
 # geom_line(aes(y=Bray_mean+Se,colour = as.factor(Depth)), size=0.3, linetype="dotted",alpha = 0.8) +
  scale_colour_manual(values = c("#CBEAF1FF", "#73C2E1FF", "#3F73B1FF", "#28407AFF", "#222E47FF", "#191819FF")) +
  scale_x_continuous(name ="Geographical distance (km)", limits=c(0,2000), breaks = c(0,500,1000,1500,2000)) +
  scale_y_continuous(name ="Bray-Curtis dissimilarity", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  theme_classic()  + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                           strip.text.x = element_text(size = 12, color = "black", face = "bold"),
                           strip.background = element_rect(fill="grey"),
                           axis.text = element_text(size=10, colour="black"),
                           axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none") 
SFig_5B

ggsave ( "~/Documents/AAASea_Science/AAA_PhD_Thesis/Photoquadrats/PhD_Diversity_Depth/Figures/SFig_5B.pdf", SFig_5B ,width = 12, height = 3)

################## Spatial beta diversity ##########################3






