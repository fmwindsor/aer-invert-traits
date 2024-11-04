# Invertebrate functional trait variation along successional gradients #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line #


## 2 - Community analysis


#### Setup ####

# Clear environment
rm(list = ls())

# Set working directory
setwd("/Users/c1513054/Library/CloudStorage/OneDrive-CardiffUniversity/Documents/Research/Papers/Published/Advances in Ecological Research (Glacier Bay)/Macroinvertebrate traits")

# Load the previous script
source("Code/1_Data-wrangling.R")


#### Community structure #### 

### 1997 wide sweep of sites

## Macroinvertebrates

# Remove sample sites with no observations 
macro_sweep_dframe_clean <- macro_sweep_dframe[-31,] # Remove NUN1

#common_taxa <- names(which(colSums(macro_sweep_dframe_clean != 0) >= 3))

# Non-metric multidimensional scaling
nmds_macro_1997 <- metaMDS(macro_sweep_dframe_clean, 
                           noshare=(engine="isoMDS"), trymax = 500, distance = "bray",
                           k = 2, autotransform = T, tidy = T)

stressplot(nmds_macro_1997)
nmds_macro_1997$stress # this is quite high, but expected based on 10 streams and 20 years

# Collate the information for plotting the results
nmds_macro_1997_sites <- as.data.frame(scores(nmds_macro_1997)$sites) 
nmds_macro_1997_sites$stream <- str_sub(rownames(macro_sweep_dframe_clean), end = -2)
nmds_macro_1997_dframe <- left_join(nmds_macro_1997_sites, chem_sweep_dframe, 
                                    by = c("stream"="Stream"))

# Collate the information on the species influences 
nmds_macro_1997_species_dframe <- as.data.frame(nmds_macro_1997$species) 

# Make a column for the labels with capitalised genera names 
nmds_macro_1997_species_dframe$genus <- str_to_sentence(rownames(nmds_macro_1997_species_dframe))

# Plot the location of the species to show relative influence
plot1a <- ggplot() + 
  geom_point(aes(x = NMDS1, y = NMDS2), pch = 21, size = 3, fill = "grey40", data = nmds_macro_1997_dframe) + 
  geom_point(aes(x = MDS1, y = MDS2), size = 1, pch = 8, colour = "darkred", data = nmds_macro_1997_species_dframe) +
  geom_text_repel(aes(x = MDS1, y = MDS2, label = genus), box.padding = 0.5, size = 2.5, data = nmds_macro_1997_species_dframe) +
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        legend.background = element_blank(), title = element_text(size=10)) + 
  ylab("NMDS2") + 
  xlab("NMDS1") + 
  ggtitle("a")
plot1a

# Recalculate convex hulls for plotting (using chull)
conv_hull_plot_macro <- nmds_macro_1997_dframe %>%
  group_by(stream) %>%
  slice(chull(NMDS1, NMDS2))

plot1b <- ggplot() + 
  geom_polygon(aes(x = NMDS1, y = NMDS2, colour = stream), fill = NA, size = 1, 
               data = conv_hull_plot_macro, alpha = 0.2) + 
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.text = element_text(size = 8), axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8, colour = "black"), legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        legend.title = element_blank(), legend.background = element_blank(), 
        title = element_text(size=10), legend.key.height = unit(2, "mm")) + 
  ggtitle("b")
plot1b

# Plot NMDS results against stream age
plot1c <- ggplot() + 
  geom_point(aes(x = NMDS1, y = NMDS2, colour = Age.y.), size = 4, 
             data = nmds_macro_1997_dframe) +
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        legend.background = element_blank(), title = element_text(size=10)) + 
  annotate(geom = "text", label = "Stress = 0.18", x = 2.5, y = -0.92, size = 3) +
  ylab("NMDS2") + 
  xlab("NMDS1") + 
  ggtitle("c")
plot1c

grid.arrange(plot1a, plot1b, plot1c, ncol = 1)

## Meiofauna 

# Non-metric multidimensional scaling
nmds_meio_1997 <- metaMDS(meio_sweep_dframe[-c(30,33),], 
                          noshare=(engine="isoMDS"), trymax = 500, distance = "bray",
                          k = 2, autotransform = T, tidy = T)

stressplot(nmds_meio_1997)
nmds_meio_1997$stress # this is quite high, but expected based on 10 streams and 20 years

# Collate the information for plotting the results
nmds_meio_1997_sites <- as.data.frame(scores(nmds_meio_1997)$sites)
nmds_meio_1997_sites$stream <- meio_1997_dframe[-c(30,33, 76:80),"stream"]
nmds_meio_1997_sites$age <- as.numeric(meio_1997_dframe[-c(30,33, 76:80),"age"])

# Collate the information on the species influences 
nmds_meio_1997_species_dframe <- as.data.frame(nmds_meio_1997$species) 

# Make a column for the labels with capitalised genera names 
nmds_meio_1997_species_dframe$genus <- str_to_sentence(rownames(nmds_meio_1997_species_dframe))

# Plot the location of the species to show relative influence
plot1d <- ggplot() + 
  geom_point(aes(x = NMDS1, y = NMDS2), pch = 21, size = 3, fill = "grey40", data = nmds_meio_1997_sites) + 
  geom_point(aes(x = MDS1, y = MDS2), size = 1, pch = 8, colour = "darkred", data = nmds_meio_1997_species_dframe) +
  geom_text_repel(aes(x = MDS1, y = MDS2, label = genus), box.padding = 0.5, size = 2.5, data = nmds_meio_1997_species_dframe) +
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        legend.background = element_blank(), title = element_text(size=10)) + 
  ylab("NMDS2") + 
  xlab("NMDS1") + 
  ggtitle("d")
plot1d

# Recalculate convex hulls for plotting (using chull)
conv_hull_plot_meio <- nmds_meio_1997_sites %>%
  group_by(stream) %>%
  slice(chull(NMDS1, NMDS2))

# Add in an extra stream (so that the colours match across the two datasets)
conv_hull_plot_meio1 <- rbind(conv_hull_plot_meio, data.frame(NMDS1 = NA,NMDS2 = NA, stream = "Car", age = 240))

plot1e <- ggplot() + 
  geom_polygon(aes(x = NMDS1, y = NMDS2, colour = stream), fill = NA, size = 1, 
               data = conv_hull_plot_meio1, alpha = 0.2) + 
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.text = element_text(size = 8), axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8, colour = "black"), legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        legend.title = element_blank(), legend.background = element_blank(), 
        title = element_text(size=10), legend.key.height = unit(2, "mm")) + 
  ggtitle("e")
plot1e

# Plot NMDS results against stream age
plot1f <- ggplot() + 
  geom_point(aes(x = NMDS1, y = NMDS2, colour = age), size = 4, data = nmds_meio_1997_sites) +
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        legend.background = element_blank(), title = element_text(size=10)) + 
  annotate(geom = "text", label = "Stress = 0.21", x = -0.7, y = -1.2, size = 3) + 
  ylab("NMDS2") + 
  xlab("NMDS1") + 
  ggtitle("f")
plot1f

grid.arrange(plot1a, plot1d, plot1b, plot1e, plot1c, plot1f, nrow = 3)

## MGLM analyses

# Macroinvertebrates

#### Multivariate GLMs #### 
spp <- dplyr::select(macro_sweep_dframe_clean, -Agathon, -Synorthocladius, 
              -Psectrocladius, -Valvatidae, -Sweltsa, -Nematoda, -Paralimnophyes)
stream <- str_sub(rownames(macro_sweep_dframe_clean), end = -2)
trial <- left_join(data.frame(stream), chem_sweep_dframe, by = c("stream"="Stream"))
age <- trial[,"Age.y."]

model1 <- manyglm(as.matrix(spp) ~ age, 
                  family = "negative.binomial", maxiter = 100)

aov.many <- anova.manyglm(model1, p.uni = "adjusted")
print(aov.many)
print.manyglm(model1)
summary.manyglm(model1, symbolic.cor = TRUE, show.est = TRUE)
residuals.manyglm(model1)
plot(model1)
predict(model1, p.uni = "adjusted")

best.r.sq(as.matrix(spp) ~ age)

# Meiofauna 

#### Multivariate GLMs #### 
spp <- meio_sweep_dframe[-c(30,33),]
age <- as.numeric(meio_1997_dframe[-c(30,33, 76:80),"age"])

model2 <- manyglm(as.matrix(spp) ~ age, # final model structure after
                  family = "negative binomial")   # backwards selection

aov.many <- anova.manyglm(model2, p.uni = "adjusted")
print(aov.many)
print.manyglm(model2)
summary.manyglm(model2, symbolic.cor = TRUE, show.est = TRUE)
residuals.manyglm(model2)
plot(model2)
predict(model2, p.uni = "adjusted")

best.r.sq(as.matrix(spp) ~ age)


### WPC long term data 

## Macroinvertebrates

# Non-metric multidimensional scaling
nmds_macro_long <- metaMDS(macro_WPC_dframe, 
                           noshare=(engine="isoMDS"), trymax = 500, distance = "bray",
                           k = 2, autotransform = T, tidy = T)

stressplot(nmds_macro_long)
nmds_macro_long$stress # this is okay!

# Collate the information for plotting the results
nmds_macro_long_sites <- as.data.frame(scores(nmds_macro_long)$sites) 
nmds_macro_long_sites$date <- rownames(macro_WPC_dframe)

# Plot the location of the species to show relative influence
plot2a <- ggplot() + 
  geom_point(aes(x = NMDS1, y = NMDS2), pch = 21, size = 3, fill = "darkred", data = nmds_macro_long_sites) + 
  geom_text_repel(aes(x = NMDS1, y = NMDS2, label = date), box.padding = 0.5, size = 2.5, data = nmds_macro_long_sites) +
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        legend.background = element_blank(), title = element_text(size=10)) + 
  ylab("NMDS2") + 
  xlab("NMDS1") + 
  ggtitle("a")
plot2a

#### Multivariate GLMs #### 
spp <- macro_WPC_dframe
year <- as.numeric(rownames(macro_WPC_dframe))

model3 <- manyglm(as.matrix(spp) ~ year, 
                  family = "negative.binomial")

aov.many <- anova.manyglm(model3, p.uni = "adjusted")
print(aov.many)
print.manyglm(model3)
summary.manyglm(model3, symbolic.cor = TRUE, show.est = TRUE)
residuals.manyglm(model3)
plot(model3)
predict(model3, p.uni = "adjusted")

best.r.sq(as.matrix(spp) ~ year)


## Meiofauna

# Non-metric multidimensional scaling
nmds_meio_long <- metaMDS(meio_WPC_dframe, 
                           noshare=(engine="isoMDS"), trymax = 500, distance = "bray",
                           k = 2, autotransform = T, tidy = T)

stressplot(nmds_meio_long)
nmds_meio_long$stress # this is okay!

# Collate the information for plotting the results
nmds_meio_long_sites <- as.data.frame(scores(nmds_meio_long)$sites) 
nmds_meio_long_sites$year <- rownames(meio_WPC_dframe)

# Plot the location of the species to show relative influence
plot2b <- ggplot() + 
  geom_point(aes(x = NMDS1, y = NMDS2), pch = 21, size = 3, fill = "darkred", data = nmds_meio_long_sites) + 
  geom_text_repel(aes(x = NMDS1, y = NMDS2, label = year), box.padding = 0.5, size = 2.5, data = nmds_meio_long_sites) +
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        legend.background = element_blank(), title = element_text(size=10)) + 
  ylab("NMDS2") + 
  xlab("NMDS1") + 
  ggtitle("b")
plot2b

grid.arrange(plot2a, plot2b, nrow = 1)

# Multivariate GLMs

spp <- meio_WPC_dframe
year <- as.numeric(rownames(meio_WPC_dframe))

model4 <- manyglm(as.matrix(spp) ~ year, # final model structure after
                  family = "negative binomial")   # backwards selection

aov.many <- anova.manyglm(model4, p.uni = "adjusted")
print(aov.many)
print.manyglm(model1)
summary.manyglm(model1, symbolic.cor = TRUE, show.est = TRUE)
residuals.manyglm(model1)
plot(model1)
predict(model1, p.uni = "adjusted")

best.r.sq(as.matrix(spp) ~ year)

