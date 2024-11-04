# Invertebrate functional trait variation along successional gradients #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line #


## 3 - Trait wrangling


#### Setup ####

# Clear environment
rm(list = ls())

# Set working directory
setwd("/Users/c1513054/Library/CloudStorage/OneDrive-CardiffUniversity/Documents/Research/Papers/Advances in Ecological Research (Glacier Bay)/Macroinvertebrate traits")

# Load the previous script
source("Code/3_Trait-wrangling.R")


#### Functional analyses #### 

#### Sweep dataset

### Macroinvertebrates

## dbFD - multivariate functional diversity calculation

# Taxonomic data
tax.mat.s <- as.matrix(macro_sweep_clean_dframe1[,order(colnames(macro_sweep_clean_dframe1))]) # Loading site x taxa matrix
traits.s <- as.matrix(macro_sweep_traits_dframe[order(rownames(macro_sweep_traits_dframe)),]) # Loading taxa x trait matrix

## Mean trait method 

# Run the mean trait analysis
FD <- dbFD(traits.s, tax.mat.s)

# Store the data in a new data frame
macro_sweep_FD_dframe <- data.frame(Stream = str_sub(rownames(macro_sweep_clean_dframe1), end = -2),
                                    FRic = FD$FRic, FEve = FD$FEve, 
                                    FDiv = FD$FDiv, FDis = FD$FDis)

# Add in the chemical data to help visualise the results
test <- left_join(macro_sweep_FD_dframe, chem_sweep_dframe, by = "Stream")

plot3a <- ggplot() + 
  stat_summary(aes(x = Age.y., y = FRic, fill = Stream), pch = 21, size = 1, data = test) + 
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
  ylab("Functional richness (FRic)") + 
  xlab("") + 
  ggtitle("b")
plot3a

plot3b <- ggplot() + 
  stat_summary(aes(x = Age.y., y = FEve, fill = Stream), pch = 21, size = 1, data = test) + 
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
  ylab("Functional evenness (FEve)") + 
  xlab("") + 
  ggtitle("(b)")
plot3b

plot3c <- ggplot() + 
  stat_smooth(method = "lm", aes(x = Age.y., y = FDiv, fill = Stream), 
              data = test, se = TRUE, fill = "lightgrey", 
              formula = y ~ poly(x, 2, raw = TRUE), 
              colour = "black", alpha = 1) + 
  stat_summary(aes(x = Age.y., y = FDiv, fill = Stream), pch = 21, size = 1, 
               data = test) + 
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
  ylab("Functional divergence (FDiv)") + 
  xlab("Stream age (years since deglaciation)") +   
  ggtitle("c")
plot3c

plot3d <- ggplot() + 
  stat_summary(aes(x = Age.y., y = FDis, fill = Stream), pch = 21, size = 1, data = test) + 
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
  ylab("Functional dispersion (FDis)") + 
  xlab("Stream age (years since deglaciation)") + 
  ggtitle("(d)")
plot3d

#grid.arrange(plot3a, plot3b, plot3c, plot3d, nrow = 2, bottom = "Site age (years since deglaciation)")

plot_grid(plot3a, plot3b, plot3c, plot3d, ncol = 2, align = "hv")


# Revised plot
gdist <- gowdis(traits.s)
PCoA1 <- dudi.pco(gdist, nf = 2, scannf = FALSE)

plot1 <- ggplot(data = data.frame(PCoA1$li), aes(x = A1, y = A2)) +
  geom_point() +
  theme_bw() + 
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
  geom_label_repel(aes(label = rownames(data.frame(PCoA1$li))), size = 3) + 
  xlab("PCoA1") + ylab("PCoA2") + 
  ggtitle("a")
plot1

grid.arrange(plot1, plot3a, plot3c, layout_matrix = rbind(c(1,1,2,2),
                                                          c(1,1,3,3)))


# Statistical analysis 

model1 <- glm(FDiv ~ Age.y., family = "gaussian" (link = "identity"), data = test)
summary.lm(model1)



## CWM across stream age

# Calculate community weighted means across sites for analysis
trial <- functcomp(traits.s, tax.mat.s)

# Trim the dataframe to traits with values 
trial_sub <- trial[, !sapply(trial, is.character)]

# Extract the streams
trial_sub$Stream <- str_sub(rownames(trial_sub), end = -2)

# Bind the chemical data
cwm_sweep <- left_join(trial_sub, chem_sweep_dframe, by = "Stream")

bb <- cwm_sweep %>% dplyr::select(order(colnames(cwm_sweep)))

# Traits
flight <- melt(dplyr::select(bb, Age.y., AdultFlyingStrength_abbrev_Strong:AdultFlyingStrength_abbrev_Weak), "Age.y.")
dispersal <- melt(dplyr::select(bb, Age.y., Female_disp_abbrev_High:Female_disp_abbrev_Low), "Age.y.")
voltinism <- melt(dplyr::select(bb, Age.y., Voltinism_abbrev_Bi_multivoltine:Voltinism_abbrev_Univoltine), "Age.y.")
synchrony <- melt(dplyr::select(bb, Age.y., Emerge_synch_abbrev_Poorly:Emerge_synch_abbrev_Well), "Age.y.")
emergence <- melt(dplyr::select(bb, Age.y., Emerge_season_1_Fall:Emerge_season_1_Winter), "Age.y.")
size <- melt(dplyr::select(bb, Age.y., Max_body_size_abbrev_Large:Max_body_size_abbrev_Small), "Age.y.")
respiration <- melt(dplyr::select(bb, Age.y., Resp_abbrev_Gills:Resp_abbrev_Tegument), "Age.y.")
rheophily <- melt(dplyr::select(bb, Age.y., Rheophily_abbrev_depo: Rheophily_abbrev_eros), "Age.y.")
thermal <- melt(dplyr::select(bb, Age.y., Thermal_pref_Cold.cool.eurythermal..0.15.C.:Thermal_pref_Warm.eurythermal..15.30.C.), "Age.y.")
habit <- melt(dplyr::select(bb, Age.y., Habit_prim_Burrower:Habit_prim_Swimmer), "Age.y.")
food <- melt(dplyr::select(bb, Age.y., Feed_prim_abbrev_CF:Feed_prim_abbrev_SH), "Age.y.")

## Extract GAM results
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = flight))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = dispersal))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = voltinism))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = synchrony))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = emergence))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = size))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = respiration))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = rheophily))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = thermal))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = habit))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = food))


# Adult flight strength
flight_labs <- c("Strong flying strength", "Weak flying strength")
names(flight_labs) <- c("AdultFlyingStrength_abbrev_Strong", "AdultFlyingStrength_abbrev_Weak")

plot4a<- ggplot(flight, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  theme_bw() +  
  guides(fill = "none") +
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = flight_labs), nrow = 1) + 
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_discrete(labels = c("Strong", "Weak"), name = "Flight") + 
  ylab("Flying strength") + 
  xlab("") + 
  ggtitle("a")
plot4a

# Female dispersal ability
dispersal_labs <- c("High female dispersal", "Low female dispersal")  
names(dispersal_labs) <- c("Female_disp_abbrev_High", "Female_disp_abbrev_Low")

plot4b<- ggplot(dispersal, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = dispersal_labs), nrow = 1) + 
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_discrete(labels = c("High", "Low"), name = "Dispersal") +
  ylab("Female dispersal") + 
  xlab("") + 
  ggtitle("g")
plot4b

# Voltinism
voltinism_labs <- c("Bivoltine", "Semivoltine", "Univoltine")
names(voltinism_labs) <- c("Voltinism_abbrev_Bi_multivoltine", "Voltinism_abbrev_Semivoltine",    
                           "Voltinism_abbrev_Univoltine")

plot4c<- ggplot(voltinism, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = voltinism_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Bivoltine", "Semivoltine", "Univoltine"), name = "Voltinism") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Voltinism") + 
  xlab("") + 
  ggtitle("b")
plot4c

# Emergence synchrony
synchrony_labs <- c("High emergence synchrony", "Low emergence synchrony")
names(synchrony_labs) <- c("Emerge_synch_abbrev_Poorly", "Emerge_synch_abbrev_Well")

plot4d <- ggplot(synchrony, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = synchrony_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("High", "Low"), name = "Emergence") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Emergence synchrony") + 
  xlab("") + 
  ggtitle("h")
plot4d

# Emergence season
emergence_labs <- c("Fall emergence", "Spring emergence", "Summer emergence", "Winter emergence")
names(emergence_labs) <- c("Emerge_season_1_Fall", "Emerge_season_1_Spring", 
                           "Emerge_season_1_Summer", "Emerge_season_1_Winter")

plot4e <- ggplot(emergence, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = emergence_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Fall", "Spring", "Summer", "Winter"), name = "Emergence season") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Emergence season") + 
  xlab("") + 
  ggtitle("c")
plot4e

# Body size
size_labs <- c("Large (>16 mm)", "Medium (9-16 mm)", "Small (<9 mm)")
names(size_labs) <- c("Max_body_size_abbrev_Large", "Max_body_size_abbrev_Medium", "Max_body_size_abbrev_Small")

plot4f <- ggplot(size, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = size_labs), nrow = 1) + 
  scale_colour_discrete(labels = c(">16 mm", "9-16 mm", "<9 mm"), name = "Body size") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Body size") + 
  xlab("") + 
  ggtitle("i")
plot4f

# Respiration
respiration_labs <- c("Gills", "Spiracle", "Tegument")
names(respiration_labs) <- c("Resp_abbrev_Gills", "Resp_abbrev_Plastron_spiracle", "Resp_abbrev_Tegument")

plot4g <- ggplot(respiration, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = respiration_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Gills", "Spiracle", "Tegument"), name = "Respiration mode") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Respiration") + 
  xlab("") + 
  ggtitle("d")
plot4g

# Rheophily
rheophily_labs <- c("Slow flow", "Mixed flow", "Fast flow")
names(rheophily_labs) <- c("Rheophily_abbrev_depo", "Rheophily_abbrev_depo_eros", "Rheophily_abbrev_eros")

plot4h <- ggplot(rheophily, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = rheophily_labs), ncol = 6) + 
  scale_colour_discrete(labels = c("Slow", "Mixed", "Fast"), name = "Flow preference") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Rheophily") + 
  xlab("") + 
  ggtitle("j")
plot4h

# Thermal
thermal_labs <- c("Cold stenothermal (<5°C)", "Cold-cool eurythermal (0-15°C)", "Cool-warm eurythermal (5-30°C)",
                  "Warm eurythermal (15-30°C)", "Hot eurythermal (>30°C)")
names(thermal_labs) <- c("Thermal_pref_Cold.cool.eurythermal..0.15.C.", "Thermal_pref_Cold.stenothermal...5.C.",
                         "Thermal_pref_Cool.warm.eurythermal..5.30.C.", "Thermal_pref_Hot.euthermal...30.C.",         
                         "Thermal_pref_Warm.eurythermal..15.30.C.")

plot4i <- ggplot(thermal, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = thermal_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Cold stenothermal (<5°C)", "Cold-cool eurythermal (0-15°C)", "Cool-warm eurythermal (5-30°C)",
                                   "Warm eurythermal (15-30°C)", "Hot eurythermal (>30°C)"), name = "Thermal preference") +
  scale_y_continuous(limits = c(0,1)) +
  #guides(colour = guide_legend(nrow = 2)) + 
  ylab("Thermal preference") + 
  xlab("Site age (years since deglaciation)") + 
  ggtitle("f")
plot4i

# Habit
habit_labs <- c("Burrower", "Climber", "Clinger", "Crawler", "Sprawler", "Swimmer")
names(habit_labs) <- c("Habit_prim_Burrower", "Habit_prim_Climber", "Habit_prim_Clinger",
                       "Habit_prim_Crawler", "Habit_prim_Sprawler", "Habit_prim_Swimmer" )

plot4j <- ggplot(habit, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = habit_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Burrower", "Climber", "Clinger", "Crawler", "Sprawler", "Swimmer"), name = "Habit") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Habit") + 
  xlab("") + 
  ggtitle("e")
plot4j

# Feeding guild
food_labs <- c("Filterers", "Gatherers", "Grazers", "Parasites", "Predators", "Shredders")
names(food_labs) <- c("Feed_prim_abbrev_CF", "Feed_prim_abbrev_CG", 
                      "Feed_prim_abbrev_HB", "Feed_prim_abbrev_PA",
                      "Feed_prim_abbrev_PR", "Feed_prim_abbrev_SH")

plot4k <- ggplot(food, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5), legend.box.background = element_rect(colour = NA, fill = NA),
        legend.key = element_rect(fill = NA)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = food_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Filterers", "Gatherers", "Grazers", "Parasites", "Predators", "Shredders"), name = "Feeding") +
  #scale_y_continuous(limits = c(0,1)) +
  ylab("Feeding guild") + 
  xlab("Site age (years since deglaciation)") + 
  ggtitle("k")
plot4k

plot_grid(plot4a, plot4b, plot4c, plot4d, plot4e, plot4f, plot4g, plot4h, 
          plot4j, plot4k, plot4i, NULL, ncol = 2, align = "hv")

## Fourth corner with environment, community and trait data

# Create the L table
test <- data.frame(Stream = str_sub(rownames(macro_sweep_clean_dframe1), end = -2), tax.mat.s)
tax.tab <- aggregate(. ~ Stream, data = test, FUN = "mean")
l_table <- data.frame(tax.tab[order(tax.tab$Stream),])
l_table$Stream <- NULL

# Physical variables (i.e., remove spot sample chemistry data)
chem_sweep <- dplyr::select(chem_sweep_dframe, -Turbidity.NTU., -Temp.C., 
                     -EC.umhos.cm., -pH, -Alkalinity.mgL., -Color.Pt., -Pctpool,
                     -Orientation.deg., -Total.P.ugL., -totalN.ugL., 
                     -NO2NO3.ugL., -SumN.ugL., -CBOM.mg., -X.wood, 
                     -totdolly.CPUE., -totfish.CPUE.)

# Create the R table
r_table <- chem_sweep[order(chem_sweep$Stream),]
r_table$Stream <- NULL

# Remove the secondary traits
traits_trimmed <- dplyr::select(data.frame(traits.s), -Habit_sec_Swimmer, -Feed_mode_sec_CG, 
                         -Feed_mode_sec_CF, -Habit_sec_Clinger, -Emerge_season_2_Spring, 
                         -Habit_sec_Climber, -Emerge_season_2_Winter, -Habit_sec_Sprawler,
                         -Feed_mode_sec_HB, -Feed_mode_sec_SH, -Feed_mode_sec_PR,
                         -Habit_sec_Planktonic, -Habit_sec_Burrower, -Feed_mode_sec_PA,
                         -Emerge_season_2_Fall, -Emerge_season_2_Summer)

# Remove traits that are invariant across sites
q_table <- data.frame(traits_trimmed)

# RLQ 
L_macro_sweep <- dudi.coa(l_table, scannf = F)
R_macro_sweep <- dudi.hillsmith(r_table, row.w = L_macro_sweep$lw, scannf = F)
Q_macro_sweep <- dudi.pca(q_table, row.w = L_macro_sweep$cw, scannf = F)
rlq_macro_sweep <- rlq(R_macro_sweep, L_macro_sweep, Q_macro_sweep, scannf = F)

summary(rlq_macro_sweep)

# Fourth corner
fourthc_macros <- fourthcorner(tabR = r_table, tabL = l_table, tabQ = q_table, 
                           nrepet = 9999,  p.adjust.method.G = "fdr", 
                           p.adjust.method.D = "fdr", modeltype = 6)
fourth_res <- summary(fourthc_macros)
plot(fourthc_macros, stat = "D")

# Combining the two methods 
test_rlq_macro_sweep <- randtest(rlq_macro_sweep, modeltype = 6, nrepet = 9999)
test_rlq_macro_sweep
plot(test_rlq_macro_sweep)

# Plot any significant individual relationships between individual traits and envs
plot(fourthc_macros, x.rlq = rlq_macro_sweep, alpha = 0.05, stat = "D2", 
     type = "biplot")

# Links between axes and other variables 
test_Qaxes_comb_macro_sweep <- fourthcorner.rlq(rlq_macro_sweep, modeltype = 6,
                                                typetest = "Q.axes", nrepet = 9999,
                                                p.adjust.method.G = "fdr", 
                                                p.adjust.method.D = "fdr")
print(test_Qaxes_comb_macro_sweep)

test_Raxes_comb_macro_sweep <- fourthcorner.rlq(rlq_macro_sweep, modeltype = 6,
                                                typetest = "R.axes", nrepet = 9999,
                                                p.adjust.method.G = "fdr", 
                                                p.adjust.method.D = "fdr")
print(test_Raxes_comb_macro_sweep)


### Meiofauna

## dbFD - multivariate functional diversity calculation

# Taxonomic data
tax.mat.s <- as.matrix(meio_taxa_sweep_dframe[,order(colnames(meio_taxa_sweep_dframe))]) # Loading site x taxa matrix
traits.s <- as.matrix(meio_sweep_fuzzy_traits_dframe[order(rownames(meio_sweep_fuzzy_traits_dframe)),]) # Loading taxa x trait matrix

sub_tax.mat.s <- tax.mat.s[-c(30,33),-c(2,3,13,19)]
sub_traits.s <- traits.s[-c(2,3,13,19),]

## Mean trait method ## 

# Run the mean trait analysis
meioFD <- dbFD(sub_traits.s, sub_tax.mat.s)

# Store the data in a new data frame
meio_sweep_FD_dframe <- data.frame(stream = str_sub(rownames(sub_tax.mat.s), end = -2),
                                    FRic = meioFD$FRic, FEve = meioFD$FEve, 
                                    FDiv = meioFD$FDiv, FDis = meioFD$FDis)

# Change the stream names to match the chem data
meio_sweep_FD_dframe$Stream <- case_match(meio_sweep_FD_dframe$stream, 
                                          "stonefly" ~ "SFC", 
                                          "Gull" ~ "Gull",
                                          "WPC" ~ "WPC", 
                                          "BBN" ~ "BBN", 
                                          "NUN" ~ "Nun", 
                                          "CarE" ~ "CarE", 
                                          "Tarr" ~ "Tar",
                                          "reid" ~ "Reid",
                                          "vivid" ~ "Viv",
                                          "ice valley" ~ "Ice", 
                                          "oyster" ~ "Oys",
                                          "tyndall" ~ "Tyn",
                                          "n fingers" ~ "NFS",
                                          "BBS" ~ "BBS", 
                                          "RPC" ~ "Rush")

# Add in the chemical data to help visualise the results
meiotest <- left_join(meio_sweep_FD_dframe, chem_sweep_dframe, by = "Stream")

plot5a <- ggplot() + 
  stat_summary(aes(x = Age.y., y = FRic, fill = Stream), pch = 21, size = 1, data = meiotest) + 
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
  ylab("Functional richness (FRic)") + 
  xlab("") + 
  ggtitle("(a)")
plot5a


plot5b <- ggplot() + 
  stat_summary(aes(x = Age.y., y = FEve, fill = Stream), pch = 21, size = 1, data = meiotest) + 
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
  ylab("Functional evenness (FEve)") + 
  xlab("") + 
  ggtitle("(b)")
plot5b

plot5c <- ggplot() + 
  stat_summary(aes(x = Age.y., y = FDiv, fill = Stream), pch = 21, size = 1, 
               data = meiotest) + 
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
  ylab("Functional diversity (FDiv)") + 
  xlab("") +   
  ggtitle("(c)")
plot5c

plot5d <- ggplot() + 
  stat_summary(aes(x = Age.y., y = FDis, fill = Stream), pch = 21, size = 1, data = meiotest) + 
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
  ylab("Functional dispersion (FDis)") + 
  xlab("") + 
  ggtitle("(d)")
plot5d

grid.arrange(plot5a, plot5b, plot5c, plot5d, nrow = 2, bottom = "Site age (years since deglaciation)")


## CWM analyses

# Fuzzy coded traits

# Calculate community weighted means across sites for analysis
meiotrial <- functcomp(sub_traits.s, sub_tax.mat.s)

# Trim the dataframe to traits with values 
meiotrial_sub <- meiotrial[, !sapply(meiotrial, is.character)]

# Extract the streams
meiotrial_sub$stream <- str_sub(rownames(meiotrial_sub), end = -2)

# Change the stream names to match the chem data
meiotrial_sub$Stream <- case_match(meiotrial_sub$stream,
                                   "stonefly" ~ "SFC", 
                                   "Gull" ~ "Gull",
                                   "WPC" ~ "WPC", 
                                   "BBN" ~ "BBN", 
                                   "NUN" ~ "Nun", 
                                   "CarE" ~ "CarE", 
                                   "Tarr" ~ "Tar",
                                   "reid" ~ "Reid",
                                   "vivid" ~ "Viv",
                                   "ice valley" ~ "Ice", 
                                   "oyster" ~ "Oys",
                                   "tyndall" ~ "Tyn",
                                   "n fingers" ~ "NFS",
                                   "BBS" ~ "BBS", 
                                   "RPC" ~ "Rush")

# Bind the chemical data
meio_cwm_sweep <- left_join(meiotrial_sub, chem_sweep_dframe, by = "Stream")

meio_bb <- meio_cwm_sweep %>% dplyr::select(order(colnames(meio_cwm_sweep)))

# Traits
body <- melt(dplyr::select(meio_bb, Age.y., Bodyform_cylindrical:Bodyform_sphaerical), "Age.y.")
fecundity <- melt(dplyr::select(meio_bb, Age.y., Fecundity_.100:Fecundity_.3000), "Age.y.")
feeding <- melt(dplyr::select(meio_bb, Age.y., Feedinghabits_absorber:Feedinghabits_shredder), "Age.y.")
food <- melt(dplyr::select(meio_bb, Age.y., Food_dead.plant.1mm:Food_livingmicrophytes), "Age.y.")
locomotion <- melt(dplyr::select(meio_bb, Age.y., Locomotionsubstrate.relation_burrower.epibenthic.:Locomotionsubstraterelation_temporaryattached), "Age.y.")
size <- melt(dplyr::select(meio_bb, Age.y., Maximal.potential.size_.0.5.1.cm:Maximalpotentialsize_.8cm), "Age.y.")
cycles <- melt(dplyr::select(meio_bb, Age.y., Nbcyclesan_1:Nbcyclesan_less1), "Age.y.")
reproduction <- melt(dplyr::select(meio_bb, Age.y., ReproductionTechnique_clutches.cementedorfixed:ReproductionTechnique_parthenogenesis), "Age.y.")
resistance <- melt(dplyr::select(meio_bb, Age.y., Resistanceforms_cocoons:Resistanceforms_none), "Age.y.")
respiration <- melt(dplyr::select(meio_bb, Age.y., Respiration_gill:Respiration_tegument), "Age.y.")

## Extract GAM results
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = body))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = fecundity))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = feeding))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = locomotion))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = size))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = cycles))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = reproduction))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = resistance))
summary(gam(value ~ s(Age.y., bs = "cs", fx = TRUE, k = -1, by = variable), data = respiration))


# Adult flight strength
body_labs <- c("Cylindrical", "Flattened", "Geometric", "Spherical")
names(body_labs) <- c("Bodyform_cylindrical", "Bodyform_flattened",
                      "Bodyform_geometric", "Bodyform_sphaerical")

plot6a<- ggplot(body, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = flight_labs), nrow = 1) + 
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_discrete(labels = c("Cylindrical", "Flattened", "Geometric", "Spherical"), name = "Body") + 
  ylab("Body form") + 
  xlab("") + 
  ggtitle("a")
plot6a

# Female dispersal ability
fecundity_labs <- c("<100", "100-1000", "1000-3000", ">3000")  
names(fecundity_labs) <- c("Fecundity_.100", "Fecundity_.100.1000", "Fecundity_.1000.3000", "Fecundity_.3000")

plot6b<- ggplot(fecundity, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = dispersal_labs), nrow = 1) + 
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_discrete(labels = c("<100", "100-1000", "1000-3000", ">3000"), name = "Fecundity") +
  ylab("Fecundity") + 
  xlab("") + 
  ggtitle("f")
plot6b

# Feeding
feeding_labs <- c("Absorber", "Deposit feeder", "Filterer", "Parasite", 
                  "Piercer", "Engulfer", "Scraper", "Shredder")
names(feeding_labs) <- c("Feedinghabits_absorber", "Feedinghabits_deposit.feeder",                     
                         "Feedinghabits_filter.feeder", "Feedinghabits_parasite",                           
                         "Feedinghabits_piercer.plantsoranimals.", 
                         "Feedinghabits_predator.carver.engulfer.swallower.",
                         "Feedinghabits_scraper", "Feedinghabits_shredder")

plot6c<- ggplot(feeding, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = voltinism_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Absorber", "Deposit feeder", "Filterer", "Parasite", 
                                   "Piercer", "Engulfer", "Scraper", "Shredder"), name = "Feeding") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Feeding guild") + 
  xlab("") + 
  ggtitle("d")
plot6c

# Food
#food_labs <- c("Dead plants", "Dead animals", "Detritus", "Microorganisms", 
#              "Macroinvertebrates", "Macrophytes", "Microinvertebrates", "Microphytes")
#names(food_labs) <- c("Food_dead.plant.1mm", "Food_deadanimal.1mm", "Food_detritus.1mm",                       
#                      "Food_finesedimentmicroorganisms", "Food_livingmacroinvertebratesvertebrates", 
#                      "Food_livingmacrophytes", "Food_livingmicroinvertebrates", "Food_livingmicrophytes")
#
#plot6d <- ggplot(food, aes(y = value, x = Age.y.)) + 
#  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
#  geom_point() + 
#  theme_bw() +  
#  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
#        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
#        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
#        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
#  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = synchrony_labs), nrow = 1) + 
#  scale_colour_discrete(labels = c("Dead plants", "Dead animals", "Detritus", "Microorganisms", 
#                                   "Macroinvertebrates", "Macrophytes", "Microinvertebrates", "Microphytes"), name = "Food") +
#  scale_y_continuous(limits = c(0,1)) +
#  ylab("") + 
#  xlab("") + 
#  ggtitle("Food")
#plot6d

# Locomotion
locomotion_labs <- c("Burrower", "Interstitial", "Crawler", "Swimmer", "Temporary attached")
names(locomotion_labs) <- c("Locomotionsubstrate.relation_burrower.epibenthic.", 
                            "Locomotionsubstrate.relation_interstitial.endobenthic.",
                            "Locomotionsubstraterelation_crawler", "Locomotionsubstraterelation_swimmer",
                            "Locomotionsubstraterelation_temporaryattached")

plot6e <- ggplot(locomotion, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = emergence_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Burrower", "Interstitial", "Crawler", "Swimmer", "Temporary attached"), name = "Locomotion") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Habit") + 
  xlab("Stream age (time since deglaciation)") + 
  ggtitle("i")
plot6e

# Size
size_labs <- c("0.5-1 mm", "0.25-0.5 mm", "<0.25 mm", "1-2 mm", "2-4 mm",
               "4-8 mm", ">8 mm") 
names(size_labs) <- c("Maximal.potential.size_.0.5.1.cm", 
                      "Maximalpotentialsize_.0.25..5.cm", 
                      "Maximalpotentialsize_.0.25.cm",
                      "Maximalpotentialsize_.1.2.cm", 
                      "Maximalpotentialsize_.2.4.cm", 
                      "Maximalpotentialsize_.4.8.cm",
                      "Maximalpotentialsize_.8cm")

## Order them by actual size classes

plot6f <- ggplot(size, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = size_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("0.5-1 mm", "0.25-0.5 mm", "<0.25 mm", "1-2 mm", "2-4 mm",
                                   "4-8 mm", ">8 mm"), name = "Body size") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Body size") + 
  xlab("") + 
  ggtitle("b")
plot6f

# Cycles
cycles_labs <- c("1", ">1", "<1")
names(cycles_labs) <- c("Nbcyclesan_1", "Nbcyclesan_greater1", "Nbcyclesan_less1")

plot6g <- ggplot(cycles, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = respiration_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("1", ">1", "<1"), name = "Number of generations") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Generations") + 
  xlab("") + 
  ggtitle("g")
plot6g

# Reproduction
reproduction_labs <- c("Attached clutches", "Free clutches", "Vegetation clutches", "Terrestrial clutches", 
                       "Attached individual", "Free individual", "Ovoviviparity", "Parthenogenesis")
names(reproduction_labs) <- c("ReproductionTechnique_clutches.cementedorfixed", "ReproductionTechnique_clutches.free",           
                              "ReproductionTechnique_clutches.invegetation", "ReproductionTechnique_clutches.terrestrial",    
                              "ReproductionTechnique_isolatedeggs.cemented", "ReproductionTechnique_isolatedeggs.free",       
                              "ReproductionTechnique_ovoviviparity", "ReproductionTechnique_parthenogenesis")

plot6h <- ggplot(reproduction, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = rheophily_labs), ncol = 6) + 
  scale_colour_discrete(labels = c("Attached clutches", "Free clutches",
                                   "Vegetation clutches", "Terrestrial clutches", 
                                   "Attached individual", "Free individual",
                                   "Ovoviviparity", "Parthenogenesis"), 
                        name = "Oviposition") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Oviposition") + 
  xlab("") + 
  ggtitle("h")
plot6h

# Resistance
resistance_labs <- c("Cocoons", "Diapause", "Statoblasts", "None")
names(resistance_labs) <- c("Resistanceforms_cocoons", 
                            "Resistanceforms_diapauseordormancy",
                            "Resistanceforms_eggs.statoblasts",  
                            "Resistanceforms_none")

plot6i <- ggplot(resistance, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = rheophily_labs), ncol = 6) + 
  scale_colour_discrete(labels = c("Cocoons", "Diapause", "Statoblasts", "None"), name = "Resistance") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Resistance forms") + 
  xlab("") + 
  ggtitle("c")
plot6i

# Respiration
respiration_labs <- c("Gills", "Spiracle", "Tegument")
names(respiration_labs) <- c("Respiration_gill", "Respiration_spiracle.aerial.",
                             "Respiration_tegument")

plot6j <- ggplot(respiration, aes(y = value, x = Age.y.)) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        legend.justification = c(0,0.5)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = rheophily_labs), ncol = 6) + 
  scale_colour_discrete(labels = c("Gills", "Spiracle", "Tegument"), name = "Respiration apparatus") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Respiration") + 
  xlab("Stream age (time since deglaciation)") + 
  ggtitle("e")
plot6j

# Multiplot
plot_grid(plot6a, plot6b, plot6f, plot6g, plot6i,
          plot6h, plot6c, plot6e, plot6j, NULL, ncol = 2, align = "hv")


## Fourth corner with environment, community and trait data

test_meio <- data.frame(Stream = str_sub(rownames(meio_taxa_sweep_dframe), end = -2), tax.mat.s)
meio.tax.tab <- aggregate(. ~ Stream, data = test_meio, FUN = "mean")
l_table <- data.frame(meio.tax.tab[order(meio.tax.tab$Stream),])
l_table$Stream <- NULL

# Trimmed environmental variables
#chem_sweep <- dplyr::select(chem_sweep_dframe, -totdolly.CPUE., -Area.km2., 
#                     -elevation.m., -Orientation.deg., -Streamlength.km., 
#                     -Gradient..., -Order, -Orientation.deg., -NO2NO3.ugL., 
#                     -SumN.ugL., -Pctpool, -Pfan_low)

# Physical variables
chem_sweep <- dplyr::select(chem_sweep_dframe, -Turbidity.NTU., -Temp.C., 
                     -EC.umhos.cm., -pH, -Alkalinity.mgL., -Color.Pt., -Pctpool,
                     -Orientation.deg., -Total.P.ugL., -totalN.ugL., 
                     -NO2NO3.ugL., -SumN.ugL., -CBOM.mg., -X.wood,
                     -totdolly.CPUE., -totfish.CPUE.)

r_table <- chem_sweep[order(chem_sweep$Stream),]
r_table$Stream <- NULL
r_table <- r_table[-16,]

# Remove traits that are invariant across sites
#traits_clean <- dplyr::select(data.frame(traits.s), -Resistanceforms_cocoons,
#                       -Bodyform_flattened, -Bodyform_geometric,
#                       -Maximalpotentialsize_.2.4.cm,
#                       -Maximalpotentialsize_.4.8.cm, -Maximalpotentialsize_.8cm,
#                       -Nbcyclesan_less1, -Resistanceforms_housingsagainstdesiccation,
#                       -Feedinghabits_shredder, -Feedinghabits_piercer.plantsoranimals.)

q_table <- data.frame(traits_clean)
q_table <- data.frame(traits.s)

# RLQ 
L_meio_sweep <- dudi.coa(l_table, scannf = F)
R_meio_sweep <- dudi.hillsmith(r_table, row.w = L_meio_sweep$lw, scannf = F)
Q_meio_sweep <- dudi.pca(q_table, row.w = L_meio_sweep$cw, scannf = F)
rlq_meio_sweep <- rlq(R_meio_sweep, L_meio_sweep, Q_meio_sweep, scannf = F)

summary(rlq_meio_sweep)

# Extract the data for plotting
#rlq_meio_sweep$

# Arrows for the stream environmental variables
rlq_meio_sweep$l1

# Arrows for the traits
rlq_meio_sweep$c1

# Labels for species 
rlq_meio_sweep$lQ


### Plot in ggplot ### 


# Fourth corner
fourthc_meios <- fourthcorner(tabR = r_table, tabL = l_table, tabQ = q_table, 
                               nrepet = 99999,  p.adjust.method.G = "fdr", 
                               p.adjust.method.D = "fdr", modeltype = 6)
fourth_res <- summary(fourthc_meios)
plot(fourthc_meios, stat = "D")

# Combining the two methods 
test_rlq_meio_sweep <- randtest(rlq_meio_sweep, modeltype = 6, nrepet = 9999)
test_rlq_meio_sweep
plot(test_rlq_meio_sweep)

# Plot any significant individual relationships between individual traits and envs
plot(fourthc_meios, x.rlq = rlq_meio_sweep, alpha = 0.05, stat = "D2", 
     type = "biplot")

# Links between axes and other variables 
test_Qaxes_comb_meio_sweep <- fourthcorner.rlq(rlq_meio_sweep, modeltype = 6,
                                                typetest = "Q.axes", nrepet = 9999,
                                                p.adjust.method.G = "fdr", 
                                                p.adjust.method.D = "fdr")
print(test_Qaxes_comb_meio_sweep)

test_Raxes_comb_meio_sweep <- fourthcorner.rlq(rlq_meio_sweep, modeltype = 6,
                                                typetest = "R.axes", nrepet = 9999,
                                                p.adjust.method.G = "fdr", 
                                                p.adjust.method.D = "fdr")
print(test_Raxes_comb_meio_sweep)



#### WPC long term dataset

### Macroinvertebrates

## dbFD - multivariate functional diversity calculation

# Taxonomic data
macro_WPC_clean_dframe1 <- dplyr::select(macro_WPC_clean_dframe, -Nyctiophylax, -Siphlonuridae)
macro_WPC_traits_dframe1 <- macro_WPC_traits_dframe[-c(14,28),]
tax.mat.s <- as.matrix(macro_WPC_clean_dframe1[,order(colnames(macro_WPC_clean_dframe1))]) # Loading site x taxa matrix
traits.s <- as.matrix(macro_WPC_traits_dframe1[order(rownames(macro_WPC_traits_dframe1)),]) # Loading taxa x trait matrix

## Mean trait method 

# Run the mean trait analysis
FD <- dbFD(traits.s, tax.mat.s)

# Store the data in a new data frame
macro_WPC_FD_dframe <- data.frame(year = str_sub(rownames(macro_WPC_clean_dframe1), end = -1),
                                    FRic = FD$FRic, FEve = FD$FEve, 
                                    FDiv = FD$FDiv, FDis = FD$FDis)

# 1978 was too small to compute so add 0s in for visualisation
#macro_WPC_FD_dframe[1,2:4] <- 0

# Statistical analysis 

model2 <- glm(FRic ~ as.numeric(year), family = Gamma (link = "log"), data = macro_WPC_FD_dframe[-1,])
summary.lm(model2)

model2_preds <- predict(model2, newdata = data.frame(year = seq(1978,2017,0.1)), se.fit = T)
FRic_preds <- data.frame(year = seq(1978,2017,0.1), preds = model2_preds$fit, se = model2_preds$se.fit)

plot7a <- ggplot() + 
  geom_line(aes(x = year, y = exp(preds)), data = FRic_preds, linetype = "solid", colour = "black") + 
  geom_line(aes(x = year, y = exp(preds-se)), data = FRic_preds, linetype = "dashed", colour = "black") + 
  geom_line(aes(x = year, y = exp(preds+se)), data = FRic_preds, linetype = "dashed", colour = "black") + 
  stat_summary(aes(x = as.numeric(year), y = FRic), geom = "point", fill = "darkred",
               pch = 21, size = 4, data = macro_WPC_FD_dframe, inherit.aes = F) + 
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10, colour = "black"), 
        legend.position = "NA", 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  ylab("Functional richness (FRic)") + 
  #coord_cartesian(ylim = c(0,500)) + 
  xlab("") + 
  ggtitle("a")
plot7a

plot7b <- ggplot() + 
  stat_summary(aes(x = as.numeric(year), y = FEve), pch = 21, size = 4, geom = "point",
               fill = "darkred", data = macro_WPC_FD_dframe) + 
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 10), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 10, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  ylab("Functional evenness (FEve)") + 
  xlab("") + 
  ggtitle("(b)")
plot7b

# Statistical analysis
model3 <- glm(FDiv ~ poly(as.numeric(year), 2, raw = T), family = "Gamma" (link = "identity"), data = macro_WPC_FD_dframe)
summary.lm(model3)

model3_preds <- predict(model3, newdata = data.frame(year = seq(1978,2017,0.1)), se.fit = T)
FDiv_preds <- data.frame(year = seq(1978,2017,0.1), preds = model3_preds$fit, se = model3_preds$se.fit)


plot7c <- ggplot() + 
  geom_line(aes(x = year, y = preds), data = FDiv_preds, linetype = "solid", colour = "black") + 
  geom_line(aes(x = year, y = preds-se), data = FDiv_preds, linetype = "dashed", colour = "black") + 
  geom_line(aes(x = year, y = preds+se), data = FDiv_preds, linetype = "dashed", colour = "black") + 
  stat_summary(aes(x = as.numeric(year), y = FDiv), pch = 21, size = 4, geom = "point",
               fill = "darkred", data = macro_WPC_FD_dframe) + 
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 10), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 10, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  ylab("Functional divergence (FDiv)") + 
  xlab("Year") +   
  ggtitle("b")
plot7c

plot7d <- ggplot() + 
  #stat_smooth(aes(x = as.numeric(year), y = FDis), method = "lm", 
  #            data = macro_WPC_FD_dframe, se = TRUE, fill = "lightgrey", 
  #            formula = y ~ poly(x, 2, raw = TRUE), 
  #            colour = "black", alpha = 1) + 
  stat_summary(aes(x = as.numeric(year), y = FDis), size = 4, pch = 21, geom = "point",
               fill = "darkred", data = macro_WPC_FD_dframe) + 
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 10), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 10, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  ylab("Functional dispersion (FDis)") + 
  xlab("") + 
  ggtitle("(d)")
plot7d

grid.arrange(plot7a, plot7c, nrow = 2)


## CWM over time in WPC

# Calculate community weighted means across sites for analysis
trial <- functcomp(traits.s, tax.mat.s)

# Convert some columns to numeric
trial[,1:2] <- sapply(trial[, 1:2], as.numeric)

# Trim the dataframe to traits with values 
trial_sub <- trial[, !sapply(trial, is.character)]

# Extract the years
trial_sub$year <- rownames(trial_sub)

bb <- trial_sub %>% dplyr::select(order(colnames(trial_sub)))

# Traits
#flight <- melt(dplyr::select(bb, year, AdultFlyingStrength_abbrev_Strong:AdultFlyingStrength_abbrev_Weak), "year")
#dispersal <- melt(dplyr::select(bb, year, Female_disp_abbrev_High:Female_disp_abbrev_Low), "year")
#voltinism <- melt(dplyr::select(bb, year, Voltinism_abbrev_Bi_multivoltine:Voltinism_abbrev_Univoltine), "year")
#synchrony <- melt(dplyr::select(bb, year, Emerge_synch_abbrev_Poorly:Emerge_synch_abbrev_Well), "year")
emergence <- melt(dplyr::select(bb, year, Emerge_season_1_Fall:Emerge_season_1_Winter), "year")
size <- melt(dplyr::select(bb, year, Max_body_size_abbrev_Large:Max_body_size_abbrev_Small), "year")
respiration <- melt(dplyr::select(bb, year, Resp_abbrev_Gills:Resp_abbrev_Tegument), "year")
rheophily <- melt(dplyr::select(bb, year, Rheophily_abbrev_depo: Rheophily_abbrev_eros), "year")
thermal <- melt(dplyr::select(bb, year, Thermal_pref_Cold.cool.eurythermal..0.15.C.:Thermal_pref_Warm.eurythermal..15.30.C.), "year")
habit <- melt(dplyr::select(bb, year, Habit_prim_Burrower:Habit_prim_Swimmer), "year")
food <- melt(dplyr::select(bb, year, Feed_prim_abbrev_CF:Feed_prim_abbrev_SH), "year")

# Extract GAM results
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = emergence))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = size))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = rheophily))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = thermal))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = habit))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = food))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = respiration))

# Adult flight strength
#flight_labs <- c("Strong flying strength", "Weak flying strength")
#names(flight_labs) <- c("AdultFlyingStrength_abbrev_Strong", "AdultFlyingStrength_abbrev_Weak")
#
#plot8a<- ggplot(flight, aes(y = value, x = year)) + 
#  geom_smooth(aes(colour = variable), fill = "lightgrey", alpha = 0) + 
#  geom_point() + 
#  theme_bw() +  
#  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
#        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
#        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
#        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
#  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = flight_labs), nrow = 1) + 
#  scale_y_continuous(limits = c(0,1)) +
#  scale_colour_discrete(labels = c("Strong", "Weak"), name = "Flight") + 
#  ylab("") + 
#  xlab("") + 
#  ggtitle("Flying strength")
#plot8a

# Emergence season
emergence_labs <- c("Fall emergence", "Spring emergence", "Summer emergence", "Winter emergence")
names(emergence_labs) <- c("Emerge_season_1_Fall", "Emerge_season_1_Spring", 
                           "Emerge_season_1_Summer", "Emerge_season_1_Winter")

plot8e <- ggplot(emergence, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = emergence_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Fall", "Spring", "Summer", "Winter"), name = "Emergence season") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Emergence season") + 
  xlab("") + 
  ggtitle("a")
plot8e

# Body size
size_labs <- c("Large (>16 mm)", "Medium (9-16 mm)", "Small (<9 mm)")
names(size_labs) <- c("Max_body_size_abbrev_Large", "Max_body_size_abbrev_Medium", "Max_body_size_abbrev_Small")

plot8f <- ggplot(size, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = size_labs), nrow = 1) + 
  scale_colour_discrete(labels = c(">16 mm", "9-16 mm", "<9 mm"), name = "Body size") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Body size") + 
  xlab("") + 
  ggtitle("d")
plot8f

# Respiration
respiration_labs <- c("Gills", "Spiracle", "Tegument")
names(respiration_labs) <- c("Resp_abbrev_Gills", "Resp_abbrev_Plastron_spiracle", "Resp_abbrev_Tegument")

plot8g <- ggplot(respiration, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = respiration_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Gills", "Spiracle", "Tegument"), name = "Respiration mode") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Respiration mode") + 
  xlab("") + 
  ggtitle("b")
plot8g

# Rheophily
rheophily_labs <- c("Slow flow", "Mixed flow", "Fast flow")
names(rheophily_labs) <- c("Rheophily_abbrev_depo", "Rheophily_abbrev_depo_eros", "Rheophily_abbrev_eros")

plot8h <- ggplot(rheophily, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", alpha = 0, method = "gam") + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = rheophily_labs), ncol = 6) + 
  scale_colour_discrete(labels = c("Slow", "Mixed", "Fast"), name = "Flow preference") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Rheophily") + 
  xlab("") + 
  ggtitle("e")
plot8h

# Thermal
thermal_labs <- c("Cold stenothermal (<5°C)", "Cold-cool eurythermal (0-15°C)", "Cool-warm eurythermal (5-30°C)",
                  "Warm eurythermal (15-30°C)", "Hot eurythermal (>30°C)")
names(thermal_labs) <- c("Thermal_pref_Cold.cool.eurythermal..0.15.C.", "Thermal_pref_Cold.stenothermal...5.C.",
                         "Thermal_pref_Cool.warm.eurythermal..5.30.C.", "Thermal_pref_Hot.euthermal...30.C.",         
                         "Thermal_pref_Warm.eurythermal..15.30.C.")

plot8i <- ggplot(thermal, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = thermal_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Cold stenothermal (<5°C)", "Cold-cool eurythermal (0-15°C)", "Cool-warm eurythermal (5-30°C)",
                                   "Warm eurythermal (15-30°C)", "Hot eurythermal (>30°C)"), name = "Thermal preference") +
  scale_y_continuous(limits = c(0,1)) +
  #guides(colour = guide_legend(nrow = 2)) + 
  ylab("Thermal preference") + 
  xlab("") + 
  ggtitle("")
plot8i

# Habit
habit_labs <- c("Burrower", "Climber", "Clinger", "Crawler", "Sprawler", "Swimmer")
names(habit_labs) <- c("Habit_prim_Burrower", "Habit_prim_Climber", "Habit_prim_Clinger",
                       "Habit_prim_Crawler", "Habit_prim_Sprawler", "Habit_prim_Swimmer" )

plot8j <- ggplot(habit, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", alpha = 0, method = "gam") + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = habit_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Burrower", "Climber", "Clinger", "Crawler", "Sprawler", "Swimmer"), name = "Habit") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Habit") + 
  xlab("Year") + 
  ggtitle("c")
plot8j

# Feeding guild
food_labs <- c("Filterers", "Gatherers", "Grazers", "Parasites", "Predators", "Shredders")
names(food_labs) <- c("Feed_prim_abbrev_CF", "Feed_prim_abbrev_CG", 
                      "Feed_prim_abbrev_HB", "Feed_prim_abbrev_PA",
                      "Feed_prim_abbrev_PR", "Feed_prim_abbrev_SH")

plot8k <- ggplot(food, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1975, 2020, 5), limits = c(1975, 2020)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = food_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Filterers", "Gatherers", "Grazers", "Parasites", "Predators", "Shredders"), name = "Feeding") +
  #scale_y_continuous(limits = c(0,1)) +
  ylab("Feeding guild") + 
  xlab("Year") + 
  ggtitle("f")
plot8k

# Multiplot of all traits
plot_grid(plot8e, plot8f, plot8g, plot8h, plot8j, plot8k, ncol = 2, align = "hv")


### Meiofauna

## dbFD - multivariate functional diversity calculation

# Taxonomic data
tax.mat.s <- as.matrix(meio_taxa_WPC_dframe[,order(colnames(meio_taxa_WPC_dframe))]) # Loading site x taxa matrix
traits.s <- as.matrix(meio_WPC_fuzzy_traits_dframe[order(rownames(meio_WPC_fuzzy_traits_dframe)),]) # Loading taxa x trait matrix

sub_tax.mat.s <- tax.mat.s[,-c(5)]
sub_traits.s <- traits.s[-c(5),]

## Mean trait method ## 

# Run the mean trait analysis
meioFD <- dbFD(sub_traits.s, sub_tax.mat.s)

# Store the data in a new data frame
meio_WPC_FD_dframe <- data.frame(year = rownames(sub_tax.mat.s),
                                 FRic = meioFD$FRic, FEve = meioFD$FEve, 
                                 FDiv = meioFD$FDiv, FDis = meioFD$FDis)

plot5a <- ggplot() + 
  stat_summary(aes(x = as.numeric(year), y = FRic), pch = 21, size = 0.5,
               fill = "darkred", data = meio_WPC_FD_dframe) + 
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 10), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 10, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  ylab("Functional richness (FRic)") + 
  xlab("") + 
  ggtitle("a")
plot5a


plot5b <- ggplot() + 
  stat_summary(aes(x = as.numeric(year), y = FEve), pch = 21, size = 1,
               fill = "darkred", data = meio_WPC_FD_dframe) + theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), ) + 
  ylab("Functional evenness (FEve)") + 
  xlab("") + 
  ggtitle("(b)")
plot5b

plot5c <- ggplot() + 
  stat_smooth(method = "lm", aes(x = as.numeric(year), y = FDiv),
              data = meio_WPC_FD_dframe, se = TRUE, 
              fill = "lightgrey", 
              formula = y ~ poly(x, 2, raw = TRUE),
              colour = "black", alpha = 1) + 
  stat_summary(aes(x = as.numeric(year), y = FDiv), pch = 21, size = 0.5,
               fill = "darkred", data = meio_WPC_FD_dframe) + 
  theme_bw() +  
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 10), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 10, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  ylab("Functional divergence (FDiv)") + 
  xlab("Year") +   
  ggtitle("b")
plot5c

plot5d <- ggplot() + 
  stat_summary(aes(x = as.numeric(year), y = FDis), pch = 21, size = 0.5,
               fill = "darkred", data = meio_WPC_FD_dframe) +
  theme_bw() +  
  scale_colour_viridis_d(option = "turbo") + 
  scale_fill_viridis_d(option = "turbo") + 
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 12), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "NA", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=10), plot.margin = unit(c(1, 2.5, 1, 1), "mm")) + 
  ylab("Functional dispersion (FDis)") + 
  xlab("") + 
  ggtitle("(d)")
plot5d

plot_grid(plot5a, plot5c, nrow = 2, align = "hv")


## CWM analyses

# Fuzzy coded traits

# Calculate community weighted means across sites for analysis
meiotrial <- functcomp(sub_traits.s, sub_tax.mat.s)

# Trim the dataframe to traits with values 
meiotrial_sub <- meiotrial[, !sapply(meiotrial, is.character)]

meiotrial_sub$year <- rownames(meiotrial_sub) 

meio_bb <- meiotrial_sub %>% dplyr::select(order(colnames(meiotrial_sub)))

# Traits
body <- melt(dplyr::select(meio_bb, year, Bodyform_cylindrical:Bodyform_sphaerical), "year")
fecundity <- melt(dplyr::select(meio_bb, year, Fecundity_.100:Fecundity_.3000), "year")
feeding <- melt(dplyr::select(meio_bb, year, Feedinghabits_absorber:Feedinghabits_shredder), "year")
food <- melt(dplyr::select(meio_bb, year, Food_dead.plant.1mm:Food_livingmicrophytes), "year")
locomotion <- melt(dplyr::select(meio_bb, year, Locomotionsubstrate.relation_burrower.epibenthic.:Locomotionsubstraterelation_temporaryattached), "year")
size <- melt(dplyr::select(meio_bb, year, Maximal.potential.size_.0.5.1.cm:Maximalpotentialsize_.8cm), "year")
cycles <- melt(dplyr::select(meio_bb, year, Nbcyclesan_less1:Nbcyclesan_greater1), "year")
reproduction <- melt(dplyr::select(meio_bb, year, ReproductionTechnique_clutches.cementedorfixed:ReproductionTechnique_parthenogenesis), "year")
resistance <- melt(dplyr::select(meio_bb, year, Resistanceforms_cocoons:Resistanceforms_none), "year")
respiration <- melt(dplyr::select(meio_bb, year, Respiration_gill:Respiration_tegument), "year")

# Extract GAM results
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = body))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = fecundity))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = feeding))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = food))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = locomotion))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = size))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = cycles))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = reproduction))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = resistance))
summary(gam(value ~ s(as.numeric(year), bs = "cs", fx = TRUE, k = -1, by = variable), data = respiration))

# Adult flight strength
body_labs <- c("Cylindrical", "Flattened", "Geometric", "Spherical")
names(body_labs) <- c("Bodyform_cylindrical", "Bodyform_flattened",
                      "Bodyform_geometric", "Bodyform_sphaerical")

plot6a<- ggplot(body, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = flight_labs), nrow = 1) + 
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_discrete(labels = c("Cylindrical", "Flattened", "Geometric", "Spherical"), name = "Body") + 
  ylab("Body form") + 
  xlab("") + 
  ggtitle("a")
plot6a

# Female dispersal ability
fecundity_labs <- c("<100", "100-1000", "1000-3000", ">3000")  
names(fecundity_labs) <- c("Fecundity_.100", "Fecundity_.100.1000", "Fecundity_.1000.3000", "Fecundity_.3000")

plot6b<- ggplot(fecundity, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) + 
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = dispersal_labs), nrow = 1) + 
  scale_y_continuous(limits = c(0,1)) +
  scale_colour_discrete(labels = c("<100", "100-1000", "1000-3000", ">3000"), name = "Fecundity") +
  ylab("Fecundity") + 
  xlab("") + 
  ggtitle("f")
plot6b

# Feeding
feeding_labs <- c("Absorber", "Deposit feeder", "Filterer", "Parasite", 
                  "Piercer", "Engulfer", "Scraper", "Shredder")
names(feeding_labs) <- c("Feedinghabits_absorber", "Feedinghabits_deposit.feeder",                     
                         "Feedinghabits_filter.feeder", "Feedinghabits_parasite",                           
                         "Feedinghabits_piercer.plantsoranimals.", 
                         "Feedinghabits_predator.carver.engulfer.swallower.",
                         "Feedinghabits_scraper", "Feedinghabits_shredder")

plot6c<- ggplot(feeding, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = voltinism_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Absorber", "Deposit feeder", "Filterer", "Parasite", 
                                   "Piercer", "Engulfer", "Scraper", "Shredder"), name = "Feeding") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Feeding") + 
  xlab("") + 
  ggtitle("h")
plot6c

# Food
food_labs <- c("Dead plants", "Dead animals", "Detritus", "Microorganisms", 
              "Macroinvertebrates", "Macrophytes", "Microinvertebrates", "Microphytes")
names(food_labs) <- c("Food_dead.plant.1mm", "Food_deadanimal.1mm", "Food_detritus.1mm",                       
                      "Food_finesedimentmicroorganisms", "Food_livingmacroinvertebratesvertebrates", 
                      "Food_livingmacrophytes", "Food_livingmicroinvertebrates", "Food_livingmicrophytes")

plot6d <- ggplot(food, aes(y = value, x = as.numeric(year))) + 
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = synchrony_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Dead plants", "Dead animals", "Detritus", "Microorganisms", 
                                   "Macroinvertebrates", "Macrophytes", "Microinvertebrates", "Microphytes"), name = "Food") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Food") + 
  xlab("") + 
  ggtitle("d")
plot6d

# Locomotion
locomotion_labs <- c("Burrower", "Interstitial", "Crawler", "Swimmer", "Temporary attached")
names(locomotion_labs) <- c("Locomotionsubstrate.relation_burrower.epibenthic.", 
                            "Locomotionsubstrate.relation_interstitial.endobenthic.",
                            "Locomotionsubstraterelation_crawler", "Locomotionsubstraterelation_swimmer",
                            "Locomotionsubstraterelation_temporaryattached")

plot6e <- ggplot(locomotion, aes(y = value, x = as.numeric(year))) +  
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = emergence_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("Burrower", "Interstitial", "Crawler", "Swimmer", "Temporary attached"), name = "Locomotion") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Habit") + 
  xlab("Year") + 
  ggtitle("i")
plot6e

# Size
size_labs <- c("0.5-1 mm", "0.25-0.5 mm", "<0.25 mm", "1-2 mm", "2-4 mm",
               "4-8 mm", ">8 mm") 
names(size_labs) <- c("Maximal.potential.size_.0.5.1.cm", 
                      "Maximalpotentialsize_.0.25..5.cm", 
                      "Maximalpotentialsize_.0.25.cm",
                      "Maximalpotentialsize_.1.2.cm", 
                      "Maximalpotentialsize_.2.4.cm", 
                      "Maximalpotentialsize_.4.8.cm",
                      "Maximalpotentialsize_.8cm")

## Order them by actual size classes

plot6f <- ggplot(size, aes(y = value, x = as.numeric(year))) +  
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +   
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = size_labs), nrow = 1) + 
  scale_colour_discrete(labels = c("0.5-1 mm", "0.25-0.5 mm", "<0.25 mm", "1-2 mm", "2-4 mm",
                                   "4-8 mm", ">8 mm"), name = "Body size") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Body size") + 
  xlab("") + 
  ggtitle("b")
plot6f

# Reproduction
reproduction_labs <- c("Attached clutches", "Free clutches", "Vegetation clutches", "Terrestrial clutches", 
                       "Attached individual", "Free individual", "Ovoviviparity", "Parthenogenesis")
names(reproduction_labs) <- c("ReproductionTechnique_clutches.cementedorfixed", "ReproductionTechnique_clutches.free",           
                              "ReproductionTechnique_clutches.invegetation", "ReproductionTechnique_clutches.terrestrial",    
                              "ReproductionTechnique_isolatedeggs.cemented", "ReproductionTechnique_isolatedeggs.free",       
                              "ReproductionTechnique_ovoviviparity", "ReproductionTechnique_parthenogenesis")

plot6h <- ggplot(reproduction, aes(y = value, x = as.numeric(year))) +  
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +  
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) +   
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = rheophily_labs), ncol = 6) + 
  scale_colour_discrete(labels = c("Attached clutches", "Free clutches",
                                   "Vegetation clutches", "Terrestrial clutches", 
                                   "Attached individual", "Free individual",
                                   "Ovoviviparity", "Parthenogenesis"), 
                        name = "Oviposition") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Oviposition") + 
  xlab("Year") + 
  ggtitle("e")
plot6h

# Resistance
resistance_labs <- c("Cocoons", "Diapause", "Statoblasts", "None")
names(resistance_labs) <- c("Resistanceforms_cocoons", 
                            "Resistanceforms_diapauseordormancy",
                            "Resistanceforms_eggs.statoblasts",  
                            "Resistanceforms_none")

plot6i <- ggplot(resistance, aes(y = value, x = as.numeric(year))) +  
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +   
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = rheophily_labs), ncol = 6) + 
  scale_colour_discrete(labels = c("Cocoons", "Diapause", "Statoblasts", "None"), name = "Resistance") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Resistance forms") + 
  xlab("") + 
  ggtitle("g")
plot6i

# Respiration
respiration_labs <- c("Gills", "Spiracle", "Tegument")
names(respiration_labs) <- c("Respiration_gill", "Respiration_spiracle.aerial.",
                             "Respiration_tegument")

plot6j <- ggplot(respiration, aes(y = value, as.numeric(year))) +  
  geom_point(aes(fill = variable), pch = 21) + 
  geom_smooth(aes(colour = variable), fill = "lightgrey", method = "gam", se = F) +
  guides(fill = "none") +
  theme_bw() +  
  theme(legend.title = element_blank(), legend.text = element_text(size = 8), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8, colour = "black"), 
        legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
        title = element_text(size=8), plot.margin = unit(c(1, 2.5, 1, 1), "mm"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification = c(0,0.5)) + 
  scale_x_continuous(breaks = seq(1994, 2008, 2), limits = c(1994, 2008)) + 
  #facet_wrap(~variable, scales = "fixed", labeller = labeller(variable = rheophily_labs), ncol = 6) + 
  scale_colour_discrete(labels = c("Gills", "Spiracle", "Tegument"), name = "Respiration apparatus") +
  scale_y_continuous(limits = c(0,1)) +
  ylab("Respiration") + 
  xlab("") + 
  ggtitle("c")
plot6j

# Multiplot
plot_grid(plot6a, plot6b, plot6f, plot6i,
          plot6j, plot6c, plot6d, plot6e, plot6h, NULL, align = "hv", ncol = 2)
