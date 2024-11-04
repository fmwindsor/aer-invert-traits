# Invertebrate functional trait variation along successional gradients #
# Code produced by Fredric M. Windsor (fmwindsor@gmail.com) #
# All code is either original or the source is cited in line #


## 1 - Data preparation


#### Setup ####

# Clear environment
rm(list = ls())

# Set working directory
setwd("/Users/c1513054/Library/CloudStorage/OneDrive-CardiffUniversity/Documents/Research/Papers/Published/Advances in Ecological Research (Glacier Bay)/Macroinvertebrate traits")

# Load the necessary libraries
library(plyr); library(dplyr); library(gridExtra); library(cowplot)
library(ggplot2); library(geometry); library(FD)
library(vegan); library(reshape2); library(stringr); library(grid)
library(tidyr); library(ade4); library(MASS); library(stringi); library(stringr)


#### Data input ####

## 1997 wide sweep samples (16 streams)

# Macroinvertebrate data
macro_1997_dframe <- read.csv("Data/1997GLB_inverts.csv")

# Meiofauna data
meio_1997_dframe <- read.csv("Data/1997GLB_meio.csv")

# Read in the 1997 sweep water chemistry data
chem_sweep_dframe <- read.csv("Data/1997GLB_water-chem.csv")


## Wolf Point Creek Long Term data (1978-2017)

# Macroinvertebrate data
macro_WPC_dframe <- read.csv("Data/WPC1978to2017_inverts.csv")

# Meiofauna data
meio_WPClong_dframe <- read.csv("Data/WPC1994to2008_meio.csv")


## Trait data 

# CONUS trait data (https://portal.edirepository.org/nis/mapbrowse?packageid=edi.481.5)
traits_macro_dframe <- read.csv("Data/Genus_Trait_Affinities.csv")

# CONUS ancillary taxonomy 
traits_macro_taxononmy_dframe <- read.csv("Data/Ancillary_Taxonomy.csv")

# Meiofauna trait data (https://doi.org/10.1890/15-1275.1)
traits_meio_dframe <- read.csv("Data/meiofauna_traits_new.csv", dec = ",")
traits_meio_dframe$taxon <- tolower(paste(traits_meio_dframe$Genus, sep = "_", traits_meio_dframe$Species))

# Additional meiofauna traits (http://dx.doi.org/10.1016/j.scitotenv.2013.06.082)
additional_traits_meio_dframe <- read.csv("Data/meiofauna_traits.csv")


#### Cleaning data #### 

## 1997 wide sweep samples (16 streams)

# Macroinvertebrate data
inverts <- macro_1997_dframe[,-1]
invert_names <- c(inverts$Genus.family)
clean_invert_names <- c("Acanthamola", "Baetis", "Acentrella", "Baetis", 
                        "Baetidae", "Barbaetis", "Cinygmula", "Drunella",
                        "Epeorus", "Heptageniidae", "Paracloeodes", 
                        "Rhithrogena", "Ephemeroptera", "Alaskaperla", "Capnia",
                        "Capniidae", "Chloroperlidae", "Haploperla", "Perlomyia",
                        "Plumiperla", "Podmosta", "Sweltsa", "Taenionema", 
                        "Zapada", "Plecoptera", "Ecclisomyia", "Moselyana", 
                        "Onocosmoecus", "Rhyacophila", "Himalopsyche",
                        "Trichoptera", "Chelifera", "Dicranota", "Hexatoma", 
                        "Hesperoconopa", "Diptera", "Molophilus", "Oreogeton",
                        "Probezzia", "Prosimulium", "Simuliidae", "Simulium",
                        "Lumbricidae", "Naididae", "Nais", "Nematoda", 
                        "Uncinais", "Oligochaeta", "Pisidium", "Planorbis", 
                        "Valvatidae", "Agathon", "Noctuidae", "Lepidoptera",
                        "Coleoptera", "Collembola","Hydracarina", 
                        "Paralimnophyes", "Orthocladiinae", "Orthocladius",
                        "Cricotopus", "Paratrichocladius", "Euorthocladius", 
                        "Synorthocladius", "Parorthocladius", "Brillia",
                        "Corynoneura", "Eukiefferiella", "Eukiefferiella",
                        "Tvetenia", "Diplocladius", "Camptocladius", 
                        "Chaetocladius", "Hydrobaenus", "Paracladius", 
                        "Parakiefferiella", "Paraphaenocladius", 
                        "Psectrocladius", "Rheocricotopus", "Thienemanniella",
                        "Telmatopelopia", "Zavrelimyia", "Thienemannimyia",
                        "Boreochlus", "Micropsectra", "Krenopsectra",
                        "Polypedilum", "Paracladopelma", "Sergentia",
                        "Zavreliella", "Endochironomus", "Endochironomus",
                        "Chironominae", "Diamesa", "Pagastia", "Potthastia", 
                        "Pseudokiefferiella", "Chironomidae")

# Create a new clean dataframe
inverts_clean <- data.frame(Taxon = clean_invert_names, 
                                   dplyr::select(inverts, -Genus.family))

# Calculate the sum for several taxa that now have multiple rows
inverts_agg <- aggregate(. ~ Taxon, data = inverts_clean, FUN = "sum")

# Transpose the dataframe so it is in a sample x species format
inverts_t <- t(inverts_agg[,-1])
colnames(inverts_t) <- inverts_agg[,1]

# Change the name back to something shorter
macro_sweep_dframe <- data.frame(inverts_t)

# Remove sample sites with no observations 
meio_sweep <- dplyr::select(meio_1997_dframe, chaetogaster:macrothricidae)
rownames(meio_sweep) <- paste(meio_1997_dframe$stream, sep = "", meio_1997_dframe$replicates)
#t_df <- as.data.frame(t(meio_sweep))
#colnames(t_df) <- rownames(meio_sweep)
#new_df <- aggregate(t_df, by=list(rownames(t_df)), sum)
#meio_agg <- data.frame(t(new_df[,-1]))
#colnames(meio_agg) <- new_df[,1]
meio_sweep_dframe <- meio_sweep[-c(76:80),] # remove Carolus river
meio_sweep_dframe[is.na((meio_sweep_dframe))] <- 0


## Wolf Point Creek Long Term data (1978-2017)

# Macroinvertebrate data 
inverts_long_t <- t(macro_WPC_dframe[,-1])
colnames(inverts_long_t) <- macro_WPC_dframe[,1]
macro_WPC_dframe <- data.frame(inverts_long_t)
macro_WPC_dframe$Time <- str_sub(rownames(macro_WPC_dframe), start = 2)
rownames(macro_WPC_dframe) <- macro_WPC_dframe$Time
macro_WPC_dframe$Time <- NULL

colnames(macro_WPC_dframe) <- c("Diamesa", "Brillia", "Pseudodiamesa", 
                                "Chaetocladius", "Cricotopus_tremulus", 
                                "Eukiefferiella_brehmi", "Eukiefferiella_claripennis",
                                "Eukiefferiella_cyanea", "Eukiefferiella_devonica",
                                "Eukiefferiella_gracei", "Eukiefferiella_rectangularis",
                                "Tokunagaia", "Eukiefferiella", "Eukiefferiella_tveta",
                                "Orthocladius", "Orthocladius_mallochi",
                                "Orthocladius_manitobensis", "Pagastia_partica",
                                "Paratrichocladius", "Potthastia", "Micropsectra",
                                "Neuropertona", "Paracladopelma", "Corynoneura", 
                                "Polypedilum", "Suwallia_forcipata", "Baetis",
                                "Cinygmula", "Seratella_tibialis", "Simuliidae",
                                "Limnophila", "Oligochaeta", "Ceratopogoniidae",
                                "Onocosmoecus", "Nyctiophylax", "Ecclisomyia",
                                "Gammaridae", "Dytiscidae", "Empididae",
                                "Hydrachnidae", "Dixidae", "Daphniidae", 
                                "Hesperophila", "Perlidae", "Leuctridae",
                                "Siphlonuridae", "Ephemerella", "Taeniopteryx",
                                "Kathroperla", "Brachycentrus", "Hydracarina") 


# Meio data
meio_WPClong_dframe$date <- as.Date(meio_WPClong_dframe$Date, format = "%d.%m.%Y")
meio_WPClong_dframe$year <- str_sub(meio_WPClong_dframe$date, end = 4)
meio_WPC_dframe <- aggregate(. ~ year, 
                                data = dplyr::select(meio_WPClong_dframe, -Stream, -Date,
                                              -replicate), FUN = "sum")
rownames(meio_WPC_dframe) <- meio_WPC_dframe$year
meio_WPC_dframe$year <- NULL
meio_WPC_dframe$date <- NULL

## Trait data

# Create a unique identifier for each of the traits
traits_macro_dframe$Trait_merged <- paste0(traits_macro_dframe$Trait_group,  
                                         sep = "_", 
                                         traits_macro_dframe$Trait)

# Spread the data into a wide format to form a species x trait dataframe
macro_traits_dframe <- pivot_wider(traits_macro_dframe, id_cols = "Genus", 
                                   names_from = "Trait_merged", 
                                   values_from = "Trait_affinity")



## Chemistry data

# Create a successional gradient metric (same as the GI)
chem_index <- data.frame(row.names = chem_sweep_dframe$Stream)

# Centre and scale variables
chem_index$turb <- scale(1/chem_sweep_dframe$Turbidity.NTU.)
chem_index$temp <- scale(chem_sweep_dframe$Temp.C.)
chem_index$pfan <- scale(1/chem_sweep_dframe$Pfan_bottom)
chem_index$cond <- scale(chem_sweep_dframe$EC)
chem_index$p <- scale(chem_sweep_dframe$Total.P.ugL.)
chem_index$n <- scale(chem_sweep_dframe$totalN.ugL.)
chem_index$wood <- scale(chem_sweep_dframe$X.wood)
chem_index$fish <- scale(chem_sweep_dframe$totfish.CPUE.)

# Rename columns
colnames(chem_index) <- c("turb", "temp", "pfan", "cond", "p", "n", "wood", 
                          "fish")

# Run a PCA 
pca1 <- princomp(chem_index, cor = F)
pca1$loadings

#### Clean the environment ####

# Remove the objects that are no longer needed
rm(clean_invert_names, invert_names, inverts, inverts_t,
   inverts_clean, inverts_agg, inverts_long_t, traits_macro_dframe, meio_sweep, 
   meio_WPClong_dframe, pca1)
