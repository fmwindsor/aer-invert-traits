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
source("Code/1_Data-wrangling.R")


#### Matching traits to datasets #### 

### Wide sweep dataset

## Macroinvertebrates 

# Those groups that are missing from the database
missing <- colnames(macro_sweep_dframe)[(colnames(macro_sweep_dframe) %in% macro_traits_dframe$Genus) == F]
missing

# Those that are not at a suitable taxonomic resolution to include
remove <- c("Diptera", "Coleoptera", "Collembola", "Oligochaeta", 
            "Ephemeroptera", "Trichoptera", "Plecoptera", "Lepidoptera", 
            "Hydracarina")

# Update list to remove the orders that won't be included
still_missing <- missing[!(missing %in% remove)]
still_missing

# Check to see whether families are represented (and take their average trait values)
families <- c("Baetidae", "Capniidae", "Chironomidae", "Chloroperlidae", 
              "Heptageniidae", "Noctuidae", "Lumbricidae", "Naididae", 
              "Simuliidae", "Valvatidae")

# Merge the ancillary taxonomy with the trait data to look for families
traits_macro_family_dframe <- left_join(macro_traits_dframe, dplyr::select(traits_macro_taxononmy_dframe, -Species), by = "Genus")

# Find the records from genera within the families
traits_macro_family_dframe_sub <- traits_macro_family_dframe[traits_macro_family_dframe$Family == families,]

# Aggregate the values (mean)
family_agg <- aggregate(. ~ Family, FUN = mean, na.rm = T, na.action = NULL,
                        data = dplyr::select(traits_macro_family_dframe_sub, 
                                      Family, AdultFlyingStrength_abbrev_Strong:Feed_mode_sec_PA))

# Which families have data been collated for?
still_still_missing <- still_missing[!(still_missing %in% family_agg$Family)]
still_still_missing

# A further list of taxa not in the database (i.e., snails, moths and worms)
further_remove <- c("Noctuidae", "Nais", "Nematoda", "Pisidium", "Planorbis",
                    "Uncinais", "Valvatidae", "Lumbricidae", "Naididae")

# Combine the lists of the taxa removed (for later subsetting of the taxonomic dataframe)
removed <- c(remove, further_remove)

# What taxa are still missing trait data
still_still_still_missing <- still_still_missing[!(still_still_missing %in% further_remove)]
still_still_still_missing

# We will need to pull other trait data across for these taxa (find better representatives based on Sandy's knowledge)
additional_taxa <- data.frame(Genus = still_still_still_missing, 
                              Match = c("Heptageniidae", "Chironomidae", 
                                        "Orthocladius", "Tanytarsus", 
                                        "Chironomidae", "Camptocladius", 
                                        "Thienemannimyia"))

additional_taxa_genera <- additional_taxa[c(3,4,6,7),]
additional_taxa_families <- additional_taxa[c(1,2,5),]

# Collate the traits for the substite taxa (i.e., the matches for those not in the dataframe)
additionalgenus_traits_dframe <- macro_traits_dframe[macro_traits_dframe$Genus %in% additional_taxa_genera$Match,]
additionalfamily_traits_dframe <- traits_macro_family_dframe[traits_macro_family_dframe$Family == additional_taxa_families$Match,]

# Aggregate the values (mean)
additionalfamily_agg <- aggregate(. ~ Family, FUN = mean, na.rm = T, na.action = NULL,
                        data = dplyr::select(additionalfamily_traits_dframe, 
                                      Family, AdultFlyingStrength_abbrev_Strong:Feed_mode_sec_PA))

# Collate the trait data and merge it together
test <- left_join(additional_taxa_genera, additionalgenus_traits_dframe, by = c("Match" = "Genus"))
test1 <- left_join(additional_taxa_families, additionalfamily_agg, by = c("Match" = "Family"))
full <- bind_rows(test, test1)


## Create taxa by trait dataset

# Subset trait data based on the list of taxa 
macro_sweep_genus_traits_dframe <- macro_traits_dframe[macro_traits_dframe$Genus %in% colnames(macro_sweep_dframe),]

# Add in the family data 
macro_sweep_genusandfamily_traits_dframe <- bind_rows(macro_sweep_genus_traits_dframe, family_agg)

# Get rid of the redundant column an rename the existing one
macro_sweep_genusandfamily_traits_dframe[65:70,"Genus"] <- macro_sweep_genusandfamily_traits_dframe[65:70,"Family"]

# Add in the additional taxa data 
macro_sweep_all_traits_dframe <- bind_rows(macro_sweep_genusandfamily_traits_dframe, full)

# Remove the final column
macro_sweep_traits_dframe <- data.frame(macro_sweep_all_traits_dframe[,-c(59,60)])
rownames(macro_sweep_traits_dframe) <- macro_sweep_traits_dframe[,1]
macro_sweep_traits_dframe[,1] <- NULL
macro_sweep_traits_dframe[is.na(macro_sweep_traits_dframe)] <- 0

## Remove the taxa that are not in the trait database from the site by taxa matrix

# Remove columns from the site x taxon dataset
macro_sweep_clean_dframe <- macro_sweep_dframe[,!(colnames(macro_sweep_dframe) %in% removed)]
macro_sweep_clean_dframe1 <- macro_sweep_clean_dframe[-31,]

# Check that the traits x taxa and the taxa x site matrices match 
sort(rownames(macro_sweep_traits_dframe)) == sort(colnames(macro_sweep_clean_dframe1))


## Meiofauna

# Generate a matching column to help bind the traits to the taxa 
meio_fuzzy_traits <- data.frame(taxon = colnames(meio_sweep_dframe),
                                match = c("Oligochaeta", "Nematoda", NA, 
                                          "Hydracarina", "Ostracoda", 
                                          "Harpaticoida", "Harpaticoida", 
                                          "Harpaticoida", "Harpaticoida",
                                          "Harpaticoida", "Harpaticoida",
                                          "Harpaticoida", "Harpaticoida",
                                          "Harpaticoida", "Harpaticoida",
                                          "Cyclopoida", "Cyclopoida",
                                          "Cyclopoida", "Cyclopoida",
                                          "Alona_quadrangularis", "Alona_guttata", 
                                          "Chydorus_sphaericus",
                                          "Chydorus_sphaericus", NA))

meio_fuzzy_traits_dframe <- left_join(meio_fuzzy_traits, additional_traits_meio_dframe, 
                                      by = c("match" = "Taxon"))

# Taxa to remove from both datasets
remove <- c("tardigrades", "macrothricidae")

# Clean the taxonomic data
meio_taxa_sweep_dframe <- dplyr::select(meio_sweep_dframe, -all_of(remove))

# Clean the trait data 
meio_sweep_fuzzy_traits_dframe <- meio_fuzzy_traits_dframe[-c(3,24),]
meio_sweep_fuzzy_traits_dframe$match <- NULL
rownames(meio_sweep_fuzzy_traits_dframe) <- meio_sweep_fuzzy_traits_dframe$taxon
meio_sweep_fuzzy_traits_dframe$taxon <- NULL

# Check that the traits x taxa and the taxa x site matrices match 
sort(rownames(meio_sweep_fuzzy_traits_dframe)) == sort(colnames(meio_taxa_sweep_dframe))


## Additional quantitative data

# Those groups that are missing from the database
missing <- colnames(meio_sweep_dframe)[(colnames(meio_sweep_dframe) %in% traits_meio_dframe$taxon) == F]
missing

# Harpacticoids are not well covered (but have a grouped set of traits)
harpacticoids <- data.frame(taxon = c("attheyella_idahoensis", "attheyella_ilinoisensis", 
                                      "epactophanes_brucei", "maraenobiotus_insignipes", 
                                      "maraenobiotus_brucei", "moraria_affinis", 
                                      "nitocra_hibernicus", "parastenocaris_sp",
                                      "bryocamptus_hiemalis", "bryocamptus_zschokkei"), 
                           Match = c("Harpacticoids", "Harpacticoids",
                                     "Harpacticoids", "Harpacticoids",
                                     "Harpacticoids", "Harpacticoids",
                                     "Harpacticoids", "Harpacticoids",
                                     "Harpacticoids", "Harpacticoids")) 

# Which taxa are still missing 
still_missing <- missing[!missing %in% harpacticoids$taxon]
still_missing

# Extract genus names from the community data
still_still_missing <- str_to_sentence(str_split_i(still_missing, pattern = "_", i = 1))
still_still_missing

# Genera that are in the database
genera <- data.frame(taxon = c("alona_sp", "chydorus_sp", "diacyclops_sp", "paracyclops_poppei"), 
                     Match = c("Alona", "Chydorus", "Diacyclops", "Paracyclops")) 

# Which taxa are still missing, even now
still_still_still_missing <- tolower(still_still_missing[!still_still_missing %in% traits_meio_dframe$Genus])
still_still_still_missing

# Additional taxa 
additional_taxa <- data.frame(taxon = c("Graptoleberis", "Macrothricidae"), 
                              Match = c("Chydorus", "Macrothrix"))

# Taxa where we cannot assign traits at a meaningful level
remove <- still_still_still_missing[!still_still_still_missing %in% additional_taxa$taxon]

## Prepare some of the dataframes for matching non-direct matches

# Select numeric columns from the trait dataframe
num_cols <- colnames(select_if(traits_meio_dframe, is.numeric))

# Aggregate the trait dataframe into genus level mean values
genus_agg <- aggregate(. ~ Genus, FUN = mean, na.rm = T, na.action = NULL,
                       data = dplyr::select(traits_meio_dframe, Genus, all_of(num_cols)))

genus_agg$taxon <- genus_agg$Genus

taxon_traits <- aggregate(. ~ taxon, FUN = mean, na.rm = T, na.action = NULL,
                          data = dplyr::select(traits_meio_dframe, taxon, all_of(num_cols)))

## Create taxa by trait dataset

# Subset trait data based on the list of taxa 
meio_sweep_taxon_traits_dframe <- taxon_traits[taxon_traits$taxon %in% colnames(meio_sweep_dframe),]
colnames(meio_sweep_taxon_traits_dframe)[1] <- "taxon"

# Add in the additional taxa from the genus
meio_sweep_genus_traits_dframe <- left_join(genera, genus_agg, 
                                            by = c("Match" = "Genus"))
colnames(meio_sweep_genus_traits_dframe)[1] <- "taxon"
meio_sweep_genus_traits_dframe$taxon.y <- NULL
meio_sweep_genus_traits_dframe$Match <- NULL

# Harpacticoids 
meio_sweep_harp_traits_dframe <- left_join(harpacticoids, genus_agg, 
                                           by = c("Match" = "Genus"))
colnames(meio_sweep_harp_traits_dframe)[1] <- "taxon"
meio_sweep_harp_traits_dframe$taxon.y <- NULL
meio_sweep_harp_traits_dframe$Match <- NULL

# Additionals 
meio_sweep_adds_traits_dframe <- left_join(additional_taxa, genus_agg, 
                                           by = c("Match" = "Genus"))
colnames(meio_sweep_adds_traits_dframe)[1] <- "taxon"
meio_sweep_adds_traits_dframe$taxon.y <- NULL
meio_sweep_adds_traits_dframe$Match <- NULL

# Add in the family data 
meio_sweep_quant_traits_dframe <- bind_rows(meio_sweep_taxon_traits_dframe, 
                                            meio_sweep_genus_traits_dframe, 
                                            meio_sweep_harp_traits_dframe, 
                                            meio_sweep_adds_traits_dframe) 

rownames(meio_sweep_quant_traits_dframe) <- meio_sweep_quant_traits_dframe$taxon

meio_sweep_quant_traits_dframe <- meio_sweep_quant_traits_dframe[,c("Body.length", "Dry.mass")]


## Clean up the environment

rm(additional_taxa, additional_taxa_families, additional_taxa_genera,
additionalfamily_agg, additionalfamily_traits_dframe, 
additionalgenus_traits_dframe, families, family_agg,
full, further_remove, genera, genus_agg, harpacticoids, 
macro_1997_dframe, macro_sweep_all_traits_dframe, macro_sweep_clean_dframe,
macro_sweep_dframe, macro_sweep_genus_traits_dframe, meio_sweep_dframe,
macro_sweep_genusandfamily_traits_dframe,
meio_1997_dframe, meio_flood_dframe, 
meio_fuzzy_traits, meio_fuzzy_traits_dframe, meio_sweep_adds_traits_dframe,
meio_sweep_genus_traits_dframe, meio_sweep_harp_traits_dframe,
meio_sweep_taxon_traits_dframe, 
missing, num_cols, remove, removed, still_missing, 
still_still_missing, still_still_still_missing, 
taxon_traits, test, test1, traits_macro_family_dframe, 
traits_macro_family_dframe_sub)

### Long term dataset

## Macroinvertebrates 

# Those groups that are missing from the database
missing <- colnames(macro_WPC_dframe)[(colnames(macro_WPC_dframe) %in% macro_traits_dframe$Genus) == F]
missing

# Those that are not at a suitable taxonomic resolution to include
remove <- c("Oligochaeta", "Hydrachnidae", "Hydracarina")

# Update list to remove the orders that won't be included
still_missing <- missing[!(missing %in% remove)]
still_missing

# Check to see whether families are represented (and take their average trait values)
families <- c("Simuliidae", "Ceratopogoniidae", "Gammaridae", "Dytiscidae", 
              "Empididae", "Dixidae", "Daphniidae", "Perlidae", "Leuctridae",
              "Siphlonuridae")

# Merge the ancillary taxonomy with the trait data to look for families
traits_macro_family_dframe <- left_join(macro_traits_dframe, dplyr::select(traits_macro_taxononmy_dframe, -Species), by = "Genus")

# Find the records from genera within the families
traits_macro_family_dframe_sub <- traits_macro_family_dframe[traits_macro_family_dframe$Family == families,]

# Aggregate the values (mean)
family_agg <- aggregate(. ~ Family, FUN = mean, na.rm = T, na.action = NULL,
                        data = dplyr::select(traits_macro_family_dframe_sub, 
                                      Family, AdultFlyingStrength_abbrev_Strong:Feed_mode_sec_PA))

# Which families have data been collated for?
still_still_missing <- still_missing[!(still_missing %in% family_agg$Family)]
still_still_missing

# A further list of taxa not in the database (i.e., snails, moths and worms)
further_remove <- c("Daphniidae", "Gammaridae", "Dixidae", "Ceratopogoniidae",
                    "Hesperophila")

# Combine the lists of the taxa removed (for later subsetting of the taxonomic dataframe)
removed <- c(remove, further_remove)

# What taxa are still missing trait data
still_still_still_missing <- still_still_missing[!(still_still_missing %in% further_remove)]
still_still_still_missing

# We will need to pull other trait data across for these taxa (find better representatives based on Sandy's knowledge)
additional_taxa <- data.frame(Genus = still_still_still_missing, 
                              Match = c("Cricotopus", "Eukiefferiella", 
                                        "Eukiefferiella", "Eukiefferiella", 
                                        "Eukiefferiella","Eukiefferiella", 
                                        "Eukiefferiella", "Brillia", 
                                        "Eukiefferiella", "Orthocladius", 
                                        "Orthocladius", "Pagastia", "Chironomidae", 
                                        "Suwallia", "Seratella"))

additional_taxa_genera <- additional_taxa[-13,]
additional_taxa_families <- additional_taxa[13,]

# Collate the traits for the substite taxa (i.e., the matches for those not in the dataframe)
additionalgenus_traits_dframe <- macro_traits_dframe[macro_traits_dframe$Genus %in% additional_taxa_genera$Match,]
additionalfamily_traits_dframe <- traits_macro_family_dframe[traits_macro_family_dframe$Family == additional_taxa_families$Match,]

# Aggregate the values (mean)
additionalfamily_agg <- aggregate(. ~ Family, FUN = mean, na.rm = T, na.action = NULL,
                                  data = dplyr::select(additionalfamily_traits_dframe, 
                                                Family, AdultFlyingStrength_abbrev_Strong:Feed_mode_sec_PA))

# Collate the trait data and merge it together
test <- left_join(additional_taxa_genera, additionalgenus_traits_dframe, by = c("Match" = "Genus"))
test1 <- left_join(additional_taxa_families, additionalfamily_agg, by = c("Match" = "Family"))
full <- bind_rows(test, test1)

## Create taxa by trait dataset

# Subset trait data based on the list of taxa 
macro_WPC_genus_traits_dframe <- macro_traits_dframe[macro_traits_dframe$Genus %in% colnames(macro_WPC_dframe),]

# Add in the family data 
macro_WPC_genusandfamily_traits_dframe <- bind_rows(macro_WPC_genus_traits_dframe, family_agg)

# Get rid of the redundant column an rename the existing one
macro_WPC_genusandfamily_traits_dframe[23:28,"Genus"] <- macro_WPC_genusandfamily_traits_dframe[23:28,"Family"]

# Add in the additional taxa data 
macro_WPC_all_traits_dframe <- bind_rows(macro_WPC_genusandfamily_traits_dframe, full)

# Remove the final column
macro_WPC_traits_dframe <- data.frame(macro_WPC_all_traits_dframe[,-c(59,60)])
rownames(macro_WPC_traits_dframe) <- macro_WPC_traits_dframe[,1]
macro_WPC_traits_dframe[,1] <- NULL
macro_WPC_traits_dframe[is.na(macro_WPC_traits_dframe)] <- 0

## Remove the taxa that are not in the trait database from the site by taxa matrix

# Remove columns from the site x taxon dataset
macro_WPC_clean_dframe <- macro_WPC_dframe[,!(colnames(macro_WPC_dframe) %in% removed)]

# Check that the traits x taxa and the taxa x site matrices match 
sort(rownames(macro_WPC_traits_dframe)) == sort(colnames(macro_WPC_clean_dframe))



## Meiofauna

# Generate a matching column to help bind the traits to the taxa 
meio_fuzzy_traits <- data.frame(taxon = colnames(meio_WPC_dframe),
                                match = c("Chironomidae", NA, NA, NA, NA,
                                          "Hydracarina", "Nematoda", 
                                          "Oligochaeta", "Oligochaeta", NA,
                                          "Ostracoda", "Harpaticoida",
                                          "Harpaticoida", "Harpaticoida", 
                                          "Harpaticoida", "Harpaticoida",
                                          "Harpaticoida", "Harpaticoida",
                                          "Cyclopoida", "Cyclopoida",
                                          "Cyclopoida", "Cyclopoida",
                                          "Chydorus_sphaericus", "Cladocera",
                                          "Cladocera"))
                                          
meio_fuzzy_traits_dframe <- left_join(meio_fuzzy_traits, additional_traits_meio_dframe, 
                                      by = c("match" = "Taxon"))

# Taxa to remove from both datasets
remove <- c("tardigrades", "simuliidae", "plecoptera", "ephemeroptera",
            "trichopter")

# Clean the taxonomic data
meio_taxa_WPC_dframe <- dplyr::select(meio_WPC_dframe, -all_of(remove))

# Clean the trait data 
meio_WPC_fuzzy_traits_dframe <- meio_fuzzy_traits_dframe[-c(2,3,4,5,10),]
meio_WPC_fuzzy_traits_dframe$match <- NULL
rownames(meio_WPC_fuzzy_traits_dframe) <- meio_WPC_fuzzy_traits_dframe$taxon
meio_WPC_fuzzy_traits_dframe$taxon <- NULL

# Check that the traits x taxa and the taxa x site matrices match 
sort(rownames(meio_WPC_fuzzy_traits_dframe)) == sort(colnames(meio_taxa_WPC_dframe))


## Additional quantitative data

# Those groups that are missing from the database
missing <- colnames(meio_WPC_dframe)[(colnames(meio_WPC_dframe) %in% traits_meio_dframe$taxon) == F]
missing

# Harpacticoids are not well covered (but have a grouped set of traits)
harpacticoids <- data.frame(taxon = c("maraenobiotus_brucei", "epactophanes_brucei",
                                      "mesochra_alaskana", "X.atheyella_sp",
                                      "moraria_affinis", "bryocamptus_zschokkei", 
                                      "bryocamptus_hiemalis"), 
                            Match = c("Harpacticoids", "Harpacticoids",
                                      "Harpacticoids", "Harpacticoids",
                                      "Harpacticoids", "Harpacticoids",
                                      "Harpacticoids")) 

# Which taxa are still missing 
still_missing <- missing[!missing %in% harpacticoids$taxon]
still_missing

# Extract genus names from the community data
still_still_missing <- str_to_sentence(str_split_i(still_missing, pattern = "_", i = 1))
still_still_missing

# Genera that are in the database
genera <- data.frame(taxon = c("paracyclops_poppei", "cyclop.s_scutifer", "chydorus_sp", "pleuroxus_sp"), 
                     Match = c("Paracyclops", "Cyclops", "Chydorus", "Pleuroxus")) 

# Which taxa are still missing, even now
still_still_still_missing <- tolower(still_still_missing[!still_still_missing %in% traits_meio_dframe$Genus])
still_still_still_missing <- still_still_still_missing[-12]
still_still_still_missing

# Taxa where we cannot assign traits at a meaningful level
remove <- still_still_still_missing[!still_still_still_missing %in% additional_taxa$taxon]

## Prepare some of the dataframes for matching non-direct matches

# Select numeric columns from the trait dataframe
num_cols <- colnames(select_if(traits_meio_dframe, is.numeric))

# Aggregate the trait dataframe into genus level mean values
genus_agg <- aggregate(. ~ Genus, FUN = mean, na.rm = T, na.action = NULL,
                       data = dplyr::select(traits_meio_dframe, Genus, all_of(num_cols)))

genus_agg$taxon <- genus_agg$Genus

taxon_traits <- aggregate(. ~ taxon, FUN = mean, na.rm = T, na.action = NULL,
                          data = dplyr::select(traits_meio_dframe, taxon, all_of(num_cols)))

## Create taxa by trait dataset

# Subset trait data based on the list of taxa 
meio_WPC_taxon_traits_dframe <- taxon_traits[taxon_traits$taxon %in% colnames(meio_WPC_dframe),]
colnames(meio_WPC_taxon_traits_dframe)[1] <- "taxon"

# Add in the additional taxa from the genus
meio_WPC_genus_traits_dframe <- left_join(genera, genus_agg, 
                                            by = c("Match" = "Genus"))
colnames(meio_WPC_genus_traits_dframe)[1] <- "taxon"
meio_WPC_genus_traits_dframe$taxon.y <- NULL
meio_WPC_genus_traits_dframe$Match <- NULL

# Harpacticoids 
meio_WPC_harp_traits_dframe <- left_join(harpacticoids, genus_agg, 
                                           by = c("Match" = "Genus"))
colnames(meio_WPC_harp_traits_dframe)[1] <- "taxon"
meio_WPC_harp_traits_dframe$taxon.y <- NULL
meio_WPC_harp_traits_dframe$Match <- NULL

# Add in the family data 
meio_WPC_quant_traits_dframe <- bind_rows(meio_WPC_taxon_traits_dframe, 
                                            meio_WPC_genus_traits_dframe, 
                                            meio_WPC_harp_traits_dframe) 

rownames(meio_WPC_quant_traits_dframe) <- meio_WPC_quant_traits_dframe$taxon

meio_WPC_quant_traits_dframe <- meio_WPC_quant_traits_dframe[,c("Body.length", "Dry.mass")]


## Clean up the environment

rm(additional_taxa, additional_taxa_families, additional_taxa_genera,
   additional_traits_meio_dframe, additionalfamily_agg, 
   additionalfamily_traits_dframe, additionalgenus_traits_dframe, families,
   family_agg, full, further_remove, genera, genus_agg, harpacticoids, 
   macro_traits_dframe, macro_WPC_dframe, meio_fuzzy_traits,
   meio_fuzzy_traits_dframe, meio_WPC_genus_traits_dframe, 
   meio_WPC_harp_traits_dframe, meio_WPC_taxon_traits_dframe,
   meio_WPC_dframe, missing, num_cols, remove, removed, still_missing, 
   still_still_missing, still_still_still_missing, 
   taxon_traits, test, test1, traits_macro_family_dframe, 
   traits_macro_family_dframe_sub, traits_macro_taxononmy_dframe, 
   traits_meio_dframe, macro_WPC_all_traits_dframe, macro_WPC_genus_traits_dframe,
   macro_WPC_genusandfamily_traits_dframe)




