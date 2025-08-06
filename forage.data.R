# Bombus forage plant analysis
# Updated 2025-07-04
# 
# Based on data downloaded from iNaturalist from the project
# Pollinator Interactions on Plants (PIP) of the NE US
# Bombus species in the state of Maine
# Query quality_grade=any&identifications=any&place_id=17&taxon_id=52775&projects[]=pollinator-interactions-on-plants-pip-of-the-ne-us
# Downloaded on July 2, 2025

# Load packages
library(tidyverse)
library(vegan)


study.species <- c("borealis","bimaculatus","griseocollis","perplexus","vagans","impatiens","ternarius","terricola")

# Filter to the research-grade observations of species in this study
forage <- read.csv("observations-590886.csv/observations-590886.csv") %>% 
  filter(
    (sub("Bombus ","",taxon_species_name) %in% study.species) &
      quality_grade == "research"
  )

names(forage)
dim(forage)


# Organize the interaction counts by bee species
forage.by.species <- forage %>% 
  select(uuid, species = taxon_species_name, plant = field.interaction..visited.flower.of) %>% 
  mutate(observed = 1) %>% 
  pivot_wider(names_from = "plant", values_from = "observed", values_fill = 0) %>% 
  group_by(species) %>% 
  select(-uuid) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  column_to_rownames(var = "species")

# dim(forage.by.species)
# rowSums(forage.by.species)
# forage.by.species[,1:4]

write_csv(forage.by.species, "forage.data.by.species.csv")


####################
# Rarefied species richness
####################

(rarefaction.depth <- min(rowSums(forage.by.species[,-1])))
# 35

# Rarefied species richness based on Hurlbert's (1971) formulation, 
# and the standard errors on Heck et al. (1975).
rarefied.richness <- rarefy(forage.by.species, sample = rarefaction.depth) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "species") %>% 
  rename_with(~ "richness", .cols = 2) 

rarefied.richness%>% 
  arrange(richness)
#               species richness
# 1 Bombus griseocollis 19.76997
# 2     Bombus borealis 20.36111
# 3    Bombus terricola 23.35433
# 4  Bombus bimaculatus 24.85353
# 5    Bombus ternarius 25.00751
# 6    Bombus perplexus 25.65634
# 7    Bombus impatiens 27.76579
# 8       Bombus vagans 28.00000


####################
# Simpson's index
####################

simpsons.index <- diversity(forage.by.species, index = "simpson") %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "species") %>% 
  rename_with(~ "Simpsons", .cols = 2) 

simpsons.index %>% 
  arrange(Simpsons)
#               species  Simpsons
# 1 Bombus griseocollis 0.9119117
# 2     Bombus borealis 0.9257726
# 3       Bombus vagans 0.9551020
# 4    Bombus terricola 0.9565354
# 5    Bombus perplexus 0.9611365
# 6    Bombus ternarius 0.9612070
# 7  Bombus bimaculatus 0.9615105
# 8    Bombus impatiens 0.9792125


####################
# NMDS on forage plant community composition
####################

# Using Bray-Curtis distances
nmds.result <- metaMDS(
  forage.by.species, 
  distance = "bray", 
  k = 2, 
  trymax = 100
)

# Check the stress (lower is better; <0.2 usually okay)
nmds.result$stress

# Plot 
plot(nmds.result, type = "t")

# Extract coordinates
nmds.result <- nmds.result$points %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "species")


####################
# Save all the diversity metrics
####################

forage.plant.diversity <- list(
  rarefied.richness,
  simpsons.index,
  nmds.result 
) %>% 
  reduce(full_join, by = "species")

forage.plant.diversity %>% 
  arrange(Simpsons)

# Save the table to file
write.csv(
  forage.plant.diversity,
  "forage.plant.diversity.csv", 
  row.names = FALSE, 
  quote = FALSE
)

# A look at some correlations
{
  with(forage.plant.diversity, plot(richness, Simpsons, type = "n"))
  with(forage.plant.diversity, text(richness, Simpsons, sub("Bombus ","",species), srt = 45))
  with(forage.plant.diversity, abline(lm(Simpsons ~ richness)))
}
{
  with(forage.plant.diversity, plot(Simpsons, MDS1, type = "n"))
  with(forage.plant.diversity, text(Simpsons, MDS1,  sub("Bombus ","",species), srt = 45))
  with(forage.plant.diversity, abline(lm(MDS1 ~ Simpsons)))
}

# End of script
# Proceed with `mouthpart.analysis.R`
 
