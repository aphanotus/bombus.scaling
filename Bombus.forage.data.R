# Bombus forage plant analysis

# Clear memory
# rm(list=ls())
# gc() 

# Load packages
x <- c("tidyverse","vegan","picante","phytools")
invisible(lapply(x, require, character.only = TRUE))

# file.choose()
forage <- read.csv("Bombus.forage.data.csv")
dim(forage)
head(forage)
tail(forage)

# Possible metrics describing bees as generalist vs specialist
# richness
# phylogenetic diversity
# morphological diversity
# dietary breadth (sensu Wood et al. 2019, see  https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2697#Rarefaction-and-dietary-breadth )
#    This metric will be impossible to replicate!
#    The authors used a vaguely-described method of "manual rarefaction"
#    based on Michigan geography unique to their study 
#    without enough details for replication.

# Filter just to the dataset from Wood et al. 2019

wood <- filter(forage, grepl("Wood",study)) %>% select(-study)

# Sample sizes obtained from Wood et al. Fig. 3
wood.sample.sizes <- c(46,52,78,37,66,86,144,41,58,62,65,76)
names(wood.sample.sizes) <- c("affinis","auricomus","bimaculatus","borealis","fervidus","griseocollis","impatiens","pensylvanicus","perplexus","ternarius","terricola","vagans")
tongue.class <- c("short","long","medium","long","long","short","medium","long","medium","medium","short","medium-long")
names(tongue.class) <- c("affinis","auricomus","bimaculatus","borealis","fervidus","griseocollis","impatiens","pensylvanicus","perplexus","ternarius","terricola","vagans")

# Import values for Wood's dietary breadth score
diet <- read.csv("Wood.et.al.2019.dietary.breadth.csv")

# PCA
wood.wide <- wood %>% 
  filter(forage.species != 'others') %>% 
  select(-common.plant.name, -family) %>% 
  pivot_wider(names_from = "forage.species", values_from = "interactions", values_fill = 0) %>% 
  column_to_rownames("Bombus.species")

wood.pca <- prcomp(wood.wide, scale = TRUE)

with(wood.pca, plot(x[,'PC1'],x[,'PC2'], type = "n"))
with(wood.pca, text(x[,'PC1'],x[,'PC2'], row.names(x), cex =0.65))

# Alpha-Diversity metrics

# Richness
S <- specnumber(wood.wide) # observed number of species
raremax <- min(rowSums(wood.wide))
Srare <- rarefy(wood.wide, raremax)
{
  par(mfrow=c(2,1))
  plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", type = 'n')
  text(S, Srare, names(S), cex=0.6)
  abline(0, 1, col='darkred')
  rarecurve(wood.wide, step = 20, sample = raremax, col = "darkred", cex = 0.6)
  par(mfrow=c(1,1))
}

wood.div <- diet %>% 
  mutate(
    richness = Srare,
    shannon = diversity(wood.wide, index = "shannon"),
    simpson = diversity(wood.wide, index = "simpson")
  ) 

alpha.plot <- function(index, ...) {
  x <- as.numeric(unlist(wood.div[,index]))
  plot(x, wood.div$dietary.breadth, type = "n", xlab = index, ...)
  text(x, wood.div$dietary.breadth, wood.div$species, cex =0.65)
}

{
  par(mfrow=c(2,2))
  alpha.plot("richness")
  alpha.plot("shannon")
  alpha.plot("simpson")
  par(mfrow=c(1,1))
}

# Calculate Faith's phylogenetic diversity for each species' forage plant assemblage
# based on family-level relationships from Qian & Zhang 2014 JSE
# https://doi.org/10.1111/jse.12086

wood.fam <- wood %>% 
  filter(!forage.species == "others") %>% 
  select(-common.plant.name, -forage.species) %>% 
  group_by(Bombus.species, family) %>% 
  summarise(interactions = sum(interactions)) %>% 
  pivot_wider(names_from = "family", values_from = "interactions", values_fill = 0) %>% 
  column_to_rownames("Bombus.species")

plant.tree <- read.newick("Qian.Zhang.2014.JSE.tre")
# plot(plant.tree)
unique(wood$family)
# Are all plant families in the Wood et al. dataset represented in the
# Qian & Zhang tree? -- Yes.
all((colnames(wood.fam) %in% plant.tree$tip.label))
# Note: the picante::pd does not require a particular order of taxa 
# in either the table or tree

wood.div <- wood.div %>% 
  mutate(faith.pd = pd(wood.fam, plant.tree)$PD) 

wood.div <- wood.div %>% 
  select(species, tongue.class, sample.size, wood.dbs = dietary.breadth, 
         richness, shannon, simpson, faith.pd)

# Quick look at correlation of Faith's PD and Wood's DBS
plot(wood.div$faith.pd, wood.div$wood.dbs, type = 'n', xlab = "Faith's phylogenetic diversity", ylab = "Wood et al.'s dietary breadth")
text(wood.div$faith.pd, wood.div$wood.dbs, wood.div$species, cex=0.75)
abline(lm(wood.div$wood.dbs ~ wood.div$faith.pd))
cor.test(wood.div$faith.pd, wood.div$wood.dbs, method = "spearman")
# rho = 0.4195804; p-value = 0.1766

# write.csv(wood.div,"Wood.et.al.2019.forage.diversity.csv", 
#           row.names = FALSE, quote = FALSE)

# Proceed into preliminary.forage.plant.analysis.R