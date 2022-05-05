# Preliminary forage plant analysis

library(tidyverse)
library(ggpubr)
library(viridis)
library(factoextra)

# Interaction matrix based on data from Wood et al. 2019 Ecology, Figure 3
wood <- read.csv("Wood.et.al.2019.Fig3.interaction.matrix.csv")
names(wood)

# Tongue-length categories for each species
tongue.class <- c("short","long","medium","long","long","short","medium","long","medium","medium","short","long-medium")
names(tongue.class) <- c("affinis","auricomus","bimaculatus","borealis","fervidus","griseocollis","impatiens","pensylvanicus","perplexus","ternarius","terricola","vagans")

# Filter the Wood et al. dataset down to Bombus species included by Lee et al.
Bombus.sp.included.by.Lee <- c("bimaculatus","borealis","fervidus","impatiens","ternarius","terricola","vagans")
x <- which(colnames(wood) %in% Bombus.sp.included.by.Lee)
wood <- wood[,c(1:3,x)]
x <- which(names(tongue.class) %in% Bombus.sp.included.by.Lee)
tongue.class <- tongue.class[x]
  
# Filter out forage plants for which fewer than 3 interactions were observed
x <- which(rowSums(wood[,-c(1:3)]) < 3)
wood <- wood[-x,]

# Filter out "other" forage plants
x <- which(wood$forage.species == "others")
wood <- wood[-x,]
dim(wood)

# Sample sizes for each Bombus species
(wood.bee.sample.sizes <- wood %>% 
  select(-forage.species, -common.name, -plant.family) %>% 
  colSums())

# Sample sizes for each forage plant species
(wood.plant.sample.sizes <- wood %>% 
  tibble() %>% 
  select(-common.name, -plant.family) %>% 
    column_to_rownames(var = "forage.species") %>% 
  rowSums())

# PCA for the interaction matrix
wood.pca <- prcomp(t(wood[,-c(1:3)]), scale = TRUE)
rownames(wood.pca$rotation) <- wood$forage.species

# Proportion of variance explained by each PC axis 
# Extract eigenvalues and variances for first 5 dimensions
get_eig(wood.pca)[1:5,]
fviz_screeplot(wood.pca, addlabels = TRUE)

# Bombus species positions in the foraging-space
Bombus.sp.foraging.space.plot <- fviz_pca_ind(wood.pca, 
                   col.ind = as.factor(tongue.class), invisible = "quali",
                   repel = TRUE, max.overlaps = 10,
                   title = "") + # PCA for Bombus species forage plants
  scale_color_viridis(discrete = TRUE, option = "magma", end = 0.75) +
  guides(color = guide_legend(title = "tongue length", override.aes = aes(label = "") ) ) +
  scale_shape_manual(name = "tongue length", values = rep(19,4)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))
Bombus.sp.foraging.space.plot$labels$x <- sub("Dim","PC",Bombus.sp.foraging.space.plot$labels$x)
Bombus.sp.foraging.space.plot$labels$y <- sub("Dim","PC",Bombus.sp.foraging.space.plot$labels$y)
Bombus.sp.foraging.space.plot

# Correlation between the interaction values for each forage plant 
# and the PC dimensions 
get_pca_var(wood.pca)$cor %>% 
  as.data.frame() %>% 
  mutate(common.name = wood$common.name,    # Add common names
         plant.family = wood$plant.family,  # Add plant family 
         n = wood.plant.sample.sizes) %>%   # Add number of interactions
  select(common.name, plant.family, n, 
         Dim.1, Dim.2) %>%                  # Just show axes 1-2
  arrange(desc(Dim.1))                      # Sort by PC1

# Plot contribution of each forage plant (loadings), colored by plant family
loading.plot <- fviz_pca_biplot(wood.pca, 
                      col.ind = "gray35",
                      col.var = wood$plant.family,
                      title = "") +
  scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.75) +
  guides(color = guide_legend(title = "plant family", override.aes = aes(label = "") ) )
loading.plot$labels$x <- sub("Dim","PC",loading.plot$labels$x)
loading.plot$labels$y <- sub("Dim","PC",loading.plot$labels$y)
loading.plot

wood.pca.plots <- ggarrange(Bombus.sp.foraging.space.plot, loading.plot + rremove("y.axis"), 
                            labels = "AUTO", widths = c(4,5))
ggsave("plots/wood.pca.plots.pdf", wood.pca.plots, width = 12, height = 5)

