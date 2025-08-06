# Bombus mouthpart morphometry analysis

# # Clear memory
# rm(list=ls())
# gc()

# Load packages
# devtools::install_github("aphanotus/borealis")
x <- c("borealis","tidyverse","magrittr","ggpubr","ggrepel","factoextra","viridis",
       "multcompView","ggbeeswarm","RRPP","phytools","nlme","hypervolume")
invisible(lapply(x, require, character.only = TRUE))

pairwise.group.comparisons <- function (pairwise.summary.output) {
  x <- pairwise.summary.output$summary.table[,4]
  names(x) <- sub(":","-",rownames(pairwise.summary.output$summary.table))
 return(multcompLetters(x))
}

####################
# Import data
####################

# # Convert file
# gnathos <- read.mmm(
#   input.filename = 'Bombus.mouthparts.rawdata.csv',
#   output.filename = 'Bombus.mouthpart.measures.csv',
#   metadata.cols = "all",
#   measurement.col = NULL,
#   apply.scale = TRUE,
#   invert.scale = TRUE
# )

# names(gnathos)
# head(gnathos$x)
# gnathos$measurement.number
# gnathos$specimen.number
# gnathos$specimens.missing.scale
# gnathos <- gnathos$x

gnathos <- read.csv('Bombus.mouthpart.measures.csv') %>%
  filter(caste == "W") %>% # workers only
  mutate(
    species = if_else(species.det!="",species.det,provisional.species),
    tongue = m01 + m03,    # Calculate tongue length
    body = body.length
  ) %>%
  select(-digitizer, -caste, -body.length, -species.det, -provisional.species, -scale, -name, -description)

# head(gnathos)

# sapply(gnathos, class)

data.cols <- grep("^m",(names(gnathos)))

# table(gnathos$species) 

# Filter out species with sample sizes (<5)
gnathos <- gnathos %>% group_by(species) %>% filter(n() > 5)
barplot(sort(c(with(gnathos, by(species, species, length)))), cex.names = 0.5)

####################
# Which body size proxy? ITS or body length?
####################

{
  gnathos %>% 
    ggplot(aes(its)) +
    facet_grid(species~.) +
    geom_histogram()
  
  gnathos %>% 
    ggplot(aes(body)) +
    facet_grid(species~.) +
    geom_histogram()
  
  aov(its ~ species, data = gnathos) %>% summary()
  aov(body ~ species, data = gnathos) %>% summary()
  
  # coefficients of variation
  gnathos %>% 
    select(species, its, body) %>% 
    pivot_longer(2:3) %>% 
    group_by(name) %>% 
    summarise(
      cv = sd(value, na.rm = TRUE)/mean(value, na.rm = TRUE)
    ) 
  
  # coefficients of variation by species
  gnathos %>% 
    select(species, its, body) %>% 
    group_by(species) %>% 
    pivot_longer(2:3) %>% 
    group_by(species, name) %>% 
    summarise(
      n = n(),
      value = sd(value, na.rm = TRUE)/mean(value, na.rm = TRUE)
    ) %>% 
    pivot_wider(names_from = "name") %>% 
    mutate(lower.cv = if_else(body < its,"body","its"))
  
  # ITS generally looks like it is more consistent by species  
}



####################
# Scaling plots
####################

tongue.scaling <- 
  scaling.plot(
    x = log10(gnathos$its), y = log10(gnathos$tongue),
    group = gnathos$species,
    xlab = "intertegular span (log10 mm)", ylab = "tongue length (log10 mm)",
    include.legend = FALSE,
    isometry.line = FALSE,
    convex.hulls = FALSE,
    groups.trendlines = TRUE
  )
tongue.scaling$plot

tongue.scaling$slopes
#          group  n    slope        p sig   ci.lo ci.hi spans.zero
# 1     borealis 20  0.81900 3.28e-02 *    0.0748 1.560      FALSE
# 2 griseocollis 14  0.59100 6.98e-02 .   -0.0559 1.240       TRUE
# 3  bimaculatus 70  0.58300 8.14e-09 ***  0.4060 0.760      FALSE
# 4       vagans 48  0.48900 1.56e-02 *    0.0968 0.881      FALSE
# 5    terricola 24  0.23100 5.17e-01     -0.4960 0.959       TRUE
# 6    impatiens 57  0.12200 3.66e-01     -0.1460 0.390       TRUE
# 7    ternarius 61  0.00919 9.40e-01     -0.2330 0.251       TRUE
# 8    perplexus  8 -0.54700 5.21e-01     -2.5100 1.420       TRUE

tongue.scaling$slopes$group <- factor(tongue.scaling$slopes$group, levels = tongue.scaling$slopes$group)

ylimits <- sort(tongue.scaling$slopes$ci.lo)[1]*1.1
ylimits <- c(ylimits, sort(tongue.scaling$slopes$ci.hi, decreasing = TRUE)[1])

tongue.scaling.slopes.plot <- tongue.scaling$slopes %>% 
  arrange(as.character(group)) %>% 
  mutate(simpson = SDI.by.species) %>% 
  mutate(group = paste0(group,"\n(n=",n,")")) %>% 
  mutate(group = fct_reorder(group, simpson)) %>%
  ggplot(aes(x=group, y=slope, fill=group)) +
  theme_bw() +
  theme(
    legend.position="none",
    axis.ticks=element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face="italic")
  ) +
  geom_hline(yintercept = c(0,1), color = "black") +
  geom_bar(stat="identity", color="grey15") +
  geom_errorbar(aes(ymin=ci.lo, ymax=ci.hi), width=.2) +
  scale_fill_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  scale_y_continuous(name="tongue\nscaling coefficent", limits=ylimits) +
  scale_x_discrete(name = NULL) +
  coord_flip()
tongue.scaling.slopes.plot


####################
# Explore which other mouthpart measurements are informative
####################

# How many specimens have each measurement?
(x <- apply(gnathos[,data.cols], 2, function(x) {sum(!(is.na(x)))}))
sort(apply(gnathos[,data.cols], 2, function(x) {sum(!(is.na(x)))}))

df <- data.frame(
  intact.specimens = x
)

# complete specimens
i <- which(!apply(gnathos[,data.cols], 1, function(x){any(is.na(x))}))
complete.gnathos <- gnathos[i,]
dim(complete.gnathos)

# PCA
gnathos.pca <- prcomp(complete.gnathos[,data.cols], center = TRUE, scale = TRUE)

shape.space(
  gnathos.pca, 
  group = complete.gnathos$species, 
  convex.hulls = TRUE,
  save.as = "plots/mouthpart.morphospace.pdf"
)

fviz_screeplot(gnathos.pca, addlabels = TRUE)
# PCA1 accounts for the vast majority of variance

# Order of measurements contributing to PC1-2
(x <- fviz_cos2(gnathos.pca, choice = "var", axes = 1:2))

# From largest to smallest contribution
# m13	Galea length
# m11	Stipes length
# m06	Labial palpomere 1 length
# m03	Glossa length
# m01	Mentum length
# 
# m14	Galea width
# m08	Labial palpomere 2 length
# m09	Labial palpomere 2 width
# m04	Glossa width, at mid-point
# m07	Labial palpomere 1 width
# 
# m15	Stipes width, across from palp
# m02	Mentum width
# m10	Distal labial palpomeres length
# m12	Maxillary palp length
# m05	Flabellum length

df$cos2 <- x$data$cos2

{
  plot(x=df$intact.specimens, y=df$cos2, type = "n")
  text(x=df$intact.specimens, y=df$cos2, labels = rownames(df))
}

x <- df$intact.specimens*df$cos2
names(x) <- rownames(df)
rev(sort(x))

informative.traits <- c("m13","m11","m06","m01","m03","m14","m08","m09","m04","m07")
informative.cols <- which(names(gnathos) %in% informative.traits)

# Repeat the PCA using only the 10 most informative characters
i.complete.gnathos <- which(!apply(gnathos[,informative.cols], 1, function(x){any(is.na(x))}))
length(i.complete.gnathos)
gnathos.pca <- prcomp(gnathos[i.complete.gnathos,informative.cols], center = TRUE, scale = TRUE)

shape.space(
  gnathos.pca, 
  group = gnathos$species[i.complete.gnathos], 
  convex.hulls = TRUE,
  save.as = "plots/mouthpart.morphospace.pdf"
)

# Record the PC1-2 values
gnathos$PC1 <- NA
gnathos$PC2 <- NA
gnathos$PC1[i.complete.gnathos] <- get_pca_ind(gnathos.pca)$coord[,1]
gnathos$PC2[i.complete.gnathos] <- get_pca_ind(gnathos.pca)$coord[,2]

# Which anatomical features load on the PC axes?
fviz_cos2(gnathos.pca, choice = "var", axes = 1)
fviz_cos2(gnathos.pca, choice = "var", axes = 2)


# A quick look at how shape (as PC1) varies by body size (ITS) and species
# Remember, this is just a univariate representation of mouthpart "shape". 
# This issue is revisited below with multivariate stats.
PC1.scaling <- scaling.plot(
  x = log10(gnathos$its), y = -gnathos$PC1,
  group = gnathos$species,
  xlab = "intertegular span (log10 mm)", 
  ylab = "PC1 for mouthpart measures",
  include.legend = TRUE,
  fixed.aspect = FALSE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  save.as = "plots/mouhtpart.PC1.v.its.pdf"
)
PC1.scaling$plot

species.order <- gnathos %>%
  group_by(species) %>%
  summarise(medianPC1 = median(PC1, na.rm = TRUE)) %>%
  arrange(medianPC1) %>% pull(species)
sp <- factor(gnathos$species, levels = species.order)

PC1.scaling$slopes
#          group  n slope        p sig   ci.lo ci.hi spans.zero
# 1       vagans 44 31.30 1.75e-07 ***  21.200 41.40      FALSE
# 2     borealis 19 29.50 2.00e-04 ***  16.300 42.60      FALSE
# 3  bimaculatus 67 26.40 3.05e-13 ***  20.600 32.10      FALSE
# 4 griseocollis 14 21.70 4.19e-02 *     0.935 42.40      FALSE
# 5    impatiens 57 21.70 1.20e-10 ***  16.200 27.20      FALSE
# 6    terricola 23  6.58 2.18e-01      -4.190 17.30       TRUE
# 7    ternarius 61  3.78 1.20e-01      -1.010  8.56       TRUE
# 8    perplexus  8 -3.92 7.21e-01     -29.600 21.70       TRUE

PC1.scaling$slopes$group <- factor(PC1.scaling$slopes$group, levels = PC1.scaling$slopes$group)

ylimits <- sort(PC1.scaling$slopes$ci.lo)[1]*1.1
ylimits <- c(ylimits, sort(PC1.scaling$slopes$ci.hi, decreasing = TRUE)[1])

scaling.PC1.slopes.plot <- PC1.scaling$slopes %>% 
  arrange(as.character(group)) %>% 
  mutate(simpson = SDI.by.species) %>% 
  mutate(group = paste0(group,"\n(n=",n,")")) %>% 
  mutate(group = fct_reorder(group, simpson)) %>%
  ggplot(aes(x=group, y=slope, fill=group, label=paste0("n=",n))) +
  theme_bw() +
  theme(
    legend.position="none",
    axis.ticks=element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face="italic")
  ) +
  geom_bar(stat="identity", color="grey15") +
  geom_errorbar(aes(ymin=ci.lo, ymax=ci.hi), width=.2) +
  scale_fill_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  scale_y_continuous(limits=ylimits) +
  xlab("species") +
  ylab("mouthpart measures\nscaling coefficient") +
  coord_flip()
scaling.PC1.slopes.plot



####################
# Calculate species morphospace hypervolumes
####################

# Cumulative variance
cumsum(gnathos.pca$sdev^2 / sum(gnathos.pca$sdev^2))
# 90% variance is covered by PC1-6

gnathos.hypervolumes <- data.frame(
  species = unique(gnathos$species),
  gnathos.hv = NA
)

{
  start.time <- Sys.time()
  for (i in 1:length(gnathos.hypervolumes$species)) {
    j <-  which(gnathos$species[i.complete.gnathos] == gnathos.hypervolumes$species[i])
    hv <- hypervolume_gaussian(gnathos.pca$x[j, 1:6])
    gnathos.hypervolumes$gnathos.hv[i] <- hv@Volume
  }
  x <- Sys.time() - start.time
  print(x)
  gnathos.hypervolumes 
}
# About 30 seconds

save(
  gnathos, gnathos.pca, gnathos.hypervolumes, i.complete.gnathos,
  file = "for.gnathos.hv.rda"
)

load("for.gnathos.hv.rda", verbose = TRUE)

# It's wise to run the following on an HOC cluster.
# It can take a while!
############

# Pairwise comparisons
trait_cols <- 1:6
iterations <- 1000
library(hypervolume)
{
  start.time <- Sys.time()
  species_list <- sort(unique(gnathos$species[i.complete.gnathos]))
  hv_list <- lapply(species_list, function (sp) {
    hypervolume_gaussian(
      gnathos.pca$x[which(gnathos$species[i.complete.gnathos] == sp), trait_cols], 
      name = sp
    )
  })
  names(hv_list) <- species_list
  
  results <- list()
  for (i in 1:(length(hv_list)-1)) {
    for (j in (i+1):length(hv_list)) {
      name_i <- names(hv_list)[i]
      name_j <- names(hv_list)[j]
      
      # Permutation test for overlap
      perm_dir <- hypervolume_permute(
        paste0("gnathos_perm_", name_i, "_", name_j),
        hv1 = hv_list[[i]],
        hv2 = hv_list[[j]],
        n = iterations
      )
      
      test <- hypervolume_overlap_test(hv_list[[i]], hv_list[[j]], path = perm_dir)
      results[[paste0(name_i, "_vs_", name_j)]] <- test$p_values
    }
  }
  gnathos.hv.pairwise <- results
  cat("\nDone!\nFinished permutation for **gnathos** hypervolume pairwise contrasts.\n")
  cat("\nThe following constrasts have significant differences:\n")
  print(names(unlist(gnathos.hv.pairwise))[which(unlist(gnathos.hv.pairwise) < 0.05)])
  save(gnathos.hv.pairwise, file = "gnathos.hv.pairwise.rda")
  cat("\nSaving results to file.\n")
  x <- Sys.time() - start.time
  print(x)
}
# Done!  Time difference of 1.877321 days

############

####################
# Model building using RRPP
####################

# Exclude specimens with missing ITS measurements
complete.gnathos <- gnathos[which(!is.na(gnathos$its)),]
# Exclude specimens with missing measurements in the informative traits
i <- which(!apply(complete.gnathos[,informative.cols], 1, function(x){any(is.na(x))}))
complete.gnathos <- complete.gnathos[i,]

gnathos.list <- list(
  mm = as.matrix(complete.gnathos[,informative.cols]),
  its = complete.gnathos$its,
  body = complete.gnathos$body,
  species = factor(complete.gnathos$species, levels = species.order)
)

# Is size (ITS) a predictor of shape (as the multivariate collection of 10 measurements)?
model.gnathos.size <- lm.rrpp(
  mm ~ log(its), 
  data = gnathos.list, 
  iter = (1e4)-1
)
anova(model.gnathos.size)
#            Df     SS      MS     Rsq      F      Z Pr(>F)    
# log(its)    1 164.76 164.758 0.25628 100.28 5.1544  1e-04 ***
# Residuals 291 478.12   1.643 0.74372                         
# Total     292 642.87          

# Is species also a predictor of shape (as the multivariate collection of 10 measurements)?
model.gnathos.species <- lm.rrpp(
  mm ~ log(its) + species, 
  data = gnathos.list, 
  iter = (1e4)-1
)
anova(model.gnathos.species)
#            Df     SS      MS     Rsq       F       Z Pr(>F)    
# log(its)    1 164.76 164.758 0.25628 177.253  5.8249  1e-04 ***
# species     7 214.13  30.591 0.33309  32.911 10.8955  1e-04 ***
# Residuals 284 263.98   0.930 0.41063                           
# Total     292 642.87 
                                     

# Pairwise species differences among species?
pw.gnathos.species <- pairwise(
  fit = model.gnathos.species,
  fit.null = model.gnathos.size,
  groups = gnathos.list$species
)
pw.gnathos.species.output <- summary(pw.gnathos.species)
(sig.letters <- pairwise.group.comparisons(pw.gnathos.species.output))
#  borealis  bimaculatus    perplexus griseocollis    impatiens    terricola       vagans    ternarius 
#       "a"          "b"         "bc"          "d"          "c"          "e"          "c"          "e" 

df <- data.frame(
  x = 1:8,
  y = rep(6.5,8),
  label = sig.letters$Letters
)

species.order <- gnathos %>%
  group_by(species) %>%
  summarise(medianPC1 = median(PC1, na.rm = TRUE)) %>%
  arrange(medianPC1) %>% pull(species)
sp <- factor(gnathos$species, levels = species.order)

PC1.by.species.plot <- gnathos %>% 
  ggplot(aes(sp, PC1, color = sp)) +
  theme_bw() +
  theme(
    legend.position="none",
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=30, vjust=1, hjust=1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face="italic")
  ) +
  geom_violin(scale = "area") +
  geom_quasirandom(varwidth = TRUE) +
  scale_color_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(data = df, aes(x = x, y = y, label = label), color = "black", size = 5 ) +
  labs(
    # caption = "Significant differences in the influence of species in the model  Y ~ intertegular span + species.\nPERMANOVA using `lm.rrpp`; contrasts using `pairwise` in the package `RRPP`."
    x = NULL, #"species",
    y = "mouthpart morphology (PC1)"
  ) 

PC1.by.species.plot


# Allometry-corrected PCA
gnathos.allo.pca <- prcomp(resid(model.gnathos.size), center = TRUE, scale = TRUE)
shape.space(
  gnathos.allo.pca, 
  group = gnathos.list$species, 
  convex.hulls = TRUE,
  main.title = "Size-corrected morphospace for mouthparts",
  save.as = "plots/allo-corrected.pca.pdf"
)

# Exclude specimens with missing ITS measurements
complete.gnathos <- gnathos[which(!is.na(gnathos$its)),]
# Exclude specimens with missing measurements in the informative traits
i <- which(!apply(complete.gnathos[,informative.cols], 1, function(x){any(is.na(x))}))
complete.gnathos <- complete.gnathos[i,]

i <- which(
  !is.na(gnathos$its) &
  !apply(gnathos[,informative.cols], 1, function(x){any(is.na(x))})
)

gnathos$alloPC1 <- NA
gnathos$alloPC1[i] <- get_pca_ind(gnathos.allo.pca)$coord[,1]

species.order <- gnathos %>%
  group_by(species) %>%
  summarise(medianAlloPC1 = median(alloPC1, na.rm = TRUE)) %>%
  arrange(medianAlloPC1) %>% pull(species)
gnathos$species <- factor(gnathos$species, levels = species.order)

# Repeat this just to get the species in the right order!
gnathos.list <- list(
  specimen.id = gnathos$specimen.id[i],
  mm = as.matrix(gnathos[i,informative.cols]),
  body = gnathos$body[i],
  its = gnathos$its[i],
  species = gnathos$species[i],
  tongue = gnathos$tongue[i],
  PC1 = gnathos$PC1[i],
  PC2 = gnathos$PC2[i],
  alloPC1 = gnathos$alloPC1[i]
)
model.gnathos.size <- lm.rrpp(mm ~ log(its), data = gnathos.list, iter = (1e4)-1)
model.gnathos.species <- lm.rrpp(mm ~ log(its) + species, data = gnathos.list, iter = (1e4)-1)
pw.gnathos.species <- pairwise(
  fit = model.gnathos.species,
  fit.null = model.gnathos.size,
  groups = gnathos.list$species
)
pw.output <- summary(pw.gnathos.species)
(sig.letters <- pairwise.group.comparisons(pw.output))
# borealis  bimaculatus griseocollis    perplexus       vagans    impatiens    ternarius    terricola 
#      "a"          "b"          "d"         "bc"          "c"          "c"          "e"          "e" 

df <- data.frame(
  x = 1:8,
  y = rep(9.5, 8),
  label = sig.letters$Letters
)

ggplot(gnathos, aes(species, alloPC1, color = species)) +
  theme_minimal() +
  theme(legend.position="none") +
  geom_violin() +
  geom_quasirandom(varwidth = TRUE, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(data = df, aes(x = x, y = y, label = label), color = "black", size = 5 ) +
  xlab("species") +
  ylab("size-corrected PC1 for mouthpart measures") +
  labs(caption = "Contrasts based on species factor in the model  Y ~ intertegular span + species.\nPERMANOVA using `RRPP::lm.rrpp`; contrasts using `pairwise`.")
ggsave("plots/size-corrected.mouthpart.PC1.by.species.png", width = 6, height = 4, scale = 1.25)


# Do species have different mouthpart "shape" allometries?
model.unique.allometries <- lm.rrpp(mm ~ log(its) * species, data = gnathos.list, iter = (1e4)-1)
anova(model.unique.allometries)
#                   Df     SS      MS     Rsq        F       Z Pr(>F)    
# log(its)           1 164.76 164.758 0.25628 188.1091  5.8971  1e-04 ***
# species            7 214.13  30.591 0.33309  34.9263 11.1214  1e-04 ***
# log(its):species   7  21.37   3.052 0.03323   3.4849  3.1396  8e-04 ***
# Residuals        277 242.61   0.876 0.37739                            
# Total            292 642.87                                               
# Yes, species have different worker mouthpart "shape" allometries

# Post hoc examination of the interaction
pw.unique.allometries <- pairwise(
  fit = model.unique.allometries,
  fit.null = model.gnathos.species,
  groups = gnathos.list$species,
  covariate = gnathos.list$its
)
pw.allometry.output <- summary(pw.unique.allometries)
(sig.letters <- pairwise.group.comparisons(pw.allometry.output))
# borealis  bimaculatus griseocollis    perplexus       vagans    impatiens    ternarius    terricola 
#      "a"         "ab"          "a"         "bc"         "ab"          "c"          "c"         "bc" 

# I confirmed with Michael Collyer that this modeling approach and use of 
# rrpp::pairwise accomplishes what we're after here.


####################
# PGLS
####################
# Phylogenetic regression using `nlme::gls`

gnathos.species.names <- unique(sort(as.character(gnathos.list$species)))

btree <- Bombus.tree
btree$tip.label <- str_split_fixed(btree$tip.label," ",2)[,2]

btree <- keep.tip(btree, btree$tip.label[which(btree$tip.label %in% gnathos.species.names)])
plot(btree)

# Create a dataframe with species means

gnathos.sp <- gnathos %>% 
  mutate(species = as.character(species)) %>% 
  group_by(species) %>% 
  summarise(
    PC1 = mean(PC1, na.rm=TRUE),
    PC2 = mean(PC2, na.rm=TRUE),
    alloPC1 = mean(alloPC1, na.rm=TRUE),
    its = mean(its, na.rm=TRUE),
    body = mean(body, na.rm=TRUE)
  ) %>% 
  mutate(simpson = SDI.by.species) %>% 
  # Reorder the dataframe to match the order of tree tips
  mutate(species = fct_reorder(species, match(species, btree$tip.label))) %>% 
  arrange(species) %>% 
  column_to_rownames(var = "species")

# Check that the names in the tree and dataset match
geiger::name.check(btree, gnathos.sp)

pgls.PC1.by.ITS <- gls(
  PC1 ~ log(its), 
  correlation = corBrownian(phy = btree), 
  method = "REML", data = gnathos.sp
)
summary(pgls.PC1.by.ITS)
# Generalized least squares fit by REML
#      AIC      BIC    logLik
# 25.01614 24.39142 -9.508069
# 
#                 Value Std.Error   t-value p-value
# (Intercept)  28.74374  7.732183  3.717416  0.0099
# log(its)    -19.38564  4.939329 -3.924752  0.0078
# Residual standard error: 1.657342
# Degrees of freedom: 8 total; 6 residual

with(gnathos.sp, plot(PC1 ~ log(its)))
abline(a = coef(pgls.PC1.by.ITS)[1], b = coef(pgls.PC1.by.ITS)[2])

# So, size is a predictor of mouthpart shape (or at least PC1 for
# these 10 measurements), even after accounting for relatedness


####################
# Correlations to forage diversity
####################

forage <- read.csv("forage.plant.diversity.csv") %>% 
  mutate(
    species = sub("Bombus ","",species),
    simpson = Simpsons
  ) %>% 
  mutate(to.rownames = species) %>% 
  column_to_rownames(var = "to.rownames") %>%
  select(species, richness, simpson, MDS1, MDS2) %T>% 
  print()

# Add forage diversity metrics to the gnathos list object
{
  gnathos.list$richness <- NA
  gnathos.list$simpson <- NA
  gnathos.list$MDS1 <- NA
  gnathos.list$MDS2 <- NA
}
for (i in 1:length(forage$species)) {
  j <- which(gnathos.list$species == forage$species[i])
  gnathos.list$richness[j] <- forage[forage$species[i],"richness"]
  gnathos.list$simpson[j] <- forage[forage$species[i],"simpson"]
  gnathos.list$MDS1[j] <- forage[forage$species[i],"MDS1"]
  gnathos.list$MDS2[j] <- forage[forage$species[i],"MDS2"]
}

# Explore potential correlations
gnathos.list %>%
  keep(names(.) %in% c("PC1", "PC2", "richness", "simpson", "MDS1", "MDS2")) %>%
  bind_rows() %>% 
  borealis::pairs()
# Both mouthparts shape PCs have reasonable correlation to the diversity metrics
# but not to the forage plant community composition dimensions (NMDS1-2) 

plot(gnathos.list$simpson, gnathos.list$PC1 )
abline(lm(gnathos.list$PC1 ~ gnathos.list$simpson))


####################
# Model the effect of foraging metrics on mouthpart shape
####################
# Comparisons among these models can be made based on Z values (effect size)
# Species should be a strong predictor, but the forage plant metric with the
# next best effect will be informative.

model.gnathos.simpson <- lm.rrpp(
  mm ~ simpson, 
  data = gnathos.list, 
  iter = 1e4-1, 
  print.progress = TRUE
)
anova(model.gnathos.simpson)
#            Df     SS     MS     Rsq      F      Z Pr(>F)    
# simpson     1  85.03 85.033 0.13227 44.358 4.1471  1e-04 ***
# Residuals 291 557.84  1.917 0.86773                         
# Total     292 642.87             

model.gnathos.richness <- lm.rrpp(
  mm ~ richness, 
  data = gnathos.list, 
  iter = 1e4-1, 
  print.progress = TRUE
)
anova(model.gnathos.richness)
#            Df     SS     MS     Rsq      F      Z Pr(>F)    
# richness    1  57.88 57.877 0.09003 28.79 3.6546  1e-04 ***
# Residuals 291 585.00  2.010 0.90997                        
# Total     292 642.87  

model.gnathos.MDS1 <- lm.rrpp(
  mm ~ MDS1, 
  data = gnathos.list, 
  iter = 1e4-1, 
  print.progress = TRUE
)
anova(model.gnathos.MDS1)
#            Df     SS     MS     Rsq      F      Z Pr(>F)
# MDS1        1   5.88 5.8822 0.00915 2.6872 1.377 0.0894 .
# Residuals 291 636.99 2.1890 0.99085                      
# Total     292 642.87

model.gnathos.MDS2 <- lm.rrpp(
  mm ~ MDS2, 
  data = gnathos.list, 
  iter = 1e4-1, 
  print.progress = TRUE
)
anova(model.gnathos.MDS2)
#            Df     SS      MS     Rsq      F      Z Pr(>F)   
# MDS2        1  13.12 13.1219 0.02041 6.0634 2.1203 0.0106 *
# Residuals 291 629.75  2.1641 0.97959                       
# Total     292 642.87                             

# How do these factors do in models with ITS?
model.gnathos.size.simpson <- lm.rrpp(
  mm ~ log(its) + simpson, 
  data = gnathos.list, 
  iter = 1e4-1, 
  print.progress = FALSE
)

model.gnathos.size.simpson.results <- 
  anova(model.gnathos.size.simpson)$table %T>%
  print() %>% 
  mutate(Rsq = signif(Rsq,3))
#            Df     SS      MS     Rsq       F      Z Pr(>F)    
# log(its)    1 164.76 164.758 0.25628 111.455 5.2770  1e-04 ***
# simpson     1  49.42  49.422 0.07688  33.433 3.9265  1e-04 ***
# Residuals 290 428.69   1.478 0.66684                          
# Total     292 642.87            

# Simpson's index still captures a large fraction of variance, 
# even when accounting for size.

# PGLS

pgls.PC1.by.SDI <- gls(
  PC1 ~ simpson, 
  correlation = corBrownian(phy = btree), 
  method = "REML", data = gnathos.sp
)
summary(pgls.PC1.by.SDI)
#      AIC     BIC    logLik
# 29.17682 28.5521 -11.58841
#                 Value Std.Error    t-value p-value
# (Intercept)  5.125239  49.23191  0.1041040  0.9205
# simpson     -6.947631  52.34936 -0.1327166  0.8988
# Residual standard error: 3.125678 
# Degrees of freedom: 8 total; 6 residual

pgls.PC1.by.ITS.SDI <- gls(
  PC1 ~ log(its) + simpson, 
  correlation = corBrownian(phy = btree), 
  method = "REML", data = gnathos.sp
)
summary(pgls.PC1.by.ITS.SDI)
#      AIC      BIC    logLik
# 18.19737 16.63512 -5.098686
#                 Value Std.Error   t-value p-value
# (Intercept)  16.29210  28.16697  0.578412  0.5881
# log(its)    -19.85468   5.39472 -3.680390  0.0143 *
# simpson      14.02374  30.31659  0.462577  0.6631
# Residual standard error: 1.777883 
# Degrees of freedom: 8 total; 5 residual

# A plot to show the relationship of Simpson's index and mouthpart morphology
df <- data.frame(
  species = gnathos.list$species,
  simpson = gnathos.list$simpson,
  PC1 = gnathos.list$PC1,
  alloPC1 = gnathos.list$alloPC1
)

species.order <- forage$species[order(forage$simpson)]
df$species <- factor(df$species, levels = species.order)

txt.df <- df %>%
  group_by(species) %>%
  summarise(simpson = median(simpson),
            .groups = "drop") %>%
  mutate(
    #         gris   bor    vag  terri   perp  tern  bimac imp
    # medians -0.98 -4.06 0.396  0.727  -1.18  2.52  -1.29 0.222
    y     = c(-3.6, -4.95,  0.3,  1.06, -1.15,  4.5, -1.20, 3.10),
    hjust = c( 0.1,  0.5,   1.2,  1.05,   1.1,  0.5, -0.10, 0.9)
  )

specialization.plot <- df %>%
  ggplot(aes(x = simpson, y = PC1)) +
  theme_bw() +
  theme(legend.position="none") +
  geom_smooth(method=lm, color = "grey40", fill = "gray85") +
  geom_point(aes(color = species), alpha = 0.65, size = 2) +
  scale_color_viridis(discrete = TRUE, end = 0.95) +
  geom_text(data = txt.df, aes(x=simpson, y=y, hjust = hjust, label = species, color = species)) +
#   annotate("text", x=0.912, y=6.00, hjust = 0, vjust = 1,
#            color="grey40", size = 2.5,
#            label=paste0("ANOVA with residual randomization
# 10 linear measurements of worker mouthparts
# Y ~ log(individual ITS) + species SDI for forage plants
#   ITS:  R^2 = ",model.gnathos.size.simpson.results$Rsq[1],", p < ",model.gnathos.size.simpson.results$`Pr(>F)`[1],"
#   SDI:  R^2 = ",model.gnathos.size.simpson.results$Rsq[2],", p < ",model.gnathos.size.simpson.results$`Pr(>F)`[2])) +
  labs(x="foraging niche breadth (SDI)",
       y="mouthpart morphology (PC1)")

specialization.plot

figure2 <- ggarrange(
  PC1.by.species.plot, specialization.plot,
  nrow = 1, labels = "AUTO"
)
ggsave("manuscript/Fig2.gnathos.by.PC1.and.SDI.png", figure2, width = 6.5, height = 3, scale = 1.65)



# Save everything
save.image("bombus.scaling.rda")
