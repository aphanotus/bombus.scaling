# Bombus leg morphometry analysis

# Clear memory
# rm(list=ls())
# gc() 

# Load packages
x <- c("borealis","tidyverse","magrittr","ggpubr","ggrepel","factoextra","viridis",
       "multcompView","ggbeeswarm","RRPP","phytools","nlme","hypervolume")
invisible(lapply(x, require, character.only = TRUE))


####################
# Import data
####################

{
  legs <- read.csv("Bombus.legs.250706.csv") %>%
    filter(caste == "W") %>% # workers only
    mutate(
      its = its.mm,
      body = body.length.mm
    ) %>%
    select(-image.initials, -caste, -its.mm, -body.length.mm, -notes)
  
  # dim(legs); names(legs)
  
  # Filter out species with sample sizes (<5)
  legs <- legs %>% group_by(species) %>% filter(n() > 5)
  barplot(sort(c(with(legs, by(species, species, length)))), cex.names = 0.5)
  # B. griseocollis is not included in this analysis
  
  data.cols <- grep("^l[123]",(names(legs)))
  
  # Average left and right sides
  long.legs <- legs %>% 
    pivot_longer(cols = all_of(data.cols), names_to = "trait") %>% 
    mutate(trait = sub("_[12]$","",trait)) %>% 
    group_by(specimen.id, species, body, its, trait) %>%
    summarise(trait.value = mean(value, na.rm=TRUE)) 
  
  # head(long.legs)
  
  # Find leg lengths (femur + tibia + basitarsus)
  legs <- long.legs %>% 
    pivot_wider(names_from = "trait", values_from = "trait.value") %>% 
    group_by(specimen.id) %>% 
    mutate(
      l1_legL = sum(c(l1_femurL, l1_tibiaL, l1_btarsL), na.rm = FALSE),
      l2_legL = sum(c(l2_femurL, l2_tibiaL, l2_btarsL), na.rm = FALSE),
      l3_legL = sum(c(l3_femurL, l3_tibiaL, l3_btarsL), na.rm = FALSE)
    ) 
  
  # head(legs)
  
  data.cols <- grep("^l[123]",(names(legs)))
  
  long.legs <- legs %>% 
    pivot_longer(cols = any_of(data.cols), names_to = "trait", values_to = "trait.value")
  tail(long.legs)
}


####################
# Scaling plots
####################

# Preliminary look at scaling of legs
long.legs %>% 
  filter(grepl("legL",trait)) %>% 
  ggplot(aes(x = log10(its), y = log10(trait.value))) +
  theme_bw() +
  facet_grid(trait ~ species) +
  geom_point(alpha = 0.35) +
  geom_smooth(method=lm, se=FALSE)
# Slopes appear to differ by species, but seem consistent among legs.

hindleg.scaling <- scaling.plot(
  x = log10(legs$its), y = log10(legs$l3_legL),
  group = legs$species,
  xlab = "intertegular span (log10 mm)", ylab = "hind leg (log10 mm)",
  include.legend = TRUE,
  isometry.line = FALSE,
  convex.hulls = FALSE,
  groups.trendlines = TRUE,
  save.as = "plots/scaling.L3leg.v.its.pdf"
)
hindleg.scaling$plot

hindleg.scaling$slopes
#         group  n    slope        p sig  ci.lo ci.hi spans.zero
# 1    borealis 21  0.95300 2.49e-05 ***  0.592 1.310      FALSE
# 2 bimaculatus 47  0.67500 3.40e-07 ***  0.448 0.903      FALSE
# 3      vagans 59  0.65800 4.27e-05 ***  0.361 0.955      FALSE
# 4   impatiens 68  0.61700 1.61e-07 ***  0.407 0.827      FALSE
# 5   terricola 18  0.56600 4.60e-03 **   0.202 0.931      FALSE
# 6   ternarius 52  0.00636 9.69e-01     -0.317 0.330       TRUE
# 7   perplexus  8 -0.20800 8.69e-01     -3.170 2.750       TRUE

hindleg.scaling$slopes$group <- factor(hindleg.scaling$slopes$group, levels = hindleg.scaling$slopes$group)

ylimits <- sort(hindleg.scaling$slopes$ci.lo)[2]*1.1
ylimits <- c(ylimits, sort(hindleg.scaling$slopes$ci.hi, decreasing = TRUE)[2])

hindleg.scaling.slopes.plot <- hindleg.scaling$slopes %>% 
  filter(group != "perplexus") %>% 
  add_row(group = "perplexus", n=0, slope = 0) %>% 
  add_row(group = "griseocollis", n=0, slope = 0) %>% 
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
  scale_y_continuous(name="hindleg\nscaling coefficent", limits=ylimits) +
  scale_x_discrete(name = NULL) +
  coord_flip()
hindleg.scaling.slopes.plot



####################
# PCA
####################

# How many specimens have each measurement?
(x <- apply(legs[,data.cols], 2, function(x) {sum(!(is.na(x)))}))
sort(apply(legs[,data.cols], 2, function(x) {sum(!(is.na(x)))}))

# Index specimens without missing data 
i.complete.legs <- which(!apply(legs[,data.cols], 1, function(x){any(is.na(x))}))
clean.legs <- as.data.frame(legs[i.complete.legs,data.cols])
rownames(clean.legs) <- legs$specimen.id[i.complete.legs]

leg.pca <- prcomp(clean.legs, center = TRUE, scale = TRUE)

shape.space(
  leg.pca, group = legs$species[i.complete.legs], 
  include.legend = TRUE,
  convex.hulls = TRUE,
  save.as = "plots/leg.morphospace.pdf"
)

legs$legPC1 <- NA
legs$legPC2 <- NA
legs$legPC1[i.complete.legs] <- get_pca_ind(leg.pca)$coord[,1]
legs$legPC2[i.complete.legs] <- get_pca_ind(leg.pca)$coord[,2]

clean.legs$body <- legs$body[i.complete.legs]
clean.legs$its <- legs$its[i.complete.legs]
clean.legs$species <- legs$species[i.complete.legs]
clean.legs$legPC1 <- legs$legPC1[i.complete.legs]
clean.legs$legPC2 <- legs$legPC2[i.complete.legs]


# A quick look at how shape (as PC1) varies by body size (ITS) and species
# Remember, this is just a univariate representation of leg "shape". 
# This issue is revisited below with multivariate stats.
leg.PC1.scaling <- scaling.plot(
  x = log10(legs$its), y = legs$legPC1,
  group = legs$species,
  xlab = "intertegular span (log10 mm)", 
  ylab = "PC1 for leg measures",
  include.legend = TRUE,
  fixed.aspect = FALSE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  save.as = "plots/leg.PC1.v.its.pdf"
)
leg.PC1.scaling$plot

ylimits <- sort(leg.PC1.scaling$slopes$ci.lo)[2]*1.1
ylimits <- c(ylimits, sort(leg.PC1.scaling$slopes$ci.hi, decreasing = TRUE)[2])

leg.scaling.PC1.slopes.plot <- leg.PC1.scaling$slopes %>% 
  filter(group != "perplexus") %>% 
  add_row(group = "perlexus", slope = 0, n=0) %>% 
  add_row(group = "grsieocollis", slope = 0, n=0) %>% 
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
  ylab("leg measures\nscaling coefficient") +
  coord_flip()
leg.scaling.PC1.slopes.plot




# A quick plot of leg "shape" (as PC1) by species
leg.sp.order <- legs %>%
  group_by(species) %>%
  summarise(medianPC1 = median(legPC1, na.rm=TRUE)) %>%
  arrange(medianPC1) %>% pull(species)
sp <- factor(legs$species, levels = leg.sp.order)
boxplot(legs$legPC1 ~ sp)


####################
# Calculate species morphospace hypervolumes
####################

# Cumulative variance
cumsum(leg.pca$sdev^2 / sum(leg.pca$sdev^2))
# 90% variance is covered by PC1-5

leg.hypervolumes <- data.frame(
  species = unique(clean.legs$species),
  leg.hv = NA
)

{
  start.time <- Sys.time()
  for (i in 1:length(leg.hypervolumes$species)) {
    j <-  which(clean.legs$species == leg.hypervolumes$species[i])
    hv <- hypervolume_gaussian(leg.pca$x[j, 1:2])
    leg.hypervolumes$leg.hv[i] <- hv@Volume
  }
  x <- Sys.time() - start.time
  print(x)
  leg.hypervolumes 
}
# Time difference of 5.482179 secs

save(
  clean.legs, leg.pca, leg.hypervolumes,
  file = "for.leg.hv.rda"
)

load("for.leg.hv.rda", verbose = TRUE)

############


# Pairwise comparisons
trait_cols <- 1:5
iterations <- 1000
{
  start.time <- Sys.time()
  species_list <- sort(unique(clean.legs$species))
  hv_list <- lapply(species_list, function (sp) {
    hypervolume_gaussian(
      leg.pca$x[which(clean.legs$species == sp), trait_cols], 
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
        paste0("leg_perm_", name_i, "_", name_j),
        hv1 = hv_list[[i]],
        hv2 = hv_list[[j]],
        n = iterations
      )
      
      test <- hypervolume_overlap_test(hv_list[[i]], hv_list[[j]], path = perm_dir)
      results[[paste0(name_i, "_vs_", name_j)]] <- test$p_values
    }
  }
  leg.hv.pairwise <- results
  
  cat("\nDone!\nFinished permutation for **leg** hypervolume pairwise contrasts.\n")
  cat("\nThe following constrasts have significant differences:\n")
  print(names(unlist(leg.hv.pairwise))[which(unlist(leg.hv.pairwise) < 0.05)])
  save(leg.hv.pairwise, file = "leg.hv.pairwise.rda")
  cat("\nSaving results to file.\n")
  x <- Sys.time() - start.time
  print(x)
}
# Done!  Time difference of 7.300066 hours



############


####################
# Model building using RRPP
####################

{
  i <- which(!is.na(clean.legs$its))
  leg.list <- list(
    specimen.id = rownames(clean.legs)[i],
    mm = as.matrix(clean.legs[i,1:12]),
    its = clean.legs$its[i],
    legPC1 = clean.legs$legPC1[i],
    legPC2 = clean.legs$legPC2[i],
    species = factor(clean.legs$species[i], levels = leg.sp.order)
  )
  i <- which(leg.list$specimen.id %in% legs$specimen.id)
  leg.list$l1_legL <- legs$l1_legL[i]
  leg.list$l2_legL <- legs$l2_legL[i]
  leg.list$l3_legL <- legs$l3_legL[i]
}

model.legs.size <- lm.rrpp(mm ~ log(its), data = leg.list, iter = 1e4-1)
anova(model.legs.size)
#            Df     SS      MS     Rsq      F      Z Pr(>F)    
# log(its)    1 129.96 129.958 0.40554 155.54 4.7335  1e-04 ***
# Residuals 228 190.50   0.836 0.59446                         
# Total     229 320.46 

model.legs.species <- lm.rrpp(mm ~ log(its) + species, data = leg.list, iter = 1e4-1)
anova(model.legs.species)
#            Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)    1 129.96 129.958 0.40554 178.2620 4.8320  1e-04 ***
# species     6  28.66   4.776 0.08943   6.5515 4.7229  1e-04 ***
# Residuals 222 161.84   0.729 0.50504                           
# Total     229 320.46       

# Pairwise species differences?
pw.legs.species <- pairwise(
  fit = model.legs.species,
  fit.null = model.legs.size,
  groups = leg.list$species
  )
pw.legs.species.output <- summary(pw.legs.species)
(sig.letters <- pairwise.group.comparisons(pw.legs.species.output))
#  ternarius   perplexus      vagans   impatiens bimaculatus   terricola    borealis 
#        "a"         "a"         "a"        "bc"        "bc"        "ab"         "c" 

df <- data.frame(
  x = 1:7,
  y = rep(11.5,7),
  label = sig.letters$Letters
)

species.order <- clean.legs %>%
  group_by(species) %>%
  summarise(medianPC1 = median(legPC1)) %>%
  arrange(medianPC1) %>% pull(species)
sp <- factor(clean.legs$species, levels = species.order)

legPC1.by.species.plot <- ggplot(clean.legs, aes(sp, legPC1, color = sp)) +
  theme_minimal() +
  theme(legend.position="none") +
  geom_violin() +
  geom_quasirandom(varwidth = TRUE) +
  scale_color_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(data = df, aes(x = x, y = y, label = label), color = "black", size = 5 ) +
  xlab("species") +
  ylab("PC1 for leg measurements") +
  labs(caption = "Significant differences in the influence of species in the model  Y ~ intertegular span + species.\nPERMANOVA using `lm.rrpp`; contrasts using `pairwise` in the package `RRPP`.")
legPC1.by.species.plot

ggsave("plots/leg.PC1.by.species.png", legPC1.by.species.plot, width = 6, height = 4, scale = 1)


# Do species have different leg "shape" allometries?
model.unique.leg.allometries <- lm.rrpp(mm ~ log(its) * species, data = leg.list, iter = (1e4)-1)
anova(model.unique.leg.allometries)
#                   Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)           1 129.96 129.958 0.40554 189.6946 4.8818  1e-04 ***
# species            6  28.66   4.776 0.08943   6.9716 4.8744  1e-04 ***
# log(its):species   6  13.86   2.311 0.04327   3.3730 3.0893  6e-04 ***
# Residuals        216 147.98   0.685 0.46177                           
# Total            229 320.46    
# Yes, species have different leg "shape" allometries

# Post hoc examination of the interaction
pw.unique.leg.allometries <- pairwise(
  fit = model.unique.leg.allometries,
  fit.null = model.legs.species,
  groups = leg.list$species,
  covariate = leg.list$its
)
pw.unique.leg.allometries.output <- summary(pw.unique.leg.allometries)
(sig.letters <- pairwise.group.comparisons(pw.unique.leg.allometries.output))
# bimaculatus   impatiens    borealis   terricola   ternarius      vagans   perplexus 
#         "b"         "b"         "b"        "ab"         "a"         "b"        "ab" 
# B. ternarius has a decidedly different allometry from the other species        

####################
# PGLS
####################

# Phylogenetic regression using `nlme::gls`

legs.species.names <- unique(sort(as.character(leg.list$species)))

btree <- Bombus.tree
btree$tip.label <- str_split_fixed(btree$tip.label," ",2)[,2]

btree <- keep.tip(btree, btree$tip.label[which(btree$tip.label %in% legs.species.names)])
plot(btree)

# Create a dataframe with species means
legs.sp <- data.frame(
  PC1 = c(by(legs$legPC1, legs$species, mean, na.rm=TRUE)),
  PC2 = c(by(legs$legPC2, legs$species, mean, na.rm=TRUE)),
  its = c(by(legs$its, legs$species, mean, na.rm=TRUE)),
  body = c(by(legs$body, legs$species, mean, na.rm=TRUE))
)
# Reorder the dataframe to match the order of tree tips
legs.sp <- legs.sp[btree$tip.label,]

# Check that the names in the tree and dataset match
geiger::name.check(btree, legs.sp)

pgls.legPC1.by.ITS <- gls(
  PC1 ~ log(its), 
  correlation = corBrownian(phy = btree), 
  method = "REML", data = legs.sp
)
summary(pgls.legPC1.by.ITS)
with(legs.sp, plot(PC1 ~ log(its)))
abline(a = coef(pgls.PC1.by.ITS)[1], b = coef(pgls.PC1.by.ITS)[2])
# Generalized least squares fit by REML
#      AIC      BIC    logLik
# 24.73747 23.56578 -9.368735
# 
#                 Value Std.Error   t-value p-value
# (Intercept) -41.04802  9.669925 -4.244916  0.0081
# log(its)     26.99528  6.053105  4.459740  0.0066
# Degrees of freedom: 7 total; 5 residual

# So, size is a predictor of leg shape (or at least PC1 for
# these measurements), even after accounting for relatedness

# With foraging plant diversity

i <- match(row.names(legs.sp),names(SDI.by.species))
legs.sp$simpson <- SDI.by.species[i]

pgls.legPC1.by.ITS.SDI <- gls(
  PC1 ~ log(its) + simpson, 
  correlation = corBrownian(phy = btree), 
  method = "REML", data = legs.sp
)
summary(pgls.legPC1.by.ITS.SDI)
#      AIC      BIC    logLik
# 13.66461 11.20979 -2.832305
#                 Value Std.Error    t-value p-value
# (Intercept) -51.58580  71.26108 -0.7238986  0.5092
# log(its)     24.09581   8.74543  2.7552465  0.0511
# simpson      15.95103  64.53135  0.2471827  0.8169
# Residual standard error: 1.884945 
# Degrees of freedom: 7 total; 4 residual



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

# Add forage diversity metrics to the legs list object
{
  leg.list$richness <- NA
  leg.list$simpson <- NA
  leg.list$MDS1 <- NA
  leg.list$MDS2 <- NA
  }
for (i in 1:length(forage$species)) {
  j <- which(leg.list$species == forage$species[i])
  leg.list$richness[j] <- forage[forage$species[i],"richness"]
  leg.list$simpson[j] <- forage[forage$species[i],"simpson"]
  leg.list$MDS1[j] <- forage[forage$species[i],"MDS1"]
  leg.list$MDS2[j] <- forage[forage$species[i],"MDS2"]
}

# Explore potential correlations
leg.list %>%
  keep(names(.) %in% c("legPC1", "legPC2", "richness", "simpson", "MDS1", "MDS2")) %>%
  bind_rows() %>% 
  borealis::pairs()
# There may be some correlation between leg shape PC1 and forage plant richness

plot(leg.list$simpson, leg.list$legPC1 )
abline(lm(leg.list$legPC1 ~ leg.list$simpson))


####################
# Model the effect of foraging metrics on mouthpart shape
####################
# Comparisons among these models can be made based on Z values (effect size)
# Species should be a strong predictor, but the forage plant metric with the
# next best effect will be informative.

model.legs.richness <- lm.rrpp(
  mm ~ richness, 
  data = leg.list, 
  iter = 1e4-1, 
  print.progress = TRUE
)
anova(model.legs.richness)
#            Df     SS      MS     Rsq      F      Z Pr(>F)    
# richness    1  15.03 15.0260 0.04689 11.217 2.6328  5e-04 ***
# Residuals 228 305.43  1.3396 0.95311                         
# Total     229 320.46   

model.legs.simpson <- lm.rrpp(
  mm ~ simpson, 
  data = leg.list, 
  iter = 1e4-1, 
  print.progress = TRUE
)
anova(model.legs.simpson)
#            Df     SS     MS     Rsq      F      Z Pr(>F)
# simpson     1   5.44 5.4445 0.01699 3.9406 1.7353 0.0359 *
# Residuals 228 315.01 1.3816 0.98301                       
# Total     229 320.46   

model.legs.MDS1 <- lm.rrpp(
  mm ~ MDS1, 
  data = leg.list, 
  iter = 1e4-1, 
  print.progress = TRUE
)
anova(model.legs.MDS1)
#            Df     SS     MS     Rsq      F     Z Pr(>F)  
# MDS1        1   3.77 3.7731 0.01177 2.7165 1.434 0.0826 .
# Residuals 228 316.69 1.3890 0.98823                      
# Total     229 320.46 

model.legs.MDS2 <- lm.rrpp(
  mm ~ MDS2, 
  data = leg.list, 
  iter = 1e4-1, 
  print.progress = TRUE
)
anova(model.legs.MDS2)
#            Df     SS      MS     Rsq      F        Z Pr(>F)
# MDS2        1   0.40 0.39761 0.00124 0.2832 -0.62402 0.7045
# Residuals 228 320.06 1.40377 0.99876                       
# Total     229 320.46                               

# How do these factors do in models with ITS?
model.legs.size.simpson <- lm.rrpp(
  mm ~ log(its) + simpson, 
  data = leg.list, 
  iter = 1e4-1, 
  print.progress = FALSE
)

model.legs.size.simpson.results <- 
  anova(model.legs.size.simpson)$table %T>%
  print() %>% 
  mutate(Rsq = signif(Rsq,3))
#            Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)    1 129.96 129.958 0.40554 155.6237 4.7330 0.0001 ***
# simpson     1   0.94   0.939 0.00293   1.1241 0.6014 0.2844    
# Residuals 227 189.56   0.835 0.59153                           
# Total     229 320.46 

# Actually, this results suggests that the forage plant richness does not 
# captures much of the shape variance, after accounting for the effect of size.


# Save everything
save.image("bombus.scaling.rda")