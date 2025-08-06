# Relative scaling of Bombus mouthpart, wing and legs 
# sensu [O’Brien et al. 2018. *Animal Behavior*](https://www.sciencedirect.com/science/article/pii/S0003347218302483)

####################
# Before running this script, run:
# - forage.data.R
# - mouthpart.analysis.R
# - leg.analysis.R
# - wing.analysis.R
# - scaling.R (this script)

# Obtain the following objects:
# - `gnathos.list`: a list containing worker mouthpart measurements
# - `gpa.fw`: geomorph dataframe with GPA-aligned forewing shape coords and metadata
# - `gpa.hw`: geomorph dataframe with GPA-aligned hindwing shape coords and metadata
# - `leg.list`: a list containing worker leg measurements
# These objects should also contain forage diversity metrics

####################
# Notes
# 
# Devin's analysis of scaling relationships describes a method to examine allometry of
# a structure by comparing the slope of a focal trait to that of a "reference trait" 
# (rather than compared to a slope of 1). The reference trait should be one that's 
# similar in many ways, but unlikely to be subject to the selective pressures relevant 
# to the focal trait (such as e.g. natural selection as a sexual signal). His analyses 
# also examined differences in variation (the CV) as it related to "hypervariability" 
# predicted for sexually selected signal traits.
# 
# Our plan here is to use bumblebee mouthparts as a test of the counter-example. Bees
# vary widely in body size, but pollinate a similar pool of available flowers. 
# 
# ## Predictions:
# - Body size (or ITS) should vary less for specialists, due to
#   the  stabilizing selection pressures to match pollinators to flowers.
# - Trait dispersion should be lower for mouthpart traits than for reference traits.
#   e.g. lower coefficient of variation (CV) for tongue length vs. legs or wings
#   e.g. lower disparity in morphospace for mouthparts vs. legs or wings
# - Scaling for mouthparts should be less than for a reference traits (wings, legs)
#   e.g. tongue length vs. length of wings or legs
#   e.g. mouthpart shape allometry vs. allometry or forewing shape, hindwing shape
#        or leg multivariate measurements
# - These effects should be more extreme for species that are specialist foragers,
#   than for more generalist pollinators.
# - These trends should be true of workers, but reproductive castes may be freed
#   from these constraints, because they do not forage to the same degree, and colony
#   success may not depend as much as on their efficiency during pollination.

####################
# Data curation
####################

load("bombus.scaling.rda", verbose = TRUE)

# Duplicate IDs? Should all be zeros.
sum(duplicated(gnathos.list$specimen.id))
sum(duplicated(leg.list$specimen.id))
sum(duplicated(gpa.fw$gdf$specimen.id))
sum(duplicated(gpa.hw$gdf$specimen.id))

# Combine the datasets
metadata.df <- read.csv("metadata.csv") %>% 
  mutate(
    species = if_else(species.det=="",provisional.species,species.det),
    its = its.mm,
    body = body.length.mm
  ) %>% 
  filter(
    !grepl("^cf", species),
    caste == "W"
  ) %>% 
  group_by(species) %>% filter(n() > 10) %>% 
  select(specimen.id, species, its, body) 

forage <- read.csv("forage.plant.diversity.csv") %>% 
  mutate(
    species = sub("Bombus ","",species),
    simpson = Simpsons
  ) %>% 
  mutate(to.rownames = species) %>% 
  column_to_rownames(var = "to.rownames") 

# Add forage diversity metrics
metadata.df$richness <- NA
metadata.df$simpson <- NA

for (i in 1:length(forage$species)) {
  j <- which(metadata.df$species == forage$species[i])
  metadata.df$richness[j] <- forage[forage$species[i],"richness"]
  metadata.df$simpson[j] <- forage[forage$species[i],"simpson"]
}

gnathos.df <- data.frame(
  specimen.id = gnathos.list$specimen.id,
  tongueL = gnathos.list$tongue,
  mpPC1 = gnathos.list$PC1,
  mpPC2 = gnathos.list$PC2
)

leg.df <- data.frame(
  specimen.id = leg.list$specimen.id,
  legPC1 = leg.list$legPC1,
  legPC2 = leg.list$legPC2,
  leg1L = leg.list$l1_legL,
  leg2L = leg.list$l2_legL,
  leg3L = leg.list$l3_legL
)

fw.df <- data.frame(
  specimen.id = gpa.fw$gdf$specimen.id,
  forewingL = gpa.fw$gdf$forewing.length.mm,
  forewingPC1 = gpa.fw$gdf$PC1,
  forewingPC2 = gpa.fw$gdf$PC2
)

hw.df <- data.frame(
  specimen.id = gpa.hw$gdf$specimen.id,
  hindwingL = gpa.hw$gdf$hindwing.length.mm,
  hindwingPC1 = gpa.hw$gdf$PC1,
  hindwingPC2 = gpa.hw$gdf$PC2
)

# Combined the datasets into one data frame
combined <- left_join(
  metadata.df, gnathos.df,  
  by = "specimen.id"
) %>% 
  left_join(., leg.df, by = "specimen.id") %>% 
  left_join(., fw.df, by = "specimen.id") %>% 
  left_join(., hw.df, by = "specimen.id") 

dim(combined)
head(combined)

####################
# Table S1. Specimen numbers by species in each dataset
####################

table(combined$species); length(combined$species)
c(with(combined, by(tongueL, species, function(x) {sum(!is.na(x))})))
c(with(combined, by(mpPC1, species, function(x) {sum(!is.na(x))})))
c(with(combined, by(leg1L, species, function(x) {sum(!is.na(x))})))
c(with(combined, by(leg2L, species, function(x) {sum(!is.na(x))})))
c(with(combined, by(leg3L, species, function(x) {sum(!is.na(x))})))
c(with(combined, by(legPC1, species, function(x) {sum(!is.na(x))})))
c(with(combined, by(forewingPC1, species, function(x) {sum(!is.na(x))})))
c(with(combined, by(hindwingPC1, species, function(x) {sum(!is.na(x))})))

combined %>% 
  filter(
    !is.na(tongueL), !is.na(mpPC1), 
    !is.na(leg1L), !is.na(leg2L), !is.na(leg3L), !is.na(legPC1), 
    !is.na(forewingPC1), !is.na(hindwingPC1)
  ) %>% 
  group_by(species) %>% 
  summarise(n=n())

####################
# Size distributions by species
####################

its.distributions.by.species.plot <- combined %>% 
  filter(!is.na(its)) %>% 
  ggplot(aes((its))) +
  theme_minimal() +
  theme(
    panel.border = element_rect(linewidth = 0.5, fill = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(face="italic", size = 12)
  ) +
  facet_grid(species~.) +
  geom_histogram() +
  labs(x = "intertegular span (mm)")
its.distributions.by.species.plot

ggsave("manuscript/FigS2.its.distributions.png", its.distributions.by.species.plot, width = 6.5, height = 6.5, scale = 1)

# Test for differences in ITS means
kruskal.test(combined$its,  combined$species)
# Kruskal-Wallis rank sum test χ2 = 145.22, df = 7, p < 2.2e-16

# Test for differences in variance
car::leveneTest(combined$its,  combined$species)
# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value    Pr(>F)    
# group   7  3.8534 0.0004231 ***
#       520               

# Post-hoc pairwise F tests
{
  species.contrasts <- combn(sort(unique(combined$species)), 2, simplify = FALSE)
  pw.its <- numeric(length(species.contrasts))
  names(pw.its) <- sapply(species.contrasts, paste, collapse = " vs ")
  x <- pw.its
  for (i in seq_along(species.contrasts)) {
    g1 <- species.contrasts[[i]][1]
    g2 <- species.contrasts[[i]][2]
    data_g1 <- combined$its[combined$species == g1]
    data_g2 <- combined$its[combined$species == g2]
    test.i <- var.test(data_g1, data_g2)
    pw.its[i] <- test.i$p.value
    x[i] <- test.i$estimate
  }
  pw.its <- data.frame(
    ratio.of.var = x,
    p = pw.its,
    # p.fdr = p.adjust(pw.its, method = "fdr"),
    p.bonferroni = p.adjust(pw.its, method = "bonferroni")
  )
  pw.its$signif <- ifelse(pw.its$p.bonferroni < 0.001,"***",
                          ifelse(pw.its$p.bonferroni < 0.01,"**",
                                 ifelse(pw.its$p.bonferroni < 0.05,"*","")))
}

pw.its
# These results become Table S4


####################
# Calculate CV for linear measurements
####################
cv <- function(x) { sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) }

trait.cv <- combined %>%
  group_by(species) %>% 
  reframe(
    cv.its = cv(its),
    cv.body = cv(body),
    cv.tongue = cv(tongueL),
    cv.leg1 = cv(leg1L),
    cv.leg2 = cv(leg2L),
    cv.leg3 = cv(leg3L),
    cv.forewing = cv(forewingL),
    cv.hindwing = cv(hindwingL)
  ) 
trait.cv

SDI.by.species <- combined %>% 
  group_by(species) %>% 
  summarise(simpson = median(simpson, na.rm=TRUE)) %>%
  pull(simpson) %>% 
  set_names(sort(unique(combined$species))) 

trait.cv.plot <- trait.cv %>% 
  mutate(simpson = SDI.by.species) %>% 
  select(-cv.its, -cv.body) %>% 
  pivot_longer(cols = -c(1,last_col()), names_to = "trait") %>% 
  mutate(type = if_else(grepl("tongue",trait),"focal","reference")) %>% 
  mutate(trait = sub("^cv\\.","",trait)) %>% 
  mutate(trait = sub("leg1"," foreleg",trait)) %>% 
  mutate(trait = sub("leg2","  midleg",trait)) %>% 
  mutate(trait = sub("leg3","  hindleg",trait)) %>% 
  mutate(trait = sub("hindwing","    hindwing",trait)) %>% 
  mutate(trait = sub("forewing","   forewing",trait)) %>% 
  mutate(species = paste0(species,"\n(SDI = ",round(simpson,3),")")) %>%
  mutate(species = fct_reorder(species, simpson)) %>%
  ggplot(aes(x=trait, y=abs(value), fill=type)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(linewidth = 0.5, fill = NA),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(face="italic", size = 10)
  ) +
  facet_wrap(.~species, nrow = 1) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("darkred","grey50")) +
  scale_y_continuous(
    name = "coefficient of variation", 
    breaks = c(0,0.1,0.2), labels = c("0","0.1","0.2")
  ) +
  coord_flip() 

trait.cv.plot


####################
# Disparity in morphospace
####################

morphospace.hypervolumes <- left_join(
  gnathos.hypervolumes, leg.hypervolumes,
  by = "species"
) %>% 
  arrange(species)

morphospace.hypervolumes

z.scores <- function(x) { (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) }


# Disparity in Procrustes variance
procrustes.variance <- data.frame(
  species = names(fw.disparity$Procrustes.var),
  forewings = fw.disparity$Procrustes.var,
  hindwings = hw.disparity$Procrustes.var
) 

procrustes.variance


multidimensional.variance.plot <- left_join(
  morphospace.hypervolumes,
  procrustes.variance,
  by = "species"
) %>% 
  mutate(simpson = SDI.by.species) %>% 
  mutate(
    gnathos.hv = z.scores(gnathos.hv),
    leg.hv = z.scores(leg.hv),
    forewings = z.scores(forewings),
    hindwings = z.scores(hindwings)
  ) %>% 
  pivot_longer(cols = 2:5, names_to = "trait") %>% 
  mutate(trait = sub("gnathos.hv","mouthparts",trait)) %>% 
  mutate(trait = sub("leg.hv","legs",trait)) %>% 
  mutate(trait = sub("hindwings"," hindwings",trait)) %>% 
  mutate(type = if_else(grepl("mouthparts",trait),"focal","reference")) %>% 
  mutate(species = paste0(species,"\n(SDI = ",round(simpson,3),")")) %>%
  mutate(species = fct_reorder(species, simpson)) %>%
  ggplot(aes(x=trait, y=abs(value), fill=type)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(linewidth = 0.5, fill = NA),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(face="italic", size = 10)
  ) +
  facet_wrap(.~species, nrow = 1) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("darkred","grey50")) +
  scale_y_continuous(
    name = "standardized variance", breaks = 0:2) +
  # ylab("standardized variance") +
  coord_flip() 

multidimensional.variance.plot

disparity.plots <- ggarrange(
  trait.cv.plot, multidimensional.variance.plot, 
  ncol = 1, labels = "AUTO", heights = c(6,4)+3
)

disparity.plots

ggsave("manuscript/Fig1.disparity.plots.png", disparity.plots, width = 6.5, height = 4, scale = 1.5)


# Scaling plots
# Figure S3

scaling.plots <- ggarrange(
  tongue.scaling$plot,
  hindleg.scaling$plot,
  forewing.scaling$plot,
  hindwing.scaling$plot,
  nrow = 2, ncol = 2, labels = "AUTO", 
  common.legend = TRUE, legend = "none"
)
scaling.plots

ggsave("manuscript/FigS3.univariate.scaling.plots.png", scaling.plots, width = 6.5, height = 6.5, scale = 1.2)


# Univariate trait scaling coefficients
# Figure S4

univariate.slope.plots <- ggarrange(
  tongue.scaling.slopes.plot,
  hindleg.scaling.slopes.plot,
  forewing.scaling.slopes.plot,
  hindwing.scaling.slopes.plot,
  nrow = 1, labels = "AUTO"
)
univariate.slope.plots

ggsave("manuscript/FigS4.univariate.slopes.plot.png", univariate.slope.plots, width = 6.5, height = 2.5, scale = 1.5)


# Multivariate trait scaling plots
# Figure S5

multivariate.scaling.plots <- ggarrange(
  PC1.scaling$plot,
  leg.PC1.scaling$plot,
  fw.shape.scaling$plot,
  hw.shape.scaling$plot,
  nrow = 2, ncol = 2, labels = "AUTO", 
  common.legend = TRUE, legend = "none"
)
multivariate.scaling.plots
ggsave("manuscript/FigS5.multivariate.scaling.plots.png", multivariate.scaling.plots, width = 6.5, height = 6.5, scale = 1.2)


# Multivariate trait scaling coefficients
# Figure S6

multivariate.slope.plots <- ggarrange(
  scaling.PC1.slopes.plot,
  leg.scaling.PC1.slopes.plot,
  scaling.forewing.PC1.slopes.plot,
  scaling.hindwing.PC1.slopes.plot,
  nrow = 1, labels = "AUTO"
)
multivariate.slope.plots
ggsave("manuscript/FigS6.multivariate.slope.plots.png", multivariate.slope.plots, width = 6.5, height = 2.5, scale = 1.5)


# 
specialization.plots <- ggarrange(
  specialization.plot,
  
  nrow = 2, ncol = 2, labels = "AUTO"
)

specialization.plots


# Save everything
save.image("bombus.scaling.rda")
