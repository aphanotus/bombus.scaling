# Bombus leg morphometry analysis

# Clear memory
# rm(list=ls())
# gc() 

# devtools::install_github("aphanotus/borealis")
x <- c("tidyverse","ggpubr","ggrepel","factoextra","viridis",
       "RRPP","multcompView","ggbeeswarm","phytools","nlme")
invisible(lapply(x, require, character.only = TRUE))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

####################
# Import data
####################

legs <- read.csv("leg.measurements.csv")
# dim(legs); names(legs)
# Create a factor combining species and caste
legs$sp.caste <- paste0(legs$species,".",legs$caste)
# Add species and caste abbreviations to the specimen IDs
legs$specimen_id <- paste0(legs$specimen_id,'.',legs$sp.caste) 
head(legs,3)

# Match leg measurements to body size measures from the mouthpart or wing datasets
legs$body <- NA
legs$its <- NA
# head(legs$specimen_id); head(m$specimen_id)
i <- match(legs$specimen_id, m$specimen_id)
legs$body <- m$body[i]
legs$its <- m$its[i]
i <- match(legs$specimen_id, sub("bimac","bim",dimnames(fw.gpa$gdf$coords)[[3]]))
x <- which(is.na(legs$body))
legs$body[x] <- fw.gpa$gdf$body[i[x]]
x <- which(is.na(legs$its))
legs$its[x] <- fw.gpa$gdf$its[i[x]]

# Manual matching
legs$body[grep("DRA180814-006.vag.W",legs$specimen_id)] <- fw.gpa$gdf$body[grep("DRA180814-006a.vag.W",dimnames(fw.gpa$gdf$coords)[[3]])]
legs$its[grep("DRA180814-006.vag.W",legs$specimen_id)] <- fw.gpa$gdf$its[grep("DRA180814-006a.vag.W",dimnames(fw.gpa$gdf$coords)[[3]])]

# There really aren't enough reproductive castes in this dataset
# So filter to just workers 
legs <- filter(legs, caste=="W")

# Filter out low sample size (<3)
barplot(sort(c(with(legs, by(species, species, length)))), cex.names = 0.5)
abline(h=c(3,10), col = "darkred")
legs <- legs %>% group_by(species) %>% filter(n() > 3)

# Average left and right sides
long.legs <- legs %>% 
  pivot_longer(cols = 6:41, names_to = "trait") %>% 
  mutate(trait = sub("_[12]$","",trait)) %>% 
  group_by(specimen_id, species, body, its, trait) %>%
  summarise(trait.value = mean(value, na.rm=TRUE)) 
# head(legs)

# Find leg lengths (femur + tibia + basitarsus)
legs <- long.legs %>% 
  pivot_wider(names_from = 5, values_from = 6) %>% 
  group_by(specimen_id) %>% 
  mutate(
    l1_legL = l1_femurL + l1_tibiaL + l1_btarsL,
    l2_legL = l2_femurL + l2_tibiaL + l2_btarsL,
    l3_legL = l3_femurL + l3_tibiaL + l3_btarsL,
    ) 

long.legs <- legs %>% 
  pivot_longer(cols = 5:25, names_to = "trait", values_to = "trait.value")
tail(long.legs)

####################
# Scaling plots
####################

# Preliminary look at scaling of legs
long.legs %>% 
  filter(grepl("legL",trait)) %>% 
  ggplot(aes(x = log10(body), y = log10(trait.value))) +
  theme_bw() +
  facet_grid(trait ~ species) +
  geom_point(alpha = 0.35) +
  geom_smooth(method=lm, se=FALSE)
# Slopes appear to differ by species, but seem consistent among legs.

scaling.plot(
  x = log10(legs$body), y = log10(legs$l3_legL),
  group = legs$species,
  xlab = "body length (log10 mm)", ylab = "hind leg (log10 mm)",
  include.legend = TRUE,
  isometry.line = TRUE,
  convex.hulls = FALSE,
  groups.trendlines = TRUE,
  save.as = "plots/scaling.L3leg.v.body.pdf"
)
#   group  n slope        p sig   ci.lo ci.hi spans.zero
# 1   bim 42 0.886 8.95e-10 ***  0.6610 1.110      FALSE
# 2   imp 52 0.874 4.34e-10 ***  0.6470 1.100      FALSE
# 3  tern 50 0.437 1.18e-02 *    0.1020 0.773      FALSE
# 4   vag 47 0.436 2.32e-02 *    0.0625 0.810      FALSE
# 5   bor 13 0.267 4.40e-01     -0.4680 1.000       TRUE
# 6 terri 14 0.220 4.28e-01     -0.3640 0.803       TRUE
# No differences in slope based on CI

k.legs <- scaling.plot(
  x = log10(legs$its), y = log10(legs$l3_legL),
  group = legs$species,
  xlab = "intertegular span (log10 mm)", ylab = "hind leg (log10 mm)",
  include.legend = TRUE,
  isometry.line = TRUE,
  convex.hulls = FALSE,
  groups.trendlines = TRUE,
  save.as = "plots/scaling.L3leg.v.its.pdf"
)
k.legs$plot
k.legs$slopes
#   group  n  slope        p sig  ci.lo ci.hi spans.zero
# 1   imp 54 0.9720 7.18e-12 ***  0.750 1.190      FALSE
# 2   bim 42 0.8590 3.79e-09 ***  0.628 1.090      FALSE
# 3   vag 49 0.6910 5.22e-05 ***  0.379 1.000      FALSE
# 4 terri 14 0.6010 2.24e-02 *    0.101 1.100      FALSE
# 5   bor 13 0.4970 3.15e-01     -0.543 1.540       TRUE
# 6  tern 50 0.0515 7.70e-01     -0.300 0.403       TRUE
# Differences in slope based on CI
# imp bim vag terri bor tern
#   a   a  ab   ab   ab    b 

k.legs$slopes$group <- factor(k.legs$slopes$group, levels = k.legs$slopes$group)
ylimits <- sort(k.legs$slopes$ci.lo)[1]
ylimits <- c(ylimits, sort(k.legs$slopes$ci.hi, decreasing = TRUE)[1])
ylimits[1] <- ylimits[1]*1.1
k.legs$slopes %>% ggplot(aes(x=group, y=slope, fill=group, label=paste0("n=",n))) +
  theme_bw() +
  theme(legend.position="none") +
  geom_bar(stat="identity", color="grey15") +
  geom_errorbar(aes(ymin=ci.lo, ymax=ci.hi), width=.2) +
  scale_fill_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(y=ylimits[1], size=3, color="grey50", vjust=1) +
  scale_y_continuous(limits=ylimits) +
  xlab("species") +
  ylab("scaling coefficent (k) for hindleg length vs. intertegular span")
ggsave("plots/scaling.hindleg.by.species.pdf", width = 6, height = 4, scale = 1)

# PCA

# Many specimens are missing the first leg, so it will be omitted
data.cols <- grep("^l[23]_[bft]",names(legs))
# Index specimens without missing data 
i <- which(apply(legs[,data.cols], 1, function(x) {!any(is.na(x))}))
clean.legs <- as.data.frame(legs[i,data.cols])
rownames(clean.legs) <- legs$specimen_id[i]

leg.pca <- prcomp(clean.legs, center = TRUE, scale = TRUE)
shape.space(leg.pca, group = legs$species[i], 
            include.legend = TRUE,
            convex.hulls = TRUE,
            save.as = "plots/leg.morphospace.pdf")
legs$legPC1 <- NA; legs$legPC2 <- NA
legs$legPC1[i] <- get_pca_ind(leg.pca)$coord[,1]
legs$legPC2[i] <- get_pca_ind(leg.pca)$coord[,2]

# A quick look at how shape (as PC1) varies by body size (ITS) and species
# Remember, this is just a univariate representation of leg "shape". 
# This issue is revisited below with multivariate stats.
scaling.plot(
  x = log10(legs$its), y = legs$legPC1,
  group = legs$species,
  xlab = "intertegular span (log10 mm)", 
  ylab = "PC1 for worker leg measures",
  include.legend = TRUE,
  fixed.aspect = FALSE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  save.as = "plots/worker.leg.PC1.v.its.pdf"
)

# A quick plot of mouthpart "shape" (as PC1) by species
leg.sp.order <- legs %>%
  group_by(species) %>%
  summarise(medianPC1 = median(legPC1, na.rm=TRUE)) %>%
  arrange(medianPC1) %>% pull(species)
sp <- factor(legs$species, levels = leg.sp.order)
boxplot(legs$legPC1 ~ sp)

####################
# Stats using RRPP
####################

# Exclude specimens with missing ITS measurements
legs <- legs[which(!is.na(legs$its)),]

clean.legs$body <- legs$body[i]
clean.legs$its <- legs$its[i]
clean.legs$species <- legs$species[i]
clean.legs$legPC1 <- legs$legPC1[i]
clean.legs$legPC2 <- legs$legPC2[i]
names(clean.legs)

leg.list <- list(
  mm = as.matrix(clean.legs[,1:12]),
  body = clean.legs$body,
  its = clean.legs$its,
  legPC1 = clean.legs$legPC1,
  legPC2 = clean.legs$legPC2,
  species = factor(clean.legs$species, levels = leg.sp.order)
)

legs.rrpp.species <- lm.rrpp(mm ~ log(its) + species, data = leg.list, iter = 1e4-1)
anova(legs.rrpp.species)
#            Df      SS      MS     Rsq       F      Z    Pr(>F)
# log(its)    1 188.27 188.270 0.43698 183.6058 4.9233  1e-04 ***
# species     5  35.44   7.087 0.08225   6.9119 4.5292  1e-04 ***
# Residuals 202 207.13   1.025 0.48076                           
# Total     208 430.84   

# Pairwise species differences?
legs.rrpp.size <- lm.rrpp(mm ~ log(its), data = leg.list, iter = 1e4-1)
legs.species.pw <- pairwise(fit = legs.rrpp.species,
                            fit.null = legs.rrpp.size,
                            groups = leg.list$species)
(legs.species.pw.output <- summary(legs.species.pw))
(sig.letters <- pairwise.group.comparisons(legs.species.pw.output))
#  bor   bim terri   imp   vag  tern 
# "a"   "a"   "a"   "a"   "b"   "b" 

df <- data.frame(
  x = 1:6,
  y = rep(11.5,6),
  label = sig.letters$Letters
)

species.order <- clean.legs %>%
  group_by(species) %>%
  summarise(medianPC1 = median(legPC1)) %>%
  arrange(medianPC1) %>% pull(species)
sp <- factor(clean.legs$species, levels = species.order)

ggplot(clean.legs, aes(sp, legPC1, color = sp)) +
  theme_minimal() +
  theme(legend.position="none") +
  geom_violin() +
  geom_quasirandom(varwidth = TRUE) +
  scale_color_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(data = df, aes(x = x, y = y, label = label), color = "black", size = 5 ) +
  xlab("species") +
  ylab("PC1 for worker leg measurements") +
  labs(caption = "Significant differences in the influence of species in the model  Y ~ intertegular span + species.\nPERMANOVA using `lm.rrpp`; contrasts using `pairwise` in the package `RRPP`.")
ggsave("plots/worker.leg.PC1.by.species.pdf", width = 6, height = 4, scale = 1)

# Allometry-corrected PCA
legs.allo.pca <- prcomp(resid(legs.rrpp.size), center = TRUE, scale = TRUE)
shape.space(legs.allo.pca, group = leg.list$species, convex.hulls = TRUE,
            main.title = "Size-corrected morphospace for worker leg measurements",
            save.as = "plots/allo-corrected.pca.worker.legs.pdf")
x <- get_pca_ind(legs.allo.pca)$coord[,1]
clean.legs$alloPC1 <- get_pca_ind(legs.allo.pca)$coord[,1]
clean.legs$alloPC2 <- get_pca_ind(legs.allo.pca)$coord[,2]

species.order <- clean.legs %>%
  group_by(species) %>%
  summarise(medianAlloPC1 = median(alloPC1)) %>%
  arrange(medianAlloPC1) %>% pull(species)
clean.legs$species <- factor(clean.legs$species, levels = species.order)

# Repeat this just to get the species in the right order!
leg.list <- list(
  specimen_id = clean.legs$specimen_id,
  species = factor(clean.legs$species, levels = species.order),
  mm = as.matrix(clean.legs[,1:12]),
  body = clean.legs$body,
  its = clean.legs$its,
  legPC1 = clean.legs$legPC1,
  legPC2 = clean.legs$legPC2,
  alloPC1 = clean.legs$alloPC1
)
legs.rrpp.size <- lm.rrpp(mm ~ log(its), data = leg.list, iter = 1e4-1)
legs.rrpp.species <- lm.rrpp(mm ~ log(its) + species, data = leg.list, iter = 1e4-1)
legs.species.pw <- pairwise(fit = legs.rrpp.species,
                            fit.null = legs.rrpp.size,
                            groups = leg.list$species)
(legs.species.pw.output <- summary(legs.species.pw))
(sig.letters <- pairwise.group.comparisons(legs.species.pw.output))
# bim   imp   bor terri  tern   vag 
# "a"   "a"   "a"   "a"   "b"   "b" 

df <- data.frame(
  x = 1:6,
  y = rep(10, 6),
  label = sig.letters$Letters
)

ggplot(clean.legs, aes(species, alloPC1, color = species)) +
  theme_minimal() +
  theme(legend.position="none") +
  geom_violin() +
  geom_quasirandom(varwidth = TRUE, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(data = df, aes(x = x, y = y, label = label), color = "black", size = 5 ) +
  xlab("species") +
  ylab("size-corrected PC1 for worker leg measures") +
  labs(caption = "Significant differences in the influence of species in the model  Y ~ intertegular span + species.\nPERMANOVA using `lm.rrpp`; contrasts using `pairwise` in the package `RRPP`.")
ggsave("plots/size-corrected.worker.leg.PC1.by.species.pdf", width = 6, height = 4, scale = 1.25)

# Do species have different worker mouthpart "shape" allometries?
legs.rrpp.unique.allometries <- lm.rrpp(mm ~ log(its) * species, data = leg.list, iter = 1e4-1)
anova(legs.rrpp.unique.allometries)
#                   Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)           1 188.27 188.270 0.43698 202.2972 4.9945  1e-04 ***
# species            5  35.44   7.087 0.08225   7.6155 4.7617  1e-04 ***
# log(its):species   5  23.79   4.758 0.05522   5.1128 3.8901  1e-04 ***
# Residuals        197 183.34   0.931 0.42554                           
# Total            208 430.84     
# Yes, species have different worker mouthpart "shape" allometries

# Post hoc examination of the interaction
legs.rrpp.unique.allometries.pw <- 
  pairwise(fit = legs.rrpp.unique.allometries,
           fit.null = legs.rrpp.species,
           groups = leg.list$species,
           covariate = leg.list$its)
(allometry.pw.legs.output <- summary(legs.rrpp.unique.allometries.pw))
#                    d UCL (95%)          Z Pr > d
# bim:imp    0.2244365 0.7147253 -0.7338118 0.7479
# bim:bor    1.0208504 1.4802817  0.9292717 0.1918
# bim:terri  0.5214521 1.3037160 -0.2984857 0.5987
# bim:tern   1.3538608 0.8125567  2.6340868 0.0009
# bim:vag    0.3952870 0.6667653  0.5589130 0.3001
# imp:bor    1.0783300 1.3350182  1.2143267 0.1258
# imp:terri  0.6485176 1.2239050  0.2894764 0.3856
# imp:tern   1.4710678 0.8786336  2.6893133 0.0002
# imp:vag    0.5127552 0.5511990  1.4978952 0.0717
# bor:terri  0.7877389 1.7603806 -0.0243850 0.5015
# bor:tern   0.8700484 1.5992171  0.4566145 0.3403
# bor:vag    0.7409478 1.3801257  0.3947277 0.3532
# terri:tern 0.9501487 1.4036385  0.8756157 0.2090
# terri:vag  0.3181454 1.2346152 -1.2974181 0.9016
# tern:vag   0.9712765 0.8075596  1.9923694 0.0160
(sig.letters <- pairwise.group.comparisons(allometry.pw.legs.output))
# bim   imp   bor terri  tern   vag 
# "a"   "a"  "ab"  "ab"   "b"   "a" 

####################
# PGLS
####################

leg.species.names <- unique(sort(as.character(leg.list$species)))
leg.species.names <- sub("bim","bimac",leg.species.names)
leg.species.names <- sub("san","sande",leg.species.names)

btree <- Bombus.tree
btree$tip.label <- btree$code.name
btree <- keep.tip(btree, btree$tip.label[which(btree$tip.label %in% leg.species.names)])
btree$tip.label <- sub("bimac","bim",btree$tip.label)
btree$tip.label <- sub("sande","san",btree$tip.label)
plot(btree)

# Create a dataframe with species means
leg.sp <- data.frame(
  PC1 = c(by(legs$legPC1, legs$species, mean, na.rm=TRUE)),
  PC2 = c(by(legs$legPC2, legs$species, mean, na.rm=TRUE)),
  alloPC1 = c(by(clean.legs$alloPC1, clean.legs$species, mean, na.rm=TRUE)),
  body = c(by(legs$body, legs$species, mean, na.rm=TRUE)),
  its = c(by(legs$its, legs$species, mean, na.rm=TRUE))
)
# Reorder the dataframe to match the order of tree tips
leg.sp <- leg.sp[match(btree$tip.label, rownames(leg.sp)),]

# Check that the names in the tree and dataset match
geiger::name.check(btree, leg.sp)

pgls.legPC1.by.its <- gls(
  PC1 ~ log(its),
  correlation = corBrownian(phy = btree),
  method = "REML", data = leg.sp)
summary(pgls.legPC1.by.its)
with(leg.sp, plot(PC1 ~ log(its)))
abline(a = coef(pgls.legPC1.by.its)[1], b = coef(pgls.legPC1.by.its)[2])
# Generalized least squares fit by REML
#      AIC      BIC    logLik
# 15.88783 14.04671 -4.943913
# 
#                 Value Std.Error   t-value p-value
# (Intercept)  31.46868  5.028289  6.258329  0.0033
# log(its)    -20.88437  3.180485 -6.566411  0.0028
# Degrees of freedom: 6 total; 4 residual

pgls.legPC2.by.its <- gls(
  PC2 ~ log(its),
  correlation = corBrownian(phy = btree),
  method = "REML", data = leg.sp)
summary(pgls.legPC2.by.its)
with(leg.sp, plot(PC2 ~ log(its)))
abline(a = coef(pgls.legPC2.by.its)[1], b = coef(pgls.legPC2.by.its)[2])
# Generalized least squares fit by REML
#      AIC      BIC     logLik
# 6.927721 5.086604 -0.4638603
# 
#                 Value Std.Error   t-value p-value
# (Intercept)  2.015683  1.640607  1.228620  0.2866
# log(its)    -1.396775  1.037714 -1.346011  0.2495
# Degrees of freedom: 6 total; 4 residual

####################
# Correlations to forage diversity
####################
( forage <- read.csv("Wood.et.al.2019.forage.diversity.csv") )

# Make the species abbreviations (code names) match
leg.list$code.name <- as.character(leg.list$species)
forage$code.name <- sub("bimac","bim",forage$code.name)

# Filter out the foraging information that covers species not in out dataset
x <- which(forage$code.name %in% unique(leg.list$code.name))
forage <- forage[x,]
# All species in our dataset are covered by the foraging information

# Add the foraging information to the leg list object
leg.list$class <- leg.list$code.name
leg.list$wood.dbs <- leg.list$code.name
leg.list$richness <- leg.list$code.name
leg.list$shannon <- leg.list$code.name
leg.list$simpson <- leg.list$code.name
leg.list$faith.pd <- leg.list$code.name
leg.list$wood.PC1 <- leg.list$code.name
leg.list$wood.PC2 <- leg.list$code.name
leg.list$wood.PC3 <- leg.list$code.name
df <- forage[,-c(1:2)]
colnames(df)[1] <- "class"
row.names(df) <- forage$code.name

for (i in 1:length(leg.list$code.name)) {
  leg.list$class[i] <- df[leg.list$code.name[i],"class"]
  leg.list$wood.dbs[i] <- df[leg.list$code.name[i],"wood.dbs"]
  leg.list$richness[i] <- df[leg.list$code.name[i],"richness"]
  leg.list$shannon[i] <- df[leg.list$code.name[i],"shannon"]
  leg.list$simpson[i] <- df[leg.list$code.name[i],"simpson"]
  leg.list$faith.pd[i] <- df[leg.list$code.name[i],"faith.pd"]
  leg.list$wood.PC1[i] <- df[leg.list$code.name[i],"wood.PC1"]
  leg.list$wood.PC2[i] <- df[leg.list$code.name[i],"wood.PC2"]
  leg.list$wood.PC3[i] <- df[leg.list$code.name[i],"wood.PC3"]
}
leg.list$class <-  as.factor(leg.list$class)
leg.list$wood.dbs <-  as.numeric(leg.list$wood.dbs)
leg.list$richness <-  as.numeric(leg.list$richness)
leg.list$shannon <-  as.numeric(leg.list$shannon)
leg.list$simpson <-  as.numeric(leg.list$simpson)
leg.list$faith.pd <-  as.numeric(leg.list$faith.pd)
leg.list$wood.PC1 <-  as.numeric(leg.list$wood.PC1)
leg.list$wood.PC2 <-  as.numeric(leg.list$wood.PC2)
leg.list$wood.PC3 <-  as.numeric(leg.list$wood.PC3)

# Explore potential correlations
x <- as.data.frame(leg.list$mm)
x$body <- unlist(leg.list$body)
x$its <- unlist(leg.list$its)
x$PC1 <- unlist(leg.list$legPC1)
x$PC2 <- unlist(leg.list$legPC2)
x$alloPC1 <- unlist(leg.list$alloPC1)
x$code.name <- unlist(leg.list$code.name)
x$class <- unlist(leg.list$class)
x$wood.dbs <- unlist(leg.list$wood.dbs)
x$richness <- unlist(leg.list$richness)
x$shannon <- unlist(leg.list$shannon)
x$simpson <- unlist(leg.list$simpson)
x$faith.pd <- unlist(leg.list$faith.pd)
x$wood.PC1 <- unlist(leg.list$wood.PC1)
x$wood.PC2 <- unlist(leg.list$wood.PC2)
x$wood.PC3 <- unlist(leg.list$wood.PC3)
pairs(x[,c("PC1","PC2","alloPC1","wood.dbs","richness","shannon","simpson",
           "faith.pd","wood.PC1","wood.PC2","wood.PC3")], 
      cor.method = "spearman")
# The strongest correlations appear to be with
# PC1 and Faith's PD (0.43), followed by Wood's PC2 (0.36), 
# richness (0.35) and Wood's DBS (0.34)

# Model the affect of foraging metrics on leg measures
# Comparisons among these models can be made based on Z values (effect size)
i <- 1e4-1
leg.forage.its.lm <- lm.rrpp(mm ~ log(its), data = leg.list, iter = i, print.progress = FALSE)
anova(leg.forage.its.lm)
#            Df      SS     MS     Rsq      F      Z Pr(>F)   
# log(its)    1 188.27 188.270 0.43698 160.66 4.8181  1e-04 ***
# Residuals 207 242.57   1.172 0.56302                         
# Total     208 430.84
leg.forage.sp.lm <- lm.rrpp(mm ~ species, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.forage.sp.lm)
#            Df     SS      MS     Rsq      F      Z   Pr(>F)    
# species     5 133.52 26.7049 0.30992 18.234 6.6769  1e-04 ***
# Residuals 203 297.31  1.4646 0.69008                         
# Total     208 430.84 
leg.forage.wood.dbs.lm <- lm.rrpp(mm ~ wood.dbs, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.forage.wood.dbs.lm)
#            Df      SS      MS    Rsq      F      Z Pr(>F)   
# wood.dbs    1  18.47 18.4679 0.04286 9.2704 2.4471 0.0014 **
# Residuals 207 412.37  1.9921 0.95714                        
# Total     208 430.84
leg.forage.richness.lm <- lm.rrpp(mm ~ richness, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.forage.richness.lm) 
#            Df     SS      MS    Rsq      F      Z Pr(>F)
# richness    1  25.25 25.2484 0.0586 12.886 2.7039  2e-04 ***
# Residuals 207 405.59  1.9594 0.9414                         
# Total     208 430.84
leg.forage.shannon.lm <- lm.rrpp(mm ~ shannon, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.forage.shannon.lm) 
#            Df     SS     MS     Rsq     F      Z Pr(>F)   
# shannon     1   2.84 2.8430 0.0066 1.375 0.83146 0.2323
# Residuals 207 427.99 2.0676 0.9934                     
# Total     208 430.84 
leg.forage.simpson.lm <- lm.rrpp(mm ~ simpson, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.forage.simpson.lm)
#            Df     SS     MS    Rsq      F      Z Pr(>F)    
# simpson     1   4.27 4.2671 0.0099 2.0707 1.1782 0.1365
# Residuals 207 426.57 2.0607 0.9901                     
# Total     208 430.84 
leg.forage.faith.lm <- lm.rrpp(mm ~ faith.pd, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.forage.faith.lm)
#            Df     SS     MS     Rsq      F      Z Pr(>F)  
# faith.pd    1  31.87 31.873 0.07398 16.537 2.8794  1e-04 ***
# Residuals 207 398.97  1.927 0.92602                         
# Total     208 430.84  
leg.forage.wood.PC1.lm <- lm.rrpp(mm ~ wood.PC1, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.forage.wood.PC1.lm)
#            Df     SS      MS     Rsq      F     Z Pr(>F)    
# wood.PC1    1  12.11 12.1149 0.02812 5.9891 2.086 0.0098 **
# Residuals 207 418.72  2.0228 0.97188                       
# Total     208 430.84  
leg.forage.wood.PC2.lm <- lm.rrpp(mm ~ wood.PC2, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.forage.wood.PC2.lm)
#            Df     SS     MS     Rsq      F      Z Pr(>F)    
# wood.PC2    1  32.41 32.410 0.07523 16.839 2.9633  1e-04 ***
# Residuals 207 398.43  1.925 0.92477                         
# Total     208 430.84
leg.forage.wood.PC3.lm <- lm.rrpp(mm ~ wood.PC3, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.forage.wood.PC3.lm)
#            Df     SS      MS     Rsq     F      Z Pr(>F)   
# wood.PC3    1  19.18 19.1808 0.04452 9.645 2.4735  8e-04 ***
# Residuals 207 411.66  1.9887 0.95548                        
# Total     208 430.84

# How do these factors do in models with ITS?
leg.size.size.sp.lm <- lm.rrpp(mm ~ log(its) + species, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.size.size.sp.lm)
#            Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)    1 188.27 188.270 0.43698 183.6058 4.9233  1e-04 ***
# species     5  35.44   7.087 0.08225   6.9119 4.5292  1e-04 ***
# Residuals 202 207.13   1.025 0.48076                           
# Total     208 430.84                                           
leg.size.size.wood.dbs.lm <- lm.rrpp(mm ~ log(its) + wood.dbs, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.size.size.wood.dbs.lm)
#            Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)    1 188.27 188.270 0.43698 162.6702 4.8295  1e-04 ***
# wood.dbs    1   4.15   4.150 0.00963   3.5857 1.7525 0.0394 *  
# Residuals 206 238.42   1.157 0.55338                           
# Total     208 430.84                                           
leg.size.size.richness.lm <- lm.rrpp(mm ~ log(its) + richness, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.size.size.richness.lm) 
#            Df     SS      MS     Rsq       F      Z Pr(>F)    
# log(its)    1 188.27 188.270 0.43698 171.810 4.8709  1e-04 ***
# richness    1  16.83  16.833 0.03907  15.361 3.0867  1e-04 ***
# Residuals 206 225.74   1.096 0.52394                          
# Total     208 430.84 
leg.size.size.shannon.lm <- lm.rrpp(mm ~ log(its) + shannon, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.size.size.shannon.lm) 
#            Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)    1 188.27 188.270 0.43698 165.2692 4.8407  1e-04 ***
# shannon     1   7.90   7.899 0.01833   6.9342 2.3919 0.0029 ** 
# Residuals 206 234.67   1.139 0.54468                           
# Total     208 430.84
leg.size.size.simpson.lm <- lm.rrpp(mm ~ log(its) + simpson, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.size.size.simpson.lm)
#            Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)    1 188.27 188.270 0.43698 162.0479 4.8251  1e-04 ***
# simpson     1   3.23   3.234 0.00751   2.7839 1.5155 0.0689 .  
# Residuals 206 239.33   1.162 0.55551                           
# Total     208 430.84
leg.size.size.faith.lm <- lm.rrpp(mm ~ log(its) + faith.pd, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.size.size.faith.lm)
#            Df     SS      MS     Rsq       F      Z Pr(>F)    
# log(its)    1 188.27 188.270 0.43698 175.861 4.8893  1e-04 ***
# faith.pd    1  22.03  22.033 0.05114  20.581 3.3215  1e-04 ***
# Residuals 206 220.54   1.071 0.51188                          
# Total     208 430.84
leg.size.size.wood.PC1.lm <- lm.rrpp(mm ~ log(its) + wood.PC1, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.size.size.wood.PC1.lm)
#            Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)    1 188.27 188.270 0.43698 162.9192 4.8300  1e-04 ***
# wood.PC1    1   4.51   4.514 0.01048   3.9064 1.8485   0.03 *  
# Residuals 206 238.05   1.156 0.55254                           
# Total     208 430.84
leg.size.size.wood.PC2.lm <- lm.rrpp(mm ~ log(its) + wood.PC2, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.size.size.wood.PC2.lm)
#            Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)    1 188.27 188.270 0.43698 162.9440 4.8306  1e-04 ***
# wood.PC2    1   4.55   4.550 0.01056   3.9383 1.8405 0.0303 *  
# Residuals 206 238.02   1.155 0.55245                           
# Total     208 430.84
leg.size.size.wood.PC3.lm <- lm.rrpp(mm ~ log(its) + wood.PC3, data = leg.list, iter = i, print.progress = FALSE)
anova(leg.size.size.wood.PC3.lm)
#            Df     SS      MS     Rsq        F      Z Pr(>F)    
# log(its)    1 188.27 188.270 0.43698 160.4258 4.8179  1e-04 ***
# wood.PC3    1   0.81   0.814 0.00189   0.6939 0.0863 0.4556    
# Residuals 206 241.75   1.174 0.56112                           
# Total     208 430.84  

# Recap of forage diversity factors (and species) by effect size
#            Df      SS     MS     Rsq        F      Z Pr(>F)    
# species     5  35.44   7.087 0.08225   6.9119 4.5292  1e-04 ***
# faith.pd    1  22.03  22.033 0.05114  20.581  3.3215  1e-04 ***
# richness    1  16.83  16.833 0.03907  15.361  3.0867  1e-04 ***
# shannon     1   7.90   7.899 0.01833   6.9342 2.3919 0.0029 ** 
# wood.PC1    1   4.51   4.514 0.01048   3.9064 1.8485   0.03 *  
# wood.PC2    1   4.55   4.550 0.01056   3.9383 1.8405 0.0303 *  
# wood.dbs    1   4.15   4.150 0.00963   3.5857 1.7525 0.0394 *  
# simpson     1   3.23   3.234 0.00751   2.7839 1.5155 0.0689 .  

# A plot to show the relationship of Faith's PD for forage data and leg morphology
plot(leg.list$faith.pd, leg.list$legPC1, col = (leg.list$class))
abline(lm(leg.list$legPC1~leg.list$faith.pd))

df <- data.frame(
  species = leg.list$species,
  class = leg.list$class,
  faith.pd = leg.list$faith.pd,
  richness = leg.list$richness,
  shannon = leg.list$shannon,
  PC1 = leg.list$legPC1
)

species.order <- forage$species[order(forage$faith.pd)]
df$species <- sub("bor","borealis",as.character(df$species))
df$species <- sub("bim","bimaculatus",as.character(df$species))
df$species <- sub("imp","impatiens",as.character(df$species))
df$species <- sub("vag","vagans",as.character(df$species))
df$species <- sub("terri","terricola",as.character(df$species))
df$species <- sub("tern","ternarius",as.character(df$species))
df$species <- factor(df$species, levels = species.order)

txt.df <- df %>%
  group_by(species) %>%
  summarise(faith.pd = median(faith.pd),
            .groups = "drop") %>%
  mutate(
    y =     c(6.9, 2,   6.4,-4.1, 4.7, 4.05),
    hjust = c(  0, 0.5, 0.5, 0, 1, 1)
  )

leg.specialization.plot <- df %>%
  ggplot(aes(x = faith.pd, y = PC1)) +
  theme_bw() +
  theme(legend.position="none") +
  geom_smooth(method=lm, color = "grey40", fill = "gray85") +
  geom_point(aes(color = species), alpha = 0.65, size = 2) +
  scale_color_viridis(discrete = TRUE, end = 0.95) +
  geom_text(data = txt.df, aes(x=faith.pd, y=y, hjust=hjust,
                               label = species, color = species)) +
  annotate("text", x=800, y=-10, hjust = 0, vjust = 1,
           color="grey40", size = 2.5,
           label="permANOVA with residual randomization
12 linear measurements of worker legs
Y ~ log(individual ITS) + species Faith's PD for forage plants
  ITS:  R^2 = 0.437, p < 10^-4
  FPD:  R^2 = 0.051, p < 10^-4") +
  labs(x="Faith's phylogenetic diversity index for forage plants\n(Data from Wood et al. 2019 Ecology)",
       y="leg morphology (PC1)")

leg.specialization.plot
ggsave("plots/leg.specialization.plot.pdf", leg.specialization.plot, width = 6.5, height = 5, scale = 1)

# Save everything
save.image("bombus.scaling.rda")
