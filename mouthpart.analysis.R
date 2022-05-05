# Bombus mouthpart morphometry analysis

# Clear memory
# rm(list=ls())
# gc() 

# devtools::install_github("aphanotus/borealis")
x <- c("borealis","tidyverse","ggpubr","ggrepel","factoextra","viridis",
       "RRPP","multcompView","ggbeeswarm","phytools","nlme")
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
# m <- read.mmm(
#   input.filename = 'Bombus.mouthparts.rawdata.210412.csv',
#   output.filename = 'Bombus.mouthpart.measures.210412.csv',
#   metadata.cols = "all",
#   measurement.col = NULL,
#   apply.scale = TRUE,
#   invert.scale = TRUE
# )
# 
# names(m)
# head(m$x)
# m$measurement.number
# m$specimen.number
# m$specimens.missing.scale
# m <- m$x

m <- read.csv('Bombus.mouthpart.measures.210412.csv')

m <- m %>%
  mutate(
    tongue = m01 + m03,    # Calculate tongue length
    body = body / scale1,  # Calculate metric body length
    its = its / scale2,    # Calculate metric intertegular span (ITS)
    sp.caste = paste(species, caste, sep=".") # Create a factor combining species and caste
  ) %>%
  mutate(specimen_id = paste(specimen_id, sp.caste, sep=".")) %>% 
  select(-scale, -scale1, -scale2, -notes, -m10)

data.cols <- 2:15

# Workers only
w <- filter(m, caste == "W")

####################
# Scaling plots
####################

# Filter out low sample size (<3)
barplot(sort(c(with(m, by(sp.caste, sp.caste, length)))), cex.names = 0.5)
abline(h=c(3,10), col = "darkred")
m <- m %>% group_by(sp.caste) %>% filter(n() >= 3)
barplot(sort(c(with(w, by(species, species, length)))), cex.names = 0.5)
abline(h=c(3,10), col = "darkred")
w <- w %>% group_by(species) %>% filter(n() >= 3)

scaling.plot(
  x = log10(m$body), y = log10(m$tongue),
  group = m$sp.caste,
  xlab = "body length (log10 mm)", ylab = "tongue length (log10 mm)",
  include.legend = TRUE,
  isometry.line = TRUE,
  convex.hulls = FALSE,
  groups.trendlines = TRUE,
  save.as = "plots/scaling.tongue.v.body.pdf"
)
#      group  n   slope        p sig   ci.lo  ci.hi spans.zero
# 1    bor.W  6  1.1700 1.62e-02 *    0.3580 1.9900      FALSE
# 2    imp.Q  3  1.1200 3.38e-01     -7.2100 9.4500       TRUE
# 3   ferv.W  4  0.9720 1.89e-01     -1.1600 3.1100       TRUE
# 4    imp.W 37  0.7380 3.67e-09 ***  0.5460 0.9310      FALSE
# 5    bim.W 56  0.6690 1.50e-12 ***  0.5220 0.8160      FALSE
# 6    vag.W 33  0.4920 2.49e-02 *    0.0664 0.9180      FALSE
# 7    san.W  4  0.4910 9.52e-02 .   -0.2120 1.1900       TRUE
# 8   tern.Q  9  0.3020 2.64e-01     -0.2860 0.8900       TRUE
# 9   tern.W 52  0.1120 2.22e-01     -0.0701 0.2950       TRUE
# 10   bim.M 14  0.0955 4.66e-01     -0.1810 0.3720       TRUE
# 11   vag.Q  3  0.0931 6.01e-01     -1.5400 1.7300       TRUE
# 12   bor.Q  3 -0.0555 6.74e-02 .   -0.1300 0.0195       TRUE
# 13 terri.W 13 -0.0705 5.54e-01     -0.3250 0.1840       TRUE
# 14   imp.M  4 -0.3940 7.43e-01     -4.9100 4.1200       TRUE
# Differences in slope based on CI
# bim.M bim.W bor.Q bor.W ferv.W imp.M imp.Q imp.W san.W tern.Q tern.W terri.W vag.Q vag.W
#     a     b     a     b     ab    ab    ab     b    ab     ab      a       a    ab    ab 
# For B. bimaculatus and B. borealis, reproductive castes have a significantly
# shallower tongue scaling allometry with body size than do workers. This trend
# also appears in B. vagans. It's definitely not seen in B. ternarius, where we
# have good sample sizes for both workers and queens.

scaling.plot(
  x = log10(m$its), y = log10(m$tongue),
  group = m$sp.caste,
  xlab = "intertegular span (log10 mm)", ylab = "tongue length (log10 mm)",
  include.legend = TRUE,
  isometry.line = TRUE,
  convex.hulls = FALSE,
  groups.trendlines = TRUE,
  save.as = "plots/scaling.tongue.v.its.pdf"
)
#      group  n   slope        p sig     ci.lo  ci.hi spans.zero
# 1    bor.W  6  1.2400 1.26e-03 **    0.81500  1.670      FALSE
# 2    imp.M  4  0.9430 1.29e-01      -0.67200  2.560       TRUE
# 3   ferv.W  4  0.7950 5.90e-01      -4.58000  6.170       TRUE
# 4    san.W  4  0.7460 8.05e-02 .    -0.22400  1.720       TRUE
# 5    vag.W 33  0.7430 8.34e-03 **    0.20500  1.280      FALSE
# 6    imp.W 39  0.7100 5.37e-09 ***   0.51900  0.900      FALSE
# 7    bim.W 56  0.6970 1.05e-12 ***   0.54600  0.848      FALSE
# 8    vag.Q  3  0.3030 5.06e-02 .    -0.00356  0.610       TRUE
# 9   tern.Q  9  0.2000 4.91e-01      -0.45200  0.853       TRUE
# 10 terri.W 13  0.1230 2.72e-01      -0.11100  0.357       TRUE
# 11  tern.W 52  0.1100 2.25e-01      -0.06980  0.290       TRUE
# 12   bim.M 13  0.0214 8.78e-01      -0.27900  0.322       TRUE
# 13   bor.Q  3 -0.0511 6.38e-01      -1.07000  0.967       TRUE
# 14   imp.Q  3 -2.9800 7.47e-01     -93.30000 87.400       TRUE
# bim.M bim.W bor.Q bor.W ferv.W imp.M imp.Q imp.W san.W tern.Q tern.W terri.W vag.Q vag.W
#     a     b    ab     b     ab    ab    ab     b    ab     ab      a       a    ab    ab 
# Same story as with body size

# Workers
scaling.plot(
  x = log10(w$body), y = log10(w$tongue),
  group = w$species,
  xlab = "body length (log10 mm)", ylab = "tongue length (log10 mm)",
  include.legend = TRUE,
  isometry.line = TRUE,
  convex.hulls = FALSE,
  groups.trendlines = TRUE,
  save.as = "plots/scaling.tongue.v.body.workers.pdf"
)
#   group  n   slope        p sig   ci.lo ci.hi spans.zero
# 1   bor  6  1.1700 1.62e-02 *    0.3580 1.990      FALSE
# 2  ferv  4  0.9720 1.89e-01     -1.1600 3.110       TRUE
# 3   imp 37  0.7380 3.67e-09 ***  0.5460 0.931      FALSE
# 4   bim 56  0.6690 1.50e-12 ***  0.5220 0.816      FALSE
# 5   vag 33  0.4920 2.49e-02 *    0.0664 0.918      FALSE
# 6   san  4  0.4910 9.52e-02 .   -0.2120 1.190       TRUE
# 7  tern 52  0.1120 2.22e-01     -0.0701 0.295       TRUE
# 8 terri 13 -0.0705 5.54e-01     -0.3250 0.184       TRUE
# Comparisons of slopes, based on CI
# bim bor ferv imp san tern terri vag
#   a   a   ab   a  ab    b     b  ab
# B. ternarius and B. terricola are significantly difference from 
# B. bimaculatus, B. borealis and B. impatiens.

k.w.tongue <- scaling.plot(
  x = log10(w$its), y = log10(w$tongue),
  group = w$species,
  xlab = "intertegular span (log10 mm)", ylab = "tongue length (log10 mm)",
  include.legend = TRUE,
  isometry.line = TRUE,
  convex.hulls = TRUE,
  save.as = "plots/scaling.tongue.v.its.workers.pdf"
)
k.w.tongue$plot
k.w.tongue$slopes
#   group  n slope        p sig   ci.lo ci.hi spans.zero
# 1   bor  6 1.240 1.26e-03 **   0.8150 1.670      FALSE
# 2  ferv  4 0.795 5.90e-01     -4.5800 6.170       TRUE
# 3   san  4 0.746 8.05e-02 .   -0.2240 1.720       TRUE
# 4   vag 33 0.743 8.34e-03 **   0.2050 1.280      FALSE
# 5   imp 39 0.710 5.37e-09 ***  0.5190 0.900      FALSE
# 6   bim 56 0.697 1.05e-12 ***  0.5460 0.848      FALSE
# 7 terri 13 0.123 2.72e-01     -0.1110 0.357       TRUE
# 8  tern 52 0.110 2.25e-01     -0.0698 0.290       TRUE

k.w.tongue$slopes$group <- factor(k.w.tongue$slopes$group, levels = k.w.tongue$slopes$group)

ylimits <- sort(k.w.tongue$slopes$ci.lo)[2]
ylimits <- c(ylimits, sort(k.w.tongue$slopes$ci.hi, decreasing = TRUE)[2])
ylimits[1] <- ylimits[1]*1.1
k.w.tongue$slopes %>% ggplot(aes(x=group, y=slope, fill=group, label=paste0("n=",n))) +
  theme_bw() +
  theme(legend.position="none") +
  geom_bar(stat="identity", color="grey15") +
  geom_errorbar(aes(ymin=ci.lo, ymax=ci.hi), width=.2) +
  geom_vline(xintercept = 2) +
  scale_fill_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(y=ylimits[1], size=3, color="grey50", vjust=1) +
  scale_y_continuous(limits=ylimits) +
  xlab("species") +
  ylab("scaling coefficent (k) for tongue vs. intertegular span")
ggsave("plots/scaling.tongue.by.species.pdf", width = 6, height = 4, scale = 1)

# PCA
m.pca <- prcomp(m[,data.cols], center = TRUE, scale = TRUE)
shape.space(m.pca, group = m$sp.caste, convex.hulls = TRUE,
               save.as = "plots/mouthpart.morphospace.all.pdf")
m$PC1 <- get_pca_ind(m.pca)$coord[,1]
m$PC2 <- get_pca_ind(m.pca)$coord[,2]

w.pca <- prcomp(w[,data.cols], center = TRUE, scale = TRUE)
shape.space(w.pca, group = w$species, convex.hulls = TRUE,
               save.as = "plots/mouthpart.morphospace.workers.pdf")
w$PC1 <- get_pca_ind(w.pca)$coord[,1]
w$PC2 <- get_pca_ind(w.pca)$coord[,2]

# A quick look at how shape (as PC1) varies by body size (ITS) and species
# Remember, this is just a univariate representation of mouthpart "shape". 
# This issue is revisited below with multivariate stats.
scaling.plot(
  x = log10(w$its), y = w$PC1,
  group = w$species,
  xlab = "intertegular span (log10 mm)", 
  ylab = "PC1 for worker mouthpart measures",
  include.legend = TRUE,
  fixed.aspect = FALSE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  save.as = "plots/worker.mouhtpart.PC1.v.its.pdf"
)

# A quick plot of mouthpart "shape" (as PC1) by species
species.order <- w %>%
  group_by(species) %>%
  summarise(medianPC1 = median(PC1)) %>%
  arrange(medianPC1) %>% pull(species)
sp <- factor(w$species, levels = species.order)
boxplot(w$PC1 ~ sp)
# The generalists (impatiens & bimaculatus) overlap everyone!

####################
# Stats using RRPP
####################

# Exclude specimens with missing ITS measurements
w <- w[which(!is.na(w$its)),]

w.list <- list(
  mm = as.matrix(w[,data.cols]),
  body = w$body,
  its = w$its,
  species = factor(w$species, levels = species.order)
)

w.rrpp.species <- lm.rrpp(mm ~ log(its) + species, data = w.list, iter = 1e4)
anova(w.rrpp.species)
#            Df      SS      MS     Rsq       F      Z    Pr(>F)
# log(its)    1 102.756 102.756 0.39144 267.297  5.4566 9.999e-05 ***
# species     7  83.638  11.948 0.31861  31.081 10.6570 9.999e-05 ***
# Residuals 198  76.117   0.384 0.28996                              
# Total     206 262.511      

# Pairwise species differences?
w.rrpp.size <- lm.rrpp(mm ~ log(its), data = w.list, iter = 1e4)
w.species.pw <- pairwise(fit = w.rrpp.species,
                         fit.null = w.rrpp.size,
                         groups = w.list$species)
(pw.output <- summary(w.species.pw))
(sig.letters <- pairwise.group.comparisons(pw.output))
# terri  tern   imp    san   vag   bim   bor  ferv 
#   "a"  "ab"   "c" "bcde"  "cd"   "e"   "e"  "de" 

df <- data.frame(
  x = 1:8,
  y = rep(11.5,8),
  label = sig.letters$Letters
)

species.order <- w %>%
  group_by(species) %>%
  summarise(medianPC1 = median(PC1)) %>%
  arrange(medianPC1) %>% pull(species)
sp <- factor(w$species, levels = species.order)

ggplot(w, aes(sp, PC1, color = sp)) +
  theme_minimal() +
  theme(legend.position="none") +
  geom_violin() +
  geom_quasirandom(varwidth = TRUE) +
  scale_color_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(data = df, aes(x = x, y = y, label = label), color = "black", size = 5 ) +
  xlab("species") +
  ylab("PC1 for worker mouthpart measurements") +
  labs(caption = "Significant differences in the influence of species in the model  Y ~ intertegular span + species.\nPERMANOVA using `lm.rrpp`; contrasts using `pairwise` in the package `RRPP`.")
ggsave("plots/worker.PC1.by.species.pdf", width = 6, height = 4, scale = 1)

# Allometry-corrected PCA
w.allo.pca <- prcomp(resid(w.rrpp.size), center = TRUE, scale = TRUE)
shape.space(w.allo.pca, group = w.list$species, convex.hulls = TRUE,
               main.title = "Size-corrected morphospace for worker mouthparts",
               save.as = "plots/allo-corrected.pca.workers.pdf")
w$alloPC1 <- get_pca_ind(w.allo.pca)$coord[,1]

species.order <- w %>%
  group_by(species) %>%
  summarise(medianAlloPC1 = median(alloPC1)) %>%
  arrange(medianAlloPC1) %>% pull(species)
w$species <- factor(w$species, levels = species.order)

# Repeat this just to get the species in the right order!
w.list <- list(
  specimen_id = w$specimen_id,
  mm = as.matrix(w[,data.cols]),
  body = w$body,
  its = w$its,
  species = w$species,
  tongue = w$tongue,
  PC1 = w$PC1,
  PC2 = w$PC2,
  alloPC1 = w$alloPC1
)
w.rrpp.size <- lm.rrpp(mm ~ log(its), data = w.list, iter = 1e4)
w.rrpp.species <- lm.rrpp(mm ~ log(its) + species, data = w.list, iter = 1e4)
w.species.pw <- pairwise(fit = w.rrpp.species,
                         fit.null = w.rrpp.size,
                         groups = w.list$species)
(pw.output <- summary(w.species.pw))
(sig.letters <- pairwise.group.comparisons(pw.output))
# terri   tern    imp    san    vag    bim    bor   ferv 
#   "a"   "ab"    "c" "bcde"   "cd"    "e"    "e"   "de"

df <- data.frame(
  x = 1:8,
  y = rep(10, 8),
  label = sig.letters$Letters
)

ggplot(w, aes(species, alloPC1, color = species)) +
  theme_minimal() +
  theme(legend.position="none") +
  geom_violin() +
  geom_quasirandom(varwidth = TRUE, alpha = 0.9) +
  scale_color_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(data = df, aes(x = x, y = y, label = label), color = "black", size = 5 ) +
  xlab("species") +
  ylab("size-corrected PC1 for worker mouthpart measures") +
  labs(caption = "Significant differences in the influence of species in the model  Y ~ intertegular span + species.\nPERMANOVA using `lm.rrpp`; contrasts using `pairwise` in the package `RRPP`.")
ggsave("plots/size-corrected.worker.PC1.by.species.pdf", width = 6, height = 4, scale = 1.25)

# Do species have different worker mouthpart "shape" allometries?
w.rrpp.unique.allometries <- lm.rrpp(mm ~ log(its) * species, data = w.list, iter = 1e4-1)
anova(w.rrpp.unique.allometries)
#                   Df      SS      MS     Rsq        F       Z    Pr(>F)    
# log(its)           1 102.756 102.756 0.39144 321.7703  5.6901  1e-04 ***
# species            7  83.638  11.948 0.31861  37.4146 10.8693  1e-04 ***
# log(its):species   7  15.121   2.160 0.05760   6.7645  3.9175  1e-04 ***
# Residuals        191  60.995   0.319 0.23235                            
# Total            206 262.511    
# Yes, species have different worker mouthpart "shape" allometries

# Post hoc examination of the interaction
w.rrpp.unique.allometries.pw <- 
  pairwise(fit = w.rrpp.unique.allometries,
           fit.null = w.rrpp.species,
           groups = w.list$species,
           covariate = w.list$its)
(allometry.pw.output <- summary(w.rrpp.unique.allometries.pw))
#                    d UCL (95%)          Z Pr > d
# terri:tern 0.1388254 0.8436683 -1.3840678 0.9117
# terri:imp  0.9475433 0.7825601  1.9271663 0.0205
# terri:san  1.0222129 1.3026127  1.2905934 0.1103
# terri:vag  1.2863968 0.9451010  2.1100826 0.0081
# terri:bim  1.0429401 0.7590406  2.1242520 0.0103
# terri:bor  1.6728389 0.8793644  2.6046661 0.0004
# terri:ferv 1.6634873 3.2464937  0.6679845 0.2704
# tern:imp   0.9470835 0.6067240  2.4265097 0.0014
# tern:san   1.0401159 1.2673791  1.3569782 0.0925
# tern:vag   1.2840508 0.7090970  2.5871643 0.0007
# tern:bim   1.0451426 0.5203583  2.7868037 0.0002
# tern:bor   1.6686259 0.7285047  2.9314483 0.0001
# tern:ferv  1.6587789 3.1953937  0.6857194 0.2666
# imp:san    0.2363699 1.1398063 -0.7311067 0.7403
# imp:vag    0.3701766 0.7334313  0.4961793 0.3389
# imp:bim    0.1413194 0.4588172 -0.3401741 0.6094
# imp:bor    0.7380433 0.6329457  1.8691713 0.0254
# imp:ferv   0.7870847 3.1972704 -0.4760302 0.6626
# san:vag    0.3850706 1.2974704 -0.2372986 0.5785
# san:bim    0.2046172 1.1697373 -1.0184933 0.8235
# san:bor    0.7240761 1.2053246  0.8872187 0.2082
# san:ferv   0.7451681 3.3095255 -0.6689656 0.7231
# vag:bim    0.2900842 0.6625513  0.2011299 0.4364
# vag:bor    0.4214576 0.8382610  0.5171682 0.3276
# vag:ferv   0.4972743 3.2214470 -1.3004579 0.8981
# bim:bor    0.6397123 0.6411345  1.6452144 0.0503
# bim:ferv   0.7041196 3.1819997 -0.6640716 0.7244
# bor:ferv   0.4858184 3.2237448 -1.3235672 0.9032
(sig.letters <- pairwise.group.comparisons(allometry.pw.output))
# terri  tern   imp   san   vag   bim   bor  ferv 
#   "a"   "a"   "b" "abc"  "bc"  "bc"   "c" "abc" 
# Interesting differences among the allometries!

# I confirmed with Michael Collyer that this modeling approach and use of 
# rrpp::pairwise accomplishes what we're after here.

####################
# PGLS
####################
# Unfortunately, geomorph::procD.pgls only takes a 3D shape array, not 
# a table of multivariate data. I can't find other implementations of PGLS
# that works with a continuous multivariate independent dataset 
# So, I'll use the simple phylogenetic regression in nlme::gls

m.species.names <- unique(sort(as.character(w.list$species)))
m.species.names <- sub("bim","bimac",m.species.names)
m.species.names <- sub("san","sande",m.species.names)

btree <- Bombus.tree
btree$tip.label <- btree$code.name
btree <- keep.tip(btree, btree$tip.label[which(btree$tip.label %in% m.species.names)])
btree$tip.label <- sub("bimac","bim",btree$tip.label)
btree$tip.label <- sub("sande","san",btree$tip.label)
plot(btree)

# Create a dataframe with species means
w.sp <- data.frame(
  PC1 = c(by(w$PC1, w$species, mean, na.rm=TRUE)),
  PC2 = c(by(w$PC2, w$species, mean, na.rm=TRUE)),
  alloPC1 = c(by(w$alloPC1, w$species, mean, na.rm=TRUE)),
  body = c(by(w$body, w$species, mean, na.rm=TRUE)),
  its = c(by(w$its, w$species, mean, na.rm=TRUE))
)
# Reorder the dataframe to match the order of tree tips
w.sp <- w.sp[btree$tip.label,]

# Check that the names in the tree and dataset match
geiger::name.check(btree, w.sp)

pgls.PC1.by.ITS <- gls(PC1 ~ log(its), 
                       correlation = corBrownian(phy = btree), 
                       method = "REML", data = w.sp)
summary(pgls.PC1.by.ITS)
with(w.sp, plot(PC1 ~ log(its)))
abline(a = coef(pgls.PC1.by.ITS)[1], b = coef(pgls.PC1.by.ITS)[2])
# Generalized least squares fit by REML
#      AIC      BIC    logLik
# 29.02208 28.39736 -11.51104
# 
#                 Value Std.Error   t-value p-value
# (Intercept) -28.64341  7.217637 -3.968531  0.0074
# log(its)     19.67571  4.555605  4.319011  0.0050
# Degrees of freedom: 8 total; 6 residual

pgls.PC2.by.ITS <- gls(PC2 ~ log(its), 
                       correlation = corBrownian(phy = btree), 
                       method = "REML", data = w.sp)
summary(pgls.PC2.by.ITS)
with(w.sp, plot(PC2 ~ log(its)))
abline(a = coef(pgls.PC2.by.ITS)[1], b = coef(pgls.PC2.by.ITS)[2])
# Generalized least squares fit by REML
#      AIC      BIC    logLik
# 21.49779 20.87307 -7.748894
# 
#                 Value Std.Error   t-value p-value
# (Intercept) -4.275946  3.855511 -1.109048  0.3099
# log(its)     2.691264  2.433509  1.105919  0.3111
# Degrees of freedom: 8 total; 6 residual

####################
# Correlations to forage diversity
####################
( forage <- read.csv("Wood.et.al.2019.forage.diversity.csv") )

# Make the species abbreviations (code names) match
w.list$code.name <- as.character(w.list$species)
forage$code.name <- sub("bimac","bim",forage$code.name)

# Filter out the foraging information that covers species not in out dataset
x <- which(forage$code.name %in% unique(w.list$code.name))
forage <- forage[x,]

# Filter out species in our dataset not covered by the foraging information
x <- which(!(w.list$code.name %in% forage$code.name))
unique(w.list$code.name[x])
wf <- w.list
wf$specimen_id <- wf$specimen_id[-x]
wf$mm <- wf$mm[-x,]
wf$body <- wf$body[-x]
wf$its <- wf$its[-x]
wf$species <- wf$species[-x]
wf$tongue <- wf$tongue[-x]
wf$PC1 <- wf$PC1[-x]
wf$PC2 <- wf$PC2[-x]
wf$alloPC1 <- wf$alloPC1[-x]
wf$code.name <- wf$code.name[-x]

# Add the foraging information to the new list object
wf$class <- wf$code.name
wf$wood.dbs <- wf$code.name
wf$richness <- wf$code.name
wf$shannon <- wf$code.name
wf$simpson <- wf$code.name
wf$faith.pd <- wf$code.name
wf$wood.PC1 <- wf$code.name
wf$wood.PC2 <- wf$code.name
wf$wood.PC3 <- wf$code.name
df <- forage[,-c(1:2)]
colnames(df)[1] <- "class"
row.names(df) <- forage$code.name

for (i in 1:length(wf$code.name)) {
  wf$class[i] <- df[wf$code.name[i],"class"]
  wf$wood.dbs[i] <- df[wf$code.name[i],"wood.dbs"]
  wf$richness[i] <- df[wf$code.name[i],"richness"]
  wf$shannon[i] <- df[wf$code.name[i],"shannon"]
  wf$simpson[i] <- df[wf$code.name[i],"simpson"]
  wf$faith.pd[i] <- df[wf$code.name[i],"faith.pd"]
  wf$wood.PC1[i] <- df[wf$code.name[i],"wood.PC1"]
  wf$wood.PC2[i] <- df[wf$code.name[i],"wood.PC2"]
  wf$wood.PC3[i] <- df[wf$code.name[i],"wood.PC3"]
}
wf$class <-  as.factor(wf$class)
wf$wood.dbs <-  as.numeric(wf$wood.dbs)
wf$richness <-  as.numeric(wf$richness)
wf$shannon <-  as.numeric(wf$shannon)
wf$simpson <-  as.numeric(wf$simpson)
wf$faith.pd <-  as.numeric(wf$faith.pd)
wf$wood.PC1 <-  as.numeric(wf$wood.PC1)
wf$wood.PC2 <-  as.numeric(wf$wood.PC2)
wf$wood.PC3 <-  as.numeric(wf$wood.PC3)

# Explore potential correlations
x <- as.data.frame(wf$mm)
x$body <- unlist(wf$body)
x$its <- unlist(wf$its)
x$PC1 <- unlist(wf$PC1)
x$PC2 <- unlist(wf$PC2)
x$alloPC1 <- unlist(wf$alloPC1)
x$code.name <- unlist(wf$code.name)
x$class <- unlist(wf$class)
x$wood.dbs <- unlist(wf$wood.dbs)
x$richness <- unlist(wf$richness)
x$shannon <- unlist(wf$shannon)
x$simpson <- unlist(wf$simpson)
x$faith.pd <- unlist(wf$faith.pd)
x$wood.PC1 <- unlist(wf$wood.PC1)
x$wood.PC2 <- unlist(wf$wood.PC2)
x$wood.PC3 <- unlist(wf$wood.PC3)
pairs(x[,c("PC1","PC2","alloPC1","wood.dbs","richness","shannon","simpson",
           "faith.pd","wood.PC1","wood.PC2","wood.PC3")], 
      cor.method = "spearman")
# The strongest correlations appear to be with
# PC1 and Simpson's diversity (0.49), Wood's forage plant PC2 (0.47),
# and Wood's DBS (0.42)

# Model the affect of foraging metrics on mouthpart shape
# Comparisons among these models can be made based on Z values (effect size)
# Species should be a strong predictor, but the forage plant metric with the
# next best effect will be informative.
i <- 1e4-1
wf.its.lm <- lm.rrpp(mm ~ log(its), data = wf, iter = i, print.progress = FALSE)
anova(wf.its.lm)
#            Df      SS     MS     Rsq      F      Z Pr(>F)   
# log(its)    1  96.316 96.316 0.37725 121.76 4.8235  1e-04 ***
# Residuals 201 158.994  0.791 0.62275                         
# Total     202 255.310 
wf.sp.lm <- lm.rrpp(mm ~ species, data = wf, iter = i, print.progress = FALSE)
anova(wf.sp.lm)
#            Df     SS      MS     Rsq     F      Z Pr(>F)    
# species     6 134.54 22.4230 0.52696 36.39 10.243  1e-04 ***
# Residuals 196 120.77  0.6162 0.47304                        
# Total     202 255.31 
wf.wood.dbs.lm <- lm.rrpp(mm ~ wood.dbs, data = wf, iter = i, print.progress = FALSE)
anova(wf.wood.dbs.lm)
#            Df      SS      MS    Rsq      F      Z Pr(>F)   
# wood.dbs    1  10.365 10.3654 0.0406 8.5058 2.2388 0.0033 **
# Residuals 201 244.945  1.2186 0.9594                        
# Total     202 255.310 
wf.richness.lm <- lm.rrpp(mm ~ richness, data = wf, iter = i, print.progress = FALSE)
anova(wf.richness.lm) 
#            Df     SS      MS     Rsq      F        Z Pr(>F)
# richness    1   0.49 0.49022 0.00192 0.3867 -0.10806 0.5463
# Residuals 201 254.82 1.26776 0.99808                       
# Total     202 255.31
wf.shannon.lm <- lm.rrpp(mm ~ shannon, data = wf, iter = i, print.progress = FALSE)
anova(wf.shannon.lm) 
#            Df      SS     MS     Rsq     F      Z Pr(>F)   
# shannon     1   9.004 9.0042 0.03527 7.348 2.1545 0.0065 **
# Residuals 201 246.306 1.2254 0.96473                       
# Total     202 255.310
wf.simpson.lm <- lm.rrpp(mm ~ simpson, data = wf, iter = i, print.progress = FALSE)
anova(wf.simpson.lm)
#            Df      SS      MS     Rsq      F      Z Pr(>F)    
# simpson     1  25.917 25.9172 0.10151 22.709 3.0802  1e-04 ***
# Residuals 201 229.393  1.1413 0.89849                         
# Total     202 255.310 
wf.faith.lm <- lm.rrpp(mm ~ faith.pd, data = wf, iter = i, print.progress = FALSE)
anova(wf.faith.lm)
#            Df      SS     MS     Rsq      F      Z Pr(>F)  
# faith.pd    1   4.424 4.4244 0.01733 3.5446 1.5266 0.0606 .
# Residuals 201 250.886 1.2482 0.98267                       
# Total     202 255.310  
wf.wood.PC1.lm <- lm.rrpp(mm ~ wood.PC1, data = wf, iter = i, print.progress = FALSE)
anova(wf.wood.PC1.lm)
#            Df     SS      MS     Rsq      F      Z Pr(>F)    
# wood.PC1    1  18.02 18.0202 0.07058 15.264 2.7627  1e-04 ***
# Residuals 201 237.29  1.1805 0.92942                         
# Total     202 255.31
wf.wood.PC2.lm <- lm.rrpp(mm ~ wood.PC2, data = wf, iter = i, print.progress = FALSE)
anova(wf.wood.PC2.lm)
#            Df     SS     MS    Rsq      F      Z Pr(>F)    
# wood.PC2    1  49.58 49.580 0.1942 48.441 3.7701  1e-04 ***
# Residuals 201 205.73  1.024 0.8058                         
# Total     202 255.31 
wf.wood.PC3.lm <- lm.rrpp(mm ~ wood.PC3, data = wf, iter = i, print.progress = FALSE)
anova(wf.wood.PC3.lm)
#            Df      SS      MS     Rsq      F      Z Pr(>F)   
# wood.PC3    1  11.252 11.2516 0.04407 9.2665 2.2965 0.0023 **
# Residuals 201 244.059  1.2142 0.95593                        
# Total     202 255.310 

# Summary of forage modeling
# - Most metrics of forage plant diversity are reasonable predictors of
#   mouthpart shape! 
# - PC2 from Wood et al. 2019's forage plant data is the best predictor 
#   of mouthpart shape, after species membership and individual ITS.
# - Among the classic diversity metrics, Simpson's, is the best predictor.
# - Only richness does a poor job of predicting mouthpart shape. This
#   could be problematic for the overall hypothesis, since richness in 
#   forage plants might be considered an indicator of "generalization" 
#   in a pollinator (contra "specialization").
# - However, Wood et al.'s dietary breadth score (DBS), which is intended
#   as a metric of generalization/specialization, is a modest predictor
#   of mouthpart shape.
#   
# The table below just repeats the key lines from the individual ANOVA tables above
# and orders them based on the highest effect size (Z).
#          Df      Rsq       F         Z     Pr(>F)  
# species   6  0.52696   36.39   10.243    1.00E-04 ***
# log(its)  1  0.37725  121.76    4.8235   1.00E-04 ***
# wood.PC2  1  0.1942    48.441   3.7701   1.00E-04 ***
# simpson   1  0.10151   22.709   3.0802   1.00E-04 ***
# wood.PC1  1  0.07058   15.264   2.7627   1.00E-04 ***
# wood.PC3  1  0.04407    9.2665  2.2965   0.0023   **
# wood.dbs  1  0.0406     8.5058  2.2388   0.0033   **
# shannon   1  0.03527    7.348   2.1545   0.0065   **
# faith.pd  1  0.01733    3.5446  1.5266   0.0606   .
# richness  1  0.00192    0.3867 -0.10806  0.5463  

# How do these factors do in models with ITS?
wf.size.sp.lm <- lm.rrpp(mm ~ log(its) + species, data = wf, iter = i, print.progress = FALSE)
anova(wf.size.sp.lm)
#            Df      SS     MS     Rsq       F       Z Pr(>F)    
# log(its)    1  96.316 96.316 0.37725 248.198  5.5726  1e-04 ***
# species     6  83.322 13.887 0.32636  35.785 10.0512  1e-04 ***
# Residuals 195  75.672  0.388 0.29639                           
# Total     202 255.310     
wf.size.wood.dbs.lm <- lm.rrpp(mm ~ log(its) + wood.dbs, data = wf, iter = i, print.progress = FALSE)
anova(wf.size.wood.dbs.lm)
#            Df      SS     MS     Rsq        F      Z Pr(>F)    
# log(its)    1  96.316 96.316 0.37725 121.8737 4.8237  1e-04 ***
# wood.dbs    1   0.935  0.935 0.00366   1.1832 0.6957 0.2764    
# Residuals 200 158.059  0.790 0.61909                           
# Total     202 255.310 
wf.size.richness.lm <- lm.rrpp(mm ~ log(its) + richness, data = wf, iter = i, print.progress = FALSE)
anova(wf.size.richness.lm) 
#            Df      SS     MS     Rsq       F      Z Pr(>F)    
# log(its)    1  96.316 96.316 0.37725 121.898 4.8257  1e-04 ***
# richness    1   0.967  0.967 0.00379   1.224 0.7241 0.2655    
# Residuals 200 158.027  0.790 0.61896                          
# Total     202 255.310    
wf.size.shannon.lm <- lm.rrpp(mm ~ log(its) + shannon, data = wf, iter = i, print.progress = FALSE)
anova(wf.size.shannon.lm) 
#            Df      SS     MS     Rsq       F      Z Pr(>F)    
# log(its)    1  96.316 96.316 0.37725 127.977 4.8758  1e-04 ***
# shannon     1   8.473  8.473 0.03319  11.259 2.5314  0.001 ** 
# Residuals 200 150.521  0.753 0.58956                          
# Total     202 255.310 
wf.size.simpson.lm <- lm.rrpp(mm ~ log(its) + simpson, data = wf, iter = i, print.progress = FALSE)
anova(wf.size.simpson.lm)
#            Df      SS     MS     Rsq       F      Z Pr(>F)    
# log(its)    1  96.316 96.316 0.37725 130.724 4.8977  1e-04 ***
# simpson     1  11.636 11.636 0.04558  15.793 2.7950  2e-04 ***
# Residuals 200 147.358  0.737 0.57717                          
# Total     202 255.310   
wf.size.faith.lm <- lm.rrpp(mm ~ log(its) + faith.pd, data = wf, iter = i, print.progress = FALSE)
anova(wf.size.faith.lm)
#            Df      SS     MS     Rsq        F      Z Pr(>F)    
# log(its)    1  96.316 96.316 0.37725 121.6362 4.8224  1e-04 ***
# faith.pd    1   0.626  0.626 0.00245   0.7911 0.3746 0.3831    
# Residuals 200 158.368  0.792 0.62029                           
# Total     202 255.310      
wf.size.wood.PC1.lm <- lm.rrpp(mm ~ log(its) + wood.PC1, data = wf, iter = i, print.progress = FALSE)
anova(wf.size.wood.PC1.lm)
#            Df      SS     MS     Rsq       F      Z Pr(>F)    
# log(its)    1  96.316 96.316 0.37725 148.610 5.0303  1e-04 ***
# wood.PC1    1  29.371 29.371 0.11504  45.318 3.7177  1e-04 ***
# Residuals 200 129.623  0.648 0.50771                          
# Total     202 255.310     
wf.size.wood.PC2.lm <- lm.rrpp(mm ~ log(its) + wood.PC2, data = wf, iter = i, print.progress = FALSE)
anova(wf.size.wood.PC2.lm)
#            Df      SS     MS     Rsq       F      Z Pr(>F)    
# log(its)    1  96.316 96.316 0.37725 137.499 4.9485  1e-04 ***
# wood.PC2    1  18.897 18.897 0.07402  26.977 3.2780  1e-04 ***
# Residuals 200 140.097  0.700 0.54873                          
# Total     202 255.310  
wf.size.wood.PC3.lm <- lm.rrpp(mm ~ log(its) + wood.PC3, data = wf, iter = i, print.progress = FALSE)
anova(wf.size.wood.PC3.lm)
#            Df      SS     MS     Rsq        F      Z Pr(>F)    
# log(its)    1  96.316 96.316 0.37725 121.5015 4.8205  1e-04 ***
# wood.PC3    1   0.451  0.451 0.00177   0.5687 0.1264  0.466    
# Residuals 200 158.543  0.793 0.62098                           
# Total     202 255.310   

# Recap of forage diversity factors (and species) by effect size
#            Df      SS     MS     Rsq       F       Z Pr(>F)    
# species     6  83.322 13.887 0.32636  35.785 10.0512  1e-04 ***
# wood.PC1    1  29.371 29.371 0.11504  45.318  3.7177  1e-04 ***
# wood.PC2    1  18.897 18.897 0.07402  26.977  3.2780  1e-04 ***
# simpson     1  11.636 11.636 0.04558  15.793  2.7950  2e-04 ***
# shannon     1   8.473  8.473 0.03319  11.259  2.5314  0.001 ** 

# These factors still capture a large fraction of variance, 
# even when accounting for size.

# A plot to show the relationship of Wood's forage data and mouthpart morphology
plot(wf$wood.PC2, wf$PC1, col = (wf$class))
abline(lm(wf$PC1~wf$wood.PC2))
# A plot to show the relationship of Simpson's index and mouthpart morphology
plot(wf$simpson, wf$PC1, col = (wf$class))
abline(lm(wf$PC1~wf$simpson))

df <- data.frame(
  species = wf$species,
  class = wf$class,
  simpson = wf$simpson,
  wood.dbs = wf$wood.dbs,
  wood.PC1 = wf$wood.PC1,
  wood.PC2 = wf$wood.PC2,
  PC1 = wf$PC1
)

species.order <- forage$species[order(forage$simpson)]
df$species <- sub("bor","borealis",as.character(df$species))
df$species <- sub("bim","bimaculatus",as.character(df$species))
df$species <- sub("imp","impatiens",as.character(df$species))
df$species <- sub("vag","vagans",as.character(df$species))
df$species <- sub("ferv","fervidus",as.character(df$species))
df$species <- sub("terri","terricola",as.character(df$species))
df$species <- sub("tern","ternarius",as.character(df$species))
df$species <- factor(df$species, levels = species.order)

txt.df <- df %>%
  group_by(species) %>%
  summarise(simpson = median(simpson),
            .groups = "drop") %>%
  mutate(
    y = c(8.35,7.25,6.4,2,7.5,6.8,-6.5),
    hjust = c(0,0.5,0.5,0.5,0.5,0.5,1)
  )

specialization.plot <- df %>%
  ggplot(aes(x = simpson, y = PC1)) +
  theme_bw() +
  theme(legend.position="none") +
  geom_smooth(method=lm, color = "grey40", fill = "gray85") +
  geom_point(aes(color = species), alpha = 0.65, size = 2) +
  scale_color_viridis(discrete = TRUE, end = 0.95) +
  geom_text(data = txt.df, aes(x=simpson, y=y, label = species, color = species)) +
  annotate("text", x=0.79, y=-2, hjust = 0, vjust = 1,
           color="grey40", size = 2.5,
           label="ANOVA with residual randomization
14 linear measurements of worker mouthparts
Y ~ log(individual ITS) + species SDI for forage plants
  ITS:  R^2 = 0.377, p < 10^-4
  SDI:  R^2 = 0.046, p < 2x10^-4") +
  labs(x="Simpson's diversity index for forage plants\n(Data from Wood et al. 2019 Ecology)",
       y="Mouthpart morphology (PC1)")

specialization.plot
ggsave("plots/specialization.plot.pdf", specialization.plot, width = 6.5, height = 5, scale = 1)

# Save everything
save.image("bombus.scaling.rda")
