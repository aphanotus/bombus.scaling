# Bombus wing GMM analysis

# Load packages
# devtools::install_github("aphanotus/borealis")
x <- c("borealis","tidyverse","magrittr","ggpubr","ggrepel","factoextra","viridis",
       "multcompView","RRPP","phytools","nlme")
invisible(lapply(x, require, character.only = TRUE))

# Custom functions
slope.comparisons <- function(scaling.plot.output) {
  scaling.plot.output <- arrange(scaling.plot.output, group)
  m <- matrix(nrow = dim(scaling.plot.output)[1], ncol = dim(scaling.plot.output)[1],
              dimnames = list(scaling.plot.output$group,scaling.plot.output$group))
  m[,] <- FALSE
  for (i in 1:(dim(m)[1])) {
    for (j in 1:(dim(m)[2])) {
      if (
        scaling.plot.output$ci.lo[i] > scaling.plot.output$ci.hi[j] |
        scaling.plot.output$ci.hi[i] < scaling.plot.output$ci.lo[j] 
      ) { m[i,j] <- TRUE }
    }
  }
  return(multcompLetters(m))
} # end slope.comparisons

pairwise.group.comparisons <- function (pairwise.summary.output) {
  x <- pairwise.summary.output$summary.table[,4]
  names(x) <- sub(":","-",rownames(pairwise.summary.output$summary.table))
  return(multcompLetters(x))
}

# Create links for landmark positions
links.forewing <- matrix(
  c(1,2, 1,5, 5,4, 4,3, 3,2, 5,6, 6,7, 7,8, 8,9, 9,4, 3,11, 11,12, 11,10, 9,10, 
    10,14, 14,15, 15,16, 16,18, 18,20, 16,17, 17,8, 12,13, 13,19, 14,13, 18,19, 
    2,12),
  ncol = 2, byrow = TRUE
)
links.hindwing <- matrix(c(1,2, 2,3, 3,4, 4,5, 5,6), ncol = 2, byrow = TRUE)

# # Convert raw data into TPS format
# create.tps(
#   input.filename = "Bombus.forewings.250706.csv",
#   output.filename = "Bombus.forewings.tps",
#   id.factors = c("species","caste","digitizer","its.mm","body.length.mm"),
#   include.scale = TRUE,
#   invert.scale =  TRUE
# )
# create.tps(
#   input.filename = "Bombus.hindwings.250706.csv",
#   output.filename = "Bombus.hindwings.tps",
#   id.factors = c("species","caste","digitizer","its.mm","body.length.mm"),
#   include.scale = TRUE,
#   invert.scale =  TRUE 
# )

# Import the TPS shape data
xy.fw <- read.tps("Bombus.forewings.tps", links = links.forewing)
xy.hw <- read.tps("Bombus.hindwings.tps", links = links.hindwing)

# filter to workers only
xy.fw <- subsetgmm(xy.fw, specimens = which(xy.fw$metadata$caste == "W"))
xy.hw <- subsetgmm(xy.hw, specimens = which(xy.hw$metadata$caste == "W"))

# Set linear measurements as numeric
xy.fw$metadata$its.mm <- as.numeric(xy.fw$metadata$its.mm)
xy.fw$metadata$body.length.mm <- as.numeric(xy.fw$metadata$body.length.mm)
xy.hw$metadata$its.mm <- as.numeric(xy.hw$metadata$its.mm)
xy.hw$metadata$body.length.mm <- as.numeric(xy.hw$metadata$body.length.mm)

# correct spelling!
xy.fw$metadata$species <- sub("terrarius","ternarius",xy.fw$metadata$species)
xy.hw$metadata$species <- sub("terrarius","ternarius",xy.hw$metadata$species)

# Remove "cf" in any species names
xy.fw$metadata$species <- sub("^cf ","",xy.fw$metadata$species)
xy.hw$metadata$species <- sub("^cf ","",xy.hw$metadata$species)

# Filter out species with sample sizes (<10)
x <- table(xy.fw$metadata$species)
included.species <- names(x)[which(x>10)]
i <- which(xy.fw$metadata$species %in% included.species)
xy.fw <- subsetgmm(xy.fw, specimens = i)
i <- which(xy.hw$metadata$species %in% included.species)
xy.hw <- subsetgmm(xy.hw, specimens = i)

# Add species names to the specimen IDs
dimnames(xy.fw$coords)[[3]] <- paste0(dimnames(xy.fw$coords)[[3]],'.',xy.fw$metadata$species) 
dimnames(xy.hw$coords)[[3]] <- paste0(dimnames(xy.hw$coords)[[3]],'.',xy.hw$metadata$species) 

# View shape data
landmark.plot(xy.fw, links = links.forewing) #flipped upside down!
landmark.plot(xy.hw, links = links.hindwing) #also flipped upside down

# Calculate wing lengths
xy.fw$metadata$forewing.length.mm <- apply(xy.fw$coords, 3, function(x) borealis::distance(x[1,],x[19,]) )
xy.hw$metadata$hindwing.length.mm <- apply(xy.hw$coords, 3, function(x) borealis::distance(x[2,],x[6,]) )



####################
# Scaling plots
####################

# Examine species-level forewing scaling (wing length vs. body length)
forewing.scaling <- scaling.plot(
  x = log10(xy.fw$metadata$its.mm), y = log10(xy.fw$metadata$forewing.length.mm),
  group = xy.fw$metadata$species,
  xlab = "intertegular span (log10 mm)", ylab = "forewing length (log10 mm)",
  include.legend = TRUE,
  isometry.line = FALSE,
  convex.hulls = FALSE,
  groups.trendlines = TRUE,
  save.as = "plots/scaling.forewing.v.its.pdf"
)
forewing.scaling$plot

forewing.scaling$slopes
#          group   n slope        p sig   ci.lo ci.hi spans.zero
# 1     borealis  41 0.646 3.24e-06 ***  0.4050 0.887      FALSE
# 2 griseocollis  15 0.642 1.06e-03 **   0.3110 0.974      FALSE
# 3  bimaculatus 125 0.612 2.95e-26 ***  0.5230 0.702      FALSE
# 4       vagans  99 0.596 6.24e-17 ***  0.4790 0.712      FALSE
# 5    impatiens  88 0.532 1.92e-12 ***  0.4030 0.660      FALSE
# 6    perplexus  11 0.321 2.03e-01     -0.2080 0.851       TRUE
# 7    terricola  33 0.308 2.91e-02 *    0.0334 0.583      FALSE
# 8    ternarius  97 0.288 8.22e-09 ***  0.1980 0.379      FALSE
# Most species some a positive wing size allometry

forewing.scaling$slopes$group <- factor(forewing.scaling$slopes$group, levels = forewing.scaling$slopes$group)

ylimits <- sort(forewing.scaling$slopes$ci.lo)[1]*1.1
ylimits <- c(ylimits, sort(forewing.scaling$slopes$ci.hi, decreasing = TRUE)[1])

forewing.scaling.slopes.plot <- forewing.scaling$slopes %>% 
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
  scale_y_continuous(name="forewing\nscaling coefficent", limits=c(ylimits[1],1)) +
  scale_x_discrete(name = NULL) +
  coord_flip()
forewing.scaling.slopes.plot


hindwing.scaling <- scaling.plot(
  x = log10(xy.hw$metadata$its.mm), y = log10(xy.hw$metadata$hindwing.length.mm),
  group = xy.hw$metadata$species,
  xlab = "intertegular span (log10 mm)", ylab = "hindwing length (log10 mm)",
  include.legend = TRUE,
  isometry.line = FALSE,
  convex.hulls = FALSE,
  groups.trendlines = TRUE,
  save.as = "plots/scaling.hindwing.v.its.pdf"
)
hindwing.scaling$plot

hindwing.scaling$slopes
#          group   n slope        p sig   ci.lo ci.hi spans.zero
# 1 griseocollis  15 0.674 6.95e-05 ***  0.4200 0.928      FALSE
# 2  bimaculatus 123 0.671 3.16e-24 ***  0.5670 0.775      FALSE
# 3     borealis  41 0.637 2.28e-06 ***  0.4040 0.870      FALSE
# 4    impatiens  87 0.632 9.43e-14 ***  0.4900 0.773      FALSE
# 5       vagans  99 0.623 6.31e-15 ***  0.4890 0.757      FALSE
# 6    perplexus  11 0.368 1.64e-01     -0.1820 0.918       TRUE
# 7    ternarius  96 0.361 1.11e-08 ***  0.2470 0.475      FALSE
# 8    terricola  33 0.315 3.56e-02 *    0.0226 0.608      FALSE
# Most species some a positive wing size allometry

hindwing.scaling$slopes$group <- factor(hindwing.scaling$slopes$group, levels = hindwing.scaling$slopes$group)

ylimits <- sort(hindwing.scaling$slopes$ci.lo)[1]*1.1
ylimits <- c(ylimits, sort(hindwing.scaling$slopes$ci.hi, decreasing = TRUE)[1])

hindwing.scaling.slopes.plot <- hindwing.scaling$slopes %>% 
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
  scale_y_continuous(name="hindwing\nscaling coefficent", limits=c(ylimits[1],1)) +
  scale_x_discrete(name = NULL) +
  coord_flip()
hindwing.scaling.slopes.plot


####################
# GPA and PCA
####################

# Reflect specimens
xy.fw <- align.reflect(xy.fw, top.pt = 2, bottom.pt = 20, left.pt = 13, links = links.forewing)
xy.hw <- align.reflect(xy.hw, top.pt = 1, bottom.pt = 4, left.pt = 6, links = links.hindwing)

# Procrustes alignment
gpa.fw <- align.procrustes(xy.fw, outlier.analysis = F) # None removed
gpa.hw <- align.procrustes(xy.hw, outlier.analysis = F) # None removed

# report data provenance
write.provenance(gpa.fw, output.filename="Bombus.forewings.provenance.md")
write.provenance(gpa.hw, output.filename="Bombus.hindwings.provenance.md")

# ORDINATION
# PCA - forewings
pca.fw <- gm.prcomp(gpa.fw$gdf$coords)

scree.plot(pca.fw)
pcvar(pca.fw, dimensions = 10)

shape.space(
  pca.fw, group = gpa.fw$gdf$species,
  group.title = 'species', convex.hulls = TRUE,
  backtransform.examples = TRUE,
  axis1 = 1, axis2 = 2,
  ref.shape = mshape(gpa.fw$gdf$coords),
  shape.method = "TPS",
  bt.shape.mag = 1,
  bt.links = links.forewing,
  main.title = "forewing shape space"
)
# - PC1 separates wings that are broad in the AP dimension (low values) 
#   from those that are narrower (high values). Low values also have a 
#   constriction of the third submarginal cell and particularly close 
#   landmarks 4 and 5. High values have the opposite, with a third
#   submarginal cell as big or larger than the second submarginal cell.
# - PC2 is roughly a distal anterior shear (low values) vs. proximal 
#   anterior shear. B. griseocollis is distinguished by low PC2 values.
# - PC3 describes the degree of "pinching in" along the proximal-to-distal
#   axis. Low-value shapes have a relative hourglass-shaped deformation grid.
#   While high-value shapes are expanded in medial PD areas. 
# - PC4 describes the relatively "fanning out" of the distal region of the wing.
#   The 3rd submarginal cell is particularly large for low PC4 shapes.

fw.shape.scaling <- scaling.plot(
  x = log10(gpa.fw$gdf$its.mm), y = -pca.fw$x[,1],
  group=gpa.fw$gdf$species, group.title = "species",
  xlab = "intertegular span (log10 mm)", 
  ylab = "PC1 for forewing shape",
  include.legend = TRUE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  # hull.alpha = 0.05,
  fixed.aspect = FALSE,
  save.as = "plots/forewing.shape.allometries.pdf",
  height = 7
)
fw.shape.scaling$plot
fw.shape.scaling$slopes
#          group   n    slope        p sig   ci.lo   ci.hi spans.zero
# 1    perplexus  13  0.09670 1.33e-01     -0.0344  0.2280       TRUE
# 2    terricola  36  0.08260 4.38e-03 **   0.0276  0.1380      FALSE
# 3 griseocollis  16 -0.00155 9.77e-01     -0.1160  0.1130       TRUE
# 4       vagans 119 -0.06440 7.48e-05 *** -0.0955 -0.0334      FALSE
# 5  bimaculatus 141 -0.07130 4.69e-05 *** -0.1050 -0.0378      FALSE
# 6    impatiens  96 -0.08560 2.54e-05 *** -0.1240 -0.0472      FALSE
# 7    ternarius 101 -0.09810 2.61e-04 *** -0.1490 -0.0467      FALSE
# 8     borealis  49 -0.10900 1.04e-03 **  -0.1710 -0.0461      FALSE

# Save the PC coordinates
gpa.fw$gdf$PC1 <- pca.fw$x[,1]
gpa.fw$gdf$PC2 <- pca.fw$x[,2]


# PCA - hindwings
pca.hw <- gm.prcomp(gpa.hw$gdf$coords)

scree.plot(pca.hw)
pcvar(pca.hw, dimensions = 10)

shape.space(
  pca.hw, 
  group = gpa.hw$gdf$species,
  group.title = 'species', convex.hulls = TRUE,
  backtransform.examples = TRUE,
  axis1 = 1, axis2 = 2,
  ref.shape = mshape(gpa.hw$gdf$coords),
  shape.method = "TPS",
  bt.shape.mag = 1,
  bt.links = links.hindwing
)
# - PC1 separates wings that are relatively thin, with large distance between 
#   Landmarks 4 & 5, (low values) from those where landmarks 4, 5, and 6 are  
#   more proximal and evenly spaced (high values)
# - PC2 separates wings that are narrow and pointed (landmarks 2 & 3 very close) 
#   (low values) from those that are wide (farther between landmarks 1 & 4),
#   more rounded (space between landmarks 2 & 3) and with landmark 4 more distal.  
# - PC3 describes the relative positions of landmarks 4 and 5. At low values, 
#   landmarks 4 and 5 are more proximal, such that landmark 5 is just anterior 
#   of landmark 6.
# - PC4 describes the relative positions of the proximal landmarks, 5 and 6. 
#   At low values they are more in-line with landmark 4. At high values, 
#   landmark 5 appears more anterior and landmark 6 is closer to 4.

hw.shape.scaling <- scaling.plot(
  x = log10(gpa.hw$gdf$Csize), y = pca.hw$x[,1],
  group=gpa.hw$gdf$species, group.title = "species",
  xlab = "intertegular span (log10 mm)", 
  ylab = "PC1 for hindwing shape",
  include.legend = TRUE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  # hull.alpha = 0.05,
  fixed.aspect = FALSE,
  save.as = "plots/hindwing.shape.allometries.pdf",
  height = 7
)
hw.shape.scaling$plot
hw.shape.scaling$slopes
#          group   n   slope        p sig   ci.lo  ci.hi spans.zero
# 1     borealis  50  0.1470 0.000126 ***  0.0764 0.2180      FALSE
# 2  bimaculatus 142  0.1400 0.000114 ***  0.0702 0.2100      FALSE
# 3    impatiens  97  0.0633 0.022100 *    0.0093 0.1170      FALSE
# 4       vagans 120  0.0621 0.018800 *    0.0105 0.1140      FALSE
# 5 griseocollis  16 -0.0165 0.753000     -0.1270 0.0936       TRUE
# 6    ternarius 101 -0.0392 0.376000     -0.1270 0.0483       TRUE
# 7    terricola  36 -0.0641 0.291000     -0.1860 0.0573       TRUE
# 8    perplexus  13 -0.1320 0.519000     -0.5690 0.3040       TRUE

# Save the PC coordinates
gpa.hw$gdf$PC1 <- pca.hw$x[,1]
gpa.hw$gdf$PC2 <- pca.hw$x[,2]

ylimits <- sort(fw.shape.scaling$slopes$ci.lo)[1]*1.1
ylimits <- c(ylimits, sort(fw.shape.scaling$slopes$ci.hi, decreasing = TRUE)[1])

scaling.forewing.PC1.slopes.plot <- fw.shape.scaling$slopes %>% 
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
  ylab("forewing shape\nscaling coefficient") +
  coord_flip()
scaling.forewing.PC1.slopes.plot

ylimits <- sort(hw.shape.scaling$slopes$ci.lo)[1]*1.1
ylimits <- c(ylimits, sort(hw.shape.scaling$slopes$ci.hi, decreasing = TRUE)[1])

scaling.hindwing.PC1.slopes.plot <- hw.shape.scaling$slopes %>% 
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
  ylab("hindwing shape\nscaling coefficient") +
  coord_flip()
scaling.hindwing.PC1.slopes.plot

####################
# Phylogenetic PCA
####################

# Find mean shape for each species
# Forewings
fw.coords.by.species <- coords.subset(
  gpa.fw$gdf$coords,
  group = gpa.fw$gdf$species
)
fw.mshape.by.species <- lapply(fw.coords.by.species, mshape)
fw.species.names <- names(fw.mshape.by.species)
fw.mshape.by.species <- array(
  data= unlist(fw.mshape.by.species),
  dim = c(dim(fw.mshape.by.species[[1]])[1],2,length(fw.species.names)),
  dimnames = list(NULL,NULL,fw.species.names)
)

# Hindwings
hw.coords.by.species <- coords.subset(
  gpa.hw$gdf$coords,
  group = gpa.hw$gdf$species
)
hw.mshape.by.species <- lapply(hw.coords.by.species, mshape)
hw.species.names <- names(hw.mshape.by.species)
hw.mshape.by.species <- array(
  data= unlist(hw.mshape.by.species),
  dim = c(dim(hw.mshape.by.species[[1]])[1],2,length(hw.species.names)),
  dimnames = list(NULL,NULL,hw.species.names)
)

# Prepare the tree
# data("Bombus.tree", package ="borealis")
btree <- Bombus.tree
btree$tip.label <- sub("^B\\. ","",btree$tip.label)

# fw.species.names %in% btree$tip.label
# hw.species.names %in% btree$tip.label
# # All species in our data are in the tree

btree <- keep.tip(btree, fw.species.names)
plot(btree)

# Overlay the phylogeny onto the PC-space
{
  pdf("plots/wing.PCA.with.phylogeny.pdf", height = 9, width = 6) 
  par(mfrow=c(2,1))
  pca.fw.w.phylo <- gm.prcomp(fw.mshape.by.species, phy = btree)
  plot(pca.fw.w.phylo, phylo = TRUE, main = "forewing shape space with phylogeny")
  pca.hw.w.phylo <- gm.prcomp(hw.mshape.by.species, phy = btree)
  plot(pca.hw.w.phylo, phylo = TRUE, main = "hindwing shape space with phylogeny")
  par(mfrow=c(1,1))
  dev.off() 
}

# Phylogenetically-aligned PCA
{
  pdf("plots/wing.phyloPCA.pdf", height = 9, width = 6) 
  par(mfrow=c(2,1))
  fw.paca <- gm.prcomp(fw.mshape.by.species, phy=btree, align.to.phy = TRUE)
  plot(fw.paca, phylo=TRUE, main = "Forewing phylogenetically-aligned PCA")
  hw.paca <- gm.prcomp(hw.mshape.by.species, phy=btree, align.to.phy = TRUE)
  plot(hw.paca, phylo=TRUE, main = "Hindwing phylogenetically-aligned PCA")
  par(mfrow=c(1,1))
  dev.off()
}
# For forewings, the phylogenetically-aligned PCA looks very similar to the 
# regular PCA, just rotated. However the phyloPCA makes placement and 
# relationships for the hindwing shapes easier to interpret.


####################
# Modeling
####################

# Forewings

# Allometric model of shape using centroid size
model.fw.size <- procD.lm(coords ~ log(Csize), data = gpa.fw$gdf, iter = 1e4-1)
anova(model.fw.size)
#             Df      SS       MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.04769 0.047693 0.07361 45.215 6.7443  1e-04 ***
# Residuals  569 0.60019 0.001055 0.92639                         
# Total      570 0.64789     

# Allometric model of shape using intertegular span
i <- which(!is.na(gpa.fw$gdf$its.mm))
x <- subsetgmm(gpa.fw$gdf, specimens = i)
model.fw.its <- procD.lm(coords ~ log(its.mm), data = x, iter = 1e4-1)
anova(model.fw.its)
#              Df      SS        MS     Rsq      F     Z Pr(>F)    
# log(its.mm)   1 0.02786 0.0278644 0.04824 25.747 5.906  1e-04 ***
# Residuals   508 0.54978 0.0010822 0.95176                        
# Total       509 0.57765          
# Perhaps not surprisingly, centroid size captures more shape variance.

# Model with size and species
model.fw.size.sp <- procD.lm(
  coords ~ log(Csize) + species, 
  data = gpa.fw$gdf, 
  iter = 1e4-1
)
anova(model.fw.size.sp)
#             Df      SS       MS     Rsq      F       Z Pr(>F)    
# log(Csize)   1 0.04769 0.047693 0.07361 72.994  7.3715  1e-04 ***
# species      7 0.23299 0.033284 0.35962 50.941 14.0585  1e-04 ***
# Residuals  562 0.36720 0.000653 0.56677                          
# Total      570 0.64789    
# Interestingly, species describes much more shape variation than size.

# Post hoc pairwise comparisons for species
pw.fw.species <- pairwise(
  fit = model.fw.size.sp,
  fit.null = model.fw.size,
  groups = gpa.fw$gdf$species,
  print.progress = TRUE
)
pw.fw.species.output <- summary(pw.fw.species)
(sig.letters <- pairwise.group.comparisons(pw.fw.species.output))
# bimaculatus     borealis griseocollis    impatiens    perplexus    ternarius 
#         "a"          "b"          "c"          "d"          "e"          "f" 
# terricola       vagans 
#       "g"          "h" 
# Every species is unique!

df <- data.frame(
  x = 1:8,
  y = rep(max(gpa.fw$gdf$PC1)*1.1,8),
  label = sig.letters$Letters
)

species.order <- data.frame(
  species = gpa.fw$gdf$species,
  PC1 = gpa.fw$gdf$PC1
) %>%
  group_by(species) %>%
  summarise(medianPC1 = median(PC1, na.rm = TRUE)) %>%
  arrange(medianPC1) %>% pull(species)
sp <- factor(gpa.fw$gdf$species, levels = species.order)

forewing.PC1.by.species.plot <- data.frame(
  species = sp,
  PC1 = gpa.fw$gdf$PC1
) %>% 
  ggplot(aes(species, PC1, color = species)) +
  theme_minimal() +
  theme(legend.position="none") +
  geom_violin(scale = "area") +
  geom_quasirandom(varwidth = TRUE) +
  scale_color_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(data = df, aes(x = x, y = y, label = label), color = "black", size = 5 ) +
  xlab("species") +
  ylab("PC1 for forewing shape") +
  labs(caption = "Significant differences in the influence of species in the model  Y ~ intertegular span + species.\nPERMANOVA using `procD.lm`; contrasts using `RRPP::pairwise`.")
forewing.PC1.by.species.plot

ggsave("plots/forewing.PC1.by.species.png", forewing.PC1.by.species.plot, width = 6, height = 4, scale = 1)


# Do species have different forewing shape allometries?
model.unique.fw.allometries <- procD.lm(
  coords ~ log(Csize) * species, 
  data = gpa.fw$gdf, 
  iter = 1e4-1
)
anova(model.unique.fw.allometries)
#                     Df      SS       MS     Rsq       F       Z Pr(>F)    
# log(Csize)           1 0.04769 0.047693 0.07361 74.6553  7.3985  1e-04 ***
# species              7 0.23299 0.033284 0.35962 52.1008 14.0654  1e-04 ***
# log(Csize):species   7 0.01264 0.001806 0.01952  2.8275  5.2892  1e-04 ***
# Residuals          555 0.35456 0.000639 0.54725                           
# Total              570 0.64789         
# Yes, species have different worker forewing shape allometries

# Post hoc examination of the interaction
pw.unique.fw.allometries <- pairwise(
  fit = model.unique.fw.allometries,
  fit.null = model.fw.size,
  groups = gpa.fw$gdf$species,
  covariate = gpa.fw$gdf$Csize,
  print.progress = TRUE
)
pw.unique.fw.allometries.output <- summary(pw.unique.fw.allometries)
(sig.letters <- pairwise.group.comparisons(pw.unique.fw.allometries.output))
# bimaculatus     borealis griseocollis    impatiens    perplexus    ternarius 
#       "abc"          "a"        "bcd"         "ab"       "abcd"       "abcd" 
# terricola       vagans 
#       "d"          "c" 


# Hindwings

# Allometric model of shape using centroid size
model.hw.size <- procD.lm(coords ~ log(Csize), data = gpa.hw$gdf, iter = 1e4-1)
anova(model.hw.size)
#             Df      SS       MS     Rsq     F      Z Pr(>F)    
# log(Csize)   1 0.04531 0.045306 0.05055 30.51 6.6902  1e-04 ***
# Residuals  573 0.85086 0.001485 0.94945                        
# Total      574 0.89617     

# Allometric model of shape using intertegular span
i <- which(!is.na(gpa.hw$gdf$its.mm))
x <- subsetgmm(gpa.hw$gdf, specimens = i)
model.hw.its <- procD.lm(coords ~ log(its.mm), data = x, iter = 1e4-1)
anova(model.hw.its)
#              Df      SS        MS   Rsq      F      Z Pr(>F)    
# log(its.mm)   1 0.01765 0.0176474 0.023 11.843 4.4512  1e-04 ***
# Residuals   503 0.74953 0.0014901 0.977                         
# Total       504 0.76717        
# Perhaps not surprisingly, centroid size captures more shape variance.

# Model with size and species
model.hw.size.sp <- procD.lm(
  coords ~ log(Csize) + species, 
  data = gpa.hw$gdf, 
  iter = 1e4-1
)
anova(model.hw.size.sp)
#             Df      SS       MS     Rsq      F       Z Pr(>F)    
# log(Csize)   1 0.04531 0.045306 0.05055 48.430  7.8157  1e-04 ***
# species      7 0.32138 0.045911 0.35861 49.078 18.1388  1e-04 ***
# Residuals  566 0.52948 0.000935 0.59083                          
# Total      574 0.89617 
# Interestingly, species describes much more shape variation than size.

# Post hoc pairwise comparisons for species
pw.hw.species <- pairwise(
  fit = model.hw.size.sp,
  fit.null = model.hw.size,
  groups = gpa.hw$gdf$species,
  print.progress = TRUE
)
pw.hw.species.output <- summary(pw.hw.species)
(sig.letters <- pairwise.group.comparisons(pw.hw.species.output))
# bimaculatus     borealis griseocollis    impatiens    perplexus    ternarius 
#         "a"          "b"          "c"          "d"          "e"          "f" 
# terricola       vagans 
#       "g"          "h" 
# Every species is unique!

df <- data.frame(
  x = 1:8,
  y = rep(max(gpa.hw$gdf$PC1)*1.1,8),
  label = sig.letters$Letters
)

species.order <- data.frame(
  species = gpa.hw$gdf$species,
  PC1 = gpa.hw$gdf$PC1
) %>%
  group_by(species) %>%
  summarise(medianPC1 = median(PC1, na.rm = TRUE)) %>%
  arrange(medianPC1) %>% pull(species)
sp <- factor(gpa.hw$gdf$species, levels = species.order)

hindwing.PC1.by.species.plot <- data.frame(
  species = sp,
  PC1 = gpa.hw$gdf$PC1
) %>% 
  ggplot(aes(species, PC1, color = species)) +
  theme_minimal() +
  theme(legend.position="none") +
  geom_violin(scale = "area") +
  geom_quasirandom(varwidth = TRUE) +
  scale_color_viridis(discrete = TRUE, option = "magma", begin = 0.15,  end = 0.85) +
  geom_text(data = df, aes(x = x, y = y, label = label), color = "black", size = 5 ) +
  xlab("species") +
  ylab("PC1 for hindwing shape") +
  labs(caption = "Significant differences in the influence of species in the model  Y ~ intertegular span + species.\nPERMANOVA using `procD.lm`; contrasts using `RRPP::pairwise`.")
hindwing.PC1.by.species.plot

ggsave("plots/hindwing.PC1.by.species.png", hindwing.PC1.by.species.plot, width = 6, height = 4, scale = 1)


# Do species have different hindwing shape allometries?
model.unique.hw.allometries <- procD.lm(
  coords ~ log(Csize) * species, 
  data = gpa.hw$gdf, 
  iter = 1e4-1
)
anova(model.unique.hw.allometries)
#                     Df      SS       MS     Rsq       F       Z Pr(>F)    
# log(Csize)           1 0.04531 0.045306 0.05055 49.7575  7.8811  1e-04 ***
# species              7 0.32138 0.045911 0.35861 50.4226 18.2703  1e-04 ***
# log(Csize):species   7 0.02050 0.002928 0.02287  3.2157  4.2601  1e-04 ***
# Residuals          559 0.50899 0.000911 0.56796                           
# Total              574 0.89617         
# Yes, species have different worker hindwing shape allometries

# Post hoc examination of the interaction
pw.unique.hw.allometries <- pairwise(
  fit = model.unique.hw.allometries,
  fit.null = model.hw.size,
  groups = gpa.hw$gdf$species,
  covariate = gpa.hw$gdf$Csize,
  print.progress = TRUE
)
pw.unique.hw.allometries.output <- summary(pw.unique.hw.allometries)
(sig.letters <- pairwise.group.comparisons(pw.unique.hw.allometries.output))
# bimaculatus     borealis griseocollis    impatiens    perplexus    ternarius 
#         "a"          "b"       "abcd"         "cd"       "abcd"        "bcd" 
#   terricola       vagans 
#        "bc"         "ad" 



####################
# Phylogenetic Generalized Least Squares (PGLS)
####################

# create custom geomorph.data.frame
fw.species.gdf <- geomorph.data.frame(
  coords = fw.mshape.by.species,
  Csize = c(by(gpa.fw$gdf$Csize, gpa.fw$gdf$species, mean, na.rm=TRUE)),
  its = c(by(gpa.fw$gdf$its, gpa.fw$gdf$species, mean, na.rm=TRUE)),
  body = c(by(gpa.fw$gdf$body, gpa.fw$gdf$species, mean, na.rm=TRUE)),
  simpson = SDI.by.species,
  tree = btree
)
hw.species.gdf <- geomorph.data.frame(
  coords = hw.mshape.by.species,
  Csize = c(by(gpa.hw$gdf$Csize, gpa.hw$gdf$species, mean, na.rm=TRUE)),
  its = c(by(gpa.hw$gdf$its, gpa.hw$gdf$species, mean, na.rm=TRUE)),
  body = c(by(gpa.hw$gdf$body, gpa.hw$gdf$species, mean, na.rm=TRUE)),
  simpson = SDI.by.species,
  tree = btree
)

fw.its.pgls <- procD.pgls(
  coords ~ log(its), 
  phy = tree, 
  data = fw.species.gdf, 
  iter = 1e4-1
)
anova(fw.its.pgls)
#            Df       SS       MS     Rsq      F      Z Pr(>F)
# log(Csize)  1 0.010477 0.010477 0.13089 0.9036 0.13345 0.4558
# Residuals  6 0.069565 0.011594 0.86911                      
# Total      7 0.080042 

hw.its.pgls <- procD.pgls(
  coords ~ log(its), 
  phy = tree,
  data = hw.species.gdf, 
  iter = 1e4-1
)
anova(hw.its.pgls)
#            Df       SS       MS     Rsq      F        Z Pr(>F)
# log(Csize)  1 0.014575 0.014575 0.13537 0.9393 0.24977  0.418
# Residuals  6 0.093094 0.015516 0.86463                      
# Total      7 0.107668  

# After correcting for relatedness, size is not a significant
# predictor of wing shapes

# With forage plant diversity
fw.its.SDI.pgls <- procD.pgls(
  coords ~ log(its) + simpson, 
  phy = tree, 
  data = fw.species.gdf, 
  iter = 1e4-1
)
anova(fw.its.SDI.pgls)
#           Df       SS       MS     Rsq      F       Z Pr(>F)
# log(its)   1 0.010477 0.010477 0.13089 0.9166 0.14689 0.4504
# simpson    1 0.012416 0.012416 0.15511 1.0862 0.33369 0.3720
# Residuals  5 0.057150 0.011430 0.71399                      
# Total      7 0.080042 

hw.its.SDI.pgls <- procD.pgls(
  coords ~ log(its) + simpson, 
  phy = tree, 
  data = hw.species.gdf, 
  iter = 1e4-1
)
anova(hw.its.SDI.pgls)
#           Df       SS       MS     Rsq      F       Z Pr(>F)
# log(its)   1 0.014575 0.014575 0.13537 1.4783 0.80437 0.2228  
# simpson    1 0.043798 0.043798 0.40679 4.4425 1.84395 0.0318 *
# Residuals  5 0.049295 0.009859 0.45784                        
# Total      7 0.107668      




####################
# Disparity comparisons
####################

fw.disparity <- morphol.disparity(
  coords ~ Csize, 
  groups = ~ species, 
  data = gpa.fw$gdf,
  iter = 1e4-1
)
fw.disparity
# Procrustes variances for defined groups
# bimaculatus     borealis griseocollis    impatiens    perplexus    ternarius    terricola       vagans 
# 0.0007511169 0.0020531930 0.0034115903 0.0010585771 0.0010932015 0.0005919834 0.0007941196 0.0011246682 
# 
# P-Values
#              bimaculatus borealis griseocollis impatiens perplexus ternarius terricola vagans
# bimaculatus       1.0000    1e-04        1e-04    0.0009    0.0960    0.0985    0.7587 0.0001
# borealis          0.0001    1e+00        1e-04    0.0001    0.0006    0.0001    0.0001 0.0001
# griseocollis      0.0001    1e-04        1e+00    0.0001    0.0001    0.0001    0.0001 0.0001
# impatiens         0.0009    1e-04        1e-04    1.0000    0.8725    0.0001    0.0695 0.5230
# perplexus         0.0960    6e-04        1e-04    0.8725    1.0000    0.0250    0.2072 0.8802
# ternarius         0.0985    1e-04        1e-04    0.0001    0.0250    1.0000    0.1643 0.0001
# terricola         0.7587    1e-04        1e-04    0.0695    0.2072    0.1643    1.0000 0.0196
# vagans            0.0001    1e-04        1e-04    0.5230    0.8802    0.0001    0.0196 1.0000

hw.disparity <- morphol.disparity(
  coords ~ Csize, 
  groups = ~ species, 
  data = gpa.hw$gdf,
  iter = 1e4-1
)
hw.disparity
# Procrustes variances for defined groups
# bimaculatus     borealis griseocollis    impatiens    perplexus    ternarius    terricola       vagans 
# 0.001376719  0.002409482  0.002410066  0.001266789  0.002359793  0.000910088  0.001585585  0.001597841 
# 
# P-Values
#              bimaculatus borealis griseocollis impatiens perplexus ternarius terricola vagans
# bimaculatus       1.0000   0.0001       0.0100    0.5329    0.0221    0.0057    0.3798 0.1774
# borealis          0.0001   1.0000       0.9979    0.0001    0.8969    0.0001    0.0061 0.0007
# griseocollis      0.0100   0.9979       1.0000    0.0065    0.9069    0.0009    0.0427 0.0285
# impatiens         0.5329   0.0001       0.0065    1.0000    0.0154    0.0534    0.1988 0.0657
# perplexus         0.0221   0.8969       0.9069    0.0154    1.0000    0.0039    0.0672 0.0479
# ternarius         0.0057   0.0001       0.0009    0.0534    0.0039    1.0000    0.0111 0.0002
# terricola         0.3798   0.0061       0.0427    0.1988    0.0672    0.0111    1.0000 0.9601
# vagans            0.1774   0.0007       0.0285    0.0657    0.0479    0.0002    0.9601 1.0000



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

# Add forage diversity metrics to the forewing shape data object
{
  gpa.fw$gdf$richness <- NA
  gpa.fw$gdf$simpson <- NA
  gpa.fw$gdf$MDS1 <- NA
  gpa.fw$gdf$MDS2 <- NA
}
for (i in 1:length(forage$species)) {
  j <- which(gpa.fw$gdf$species == forage$species[i])
  gpa.fw$gdf$richness[j] <- forage[forage$species[i],"richness"]
  gpa.fw$gdf$simpson[j] <- forage[forage$species[i],"simpson"]
  gpa.fw$gdf$MDS1[j] <- forage[forage$species[i],"MDS1"]
  gpa.fw$gdf$MDS2[j] <- forage[forage$species[i],"MDS2"]
}

# Explore potential correlations
gpa.fw$gdf %>%
  keep(names(.) %in% c("PC1", "PC2", "richness", "simpson", "MDS1", "MDS2")) %>%
  bind_rows() %>% 
  borealis::pairs()
# Wow, really strong correlations between wing shape PCs and forage plant diversity metrics!

plot(gpa.fw$gdf$simpson, gpa.fw$gdf$PC1 )
abline(lm(gpa.fw$gdf$PC1 ~ gpa.fw$gdf$simpson))
plot(gpa.fw$gdf$simpson, gpa.fw$gdf$PC2 )
abline(lm(gpa.fw$gdf$PC2 ~ gpa.fw$gdf$simpson))


####################
# Model the effect of foraging metrics on forewing shape
####################

model.fw.simpson <- procD.lm(
  coords ~ simpson, 
  data = gpa.fw$gdf, 
  iter = 1e4-1, 
  print.progress = TRUE
)
anova(model.fw.simpson)
#            Df      SS       MS     Rsq      F      Z Pr(>F)
# simpson     1 0.06845 0.068451 0.10565 67.218 6.9721  1e-04 ***
# Residuals 569 0.57943 0.001018 0.89435                         
# Total     570 0.64789   

# How do these factors do in models with ITS?

i <- !is.na(gpa.fw$gdf$its.mm)
x <- subsetgmm(gpa.fw$gdf,i)

model.fw.size.simpson <- procD.lm(
  coords ~ log(its.mm) + simpson, 
  data = x, 
  iter = 1e4-1 
)
model.fw.size.simpson.results <- 
  anova(model.fw.size.simpson)$table %T>%
  print() %>% 
  mutate(Rsq = signif(Rsq,3))
#              Df      SS       MS     Rsq      F      Z Pr(>F)    
# log(Csize)    1 0.02770 0.027699 0.04805 28.077 6.1921  1e-04 ***
# simpson       1 0.04952 0.049523 0.08592 50.199 6.8811  1e-04 ***
# Residuals   506 0.49918 0.000987 0.86603                         
# Total       508 0.57641            

# Amazingly, forage plant diversity does almost as well as size in
# predicting forewing shape.

# A plot to show the relationship of Simpson's index and mouthpart morphology
plot(gpa.fw$gdf$simpson, gpa.fw$gdf$PC1, col = as.factor(gpa.fw$gdf$species))
abline(lm(gpa.fw$gdf$PC1~gpa.fw$gdf$simpson))

df <- data.frame(
  species = gpa.fw$gdf$species,
  simpson = gpa.fw$gdf$simpson,
  PC1 = gpa.fw$gdf$PC1,
  PC2 = gpa.fw$gdf$PC2
)

species.order <- forage$species[order(forage$simpson)]
df$species <- factor(df$species, levels = species.order)

txt.df <- df %>%
  group_by(species) %>%
  summarise(simpson = median(simpson),
            .groups = "drop") %>%
  mutate(
    #          gris    bor    vag    terri   perp  tern   bimac  imp
    # medians  0.03   -0.03 -0.01    0.002  0.004  0.002 -0.002  0.022
    y     = c(-0.010, -0.05,  0.002, 0.005, 0.02,  0.027, 0.020, 0.002),
    hjust = c( 0.1,    0.5,   1.1,    1.05,   1.1,  0.5, -0.10,  1.1)
  )

fw.specialization.plot <- df %>%
  ggplot(aes(x = simpson, y = PC1)) +
  theme_bw() +
  theme(legend.position="none") +
  geom_smooth(method=lm, color = "grey40", fill = "gray85") +
  geom_point(aes(color = species), alpha = 0.65, size = 2) +
  scale_color_viridis(discrete = TRUE, end = 0.95) +
  geom_text(data = txt.df, aes(x=simpson, y=y, hjust = hjust, label = species, color = species)) +
  annotate("text", x=0.912, y=-0.060, hjust = 0, vjust = 1,
           color="grey40", size = 2.5,
           label=paste0("ANOVA with residual randomization
10 linear measurements of worker mouthparts
Y ~ log(individual ITS) + species SDI for forage plants
  ITS:  R^2 = ",model.gnathos.size.simpson.results$Rsq[1],", p < ",model.gnathos.size.simpson.results$`Pr(>F)`[1],"
  SDI:  R^2 = ",model.gnathos.size.simpson.results$Rsq[2],", p < ",model.gnathos.size.simpson.results$`Pr(>F)`[2])) +
  labs(x="Simpson's diversity index for forage plants",
       y="Mouthpart morphology (PC1)")

fw.specialization.plot
ggsave("plots/forewing.specialization.plot.png", fw.specialization.plot, width = 5, height = 5, scale = 1)



# Hindwings

# Add forage diversity metrics to the hindwing shape data object
{
  gpa.hw$gdf$richness <- NA
  gpa.hw$gdf$simpson <- NA
  gpa.hw$gdf$MDS1 <- NA
  gpa.hw$gdf$MDS2 <- NA
}
for (i in 1:length(forage$species)) {
  j <- which(gpa.hw$gdf$species == forage$species[i])
  gpa.hw$gdf$richness[j] <- forage[forage$species[i],"richness"]
  gpa.hw$gdf$simpson[j] <- forage[forage$species[i],"simpson"]
  gpa.hw$gdf$MDS1[j] <- forage[forage$species[i],"MDS1"]
  gpa.hw$gdf$MDS2[j] <- forage[forage$species[i],"MDS2"]
}

# Explore potential correlations
gpa.hw$gdf %>%
  keep(names(.) %in% c("PC1", "PC2", "richness", "simpson", "MDS1", "MDS2")) %>%
  bind_rows() %>% 
  borealis::pairs()
# Strong correlations between wing shape PCs and forage plant diversity metrics!

plot(gpa.hw$gdf$richness, gpa.hw$gdf$PC1 )
abline(lm(gpa.hw$gdf$PC1 ~ gpa.hw$gdf$richness))


####################
# Model the effect of foraging metrics on hindwing shape
####################

model.hw.simpson <- procD.lm(
  coords ~ simpson, 
  data = gpa.hw$gdf, 
  iter = 1e4-1, 
  print.progress = TRUE
)
anova(model.hw.simpson)
#            Df      SS       MS     Rsq      F      Z Pr(>F)
# simpson     1 0.11012 0.110123 0.12288 80.276 9.2089  1e-04 ***
# Residuals 573 0.78604 0.001372 0.87712                         
# Total     574 0.89617      

# How do these factors do in models with ITS?

i <- !is.na(gpa.hw$gdf$its.mm)
x <- subsetgmm(gpa.hw$gdf,i)

model.hw.size.simpson <- procD.lm(
  coords ~ log(its.mm) + simpson, 
  data = x, 
  iter = 1e4-1 
)
model.hw.size.simpson.results <- 
  anova(model.hw.size.simpson)$table %T>%
  print() %>% 
  mutate(Rsq = signif(Rsq,3))
#             Df      SS       MS     Rsq      F      Z Pr(>F)    
# log(its.mm)   1 0.01765 0.017647 0.02300 13.401 4.6895  1e-04 ***
# simpson       1 0.08846 0.088462 0.11531 67.176 8.3822  1e-04 ***
# Residuals   502 0.66107 0.001317 0.86169                         
# Total       504 0.76717          

# Amazingly, forage plant diversity does better than size in
# predicting hindwing shape.

# A plot to show the relationship of Simpson's index and mouthpart morphology

df <- data.frame(
  species = gpa.hw$gdf$species,
  simpson = gpa.hw$gdf$simpson,
  PC1 = gpa.hw$gdf$PC1,
  PC2 = gpa.hw$gdf$PC2
)

species.order <- forage$species[order(forage$simpson)]
df$species <- factor(df$species, levels = species.order)

txt.df <- df %>%
  group_by(species) %>%
  summarise(simpson = median(simpson),
            .groups = "drop") %>%
  mutate(
    #          gris    bor    vag   terri   perp  tern   bimac  imp
    y     = c( 0.063,  0.072, -0.01,  0.01, 0.02,  0.03, 0.020, 0.002),
    hjust = c( 0.1,    0.5,   1.1,   0.5,   1.1,  1.1, -0.10,  1.1)
  )

hw.specialization.plot <- df %>%
  ggplot(aes(x = simpson, y = PC1)) +
  theme_bw() +
  theme(legend.position="none") +
  geom_smooth(method=lm, color = "grey40", fill = "gray85") +
  geom_point(aes(color = species), alpha = 0.65, size = 2) +
  scale_color_viridis(discrete = TRUE, end = 0.95) +
  geom_text(data = txt.df, aes(x=simpson, y=y, hjust = hjust, label = species, color = species)) +
  annotate("text", x=0.912, y=-0.060, hjust = 0, vjust = 1,
           color="grey40", size = 2.5,
           label=paste0("ANOVA with residual randomization
10 linear measurements of worker mouthparts
Y ~ log(individual ITS) + species SDI for forage plants
  ITS:  R^2 = ",model.gnathos.size.simpson.results$Rsq[1],", p < ",model.gnathos.size.simpson.results$`Pr(>F)`[1],"
  SDI:  R^2 = ",model.gnathos.size.simpson.results$Rsq[2],", p < ",model.gnathos.size.simpson.results$`Pr(>F)`[2])) +
  labs(x="Simpson's diversity index for forage plants",
       y="Mouthpart morphology (PC1)")

hw.specialization.plot
ggsave("plots/hindwing.specialization.plot.png", hw.specialization.plot, width = 5, height = 5, scale = 1)



# Save everything
save.image("bombus.scaling.rda")
