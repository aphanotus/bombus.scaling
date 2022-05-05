# Bombus wing GMM analysis

# Clear memory
# rm(list=ls())
# gc() 

# Load packages
# devtools::install_github("aphanotus/borealis")
x <- c("borealis","geomorph","tidyverse","ggpubr","phytools",
       "RRPP","multcompView","rgl")
invisible(lapply(x, require, character.only = TRUE))

library(borealis)
library(tidyverse)
library(ggpubr)
library(geomorph)
library(phytools)
library(rgl)
library(RRPP)
library(multcompView)

pairwise.group.comparisons <- function (pairwise.summary.output) {
  x <- pairwise.summary.output$summary.table[,4]
  names(x) <- sub(":","-",rownames(pairwise.summary.output$summary.table))
  return(multcompLetters(x))
}

####################
# Import data
####################
# # Convert raw XY coordinates to TPS format
# ?create.tps
# 
# create.tps(
#   input.filename = 'Bombus.forewings.200711.curated.csv',
#   output.filename = 'Bombus.forewings.200711.curated.tps',
#   id.factors = c('species', 'caste', 'digitizer', 'scale1', 'body', 'scale2', 'its'),
#   include.scale = TRUE,
#   invert.scale = TRUE
# )
# create.tps(
#   input.filename = 'Bombus.hindwings.200711.curated.csv',
#   output.filename = 'Bombus.hindwings.200711.curated.tps',
#   id.factors = c('species', 'caste', 'digitizer', 'scale1', 'body', 'scale2', 'its'),
#   include.scale = TRUE,
#   invert.scale = TRUE
# )

# Create links
fw.links <- matrix(c(1,2, 1,5, 5,4, 4,3, 3,2, 5,6, 6,7, 7,8, 8,9, 9,4, 3,11, 11,12, 11,10, 9,10, 10,14, 14,15, 15,16, 16,18, 18,20, 16,17, 17,8, 12,13, 13,19, 14,13, 18,19, 2,12),
                   ncol = 2, byrow = TRUE)
hw.links <- matrix(c(1,2, 2,3, 3,4, 4,5, 5,6), ncol = 2, byrow = TRUE)

# Import the TPS file data
fw.xy <- read.tps("Bombus.forewings.200711.curated.tps", links = fw.links)
hw.xy <- read.tps("Bombus.hindwings.200711.curated.tps", links = hw.links)

# BODY SIZE AND ITS
# Scale body and ITS measurements from pixel values to metric
fw.xy$metadata$body <- as.numeric(fw.xy$metadata$body) / as.numeric(fw.xy$metadata$scale1)
fw.xy$metadata$its <- as.numeric(fw.xy$metadata$its) / as.numeric(fw.xy$metadata$scale2)
# Remove the scaling factors
fw.xy$metadata <- fw.xy$metadata[,-c(5,7)]
# Create a factor combining species and caste
fw.xy$metadata$sp.caste <- paste(fw.xy$metadata$species,fw.xy$metadata$caste,sep = ".")
# Add species and caste abbreviations to the specimen IDs
dimnames(fw.xy$coords)[[3]] <- paste0(dimnames(fw.xy$coords)[[3]],'.',fw.xy$metadata$sp.caste) 

# Repeat these steps for the hindwing
hw.xy$metadata$body <- as.numeric(hw.xy$metadata$body) / as.numeric(hw.xy$metadata$scale1)
hw.xy$metadata$its <- as.numeric(hw.xy$metadata$its) / as.numeric(hw.xy$metadata$scale2)
hw.xy$metadata <- hw.xy$metadata[,-c(5,7)]
hw.xy$metadata$sp.caste <- paste(hw.xy$metadata$species,hw.xy$metadata$caste,sep = ".")
dimnames(hw.xy$coords)[[3]] <- paste0(dimnames(hw.xy$coords)[[3]],'.',hw.xy$metadata$sp.caste) 

# View shape data
landmark.plot(fw.xy, links = fw.links) #flipped upside down!
landmark.plot(hw.xy, links = hw.links) #also flipped upside down

# Reflect specimens
fw.xy <- align.reflect(fw.xy, top.pt = 1, links = fw.links)
hw.xy <- align.reflect(hw.xy, top.pt = 1, links = hw.links)

# Calculate wing lengths
fw.xy$metadata$fw.length <- apply(fw.xy$coords, 3, function(x) borealis::distance(x[1,],x[19,]) )
hw.xy$metadata$hw.length <- apply(hw.xy$coords, 3, function(x) borealis::distance(x[2,],x[6,]) )

####################
# Scaling
####################

# Examine species-level forewing scaling (wing length vs. body length)
x <- which(is.na(fw.xy$metadata$body) | (fw.xy$metadata$body < 10))
scaling.plot(
  x = log10(fw.xy$metadata$body[-x]), y = log10(fw.xy$metadata$fw.length[-x]),
  group=fw.xy$metadata$species[-x], group.title = "species",
  xlab = "body length (log10 mm)", ylab = "forewing length (log10 mm)",
  include.legend = TRUE,
  isometry.line = TRUE,
  convex.hulls = FALSE
)
#        group   n slope        p sig  ci.lo ci.hi spans.zero
# 1       perp   5 1.400 1.29e-01     -0.741 3.550       TRUE
# 2        vag  97 0.844 2.50e-21 ***  0.708 0.981      FALSE
# 3       rufo   9 0.800 2.06e-04 ***  0.531 1.070      FALSE
# 4        imp  77 0.777 6.30e-35 ***  0.708 0.847      FALSE
# 5      bimac 152 0.741 3.66e-37 ***  0.655 0.826      FALSE
# 6        bor  35 0.720 4.09e-07 ***  0.487 0.952      FALSE
# 7       tern 102 0.688 3.45e-33 ***  0.613 0.764      FALSE
# 8      terri  29 0.632 3.22e-07 ***  0.439 0.825      FALSE
# 9      sande   8 0.609 1.07e-02 *    0.200 1.020      FALSE
# 10      ferv  13 0.452 8.44e-03 **   0.141 0.763      FALSE
# The CIs on all these slopes overlap

# Examine species-level forewing scaling (wing length vs. ITS)
x <- which((is.na(fw.xy$metadata$its) & fw.xy$metadata$species!="rufo-perp"))
scaling.plot(
  x = log10(fw.xy$metadata$its[-x]), y = log10(fw.xy$metadata$fw.length[-x]),
  group=fw.xy$metadata$species[-x], group.title = "species",
  xlab = "intertegular span (log10 mm)", ylab = "forewing length (log10 mm)",
  include.legend = TRUE,
  isometry.line = TRUE,
  convex.hulls = FALSE,
  groups.trendlines = TRUE,
  save.as = "plots/scaling.forewing.v.its.pdf"
)
#        group   n slope        p sig  ci.lo ci.hi spans.zero
# 1       perp   5 0.942 9.21e-02 .   -0.284 2.170       TRUE
# 2        imp  92 0.926 1.52e-43 ***  0.855 0.997      FALSE
# 3        bor  36 0.724 6.81e-09 ***  0.532 0.916      FALSE
# 4      bimac 153 0.708 3.63e-43 ***  0.636 0.779      FALSE
# 5       rufo   9 0.707 1.24e-02 *    0.206 1.210      FALSE
# 6       tern 101 0.676 7.65e-29 ***  0.591 0.761      FALSE
# 7        vag 103 0.648 6.18e-19 ***  0.531 0.765      FALSE
# 8      sande   8 0.581 2.05e-03 **   0.307 0.855      FALSE
# 9       ferv  13 0.533 1.14e-02 *    0.147 0.920      FALSE
# 10     terri  30 0.509 2.05e-05 ***  0.305 0.712      FALSE
# Differences in slope based on overlap of CI
# bimac bor ferv imp perp rufo sande tern terri vag
#     a  ab   ab   b   ab   ab     a    a     a   a

# Examine species-level hindwing scaling (wing length vs. body length)
x <- which(is.na(hw.xy$metadata$body) | (hw.xy$metadata$body < 10))
scaling.plot(
  x = log10(hw.xy$metadata$body[-x]), y = log10(hw.xy$metadata$hw.length[-x]),
  group=hw.xy$metadata$species[-x], group.title = "species",
  xlab = "body length (log10 mm)", ylab = "hindwing length (log10 mm)",
  include.legend = TRUE,
  isometry.line = TRUE,
  convex.hulls = FALSE
)
#        group   n slope        p sig  ci.lo ci.hi spans.zero
# 1       perp   4 1.480 2.48e-01     -2.470 5.420       TRUE
# 2        vag  98 0.948 2.23e-20 ***  0.788 1.110      FALSE
# 3        imp  75 0.866 4.30e-36 ***  0.793 0.939      FALSE
# 4      bimac 147 0.842 1.57e-39 ***  0.751 0.932      FALSE
# 5       rufo   9 0.810 1.82e-03 **   0.417 1.200      FALSE
# 6       tern 101 0.777 2.63e-33 ***  0.693 0.862      FALSE
# 7        bor  34 0.721 9.59e-08 ***  0.507 0.936      FALSE
# 8      terri  29 0.704 1.96e-06 ***  0.465 0.944      FALSE
# 9      sande   8 0.676 9.54e-03 **   0.234 1.120      FALSE
# 10      ferv  13 0.447 4.16e-03 **   0.174 0.721      FALSE
# bimac bor ferv imp perp rufo sande tern terri vag
#     a  ab    b   a   ab   ab    ab    a    ab   b

# Examine species-level hindwing scaling (wing length vs. ITS)
x <- which(is.na(hw.xy$metadata$its))
scaling.plot(
  x = log10(hw.xy$metadata$its[-x]), y = log10(hw.xy$metadata$hw.length[-x]),
  group=hw.xy$metadata$species[-x], group.title = "species",
  xlab = "intertegular span (log10 mm)", ylab = "hindwing length (log10 mm)",
  include.legend = TRUE,
  isometry.line = TRUE,
  convex.hulls = FALSE
)
#        group   n slope        p sig ci.lo ci.hi spans.zero
# 1       perp   4 1.120 4.17e-02 *   0.104 2.140      FALSE
# 2        imp  90 1.020 1.57e-39 *** 0.934 1.110      FALSE
# 3       rufo   9 0.807 5.65e-03 **  0.322 1.290      FALSE
# 4      bimac 149 0.796 3.50e-47 *** 0.723 0.869      FALSE
# 5       tern 100 0.767 2.11e-29 *** 0.673 0.861      FALSE
# 6        bor  35 0.737 3.23e-08 *** 0.528 0.946      FALSE
# 7        vag 104 0.712 1.22e-17 *** 0.576 0.848      FALSE
# 8      sande   8 0.643 1.77e-03 **  0.348 0.937      FALSE
# 9      terri  30 0.577 5.66e-05 *** 0.327 0.826      FALSE
# 10      ferv  13 0.563 2.41e-03 **  0.246 0.879      FALSE
# bimac bor ferv imp perp rufo sande tern terri vag
#     a  ab    a   b   ab   ab     a    a     a   a

####################
# GPA and PCA
####################

# Procrustes alignment - GPA with outlier detection
fw.gpa <- align.procrustes(fw.xy, outlier.analysis = F) # None removed
hw.gpa <- align.procrustes(hw.xy, outlier.analysis = F) # None removed

# report data provenance
write.provenance(fw.gpa, output.filename="Bombus.forewings.provenance.200711.md")
write.provenance(hw.gpa, output.filename="Bombus.hindwings.provenance.200711.md")

# ORDINATION
# PCA - forewings
fw.pca <- gm.prcomp(fw.gpa$gdf$coords)
shape.space(fw.pca, group=fw.gpa$gdf$species,
            group.title = 'species',
            convex.hulls = TRUE, include.legend=TRUE,
            # fixed.aspect = TRUE,
            main.title = "Forewing PCA")

shape.space(fw.pca, group = fw.gpa$gdf$species,
            group.title = 'species', convex.hulls = TRUE,
            backtransform.examples = TRUE,
            # axis1 = 1, axis2 = 3,
            ref.shape = mshape(fw.gpa$gdf$coords),
            shape.method = "TPS",
            bt.shape.mag = 3,
            bt.links = fw.links,
            main.title = "Forewing PCA")
# PC1: 29.8%, PC2: 17.5%, PC3: 8.04%
# - PC1 separates wings that are broad (low values) from those that
#   are narrow (high values)
# - PC2 is roughly a distal anterior skew (low values) vs. proximal 
#   anterior skew. 
# - PC3 is a proximal-to-distal expansion of the anterior/costal edge 
#   (low values) vs the opposite and a anterior-posterior stretch. 
#   PC3 shows the most separation by species.

# PCA - hindwings
hw.pca <- gm.prcomp(hw.gpa$gdf$coords)
shape.space(hw.pca, group=hw.gpa$gdf$species,
            group.title = 'species',
            convex.hulls = TRUE, include.legend=TRUE)

shape.space(hw.pca, group = hw.gpa$gdf$species,
            group.title = 'species', convex.hulls = TRUE,
            backtransform.examples = TRUE,
            # axis1 = 1, axis2 = 3,
            ref.shape = mshape(hw.gpa$gdf$coords),
            shape.method = "TPS",
            bt.shape.mag = 3,
            bt.links = hw.links)
# PC1: 43.45%, PC2: 23.24%, PC3 15.2%
# - PC1 separates wings that are relatively thin, with large distance between 
#   Landmarks 4 &5, (low values) from those where landmarks 4,5,6 are more 
#   proximal and evenly spaced (high values)
# - PC2 separates wings that are narrow and pointed (landmarks 2&3 very close) 
#   (low values) from those that are wide (distance between landmarks 1&4),
#   more rounded (space between landmarks 2&3) and with landmark 4 more distal.  

# MORE PCA
# JUST WORKERS - forewings
fw.workers <- subsetgmm(fw.gpa, fw.gpa$gdf$caste == "W")
fw.w.pca <- gm.prcomp(fw.workers$gdf$coords)
shape.space(fw.w.pca, group = fw.workers$gdf$species,
          group.title = 'species', convex.hulls = TRUE,
          backtransform.examples = TRUE,
          # axis1 = 1, axis2 = 3,
          ref.shape = mshape(fw.workers$gdf$coords),
          shape.method = "TPS",
          bt.shape.mag = 3,
          bt.links = fw.links,
          main.title ="Worker forewing PCA",
          # save.as = "plots/worker.forewing.PCA.pdf", 
          height = 8, width = 10)
# PC1: 26.4%, PC2: 18.2%, PC3: 9.25%, PC4: 7.15%
# - B. borealis and B. fervidus have notably lower values of PC1, PC2 and PC3 
# - PC1 separates wings that are broad (low values) from those that
#   are narrow (high values)
# - PC2 is roughly a distal anterior skew (low values) vs. proximal 
#   anterior skew. 
# - PC3 is a proximal-to-distal expansion of the anterior/costal edge 
#   (low values) vs the opposite and an anterior-posterior stretch.
# - PC4 is an anterior/costal bend to the distal tip (low values) 

# JUST WORKERS - hindwings
hw.workers <- subsetgmm(hw.gpa, hw.gpa$gdf$caste == "W")
hw.w.pca <- gm.prcomp(hw.workers$gdf$coords)
shape.space(hw.w.pca, group = hw.workers$gdf$species,
          group.title = 'species', convex.hulls = TRUE,
          backtransform.examples = TRUE,
          axis1 = 1, axis2 = 5,
          ref.shape = mshape(hw.workers$gdf$coords),
          shape.method = "TPS",
          bt.shape.mag = 3,
          bt.links = hw.links,
          main.title ="Worker hindwing PCA")
# PC1: 42.1%, PC2: 22.7%, PC3: 15.6%, PC4: 8.9%
# - B. borealis and B. fervidus have notably higher values of PC1, and PC2 
# - PC1 separates wings that are relatively thin with large distance between 
#   Landmarks 4 & 5 (low values) from those where landmarks 4,5,6 are more 
#   proximal and evenly spaced (high values)
# - PC2 separates wings that are narrow and pointed (landmarks 2 & 3 are very 
#   close) (low values) from those that are wide (large distance between 
#   landmarks 1 & 4), more rounded (space between landmarks 2 & 3) and with 
#   landmark 4 more distal.  
# - PC3 separates wings where landmark 6 is distal of landmark 5 (low values)
#   from those where it is more proximal

####################
# Phylogenetic PCA
####################
# Find mean shape for each species
# Forewings
fw.coords.by.species <- coords.subset(fw.workers$gdf$coords,
                                      group = fw.workers$gdf$species)
fw.mshape.by.species <- lapply(fw.coords.by.species, mshape)
fw.species.names <- names(fw.mshape.by.species)
fw.mshape.by.species <- array(
  data= unlist(fw.mshape.by.species),
  dim = c(dim(fw.mshape.by.species[[1]])[1],2,length(fw.species.names)),
  dimnames = list(NULL,NULL,fw.species.names)
)
# Hindwings
hw.coords.by.species <- coords.subset(hw.workers$gdf$coords,
                                      group = hw.workers$gdf$species)
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
btree$tip.label <- btree$code.name
btree$tip.label

fw.species.names %in% btree$tip.label
hw.species.names %in% btree$tip.label
# All species in our data are in the tree

fw.btree <- keep.tip(btree, fw.species.names)
hw.btree <- keep.tip(btree, hw.species.names)
plot(fw.btree)
plot(hw.btree)

# Overlay the phylogeny onto the PC-space
{
  pdf("plots/wing.PCA.with.phylogeny.pdf", height = 9, width = 6) 
  par(mfrow=c(2,1))
  fw.pca.w.phylo <- gm.prcomp(fw.mshape.by.species, phy = fw.btree)
  plot(fw.pca.w.phylo, phylo = TRUE, main = "Forewing PCA with phylogeny")
  # PC1: 60.45%, PC2: 23.07%
  hw.pca.w.phylo <- gm.prcomp(hw.mshape.by.species, phy = hw.btree)
  plot(hw.pca.w.phylo, phylo = TRUE, main = "Hindwing PCA with phylogeny")
  #PC1: 71.57%, PC2: 14.95%
  par(mfrow=c(1,1))
  dev.off() 
}

# Phylogenetically-aligned PCA
{
  pdf("plots/wing.phyloPCA.pdf", height = 9, width = 6) 
  par(mfrow=c(2,1))
  fw.paca <- gm.prcomp(fw.mshape.by.species, phy=fw.btree, align.to.phy = TRUE)
  plot(fw.paca, phylo=TRUE, main = "Forewing phylogenetically-aligned PCA")
  hw.paca <- gm.prcomp(hw.mshape.by.species, phy=hw.btree, align.to.phy = TRUE)
  plot(hw.paca, phylo=TRUE, main = "Hindwing phylogenetically-aligned PCA")
  par(mfrow=c(1,1))
  dev.off()
  # The phylogenetically-aligned PCA looks very similar to the regular PCA.
  # The long-faced bees, B. borealis and B. fervidus, have a distinct space
  # and that corresponds to their divergence from the other species.
}

####################
# Modeling
####################
i <- 1e3-1

# Examine sample sizes
barplot(sort(c(with(fw.gpa$gdf, by(species, species, length)))), cex.names = 0.5)
abline(h=c(3,10), col="darkred")
barplot(sort(c(with(fw.gpa$gdf, by(sp.caste, sp.caste, length)))), cex.names = 0.5)
abline(h=c(3,10), col="darkred")
(x <- sort(c(with(fw.gpa$gdf, by(sp.caste, sp.caste, length)))))
# Filter to species/caste combinations with at least 10 individuals
(x <- x[which(x>10)])
fw.sub <- subsetgmm(fw.gpa, specimens = which(fw.gpa$gdf$sp.caste %in% names(x)))
fw.gpa$specimen.number
fw.sub$specimen.number

# Repeat for hindwings
barplot(sort(c(with(hw.gpa$gdf, by(sp.caste, sp.caste, length)))), cex.names = 0.5)
abline(h=c(3,10), col="darkred")
(x <- sort(c(with(hw.gpa$gdf, by(sp.caste, sp.caste, length)))))
(x <- x[which(x>10)])
hw.sub <- subsetgmm(hw.gpa, specimens = which(hw.gpa$gdf$sp.caste %in% names(x)))
hw.gpa$specimen.number
hw.sub$specimen.number

# Allometric model of shape
fw.size.model <- procD.lm(coords ~ log(Csize), data = fw.sub$gdf, iter =i)
anova(fw.size.model)
#             Df      SS       MS     Rsq      F      Z Pr(>F)
# log(Csize)   1 0.05378 0.053783 0.09598 61.257 9.0624  0.001 **
# Residuals  577 0.50659 0.000878 0.90402                        
# Total      578 0.56038    
hw.size.model <- procD.lm(coords ~ log(Csize), data = hw.sub$gdf, iter =i)
anova(hw.size.model)
#             Df      SS       MS     Rsq      F      Z Pr(>F)
# log(Csize)   1 0.06491 0.064910 0.08183 51.422 7.8584  0.001 **
# Residuals  577 0.72834 0.001262 0.91817                        
# Total      578 0.79325  

# Model with size and species
fw.species.model <- procD.lm(coords ~ log(Csize) + species, 
                             data=fw.sub$gdf, iter=i)
anova(fw.species.model)
#             Df      SS       MS     Rsq      F       Z Pr(>F)
# log(Csize)   1 0.05378 0.053783 0.09598 88.315  9.8149  0.001 **
# species      6 0.15886 0.026477 0.28349 43.477 18.6085  0.001 **
# Residuals  571 0.34773 0.000609 0.62053                         
# Total      578 0.56038     

hw.species.model <- procD.lm(coords ~ log(Csize) + species, 
                             data=hw.sub$gdf, iter=i)
anova(hw.species.model)
#             Df      SS       MS    Rsq      F        Z Pr(>F)
# log(Csize)   1 0.06491 0.064910 0.08183 74.231  8.7652  0.001 **
# species      6 0.22904 0.038173 0.28873 43.655 16.0488  0.001 **
# Residuals  571 0.49930 0.000874 0.62944                         
# Total      578 0.79325    

# Model with size, species, caste
fw.sc.model <- procD.lm(coords ~ log(Csize) + species + caste,
                        data=fw.sub$gdf, iter=i)
anova(fw.sc.model)
#             Df      SS       MS     Rsq      F       Z Pr(>F)
# log(Csize)   1 0.05378 0.053783 0.09598 91.890  9.8915  0.001 **
# species      6 0.15886 0.026477 0.28349 45.237 18.7487  0.001 **
# caste        2 0.01470 0.007350 0.02623 12.557  9.0105  0.001 **
# Residuals  569 0.33303 0.000585 0.59430                         
# Total      578 0.56038

hw.sc.model <- procD.lm(coords ~ log(Csize) + species + caste,
                        data=hw.sub$gdf, iter=i)
anova(hw.sc.model)
#             Df      SS       MS     Rsq       F      Z Pr(>F)
# log(Csize)   1 0.06491 0.064910 0.08183 75.3495  8.804  0.001 **
# species      6 0.22904 0.038173 0.28873 44.3126 16.112  0.001 **
# caste        2 0.00913 0.004567 0.01152  5.3019  4.246  0.001 **
# Residuals  569 0.49017 0.000861 0.61792                         
# Total      578 0.79325 

# "Species" as a factor has by far the largest effect on wing shapes (fore and hind),
# based on Z-scores. Size is also a strong influence on shape. Caste has
# the weakest effect, though it is still highly significant. 

####################
# Post hoc pairwise comparisons for species
####################

fw.sp.pw <- pairwise(fit = fw.species.model,
                     fit.null = fw.size.model,
                     groups = fw.sub$gdf$species)
summary(fw.sp.pw)
#                      d   UCL (95%)        Z Pr > d
# bimac:bor   0.03808317 0.007127702 8.578943  0.001
# bimac:imp   0.02364539 0.005099359 7.017688  0.001
# bimac:rufo  0.02128584 0.012982237 3.743490  0.001
# bimac:tern  0.02111468 0.005431960 6.010731  0.001
# bimac:terri 0.01867264 0.008417912 4.926807  0.001
# bimac:vag   0.02002867 0.005088553 7.475574  0.001
# bor:imp     0.04760047 0.007642615 8.203169  0.001
# bor:rufo    0.02914266 0.013877540 4.649655  0.001
# bor:tern    0.03072756 0.007898087 6.538529  0.001
# bor:terri   0.04068906 0.010096281 6.749241  0.001
# bor:vag     0.03673413 0.007752510 6.881735  0.001
# imp:rufo    0.03656305 0.013371497 5.751588  0.001
# imp:tern    0.02637340 0.006071257 7.003233  0.001
# imp:terri   0.01994303 0.008861059 5.125002  0.001
# imp:vag     0.03264545 0.005683793 8.977442  0.001
# rufo:tern   0.02426745 0.012967226 4.166129  0.001
# rufo:terri  0.03133875 0.014485257 4.865289  0.001
# rufo:vag    0.01219410 0.012882087 1.381178  0.086
# tern:terri  0.01756777 0.008720666 4.242021  0.001
# tern:vag    0.02316485 0.005500237 7.000131  0.001
# terri:vag   0.02837955 0.008630841 6.004502  0.001
# Wow, they're all different! -- Except rufo:vag

hw.sp.pw <- pairwise(fit = hw.species.model,
                     fit.null = hw.size.model,
                     groups = hw.sub$gdf$species)
summary(hw.sp.pw)
#                       d   UCL (95%)          Z Pr > d
# bimac:bor   0.038758900 0.009566949  7.12751357  0.001
# bimac:imp   0.034370459 0.006877851  7.20740110  0.001
# bimac:rufo  0.028864115 0.016915817  3.44109930  0.001
# bimac:tern  0.017157111 0.007002143  5.18139293  0.001
# bimac:terri 0.038763759 0.011428300  5.54665757  0.001
# bimac:vag   0.031059604 0.006848305  7.50460673  0.001
# bor:imp     0.045953344 0.010109190  8.30033207  0.001
# bor:rufo    0.049997774 0.019138077  5.12299775  0.001
# bor:tern    0.040128455 0.010734639  6.49424706  0.001
# bor:terri   0.057205397 0.013850478  7.60830071  0.001
# bor:vag     0.052078732 0.010192888  8.00048929  0.001
# imp:rufo    0.017257356 0.017518740  1.63083583  0.057
# imp:tern    0.033575395 0.008181679  8.13384484  0.001
# imp:terri   0.025917449 0.011909326  4.55006715  0.001
# imp:vag     0.023283299 0.007640533  5.29462621  0.001
# rufo:tern   0.025465872 0.016793899  2.88395245  0.005
# rufo:terri  0.020672201 0.019621891  1.77727910  0.037
# rufo:vag    0.009874854 0.017122205 -0.01791209  0.509
# tern:terri  0.041089109 0.012193672  5.69552168  0.001
# tern:vag    0.028413446 0.007319128  7.18478111  0.001
# terri:vag   0.016427260 0.011536757  2.83956056  0.004
# Everything's different! -- Except comparisons involving B. rufocinctus.
# That species is hard to ID, but it's not easily confused with any of
# the common species in this group. Instead, it's likely low sample 
# size that leads to the inability to distinguish B. rufocinctus.

####################
# Test for common allometry among workers
####################

fw.w <- subsetgmm(fw.sub, specimens = which(fw.sub$gdf$caste == "W"))
fw.sub$specimen.number
fw.w$specimen.number
hw.w <- subsetgmm(hw.sub, specimens = which(hw.sub$gdf$caste == "W"))
hw.sub$specimen.number
hw.w$specimen.number

fw.w.pca <- gm.prcomp(fw.w$gdf$coords)
scaling.plot(
  x = log10(fw.w$gdf$Csize), y = fw.w.pca$x[,1],
  group=fw.w$gdf$species, group.title = "species",
  xlab = "log10 worker forewing centroid size", ylab = "PC1 for worker forewing shape",
  include.legend = TRUE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  # hull.alpha = 0.05,
  fixed.aspect = FALSE,
  save.as = "plots/worker.forewing.shape.allometries.pdf",
  height = 7
)
#   group   n   slope        p sig   ci.lo     ci.hi spans.zero
# 1 terri  28  0.0562 1.02e-01     -0.0119  0.12400       TRUE
# 2   vag  98 -0.0381 4.39e-02 *   -0.0751 -0.00107      FALSE
# 3  rufo  11 -0.0762 6.26e-02 .   -0.1570  0.00495       TRUE
# 4 bimac 133 -0.0770 5.54e-06 *** -0.1090 -0.04480      FALSE
# 5   imp  78 -0.0943 3.17e-05 *** -0.1370 -0.05190      FALSE
# 6  tern  97 -0.1100 1.95e-04 *** -0.1670 -0.05390      FALSE
# 7   bor  28 -0.1260 4.08e-05 *** -0.1780 -0.07330      FALSE
# B. terricola appears to be the only species with a flat slope.

fw.w.sp.model <- procD.lm(coords ~ log(Csize) + species, data=fw.w$gdf, iter=i)
fw.w.sp.unique.model <- procD.lm(coords ~ log(Csize) * species,
                                    data=fw.w$gdf, iter=i)
anova(fw.w.sp.unique.model) 
#                     Df      SS        MS     Rsq       F       Z Pr(>F)   
# log(Csize)           1 0.02696 0.0269637 0.06233 48.3671  7.5937  0.001 **
# species              6 0.14162 0.0236030 0.32738 42.3387 14.4612  0.001 **
# log(Csize):species   6 0.00811 0.0013514 0.01875  2.4242  5.0135  0.001 **
# Residuals          459 0.25588 0.0005575 0.59154                          
# Total              472 0.43257
fw.w.sp.unique.model.pw <- pairwise(fit = fw.w.sp.unique.model,
                                    fit.null = fw.w.sp.model,
                                    groups = fw.w$gdf$species,
                                    covariate = log(fw.w$gdf$Csize))
(fw.w.sp.unique.model.pw.results <- summary(fw.w.sp.unique.model.pw))
#                      d  UCL (95%)         Z Pr > d
# bimac:bor   0.07081041 0.05919111 2.3954845  0.007
# bimac:imp   0.04383832 0.04104749 1.9921931  0.025
# bimac:rufo  0.05539796 0.06342110 0.9557489  0.167
# bimac:tern  0.05423949 0.05881664 1.2233731  0.121
# bimac:terri 0.08262553 0.06733798 2.5356581  0.008
# bimac:vag   0.04692101 0.03554261 3.0160059  0.004
# bor:imp     0.07414501 0.05968024 2.7375723  0.003
# bor:rufo    0.08461304 0.07612660 2.1313508  0.015
# bor:tern    0.07319213 0.07276249 1.6931612  0.047
# bor:terri   0.11455664 0.08015882 3.2449286  0.001
# bor:vag     0.06728377 0.05702109 2.3828598  0.007
# imp:rufo    0.06176112 0.06362769 1.4689798  0.078
# imp:tern    0.05050523 0.05956746 0.7912728  0.215
# imp:terri   0.08856363 0.06833724 2.7917096  0.003
# imp:vag     0.05621693 0.03777757 3.8794946  0.001
# rufo:tern   0.06312868 0.07642624 0.6878006  0.253
# rufo:terri  0.10971191 0.08298464 2.9602111  0.003
# rufo:vag    0.05004292 0.06176917 0.6798831  0.248
# tern:terri  0.09992560 0.07951620 2.7082826  0.003
# tern:vag    0.06693147 0.05672955 2.5256596  0.007
# terri:vag   0.08530446 0.06591727 2.7880006  0.004
pairwise.group.comparisons(fw.w.sp.unique.model.pw.results)
# bimac   bor   imp  rufo  tern terri   vag 
#   "a"   "b"   "c" "acd" "ac"   "e"   "d" 

hw.w.pca <- gm.prcomp(hw.w$gdf$coords)
scaling.plot(
  x = log10(hw.w$gdf$Csize), y = hw.w.pca$x[,1],
  group=hw.w$gdf$species, group.title = "species",
  xlab = "log10 worker hindwing centroid size", ylab = "PC1 for worker hindwing shape",
  include.legend = TRUE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  # hull.alpha = 0.05,
  fixed.aspect = FALSE,
  save.as = "plots/worker.hindwing.shape.allometries.pdf",
  height = 7
)
#   group   n   slope        p sig   ci.lo  ci.hi spans.zero
# 1  rufo  11  0.23100 2.88e-02 *    0.0300 0.4320      FALSE
# 2   bor  29  0.19300 2.09e-04 ***  0.1010 0.2860      FALSE
# 3 bimac 134  0.16400 9.04e-06 ***  0.0940 0.2350      FALSE
# 4   vag  99  0.10700 5.55e-04 ***  0.0475 0.1660      FALSE
# 5   imp  79  0.10300 1.56e-03 **   0.0404 0.1650      FALSE
# 6  tern  97  0.00934 8.33e-01     -0.0782 0.0969       TRUE
# 7 terri  28 -0.07230 2.35e-01     -0.1950 0.0500       TRUE

hw.w.sp.model <- procD.lm(coords ~ log(Csize) + species, data=hw.w$gdf, iter=i)
hw.w.sp.unique.model <- procD.lm(coords ~ log(Csize) * species,
                                 data=hw.w$gdf, iter=i)
anova(hw.w.sp.unique.model) 
#                     Df      SS       MS     Rsq       F       Z Pr(>F)   
# log(Csize)           1 0.02612 0.026117 0.04126 31.1716  6.7972  0.001 **
# species              6 0.20416 0.034027 0.32251 40.6123 16.5804  0.001 **
# log(Csize):species   6 0.01483 0.002472 0.02343  2.9506  3.8719  0.001 **
# Residuals          463 0.38792 0.000838 0.61280                          
# Total              476 0.63303                    
hw.w.sp.unique.model.pw <- pairwise(fit = hw.w.sp.unique.model,
                                    fit.null = hw.w.sp.model,
                                    groups = hw.w$gdf$species,
                                    covariate = log(hw.w$gdf$Csize))
(hw.w.sp.unique.model.pw.results <- summary(hw.w.sp.unique.model.pw))
#                      d  UCL (95%)         Z Pr > d
# bimac:bor   0.10005944 0.07670320 2.5114070  0.006
# bimac:imp   0.05268840 0.05519898 1.5172490  0.065
# bimac:rufo  0.07742065 0.08573557 1.3424043  0.090
# bimac:tern  0.09957748 0.07025800 2.9116634  0.004
# bimac:terri 0.12989655 0.08600443 3.0010154  0.002
# bimac:vag   0.04857430 0.04616383 1.7739058  0.034
# bor:imp     0.08851835 0.07631038 2.0587185  0.020
# bor:rufo    0.06479397 0.10206385 0.1993803  0.437
# bor:tern    0.10189322 0.09228511 1.9762921  0.027
# bor:terri   0.12802701 0.10324085 2.3872034  0.010
# bor:vag     0.10523622 0.07380071 2.8007290  0.005
# imp:rufo    0.06865347 0.08696461 0.9404096  0.178
# imp:tern    0.05485756 0.07261387 0.7847684  0.214
# imp:terri   0.09960384 0.08785594 2.1001763  0.016
# imp:vag     0.04209865 0.05021517 1.1187883  0.124
# rufo:tern   0.09839152 0.09525777 1.6956846  0.043
# rufo:terri  0.14760502 0.10990296 2.5274965  0.010
# rufo:vag    0.08187592 0.08157363 1.6264895  0.047
# tern:terri  0.07097688 0.10104592 0.5658313  0.298
# tern:vag    0.07446403 0.06625492 2.0288951  0.020
# terri:vag   0.11835516 0.08499872 2.8474496  0.003
pairwise.group.comparisons(hw.w.sp.unique.model.pw.results)
# bimac   bor   imp  rufo  tern terri   vag 
#   "a"   "b" "acd"  "ab"  "ce"   "e"   "d" 
# Similar pattern to forewings, although that also distinguished
# B. impatiens

####################
# Phylogenetic Generalized Least Squares (PGLS)
####################
# create custom geomorph.data.frame
fw.species.gdf <- geomorph.data.frame(
  coords = fw.mshape.by.species,
  Csize = c(by(fw.workers$gdf$Csize, fw.workers$gdf$species, mean, na.rm=TRUE)),
  body = c(by(fw.workers$gdf$body, fw.workers$gdf$species, mean, na.rm=TRUE)),
  its = c(by(fw.workers$gdf$its, fw.workers$gdf$species, mean, na.rm=TRUE)),
  tree = fw.btree
)
hw.species.gdf <- geomorph.data.frame(
  coords = hw.mshape.by.species,
  Csize = c(by(hw.workers$gdf$Csize, hw.workers$gdf$species, mean, na.rm=TRUE)),
  body = c(by(hw.workers$gdf$body, hw.workers$gdf$species, mean, na.rm=TRUE)),
  its = c(by(hw.workers$gdf$its, hw.workers$gdf$species, mean, na.rm=TRUE)),
  tree = hw.btree
)
i <- 1e4-1
fw.Csize.pgls <- procD.pgls(coords ~ log(Csize), phy = tree, 
                            data = fw.species.gdf, iter = i)
anova(fw.Csize.pgls)
#            Df     SS      MS     Rsq     F      Z Pr(>F)
# log(Csize)  1 2.7686 2.76860 0.38631 5.036 1.5683 0.0535 .
# Residuals   8 4.3981 0.54977 0.61369                      
# Total       9 7.1667 
# Marginal allometric effect, after correcting for relatedness

hw.Csize.pgls <- procD.pgls(coords ~ log(Csize), phy = tree,
                            data = hw.species.gdf, iter = i)
anova(hw.Csize.pgls)
#            Df     SS      MS     Rsq      F      Z Pr(>F)
# log(Csize)  1 2.9900 2.98999 0.41586 5.6953 1.6486 0.0419 *
# Residuals   8 4.1999 0.52499 0.58414                       
# Total       9 7.1899    
# Marginal allometric effect, after correcting for relatedness

fw.body.pgls <- procD.pgls(coords ~ log(body), phy = tree, 
                           data = fw.species.gdf, iter = i)
anova(fw.body.pgls)
#           Df     SS      MS     Rsq      F      Z Pr(>F)
# log(body)  1 2.6950 2.69496 0.37604 4.8213 1.5151 0.0592 .
# Residuals  8 4.4718 0.55897 0.62396                       
# Total      9 7.1667  

hw.body.pgls <- procD.pgls(coords ~ log(body), phy = tree,
                           data = hw.species.gdf, iter = i)
anova(hw.body.pgls)
#           Df     SS      MS     Rsq     F      Z Pr(>F)
# log(body)  1 2.7088 2.70881 0.37675 4.836 1.5251 0.0579 .
# Residuals  8 4.4811 0.56013 0.62325                      
# Total      9 7.1899  

fw.its.pgls <- procD.pgls(coords ~ log(its), phy = tree, 
                          data = fw.species.gdf, iter = i)
anova(fw.its.pgls)
#           Df     SS      MS     Rsq      F       Z Pr(>F)
# log(its)   1 0.7748 0.77478 0.10811 0.9697 0.44396  0.354
# Residuals  8 6.3919 0.79899 0.89189                      
# Total      9 7.1667

hw.its.pgls <- procD.pgls(coords ~ log(its), phy = tree,
                          data = hw.species.gdf, iter = i)
anova(hw.its.pgls)
#           Df     SS      MS     Rsq      F       Z Pr(>F)
# log(its)   1 0.7927 0.79271 0.11025 0.9913 0.45616 0.3494
# Residuals  8 6.3972 0.79965 0.88975                      
# Total      9 7.1899 

####################
# Disparity comparisons
####################
# Disparity differences among workers for species with n>10 
morphol.disparity(coords ~ Csize, groups = ~ species, data = fw.w$gdf)
# Procrustes variances for defined groups
#        bimac          bor          imp         rufo         tern        terri          vag 
# 0.0006465557 0.0019382406 0.0010073911 0.0011280813 0.0005731222 0.0007029508 0.0010121767
# P-Values
#       bimac   bor   imp  rufo  tern terri   vag
# bimac 1.000 0.001 0.001 0.008 0.265 0.585 0.001
# bor   0.001 1.000 0.001 0.001 0.001 0.001 0.001
# imp   0.001 0.001 1.000 0.435 0.001 0.007 0.954
# rufo  0.008 0.001 0.435 1.000 0.002 0.017 0.468
# tern  0.265 0.001 0.001 0.002 1.000 0.213 0.001
# terri 0.585 0.001 0.007 0.017 0.213 1.000 0.003
# vag   0.001 0.001 0.954 0.468 0.001 0.003 1.000

morphol.disparity(coords ~ Csize, groups = ~ species, data = hw.w$gdf)
# Procrustes variances for defined groups
#        bimac          bor          imp         rufo         tern        terri          vag 
# 0.0012653944 0.0023214070 0.0011230135 0.0015552815 0.0008657069 0.0012499736 0.0014713723
# P-Values
#       bimac   bor   imp  rufo  tern terri   vag
# bimac 1.000 0.004 0.398 0.320 0.011 0.942 0.182
# bor   0.004 1.000 0.001 0.056 0.001 0.002 0.005
# imp   0.398 0.001 1.000 0.167 0.144 0.567 0.040
# rufo  0.320 0.056 0.167 1.000 0.048 0.363 0.774
# tern  0.011 0.001 0.144 0.048 1.000 0.102 0.001
# terri 0.942 0.002 0.567 0.363 0.102 1.000 0.329
# vag   0.182 0.005 0.040 0.774 0.001 0.329 1.000

####################################
# Scaling using BODY SIZE and ITS
####################################
# Filter to specimens with body size
x <- which(fw.gpa$gdf$caste == "W" & !is.na(fw.gpa$gdf$body) & (fw.gpa$gdf$body > 10))
fw.worker.bsz <- subsetgmm(fw.gpa, specimens = x)
fw.worker.bsz$specimen.number
x <- which(hw.gpa$gdf$caste == "W" & !is.na(hw.gpa$gdf$body) & (hw.gpa$gdf$body > 10))
hw.worker.bsz <- subsetgmm(hw.gpa, specimens = x)
hw.worker.bsz$specimen.number
# Filter to specimens with ITS
x <- which(fw.gpa$gdf$caste == "W" & !is.na(fw.gpa$gdf$its))
fw.worker.its <- subsetgmm(fw.gpa, specimens = x)
fw.worker.its$specimen.number
x <- which(hw.gpa$gdf$caste == "W" & !is.na(hw.gpa$gdf$its))
hw.worker.its <- subsetgmm(hw.gpa, specimens = x)
hw.worker.its$specimen.number

x <- gm.prcomp(fw.worker.bsz$gdf$coords)
fw.worker.bsz$gdf$PC1 <- x$x[,1]
fw.worker.bsz$gdf$PC2 <- x$x[,2]
x <- gm.prcomp(hw.worker.bsz$gdf$coords)
hw.worker.bsz$gdf$PC1 <- x$x[,1]
hw.worker.bsz$gdf$PC2 <- x$x[,2]
x <- gm.prcomp(fw.worker.its$gdf$coords)
fw.worker.its$gdf$PC1 <- x$x[,1]
fw.worker.its$gdf$PC2 <- x$x[,2]
x <- gm.prcomp(hw.worker.its$gdf$coords)
hw.worker.its$gdf$PC1 <- x$x[,1]
hw.worker.its$gdf$PC2 <- x$x[,2]

scaling.plot(
  x = log10(fw.worker.bsz$gdf$body), y = fw.worker.bsz$gdf$PC1,
  group=fw.worker.bsz$gdf$species, group.title = "species",
  xlab = "log10 body size (mm)", ylab = "PC1 for forewing shape",
  include.legend = TRUE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  fixed.aspect = FALSE,
  save.as = "plots/worker.forewing.PC1.v.body.size.pdf",
  height = 7
)
#    group   n   slope        p sig   ci.lo   ci.hi spans.zero
# 1   ferv   5  0.2150 1.14e-01     -0.0941  0.5240       TRUE
# 2  terri  24  0.0523 1.83e-01     -0.0266  0.1310       TRUE
# 3   perp   2 -0.0214 1.00e+00         NaN     NaN      FALSE
# 4    vag  76 -0.0376 1.50e-01     -0.0892  0.0139       TRUE
# 5   tern  93 -0.0641 1.43e-03 **  -0.1030 -0.0254      FALSE
# 6    imp  61 -0.0752 3.81e-05 *** -0.1090 -0.0414      FALSE
# 7  bimac 118 -0.0889 6.67e-09 *** -0.1170 -0.0607      FALSE
# 8   rufo   6 -0.1090 1.46e-01     -0.2770  0.0589       TRUE
# 9    bor  22 -0.1100 7.14e-03 **  -0.1870 -0.0335      FALSE
# 10 sande   6 -0.1850 6.18e-02 .   -0.3850  0.0147       TRUE
# See below for statistic comparisons of slopes via rrpp::pairwise

scaling.plot(
  x = log10(fw.worker.its$gdf$its), y = fw.worker.its$gdf$PC1,
  group=fw.worker.its$gdf$species, group.title = "species",
  xlab = "log10 intertegular span (mm)", ylab = "PC1 for forewing shape",
  include.legend = TRUE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  fixed.aspect = FALSE,
  save.as = "plots/worker.forewing.PC1.v.ITS.pdf",
  height = 7
)
#    group   n   slope        p sig   ci.lo    ci.hi spans.zero
# 1   ferv   5  0.1070 4.70e-01     -0.3060  0.52000       TRUE
# 2  terri  25  0.0256 5.02e-01     -0.0521  0.10300       TRUE
# 3   perp   2 -0.0152 1.00e+00         NaN      NaN      FALSE
# 4    vag  81 -0.0311 8.76e-02 .   -0.0668  0.00468       TRUE
# 5   tern  92 -0.0497 5.10e-03 **  -0.0842 -0.01530      FALSE
# 6   rufo   6 -0.0811 3.53e-01     -0.2950  0.13300       TRUE
# 7    imp  67 -0.0844 2.01e-05 *** -0.1210 -0.04780      FALSE
# 8  bimac 118 -0.0867 2.48e-08 *** -0.1150 -0.05800      FALSE
# 9    bor  22 -0.1120 8.55e-03 **  -0.1910 -0.03170      FALSE
# 10 sande   6 -0.1160 2.27e-01     -0.3410  0.11000       TRUE

scaling.plot(
  x = log10(hw.worker.bsz$gdf$body), y = hw.worker.bsz$gdf$PC1,
  group=hw.worker.bsz$gdf$species, group.title = "species",
  xlab = "log10 body size (mm)", ylab = "PC1 for hindwing shape",
  include.legend = TRUE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  fixed.aspect = FALSE,
  save.as = "plots/worker.hindwing.PC1.v.body.size.pdf",
  height = 7
)
#    group   n   slope        p sig   ci.lo   ci.hi spans.zero
# 1   rufo   6  0.2890 7.44e-02 .   -0.0454 0.6240       TRUE
# 2  bimac 116  0.1670 1.49e-06 ***  0.1020 0.2320      FALSE
# 3    bor  22  0.1540 1.41e-02 *    0.0345 0.2730      FALSE
# 4  sande   6  0.1530 9.33e-02 .   -0.0405 0.3460       TRUE
# 5    imp  60  0.1080 2.43e-03 **   0.0397 0.1760      FALSE
# 6   tern  92  0.0391 2.47e-01     -0.0275 0.1060       TRUE
# 7    vag  76  0.0368 4.29e-01     -0.0555 0.1290       TRUE
# 8  terri  24 -0.0526 4.27e-01     -0.1870 0.0822       TRUE
# 9   perp   2 -0.2300 1.00e+00         NaN    NaN      FALSE
# 10  ferv   5 -0.3350 2.59e-01     -1.1000 0.4340       TRUE

scaling.plot(
  x = log10(hw.worker.its$gdf$its), y = hw.worker.its$gdf$PC1,
  group=hw.worker.its$gdf$species, group.title = "species",
  xlab = "log10 intertegular span (mm)", ylab = "PC1 for hindwing shape",
  include.legend = TRUE,
  groups.trendlines = TRUE,
  convex.hulls = FALSE,
  fixed.aspect = FALSE,
  save.as = "plots/worker.hindwing.PC1.v.ITS.pdf",
  height = 7
)
#    group   n   slope        p sig   ci.lo    ci.hi spans.zero
# 1   rufo   6  0.3000 1.10e-01     -0.106000 0.7060       TRUE
# 2  bimac 116  0.1660 2.27e-06 ***  0.099800 0.2320      FALSE
# 3    bor  22  0.1340 4.85e-02 *    0.000994 0.2670      FALSE
# 4    imp  66  0.0786 2.49e-02 *    0.010200 0.1470      FALSE
# 5  sande   6  0.0509 5.80e-01     -0.184000 0.2860       TRUE
# 6    vag  81  0.0307 3.27e-01     -0.031200 0.0925       TRUE
# 7   tern  91  0.0158 5.98e-01     -0.043600 0.0753       TRUE
# 8  terri  25 -0.0849 1.93e-01     -0.216000 0.0461       TRUE
# 9   perp   2 -0.1470 1.00e+00           NaN    NaN      FALSE
# 10  ferv   5 -0.2840 2.99e-01     -1.010000 0.4380       TRUE

####################
# Modeling with worker wing shapes with body size and ITS
####################
i <- 1e4-1 # number of iterations
fw.body.lm <- procD.lm(coords ~ log(body), data = fw.worker.bsz$gdf, iter = i)
anova(fw.body.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(body)   1 0.02262 0.0226250 0.05942 25.964 6.6141  1e-04 ***
# Residuals 411 0.35814 0.0008714 0.94058                         
# Total     412 0.38077
fw.body.sp.lm <- procD.lm(coords ~ log(body) + species, data = fw.worker.bsz$gdf, iter = i)
anova(fw.body.sp.lm) 
#            Df      SS        MS     Rsq      F       Z Pr(>F)    
# log(body)   1 0.02262 0.0226250 0.05942 39.048  7.3187  1e-04 ***
# species     9 0.12521 0.0139127 0.32885 24.011 16.6803  1e-04 ***
# Residuals 402 0.23293 0.0005794 0.61173                          
# Total     412 0.38077 

hw.body.lm <- procD.lm(coords ~ log(body), data = hw.worker.bsz$gdf, iter = i)
anova(hw.body.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(body)   1 0.02279 0.0227908 0.04295 18.266 5.4938  1e-04 ***
# Residuals 407 0.50783 0.0012477 0.95705                         
# Total     408 0.53062
hw.body.sp.lm <- procD.lm(coords ~ log(body) + species, data = hw.worker.bsz$gdf, iter = i)
anova(hw.body.sp.lm) 
#            Df      SS        MS     Rsq      F       Z Pr(>F)    
# log(body)   1 0.02279 0.0227908 0.04295 27.372  6.4034  1e-04 ***
# species     9 0.17644 0.0196046 0.33252 23.545 13.2631  1e-04 ***
# Residuals 398 0.33139 0.0008326 0.62453                          
# Total     408 0.53062  

fw.its.lm <- procD.lm(coords ~ log(its), data = fw.worker.its$gdf, iter = i)
anova(fw.its.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(its)    1 0.02050 0.0204992 0.05217 23.227 6.2291  1e-04 ***
# Residuals 422 0.37244 0.0008826 0.94783                         
# Total     423 0.39294  
fw.its.sp.lm <- procD.lm(coords ~ log(its) + species, data = fw.worker.its$gdf, iter = i)
anova(fw.its.sp.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(its)    1 0.02050 0.0204992 0.05217 35.120  6.894  1e-04 ***
# species     9 0.13138 0.0145974 0.33434 25.009 15.840  1e-04 ***
# Residuals 413 0.24107 0.0005837 0.61349                         
# Total     423 0.39294 

hw.its.lm <- procD.lm(coords ~ log(its), data = hw.worker.its$gdf, iter = i)
anova(hw.its.lm)
#            Df      SS       MS     Rsq      F      Z Pr(>F)    
# log(its)    1 0.01130 0.011296 0.02079 8.8734 4.0245  1e-04 ***
# Residuals 418 0.53213 0.001273 0.97921                         
# Total     419 0.54343   

hw.its.sp.lm <- procD.lm(coords ~ log(its) + species, data = hw.worker.its$gdf, iter = i)
anova(hw.its.sp.lm)
#            Df      SS        MS     Rsq      F       Z Pr(>F)    
# log(its)    1 0.01130 0.0112962 0.02079 13.425  4.8724  1e-04 ***
# species     9 0.18797 0.0208860 0.34591 24.821 15.0819  1e-04 ***
# Residuals 409 0.34416 0.0008415 0.63331                          
# Total     419 0.54343 

# Testing common allometry
fw.bodyXsp.lm <- procD.lm(coords ~ log(body) * species, data = fw.worker.bsz$gdf, iter = i)
fw.itsXsp.lm <- procD.lm(coords ~ log(its) * species, data = fw.worker.its$gdf, iter = i)
hw.bodyXsp.lm <- procD.lm(coords ~ log(body) * species, data = hw.worker.bsz$gdf, iter = i)
hw.itsXsp.lm <- procD.lm(coords ~ log(its) * species, data = hw.worker.its$gdf, iter = i)

anova(fw.body.sp.lm, fw.bodyXsp.lm) 
# ResDf Df     RSS       SS        MS     Rsq      F      Z     P Pr(>F)
#   393  9 0.22291 0.010022 0.0011136 0.02632 1.9633 3.7649 3e-04   
anova(fw.its.sp.lm, fw.itsXsp.lm) 
#   404  9 0.23224 0.0088244 0.00098049 0.022457 1.7056 3.0805 0.0023
anova(hw.body.sp.lm, hw.bodyXsp.lm)
#   389  9 0.31566 0.015725 0.0017472 0.029635 2.1531 2.5417 0.0115
anova(hw.its.sp.lm, hw.itsXsp.lm)
#   400  9 0.32914 0.01502 0.0016689 0.027639 2.0282 2.4082 0.0134 

# Pairwise comparisons species-specific allometric slopes
fw.bodyXsp.pw <- pairwise(fit = fw.bodyXsp.lm, fit.null = fw.body.sp.lm, 
                          groups = fw.worker.bsz$gdf$species, 
                          covariate = fw.worker.bsz$gdf$body)
(fw.bodyXsp.pw.output <- summary(fw.bodyXsp.pw))
pairwise.group.comparisons(fw.bodyXsp.pw.output)
# bimac     bor    ferv     imp    perp    rufo   sande    tern   terri     vag 
#  "ab"     "c" "abcde"    "ac" "abcde"  "abcd"     "b"     "c"     "e"     "d" 

fw.itsXsp.pw <- pairwise(fit = fw.itsXsp.lm, fit.null = fw.its.sp.lm, 
                          groups = fw.worker.its$gdf$species, 
                          covariate = fw.worker.its$gdf$its)
(fw.itsXsp.pw.output <- summary(fw.itsXsp.pw))
pairwise.group.comparisons(fw.bodyXsp.pw.output)
# bimac     bor    ferv     imp    perp    rufo   sande    tern   terri     vag 
#  "ab"     "c" "abcde"    "ac" "abcde"  "abcd"     "b"     "c"     "e"     "d" 

hw.bodyXsp.pw <- pairwise(fit = hw.bodyXsp.lm, fit.null = hw.body.sp.lm, 
                          groups = hw.worker.bsz$gdf$species, 
                          covariate = hw.worker.bsz$gdf$body)
(hw.bodyXsp.pw.output <- summary(hw.bodyXsp.pw))
pairwise.group.comparisons(hw.bodyXsp.pw.output)
# bimac    bor   ferv    imp   perp   rufo  sande   tern  terri    vag 
#   "a"   "bc" "abcd" "abcd" "abcd"    "b"   "ab"   "bc"   "cd"    "d" 

hw.itsXsp.pw <- pairwise(fit = hw.itsXsp.lm, fit.null = hw.its.sp.lm, 
                         groups = hw.worker.its$gdf$species, 
                         covariate = hw.worker.its$gdf$its)
(hw.itsXsp.pw.output <- summary(hw.itsXsp.pw))
pairwise.group.comparisons(hw.bodyXsp.pw.output)
# bimac    bor   ferv    imp   perp   rufo  sande   tern  terri    vag 
#   "a"   "bc" "abcd" "abcd" "abcd"    "b"   "ab"   "bc"   "cd"    "d" 
# Body size and ITS lead to the same pattern

# Allometry-corrected PCA
fw.body.lm.pca <- gm.prcomp(fw.body.lm$residuals)
shape.space(fw.body.lm.pca, group = fw.worker.bsz$gdf$species,
          main.title = "Forewing morphospace with body size correction",
          group.title = 'species', convex.hulls = TRUE,
          include.legend = TRUE)
fw.its.lm.pca <- gm.prcomp(fw.its.lm$residuals)
shape.space(fw.its.lm.pca, group = fw.worker.its$gdf$species,
          main.title = "Forewing morphospace with its correction",
          group.title = 'species', convex.hulls = TRUE,
          include.legend = TRUE)
hw.body.lm.pca <- gm.prcomp(hw.body.lm$residuals)
shape.space(hw.body.lm.pca, group = hw.worker.bsz$gdf$species,
          main.title = "Hindwing morphospace with body size correction",
          group.title = 'species', convex.hulls = TRUE,
          include.legend = TRUE)
hw.its.lm.pca <- gm.prcomp(hw.its.lm$residuals)
shape.space(hw.its.lm.pca, group = hw.worker.its$gdf$species,
          main.title = "Hindwing morphospace with its correction",
          group.title = 'species', convex.hulls = TRUE,
          include.legend = TRUE)
# These plots tell a similar story, compared to those before the size-correction.

fw.worker.bsz$gdf$alloPC1 <- fw.body.lm.pca$x[,1]
fw.worker.its$gdf$alloPC1 <- fw.its.lm.pca$x[,1]
hw.worker.bsz$gdf$alloPC1 <- hw.body.lm.pca$x[,1]
hw.worker.its$gdf$alloPC1 <- hw.its.lm.pca$x[,1]

####################
# Correlations to forage diversity
####################
forage <- read.csv("Wood.et.al.2019.forage.diversity.csv") 

# Make the species abbreviations (code names) match
fw.worker.its$gdf$code.name <- as.character(fw.worker.its$gdf$species)
# Both the forage and wing shape datasets use the same species abbreviations
# e.g. "bimac" for B. bimaculatus

# Filter out the foraging information that covers species not in our dataset
x <- which(forage$code.name %in% unique(fw.worker.its$gdf$code.name))
forage <- forage[x,]
# All species in our dataset are covered by the foraging information
# Filter out species in our dataset not covered by the foraging information
x <- which(!(fw.worker.its$gdf$code.name %in% forage$code.name))
unique(fw.worker.its$gdf$code.name[x])
fw.frgls <- fw.worker.its$gdf # This object's name stands for "forewing forage list"
fw.frgls$specimen_id <- fw.frgls$specimen_id[-x]
fw.frgls$coords <- fw.frgls$coords[,,-x]
fw.frgls$body <- fw.frgls$body[-x]
fw.frgls$its <- fw.frgls$its[-x]
fw.frgls$species <- fw.frgls$species[-x]
fw.frgls$fw.length <- fw.frgls$fw.length[-x]
fw.frgls$PC1 <- fw.frgls$PC1[-x]
fw.frgls$PC2 <- fw.frgls$PC2[-x]
fw.frgls$alloPC1 <- fw.frgls$alloPC1[-x]
fw.frgls$code.name <- fw.frgls$code.name[-x]

# Add the foraging information to the wing GDF
fw.worker.its$gdf$class <- fw.worker.its$gdf$code.name
fw.worker.its$gdf$wood.dbs <- fw.worker.its$gdf$code.name
fw.worker.its$gdf$richness <- fw.worker.its$gdf$code.name
fw.worker.its$gdf$shannon <- fw.worker.its$gdf$code.name
fw.worker.its$gdf$simpson <- fw.worker.its$gdf$code.name
fw.worker.its$gdf$faith.pd <- fw.worker.its$gdf$code.name
fw.worker.its$gdf$wood.PC1 <- fw.worker.its$gdf$code.name
fw.worker.its$gdf$wood.PC2 <- fw.worker.its$gdf$code.name
fw.worker.its$gdf$wood.PC3 <- fw.worker.its$gdf$code.name
df <- forage[,-c(1:2)]
colnames(df)[1] <- "class"
row.names(df) <- forage$code.name

for (i in 1:length(fw.worker.its$gdf$code.name)) {
  fw.worker.its$gdf$class[i] <- df[fw.worker.its$gdf$code.name[i],"class"]
  fw.worker.its$gdf$wood.dbs[i] <- df[fw.worker.its$gdf$code.name[i],"wood.dbs"]
  fw.worker.its$gdf$richness[i] <- df[fw.worker.its$gdf$code.name[i],"richness"]
  fw.worker.its$gdf$shannon[i] <- df[fw.worker.its$gdf$code.name[i],"shannon"]
  fw.worker.its$gdf$simpson[i] <- df[fw.worker.its$gdf$code.name[i],"simpson"]
  fw.worker.its$gdf$faith.pd[i] <- df[fw.worker.its$gdf$code.name[i],"faith.pd"]
  fw.worker.its$gdf$wood.PC1[i] <- df[fw.worker.its$gdf$code.name[i],"wood.PC1"]
  fw.worker.its$gdf$wood.PC2[i] <- df[fw.worker.its$gdf$code.name[i],"wood.PC2"]
  fw.worker.its$gdf$wood.PC3[i] <- df[fw.worker.its$gdf$code.name[i],"wood.PC3"]
}
fw.worker.its$gdf$class <-  as.factor(fw.worker.its$gdf$class)
fw.worker.its$gdf$wood.dbs <-  as.numeric(fw.worker.its$gdf$wood.dbs)
fw.worker.its$gdf$richness <-  as.numeric(fw.worker.its$gdf$richness)
fw.worker.its$gdf$shannon <-  as.numeric(fw.worker.its$gdf$shannon)
fw.worker.its$gdf$simpson <-  as.numeric(fw.worker.its$gdf$simpson)
fw.worker.its$gdf$faith.pd <-  as.numeric(fw.worker.its$gdf$faith.pd)
fw.worker.its$gdf$wood.PC1 <-  as.numeric(fw.worker.its$gdf$wood.PC1)
fw.worker.its$gdf$wood.PC2 <-  as.numeric(fw.worker.its$gdf$wood.PC2)
fw.worker.its$gdf$wood.PC3 <-  as.numeric(fw.worker.its$gdf$wood.PC3)

# Explore potential correlations
x <- data.frame(
  PC1 = unlist(fw.worker.its$gdf$PC1),
  PC2 = unlist(fw.worker.its$gdf$PC2),
  alloPC1 = unlist(fw.worker.its$gdf$alloPC1),
  wood.dbs = unlist(fw.worker.its$gdf$wood.dbs),
  richness = unlist(fw.worker.its$gdf$richness),
  shannon = unlist(fw.worker.its$gdf$shannon),
  simpson = unlist(fw.worker.its$gdf$simpson),
  faith.pd = unlist(fw.worker.its$gdf$faith.pd),
  wood.PC1 = unlist(fw.worker.its$gdf$wood.PC1),
  wood.PC2 = unlist(fw.worker.its$gdf$wood.PC2),
  wood.PC3 = unlist(fw.worker.its$gdf$wood.PC3)
)

pairs(x, cor.method = "spearman")
# There are strong correlations with all diversity metrics!
# The strongest correlations appear to be with
# PC1 and Wood's PC1 (0.62), followed by Shannon's Index (0.59), 
# richness (0.48) and Simpson's Index (0.46)

# Model the effect of foraging metrics on wing shape
# Comparisons among these models can be made based on Z values (effect size)
i <- 1e4-1
fw.forage.its.lm <- procD.lm(coords ~ log(Csize), data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.forage.its.lm)
#             Df      SS        MS     Rsq      F      Z Pr(>F)   
# log(Csize)   1 0.02489 0.0248934 0.06335 28.542 6.5584  1e-04 ***
# Residuals  422 0.36805 0.0008722 0.93665                         
# Total      423 0.39294 
fw.forage.sp.lm <- procD.lm(coords ~ species, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.forage.sp.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# species     9 0.14015 0.0155727 0.35668 25.504 14.774  1e-04 ***
# Residuals 414 0.25279 0.0006106 0.64332                         
# Total     423 0.39294                                           
fw.forage.wood.dbs.lm <- procD.lm(coords ~ wood.dbs, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.forage.wood.dbs.lm)
#            Df      SS       MS     Rsq     F      Z Pr(>F)    
# wood.dbs    1 0.04451 0.044513 0.11813 54.92 7.9101  1e-04 ***
# Residuals 410 0.33231 0.000811 0.88187                        
# Total     411 0.37682                                         
fw.forage.richness.lm <- procD.lm(coords ~ richness, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.forage.richness.lm) 
#            Df      SS       MS     Rsq      F      Z Pr(>F)    
# richness    1 0.04256 0.042561 0.11295 52.204 7.6394  1e-04 ***
# Residuals 410 0.33426 0.000815 0.88705                         
# Total     411 0.37682                                          
fw.forage.shannon.lm <- procD.lm(coords ~ shannon, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.forage.shannon.lm) 
#            Df      SS       MS     Rsq      F      Z Pr(>F)    
# shannon     1 0.04802 0.048017 0.12743 59.874 7.7826  1e-04 ***
# Residuals 410 0.32880 0.000802 0.87257                         
# Total     411 0.37682                                          
fw.forage.simpson.lm <- procD.lm(coords ~ simpson, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.forage.simpson.lm)
#            Df      SS       MS     Rsq      F      Z Pr(>F)    
# simpson     1 0.04342 0.043418 0.11522 53.393 7.8419  1e-04 ***
# Residuals 410 0.33340 0.000813 0.88478                         
# Total     411 0.37682                                          
fw.forage.faith.lm <- procD.lm(coords ~ faith.pd, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.forage.faith.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# faith.pd    1 0.02037 0.0203679 0.05405 23.428 6.4576  1e-04 ***
# Residuals 410 0.35645 0.0008694 0.94595                         
# Total     411 0.37682                                           
fw.forage.wood.PC1.lm <- procD.lm(coords ~ wood.PC1, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.forage.wood.PC1.lm)
#            Df      SS       MS     Rsq      F      Z Pr(>F)    
# wood.PC1    1 0.05239 0.052390 0.13903 66.209 7.6791  1e-04 ***
# Residuals 410 0.32443 0.000791 0.86097                         
# Total     411 0.37682                                          
fw.forage.wood.PC2.lm <- procD.lm(coords ~ wood.PC2, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.forage.wood.PC2.lm)
#            Df      SS       MS     Rsq      F     Z Pr(>F)    
# wood.PC2    1 0.02215 0.022150 0.05878 25.605 6.673  1e-04 ***
# Residuals 410 0.35467 0.000865 0.94122                        
# Total     411 0.37682                                         
fw.forage.wood.PC3.lm <- procD.lm(coords ~ wood.PC3, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.forage.wood.PC3.lm)
#            Df      SS       MS     Rsq      F      Z Pr(>F)    
# wood.PC3    1 0.01354 0.013543 0.03594 15.285 5.6004  1e-04 ***
# Residuals 410 0.36328 0.000886 0.96406                         
# Total     411 0.37682                                          

# How do these factors do in models with ITS?
fw.size.size.sp.lm <- procD.lm(coords ~ log(Csize) + species, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.size.size.sp.lm)
#             Df      SS       MS     Rsq      F       Z Pr(>F)    
# log(Csize)   1 0.02489 0.024893 0.06335 43.517  7.2327  1e-04 ***
# species      9 0.13180 0.014644 0.33541 25.600 15.8743  1e-04 ***
# Residuals  413 0.23625 0.000572 0.60124                          
# Total      423 0.39294                                           
fw.size.size.wood.dbs.lm <- procD.lm(coords ~ log(Csize) + wood.dbs, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.size.size.wood.dbs.lm)
#             Df      SS       MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.02305 0.023051 0.06117 30.859 6.8262  1e-04 ***
# wood.dbs     1 0.04826 0.048255 0.12806 64.601 8.3307  1e-04 ***
# Residuals  409 0.30551 0.000747 0.81077                         
# Total      411 0.37682                                          
fw.size.size.richness.lm <- procD.lm(coords ~ log(Csize) + richness, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.size.size.richness.lm) 
#             Df      SS       MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.02305 0.023051 0.06117 30.401 6.8005  1e-04 ***
# richness     1 0.04365 0.043653 0.11585 57.572 7.5984  1e-04 ***
# Residuals  409 0.31012 0.000758 0.82298                         
# Total      411 0.37682                                          
fw.size.size.shannon.lm <- procD.lm(coords ~ log(Csize) + shannon, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.size.size.shannon.lm) 
#             Df      SS       MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.02305 0.023051 0.06117 30.622 6.8131  1e-04 ***
# shannon      1 0.04589 0.045893 0.12179 60.967 7.5886  1e-04 ***
# Residuals  409 0.30788 0.000753 0.81704                         
# Total      411 0.37682
fw.size.size.simpson.lm <- procD.lm(coords ~ log(Csize) + simpson, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.size.size.simpson.lm)
#             Df      SS       MS     Rsq      F     Z Pr(>F)    
# log(Csize)   1 0.02305 0.023051 0.06117 29.931 6.775  1e-04 ***
# simpson      1 0.03878 0.038782 0.10292 50.357 7.364  1e-04 ***
# Residuals  409 0.31499 0.000770 0.83591                        
# Total      411 0.37682 
fw.size.size.faith.lm <- procD.lm(coords ~ log(Csize) + faith.pd, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.size.size.faith.lm)
#             Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.02305 0.0230511 0.06117 28.312 6.6795  1e-04 ***
# faith.pd     1 0.02077 0.0207674 0.05511 25.507 6.5997  1e-04 ***
# Residuals  409 0.33300 0.0008142 0.88372                         
# Total      411 0.37682
fw.size.size.wood.PC1.lm <- procD.lm(coords ~ log(Csize) + wood.PC1, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.size.size.wood.PC1.lm)
#             Df      SS       MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.02305 0.023051 0.06117 31.139 6.8413  1e-04 ***
# wood.PC1     1 0.05100 0.051000 0.13534 68.894 8.0539  1e-04 ***
# Residuals  409 0.30277 0.000740 0.80348                         
# Total      411 0.37682 
fw.size.size.wood.PC2.lm <- procD.lm(coords ~ log(Csize) + wood.PC2, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.size.size.wood.PC2.lm)
#             Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.02305 0.0230511 0.06117 28.719 6.7054  1e-04 ***
# wood.PC2     1 0.02549 0.0254889 0.06764 31.756 6.8170  1e-04 ***
# Residuals  409 0.32828 0.0008026 0.87119                         
# Total      411 0.37682 
fw.size.size.wood.PC3.lm <- procD.lm(coords ~ log(Csize) + wood.PC3, data = fw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(fw.size.size.wood.PC3.lm)
#             Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.02305 0.0230511 0.06117 27.708 6.6438  1e-04 ***
# wood.PC3     1 0.01351 0.0135055 0.03584 16.234 5.6432  1e-04 ***
# Residuals  409 0.34026 0.0008319 0.90299                         
# Total      411 0.37682 

# Recap of forage diversity factors (and species) by effect size
#             Df      SS       MS      Rsq      F       Z Pr(>F)    
# species      9 0.13180 0.014644  0.33541 25.600 15.8743  1e-04 ***
# wood.dbs     1 0.04826 0.048255  0.12806 64.601  8.3307  1e-04 ***
# wood.PC1     1 0.05100 0.051000  0.13534 68.894  8.0539  1e-04 ***
# richness     1 0.04365 0.043653  0.11585 57.572  7.5984  1e-04 ***
# shannon      1 0.04589 0.045893  0.12179 60.967  7.5886  1e-04 ***
# simpson      1 0.03878 0.038782  0.10292 50.357  7.364   1e-04 ***
# wood.PC2     1 0.02549 0.0254889 0.06764 31.756  6.8170  1e-04 ***
# faith.pd     1 0.02077 0.0207674 0.05511 25.507  6.5997  1e-04 ***
# wood.PC3     1 0.01351 0.0135055 0.03584 16.234  5.6432  1e-04 ***

# A plot to show the relationship of Wood's DBS for forage data and wing shape
plot(fw.worker.its$gdf$wood.dbs, fw.worker.its$gdf$PC1, col = (fw.worker.its$gdf$class))
abline(lm(fw.worker.its$gdf$PC1 ~ fw.worker.its$gdf$wood.dbs))

df <- data.frame(
  species = fw.worker.its$gdf$species,
  class = fw.worker.its$gdf$class,
  wood.dbs = fw.worker.its$gdf$wood.dbs,
  richness = fw.worker.its$gdf$richness,
  shannon = fw.worker.its$gdf$shannon,
  PC1 = fw.worker.its$gdf$PC1
) %>% 
  filter(!(species %in% c("rufo","san"))) %>% 
  filter(!is.na(wood.dbs))

{
  species.order <- forage$species[order(forage$wood.dbs)]
  df$species <- sub("bimac","bimaculatus",as.character(df$species))
  df$species <- sub("bor","borealis",as.character(df$species))
  df$species <- sub("ferv","fervidus",as.character(df$species))
  df$species <- sub("imp","impatiens",as.character(df$species))
  df$species <- sub("perp","perplexus",as.character(df$species))
  df$species <- sub("tern","ternarius",as.character(df$species))
  df$species <- sub("terri","terricola",as.character(df$species))
  df$species <- sub("vag","vagans",as.character(df$species))
  df$species <- factor(df$species, levels = species.order)
  }

txt.df <- df %>%
  group_by(species) %>%
  summarise(wood.dbs = median(wood.dbs),
            .groups = "drop") %>%
  mutate(
    y =     c(0.006, -0.027, 0.023, 0.018, 0.004, 0.023, 0.022, -0.0075),
    hjust = c(0,      1,     1,     0,     0,     1,     0,      1)
  )

fw.specialization.plot <- df %>%
  ggplot(aes(x = wood.dbs, y = PC1)) +
  theme_bw() +
  theme(legend.position="none") +
  geom_smooth(method=lm, color = "grey40", fill = "gray85") +
  geom_point(aes(color = species), alpha = 0.65, size = 2) +
  scale_color_viridis(discrete = TRUE, end = 0.95) +
  geom_text(data = txt.df, aes(x=wood.dbs, y=y, hjust=hjust,
                               label = species, color = species)) +
  annotate("text", x=31.5, y=-0.04, hjust = 0, vjust = 1,
           color="grey40", size = 2.5,
           label="permANOVA with residual randomization
12 linear measurements of worker legs
Y ~ log(Csize) + Wood's dietary breadth score (DBS)
  Csize: R^2 = 0.061, p < 10^-4
  DBS:   R^2 = 0.128, p < 10^-4") +
  labs(x="dietary breadth score (Wood et al. 2019 Ecology)",
       y="forewing shape (PC1)")

fw.specialization.plot
ggsave("plots/fw.specialization.plot.pdf", fw.specialization.plot, width = 6.5, height = 5, scale = 1)

# Repeat for hindwings
forage <- read.csv("Wood.et.al.2019.forage.diversity.csv")

# Make the species abbreviations (code names) match
hw.worker.its$gdf$code.name <- as.character(hw.worker.its$gdf$species)
# Both the forage and wing shape datasets use the same species abbreviations
# e.g. "bimac" for B. bimaculatus

# Filter out the foraging information that covers species not in our dataset
x <- which(forage$code.name %in% unique(hw.worker.its$gdf$code.name))
forage <- forage[x,]
# All species in our dataset are covered by the foraging information
# Filter out species in our dataset not covered by the foraging information
x <- which(!(hw.worker.its$gdf$code.name %in% forage$code.name))
unique(hw.worker.its$gdf$code.name[x])
hw.frgls <- hw.worker.its$gdf # This object's name stands for "hindwing forage list"
hw.frgls$specimen_id <- hw.frgls$specimen_id[-x]
hw.frgls$coords <- hw.frgls$coords[,,-x]
hw.frgls$body <- hw.frgls$body[-x]
hw.frgls$its <- hw.frgls$its[-x]
hw.frgls$species <- hw.frgls$species[-x]
hw.frgls$hw.length <- hw.frgls$hw.length[-x]
hw.frgls$PC1 <- hw.frgls$PC1[-x]
hw.frgls$PC2 <- hw.frgls$PC2[-x]
hw.frgls$alloPC1 <- hw.frgls$alloPC1[-x]
hw.frgls$code.name <- hw.frgls$code.name[-x]

# Add the foraging information to the wing GDF
hw.worker.its$gdf$class <- hw.worker.its$gdf$code.name
hw.worker.its$gdf$wood.dbs <- hw.worker.its$gdf$code.name
hw.worker.its$gdf$richness <- hw.worker.its$gdf$code.name
hw.worker.its$gdf$shannon <- hw.worker.its$gdf$code.name
hw.worker.its$gdf$simpson <- hw.worker.its$gdf$code.name
hw.worker.its$gdf$faith.pd <- hw.worker.its$gdf$code.name
hw.worker.its$gdf$wood.PC1 <- hw.worker.its$gdf$code.name
hw.worker.its$gdf$wood.PC2 <- hw.worker.its$gdf$code.name
hw.worker.its$gdf$wood.PC3 <- hw.worker.its$gdf$code.name
df <- forage[,-c(1:2)]
colnames(df)[1] <- "class"
row.names(df) <- forage$code.name

for (i in 1:length(hw.worker.its$gdf$code.name)) {
  hw.worker.its$gdf$class[i] <- df[hw.worker.its$gdf$code.name[i],"class"]
  hw.worker.its$gdf$wood.dbs[i] <- df[hw.worker.its$gdf$code.name[i],"wood.dbs"]
  hw.worker.its$gdf$richness[i] <- df[hw.worker.its$gdf$code.name[i],"richness"]
  hw.worker.its$gdf$shannon[i] <- df[hw.worker.its$gdf$code.name[i],"shannon"]
  hw.worker.its$gdf$simpson[i] <- df[hw.worker.its$gdf$code.name[i],"simpson"]
  hw.worker.its$gdf$faith.pd[i] <- df[hw.worker.its$gdf$code.name[i],"faith.pd"]
  hw.worker.its$gdf$wood.PC1[i] <- df[hw.worker.its$gdf$code.name[i],"wood.PC1"]
  hw.worker.its$gdf$wood.PC2[i] <- df[hw.worker.its$gdf$code.name[i],"wood.PC2"]
  hw.worker.its$gdf$wood.PC3[i] <- df[hw.worker.its$gdf$code.name[i],"wood.PC3"]
}
hw.worker.its$gdf$class <-  as.factor(hw.worker.its$gdf$class)
hw.worker.its$gdf$wood.dbs <-  as.numeric(hw.worker.its$gdf$wood.dbs)
hw.worker.its$gdf$richness <-  as.numeric(hw.worker.its$gdf$richness)
hw.worker.its$gdf$shannon <-  as.numeric(hw.worker.its$gdf$shannon)
hw.worker.its$gdf$simpson <-  as.numeric(hw.worker.its$gdf$simpson)
hw.worker.its$gdf$faith.pd <-  as.numeric(hw.worker.its$gdf$faith.pd)
hw.worker.its$gdf$wood.PC1 <-  as.numeric(hw.worker.its$gdf$wood.PC1)
hw.worker.its$gdf$wood.PC2 <-  as.numeric(hw.worker.its$gdf$wood.PC2)
hw.worker.its$gdf$wood.PC3 <-  as.numeric(hw.worker.its$gdf$wood.PC3)

# Explore potential correlations
x <- data.frame(
  PC1 = unlist(hw.worker.its$gdf$PC1),
  PC2 = unlist(hw.worker.its$gdf$PC2),
  alloPC1 = unlist(hw.worker.its$gdf$alloPC1),
  wood.dbs = unlist(hw.worker.its$gdf$wood.dbs),
  richness = unlist(hw.worker.its$gdf$richness),
  shannon = unlist(hw.worker.its$gdf$shannon),
  simpson = unlist(hw.worker.its$gdf$simpson),
  faith.pd = unlist(hw.worker.its$gdf$faith.pd),
  wood.PC1 = unlist(hw.worker.its$gdf$wood.PC1),
  wood.PC2 = unlist(hw.worker.its$gdf$wood.PC2),
  wood.PC3 = unlist(hw.worker.its$gdf$wood.PC3)
)

pairs(x, cor.method = "spearman")
# Correlations are generally weak to the diversity metrics.
# The strongest correlations appear to be with PC1 and 
# Wood's PC3 (0.39), followed by Faith's PD (0.20) 

# Model the effect of foraging metrics on wing shape
# Comparisons among these models can be made based on Z values (effect size)
i <- 1e4-1
hw.forage.its.lm <- procD.lm(coords ~ log(Csize), data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.forage.its.lm)
#             Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.01988 0.0198831 0.03659 15.875 5.2571  1e-04 ***
# Residuals  418 0.52354 0.0012525 0.96341                         
# Total      419 0.54343 
hw.forage.sp.lm <- procD.lm(coords ~ species, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.forage.sp.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# species     9 0.18926 0.0210292 0.34828 24.345 15.469  1e-04 ***
# Residuals 410 0.35416 0.0008638 0.65172                         
# Total     419 0.54343  
hw.forage.wood.dbs.lm <- procD.lm(coords ~ wood.dbs, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.forage.wood.dbs.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# wood.dbs    1 0.02765 0.0276532 0.05229 22.403 5.8173  1e-04 ***
# Residuals 406 0.50114 0.0012343 0.94771                         
# Total     407 0.52880    
hw.forage.richness.lm <- procD.lm(coords ~ richness, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.forage.richness.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# richness    1 0.01558 0.0155792 0.02946 12.325 4.6187  1e-04 ***
# Residuals 406 0.51322 0.0012641 0.97054                         
# Total     407 0.52880 
hw.forage.shannon.lm <- procD.lm(coords ~ shannon, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.forage.shannon.lm) 
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# shannon     1 0.01615 0.0161524 0.03055 12.792 4.6824  1e-04 ***
# Residuals 406 0.51264 0.0012627 0.96945                         
# Total     407 0.52880      
hw.forage.simpson.lm <- procD.lm(coords ~ simpson, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.forage.simpson.lm)
#            Df     SS        MS     Rsq      F      Z Pr(>F)    
# simpson     1 0.0216 0.0215969 0.04084 17.288 5.2323  1e-04 ***
# Residuals 406 0.5072 0.0012493 0.95916                         
# Total     407 0.5288 
hw.forage.faith.lm <- procD.lm(coords ~ faith.pd, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.forage.faith.lm)
#            Df      SS       MS     Rsq      F      Z Pr(>F)    
# faith.pd    1 0.03366 0.033656 0.06365 27.597 6.2965  1e-04 ***
# Residuals 406 0.49514 0.001220 0.93635                         
# Total     407 0.52880
hw.forage.wood.PC1.lm <- procD.lm(coords ~ wood.PC1, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.forage.wood.PC1.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# wood.PC1    1 0.02431 0.0243094 0.04597 19.564 5.5574  1e-04 ***
# Residuals 406 0.50449 0.0012426 0.95403                         
# Total     407 0.52880 
hw.forage.wood.PC2.lm <- procD.lm(coords ~ wood.PC2, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.forage.wood.PC2.lm)
#            Df      SS        MS     Rsq      F      Z Pr(>F)    
# wood.PC2    1 0.01171 0.0117139 0.02215 9.1975 4.0783  1e-04 ***
# Residuals 406 0.51708 0.0012736 0.97785                         
# Total     407 0.52880 
hw.forage.wood.PC3.lm <- procD.lm(coords ~ wood.PC3, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.forage.wood.PC3.lm)
#            Df      SS       MS     Rsq      F      Z Pr(>F)    
# wood.PC3    1 0.03763 0.037627 0.07116 31.102 6.5224  1e-04 ***
# Residuals 406 0.49117 0.001210 0.92884                         
# Total     407 0.52880 

# How do these factors do in models with ITS?
hw.size.size.sp.lm <- procD.lm(coords ~ log(Csize) + species, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.size.size.sp.lm)
#             Df      SS       MS     Rsq      F       Z Pr(>F)    
# log(Csize)   1 0.01988 0.019883 0.03659 24.486  6.2281  1e-04 ***
# species      9 0.19143 0.021270 0.35226 26.193 15.4562  1e-04 ***
# Residuals  409 0.33212 0.000812 0.61115                          
# Total      419 0.54343  
hw.size.size.wood.dbs.lm <- procD.lm(coords ~ log(Csize) + wood.dbs, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.size.size.wood.dbs.lm)
#             Df      SS       MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.01923 0.019234 0.03637 16.313 5.1110  1e-04 ***
# wood.dbs     1 0.03204 0.032043 0.06060 27.177 6.2588  1e-04 ***
# Residuals  405 0.47752 0.001179 0.90303                         
# Total      407 0.52880 
hw.size.size.richness.lm <- procD.lm(coords ~ log(Csize) + richness, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.size.size.richness.lm) 
#             Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.01923 0.0192338 0.03637 15.738 5.0388  1e-04 ***
# richness     1 0.01461 0.0146146 0.02764 11.959 4.5756  1e-04 ***
# Residuals  405 0.49495 0.0012221 0.93599                         
# Total      407 0.52880 
hw.size.size.shannon.lm <- procD.lm(coords ~ log(Csize) + shannon, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.size.size.shannon.lm) 
#             Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.01923 0.0192338 0.03637 15.787 5.0449  1e-04 ***
# shannon      1 0.01612 0.0161234 0.03049 13.234 4.7794  1e-04 ***
# Residuals  405 0.49344 0.0012184 0.93314                         
# Total      407 0.52880  
hw.size.size.simpson.lm <- procD.lm(coords ~ log(Csize) + simpson, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.size.size.simpson.lm)
#             Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.01923 0.0192338 0.03637 15.973 5.0686  1e-04 ***
# simpson      1 0.02188 0.0218805 0.04138 18.171 5.3874  1e-04 ***
# Residuals  405 0.48768 0.0012042 0.92225                         
# Total      407 0.52880 
hw.size.size.faith.lm <- procD.lm(coords ~ log(Csize) + faith.pd, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.size.size.faith.lm)
#             Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.01923 0.0192338 0.03637 16.112 5.0865  1e-04 ***
# faith.pd     1 0.02609 0.0260938 0.04935 21.859 5.8533  1e-04 ***
# Residuals  405 0.48347 0.0011937 0.91428                         
# Total      407 0.52880 
hw.size.size.wood.PC1.lm <- procD.lm(coords ~ log(Csize) + wood.PC1, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.size.size.wood.PC1.lm)
#             Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.01923 0.0192338 0.03637 15.997 5.0718  1e-04 ***
# wood.PC1     1 0.02262 0.0226187 0.04277 18.812 5.5834  1e-04 ***
# Residuals  405 0.48694 0.0012023 0.92085                         
# Total      407 0.52880    
hw.size.size.wood.PC2.lm <- procD.lm(coords ~ log(Csize) + wood.PC2, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.size.size.wood.PC2.lm)
#             Df      SS        MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.01923 0.0192338 0.03637 15.783 5.0436  1e-04 ***
# wood.PC2     1 0.01601 0.0160115 0.03028 13.139 4.7465  1e-04 ***
# Residuals  405 0.49355 0.0012186 0.93335                         
# Total      407 0.52880  
hw.size.size.wood.PC3.lm <- procD.lm(coords ~ log(Csize) + wood.PC3, data = hw.worker.its$gdf, iter = i, print.progress = FALSE)
anova(hw.size.size.wood.PC3.lm)
#             Df      SS       MS     Rsq      F      Z Pr(>F)    
# log(Csize)   1 0.01923 0.019234 0.03637 16.787 5.1696  1e-04 ***
# wood.PC3     1 0.04554 0.045539 0.08612 39.747 7.0856  1e-04 ***
# Residuals  405 0.46402 0.001146 0.87751                         
# Total      407 0.52880  

# Recap of forage diversity factors (and species) by effect size
#             Df      SS       MS      Rsq      F       Z Pr(>F)    
# species      9 0.19143 0.021270  0.35226 26.193 15.4562  1e-04 ***
# wood.PC3     1 0.04554 0.045539  0.08612 39.747  7.0856  1e-04 ***
# wood.dbs     1 0.03204 0.032043  0.06060 27.177  6.2588  1e-04 ***
# faith.pd     1 0.02609 0.0260938 0.04935 21.859  5.8533  1e-04 ***
# wood.PC1     1 0.02262 0.0226187 0.04277 18.812  5.5834  1e-04 ***
# simpson      1 0.02188 0.0218805 0.04138 18.171  5.3874  1e-04 ***
# shannon      1 0.01612 0.0161234 0.03049 13.234  4.7794  1e-04 ***
# wood.PC2     1 0.01601 0.0160115 0.03028 13.139  4.7465  1e-04 ***
# richness     1 0.01461 0.0146146 0.02764 11.959  4.5756  1e-04 ***

# A plot to show the relationship of Wood's DBS for forage data and wing shape
plot(hw.worker.its$gdf$wood.dbs, hw.worker.its$gdf$PC1, col = (hw.worker.its$gdf$class))
abline(lm(hw.worker.its$gdf$PC1 ~ hw.worker.its$gdf$wood.dbs))

df <- data.frame(
  species = hw.worker.its$gdf$species,
  class = hw.worker.its$gdf$class,
  wood.dbs = hw.worker.its$gdf$wood.dbs,
  richness = hw.worker.its$gdf$richness,
  shannon = hw.worker.its$gdf$shannon,
  PC1 = hw.worker.its$gdf$PC1
) %>% 
  filter(!(species %in% c("rufo","san"))) %>% 
  filter(!is.na(wood.dbs))

{
  species.order <- forage$species[order(forage$wood.dbs)]
  df$species <- sub("bimac","bimaculatus",as.character(df$species))
  df$species <- sub("bor","borealis",as.character(df$species))
  df$species <- sub("ferv","fervidus",as.character(df$species))
  df$species <- sub("imp","impatiens",as.character(df$species))
  df$species <- sub("perp","perplexus",as.character(df$species))
  df$species <- sub("tern","ternarius",as.character(df$species))
  df$species <- sub("terri","terricola",as.character(df$species))
  df$species <- sub("vag","vagans",as.character(df$species))
  df$species <- factor(df$species, levels = species.order)
  }

txt.df <- df %>%
  group_by(species) %>%
  summarise(wood.dbs = median(wood.dbs),
            .groups = "drop") %>%
  mutate(
    y =     c(0.053, 0.03, 0.04, 0.01, -0.01, 0.017, 0.058, 0.02),
    hjust = c(0,     1,    0.5,  0,     0,    1,     0,     1)
  )

hw.specialization.plot <- df %>%
  ggplot(aes(x = wood.dbs, y = PC1)) +
  theme_bw() +
  theme(legend.position="none") +
  geom_smooth(method=lm, color = "grey40", fill = "gray85") +
  geom_point(aes(color = species), alpha = 0.65, size = 2) +
  scale_color_viridis(discrete = TRUE, end = 0.95) +
  geom_text(data = txt.df, aes(x=wood.dbs, y=y, hjust=hjust,
                               label = species, color = species)) +
  annotate("text", x=15, y=-0.042, hjust = 0, vjust = 1,
           color="grey40", size = 2.5,
           label="permANOVA with residual randomization
12 linear measurements of worker legs
Y ~ log(Csize) Wood's dietary breadth score (DBS)
  Csize: R^2 = 0.036, p < 10^-4
  DBS:   R^2 = 0.061, p < 10^-4") +
  labs(x="dietary breadth score (Wood et al. 2019 Ecology)",
       y="hindwing shape (PC1)")

hw.specialization.plot
ggsave("plots/hw.specialization.plot.pdf", hw.specialization.plot, width = 6.5, height = 5, scale = 1)

# Save everything
save.image("bombus.scaling.rda")
