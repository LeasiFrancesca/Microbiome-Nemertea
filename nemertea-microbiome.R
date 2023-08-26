##### load libraries #####

library(ape)
library(car)
library(ecodist)
library(geosphere)
library(ggplot2)
library(ggsci)
library(marmap)
library(topicmodels)
library(vegan)

##### load data anche check it #####

OTUs <- read.csv("nemertea_ASV.csv", header=T, as.is=F)
names(OTUs[,1:15])
summary(OTUs[,1:10])
head(OTUs[,1:10])
table(OTUs$species)
dim(OTUs)

tree <- ape::read.tree('nemertea_18s.tree')
tree
plot(tree)

# obtain a tree with only the 55 samples for which a microbiome is available
tree55 <- ape::drop.tip(tree, setdiff(tree$tip.label, OTUs$species));
plot(tree55)
tree55$tip.label

##### microbiome distances between samples #####

#try with both Hellinger & Bray-Curtis
OTUs_dist_H <- topicmodels::distHellinger(as.matrix(OTUs[,-1]))
OTUs_dist_BC <- ecodist::bcdist(as.matrix(OTUs[,-1]))

par(mfrow=c(2,1))
plot(hclust(as.dist(OTUs_dist_H), method="average"), labels=OTUs$species, hang=-1)
plot(hclust(as.dist(OTUs_dist_BC), method="average"), labels=OTUs$species, hang=-1)
dev.off()
# to choose Bray-Curtis

##### correlation phylogeny-microbiome #####

DNA_dist <- ape::cophenetic.phylo(tree55)

# check that names of the samples are ordered in the same way

all(rownames(DNA_dist)==OTUs$species)
which(rownames(DNA_dist)!=OTUs$species)

# they are not, thus, reorder and check that orders are fine
OTUs_ord <- OTUs[match(rownames(DNA_dist), OTUs$species),]
all(rownames(DNA_dist)==OTUs_ord$species)
which(rownames(DNA_dist)!=OTUs_ord$species)

# obtain matrix of Bray-Curtis distances, ordered as in the matrix of phylogenetic distances

OTUs_ord_BC <- as.matrix(ecodist::bcdist(as.matrix(OTUs_ord[,-1])))

dim(OTUs_ord_BC)
dim(DNA_dist)

vegan::mantel(OTUs_ord_BC, DNA_dist)

plot(OTUs_ord_BC, DNA_dist)

plot(as.vector(DNA_dist),as.vector(1-OTUs_ord_BC), xlab="phylogenetic distance", ylab="community similarity")
plot(vegan::mantel.correlog((OTUs_ord_BC), DNA_dist, cutoff=F, n.class=7))

# tanglegram

t1 <- ape::chronos(tree55)
t2 <- ape::chronoMPL(tree55)
t3 <- ape::chronopl(tree55,1)

par(mfrow=c(2,2))
plot(tree55)
plot(t1)
plot(t2)
plot(t3)
dev.off()

plot(t1, show.tip.label=F)
axisPhylo()

write.tree(t1,"nemertea_18s_ULTRA.tree")

plot(hclust(as.dist(OTUs_dist_BC), method="average"), labels=OTUs$species, hang=-1)

##### effect of host genus or geography #####

nemert <- read.csv("nemertea_taxonomy 230812.csv", header=T, as.is=F)

nemert_reduced <- droplevels(nemert[nemert$rank!="none",])
OTUs_reduced <- droplevels(OTUs[nemert$rank!="none",])

dim(nemert)
dim(nemert_reduced)
dim(OTUs)
dim(OTUs_reduced)

OTUs_dist_reduced <- ecodist::bcdist(as.matrix(OTUs_reduced[,-1]))

effects <- vegan::adonis2(as.matrix(OTUs_dist_reduced) ~ nemert_reduced$Genus + nemert_reduced$Locality)
effects

##### analyses on only the genus Ototyphlonemertes #####

# read the data with Nemertea information
names(nemert)
summary(nemert)
head(nemert)
table(nemert$taxonomy)
table(nemert$Genus)
dim(nemert)

# check that the names match and are also in the same order
all(nemert$SampleID==OTUs$species)
which(nemert$SampleID!=OTUs$species)

# subset only for the genus Ototyphlonemertes
oto <- na.omit(droplevels(nemert[nemert$Genus=="Ototyphlonemertes",]))
names(oto)
dim(oto)

oto_m <- OTUs[OTUs$species %in% levels(oto$SampleID),]
dim(oto_m)
str(oto_m)
oto_m2 <- oto_m[, colSums(oto_m[-1])>0]
dim(oto_m2)
str(oto_m2)

# plot nMDS

pchs <- c(20,10,8)
col.gr <- c(pal_jco()(7))
gr.species <- factor(oto$taxonomy)
gr.area <- factor(oto$Locality)

oto_nMDS <- vegan::metaMDS(oto_m2[-1],"bray")
plot(oto_nMDS, type = "n", display = "sites")
points(oto_nMDS, display = "sites", pch=pchs[gr.area], col=col.gr[gr.species], cex=2)
legend("topright", legend=levels(gr.area), bty="n", col=c("black"), pch=pchs)
legend("bottomright", legend=levels(gr.species), bty="n", col=col.gr, pch=c(20), cex=0.8, pt.cex=1.5)
for(i in unique(gr.species)) {
  ordihull(oto_nMDS$point[grep(i,gr.species),], draw="polygon",
   groups=gr.species[gr.species==i], label=F) } 

# effect of host species or geography

dim(oto)
dim(oto_m2)

oto_m2_dist <- ecodist::bcdist(as.matrix(oto_m2[,-1]))

effects_oto <- vegan::adonis2(as.matrix(oto_m2_dist) ~ oto$taxonomy + oto$Locality)
effects_oto

##### correlation phylogeny-microbiome #####

DNA_dist <- ape::cophenetic.phylo(tree55)
dim(DNA_dist)

oto_34 <- droplevels(oto_m$species)

DNA_dist_oto1 <- DNA_dist[rownames(DNA_dist) %in% oto_34,]
DNA_dist_oto2 <- t(DNA_dist_oto1)
DNA_dist_oto <- DNA_dist_oto2[rownames(DNA_dist_oto2) %in% oto_34,]

dim(DNA_dist_oto)

rownames(DNA_dist_oto)

# check that names of the samples are ordered in the same way

all(rownames(DNA_dist_oto)==oto_m$species)
which(rownames(DNA_dist_oto)!=oto_m$species)

# they are not, thus, reorder and check that orders are fine
oto_m_ord <- oto_m[match(rownames(DNA_dist_oto), oto_m$species),]
dim(oto_m_ord)
all(rownames(DNA_dist_oto)==oto_m_ord$species)
which(rownames(DNA_dist_oto)!=oto_m_ord$species)


# obtain matrix of Bray-Curtis distances, ordered as in the matrix of phylogenetic distances

oto_m_ord_BC <- as.matrix(ecodist::bcdist(as.matrix(oto_m_ord[,-1])))

dim(oto_m_ord_BC)
dim(DNA_dist_oto)

vegan::mantel(oto_m_ord_BC, DNA_dist_oto)

plot(oto_m_ord_BC, DNA_dist_oto)

plot(as.vector(DNA_dist_oto),as.vector(1-oto_m_ord_BC), xlab="phylogenetic distance", ylab="community similarity")
plot(vegan::mantel.correlog((oto_m_ord_BC), DNA_dist_oto, cutoff=F, n.class=7))

##### correlation phylogeny-microbiome checking for geography ######

# check that names of the samples are ordered in the same way
all(rownames(DNA_dist_oto)==oto$sampleID)
which(rownames(DNA_dist_oto)!=oto$sampleID)
all(oto_m_ord$species==oto$sampleID)
which(oto_m_ord$species!=oto$sampleID)

# map of the sampling area
centrAmer <- marmap::getNOAA.bathy(lon1=-98,lon2=-72, lat1=4,lat2=22.9,resolution=1)  

# check that dephts make sense
depths <- get.depth(centrAmer,oto$longitude, oto$latitude, locator = F)
depths$depth > -1 ## all good
rm(depths)

# Euclidean distances between beaches  
d.euclidean <- as.matrix(geodist::geodist(oto[c("longitude", "latitude")], measure = "geodesic"))

oto_geogr_dist <- geosphere::distm(cbind(oto$longitude, oto$latitude), fun=distGeo) # an alternative

# Least cost distance between beaches, considering 1,100 m depth
dist.lc <- marmap::lc.dist(trans.mat(centrAmer, min.depth=-1,  max.depth=-200), oto[c("longitude","latitude")], res="dist") 

# partial Mantel tests
vegan::mantel.partial(oto_m_ord_BC, DNA_dist_oto, d.euclidean)
vegan::mantel.partial(oto_m_ord_BC, DNA_dist_oto, oto_geogr_dist) # alternative
vegan::mantel.partial(oto_m_ord_BC, DNA_dist_oto, dist.lc)



