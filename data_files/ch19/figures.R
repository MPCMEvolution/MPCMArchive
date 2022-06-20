#Code to reproduce figure 2 in Metrics and Models of Community Phylogenetics; Pearse et al. 2014
#Code mostly by Matthew R. Helmus; some (cosmetic) commenting/restructuring/plotting by Will Pearse
# - greatly appreciated code contributions from Gustavo Carvalho and Marc Cadotte (as described in chapter)
#*Important* - at the time of writing, Helmus and Pearse are finishing an R package called 'pez' that implements all of this and more, and is probably (by now) easier to use. Have a look for it!
#Notes: because of R's plotting order, the plotting order is A-D-B-E-C-F. This script takes *hours* to run; consider writing RDS files (see comments to see which variabels should be saved)

#################
#HEADERS#########
#################
source("functions.R") #...which loads 'metrics.R'

#################
#SIMULATE DATA###
#################
tree <- stree.b(8)
psv.f.dat <- psv.feasible(tree,pasave=TRUE)#NOTE this takes a very long time to compute for a large tree since it makes all species combinations
f.pa <- psv.f.dat[[2]]
psv.f.dat <- psv.feasible(tree)[[1]]
num <- 5000#number of communities
sc<-10#value of c (see text of chapter)

#################
#SETUP FIGURE####
#################
graphics.off()
pdf("figure_2.pdf", width=10,height=8)
par(las=1,mfcol=c(2,3),mar=c(5, 4, 3,0),cex=1)

#################
#2A##############
#################
pt.cex <- 10*psv.f.dat$N/max(psv.f.dat$N)#size of pixels
plot(psv.f.dat$SR, psv.f.dat$PSV, pch = 21, cex = pt.cex, col = "black", lwd = 1,ylim=c(0,1) ,ylab="phylogenetic species variability",xlab="",main="feasible")
points(psv.f.dat$SR, psv.f.dat$PSV, pch = 20, cex=0.3, col="darkgrey")

#################
#2D##############
#################
#shape.divs<-readRDS("f.shape.rds")
#disp.divs<-readRDS("f.disp.rds")
#dis.divs<-readRDS("f.dis.rds")
shape.divs <- shape.metrics(f.pa,tree)
disp.divs <- dispersion.metrics(f.pa,tree)
dis.divs <- evenness.metrics(f.pa,tree,reps=2000)
divs <- data.frame(shape.divs,disp.divs)
divs <- divs[!apply(is.na(divs),1,any),]
divs <- divs[colnames(divs)!="SR"]
names(divs)[c(8:10,16)] <- c("D+", "SRD+", "PE", "SESpd")
hclust.pdivs.f <- hclust(dist(t(scale(divs))),method="complete")
par(mar=c(0, 0, 0, 1))
hclust.pdivs.f <- as.phylo(hclust.pdivs.f)
plot(hclust.pdivs.f,edge.width = 2, underscore=TRUE, label.offset=1)
tiplabels(tip=13:15, pch=21, adj=1, cex=1.5, col="black", bg="black")
tiplabels(tip=seq_along(hclust.pdivs.f$tip.label)[-13:-15], pch=21, adj=1, cex=1.5, bg="white", lwd=2)

#################
#2B##############
#################
u.ab <- u.pa <- pa <- ab.struct(num=num,tree,sc=sc)
u.pa[u.pa>0] <- 1
psv.u.dat <- psv.unique(u.pa,tree)
pt.cex <- 10*psv.u.dat$N/max(psv.u.dat$N)#size of pixels
par(mar=c(5, 2, 3, 2))
plot(psv.u.dat$SR, psv.u.dat$PSV, pch = 21, cex = pt.cex, col = "black", lwd = 1,xlim=c(1.8,8),ylim=c(0,1),ylab="",xlab="species richness",yaxt="n",main="attraction")
points(psv.u.dat$SR, psv.u.dat$PSV, pch = 20, cex = 0.3, col = "darkgrey")

#################
#2E##############
#################
#shape.divs.u<-readRDS("u.shape.rds")
#disp.divs.u<-readRDS("u.disp.rds")
#evenness.divs.u<-readRDS("u.diversity.rds")
shape.divs.u <- shape.metrics(u.pa,tree)
disp.divs.u <- dispersion.metrics(u.pa,tree)
evenness.divs.u <- evenness.metrics(u.ab,tree)
divs.u <- data.frame(shape.divs.u,disp.divs.u,evenness.divs.u)
divs.u <- divs.u[!apply(is.na(divs.u),1,any),]
divs.u <- divs.u[colnames(divs.u)!="SR"]
divs.u <- divs.u[colnames(divs.u)!="SR.1"]
names(divs.u)[c(8:10, 16, 19:21)] <- c("D+", "SRD+", "PE", "SESpd", "TD_D", "Dstar", "TD_L")
par(mar=c(0, 0, 0, 1))
hclust.pdivs.u <- hclust(dist(t(scale(divs.u))))
hclust.pdivs.u <- as.phylo(hclust.pdivs.u)
plot(hclust.pdivs.u, edge.width = 2, underscore=TRUE, label.offset=1.5, x.lim=c(-1,59))
tiplabels(tip=13:15, pch=21, adj=1, cex=1.5, bg="black")
tiplabels(tip=17:27, pch=21, adj=1, cex=1.5, bg="darkgrey", col="darkgrey")
tiplabels(tip=seq_along(hclust.pdivs.u$tip.label)[c(-13:-15,-17:-27)], adj=1, cex=1.5, bg="white", pch=21, lwd=2)

#################
#2C##############
#################
o.ab <- o.pa <- pa <- ab.struct(num=num,tree,sc=sc,repulse=TRUE)
o.pa[o.pa>0] <- 1
psv.o.dat <- psv.unique(o.pa,tree)
pt.cex <- 10*psv.o.dat$N/max(psv.o.dat$N)#point size
par(mar=c(5, 2, 3, 2))
plot(psv.o.dat$SR, psv.o.dat$PSV, pch = 21, cex = pt.cex, col = "black", lwd = 1,xlim=c(1.8,8),ylim=c(0,1),ylab="",xlab="",yaxt="n",main="repulsion")
points(psv.o.dat$SR, psv.o.dat$PSV, pch = 20, cex = 0.3, col = "darkgrey")

#################
#2F##############
#################
#shape.divs.o<-readRDS("o.shape.rds")
#disp.divs.o<-readRDS("o.disp.rds")
#evenness.divs.o<-readRDS("o.diversity.rds")
shape.divs.o<-shape.metrics(o.pa,tree)
disp.divs.o<-dispersion.metrics(o.pa,tree)
evenness.divs.o<-evenness.metrics(o.ab,tree)
divs.o <- data.frame(shape.divs.o,disp.divs.o,evenness.divs.o)
divs.o <- divs.o[!apply(is.na(divs.o),1,any),]
divs.o <- divs.o[colnames(divs.o)!="SR"]
divs.o <- divs.o[colnames(divs.o)!="SR.1"]
names(divs.o)[c(8:10, 16, 19:21)] <- c("D+", "SRD+", "PE", "SESpd", "TD_D", "Dstar", "TD_L")
par(mar=c(0, 0, 0, 1))
hclust.pdivs.o <- hclust(dist(t(scale(divs.o))))
hclust.pdivs.o <- as.phylo(hclust.pdivs.o)
plot(hclust.pdivs.o,edge.width = 2, underscore=TRUE, label.offset=1.5)
tiplabels(tip=13:15, pch=21, adj=1, cex=1.5, bg="black")
tiplabels(tip=17:27, pch=21, adj=1, cex=1.5, bg="darkgrey", col="darkgrey")
tiplabels(tip=seq_along(hclust.pdivs.o$tip.label)[c(-13:-15,-17:-27)], adj=1, cex=1.5, bg="white", pch=21, lwd=2)

#################
#FINISH FIGURE###
#################
dev.off()
