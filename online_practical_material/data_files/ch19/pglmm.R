#Code to reproduce PGLMMs in Metrics and Models of Community Phylogenetics; Pearse et al. 2014
#Code mostly by Matthew R. Helmus; some (cosmetic) commenting/restructuring by Will Pearse
#*Important* - at the time of writing, Helmus and Pearse are finishing an R package called 'pez' that implements much of this. Do a Google search for that package now!
#*Important* - these models take a *very* long time to fit. Time for a cup of tea!
#2014-1-22

#################
#HEADERS#########
#################
library(Biobase)
source("functions.R") #...which loads 'metrics.R'

#################
#SIMULATE DATA###
#################
tree<-stree.b(8)
psv.f.dat<-psv.feasible(tree,pasave=TRUE)#NOTE this takes a very long time to compute for a large tree since it makes all species combinations
f.pa<-psv.f.dat[[2]]
psv.f.dat<-psv.feasible(tree)[[1]]
num <- 5000#number of communities
sc<-10#value of c (see text of chapter)
dat <- pglmm.data.(samp=f.pa,tree=tree)

#################
#PGLMM 1#########
#################
#underdispersed set
pa <- u.pa[rowSums(u.pa)>1,]
dat <- pglmm.data.(samp=pa[1:50,],tree=tree)
u.fit.a <- pglmm.fit(dat,maxit=50,exitcountermax=100)

#overdispersed set
pa <- o.pa[rowSums(o.pa)>1,]
dat <- pglmm.data.(samp=pa[1:50,],tree=tree,modelflag=1)
o.fit.a <- pglmm.fit(dat,maxit=50,exitcountermax=100)


#################
#PGLMM 2#########
#################
#Note: change in 'modelflag'
#underdispersed set
pa <- u.pa[rowSums(u.pa)>1,]
dat <- pglmm.data.(samp=pa[1:50,],tree=tree,modelflag=6)
u.fit.r <- pglmm.fit(dat,maxit=50,exitcountermax=100)

#overdispersed set
pa <- o.pa[rowSums(o.pa)>1,]
dat <- pglmm.data.(samp=pa[1:50,],tree=tree,modelflag=6)
o.fit.r <- pglmm.fit(dat,maxit=50,exitcountermax=100)