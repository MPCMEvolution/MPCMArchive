#Wrapper functions to calculate community phylogenetic metrics
#Taken from Metrics and Models of Community Phylogenetics; Pearse et al. 2014
#Code mostly by Matthew R. Helmus; some (cosmetic) commenting/restructuring/plotting by Will Pearse
# - greatly appreciated code contributions from Gustavo Carvalho and Marc Cadotte (as described in chapter)
#*Important* - at the time of writing, Helmus and Pearse are finishing an R package called 'pez' that implements all of this and more, and is probably (by now) easier to use. Have a look for it!

#################
#HEADERS#########
#################
require(MASS)
require(ade4)
require(picante)
require(apTreeshape)
require(spacodiR)
require(PVR)
require(phylobase)
require(phyloseq)
source("metrics.R")

#################
#DATASET#########
#SIMULATION######
#################

#Tree construction functions
# - lambda=0.63093 for 0.5 cov for 4 species tree
stree.b <- function(n=32,lambda=0.63093,method="Grafen") {compute.brlen(stree(n,type="b"),method=method,power=lambda)}
stree.l <- function(n=32,lambda=0.63093,method="Grafen") {compute.brlen(stree(n,type="l"),method=method,power=lambda)}

#Simulates the pa of a community
sim.phylostruct<-function(tree=stree.b(),sc=10,repulse=FALSE,random=FALSE)
{
  #V=stree.b(4)
  #sc=1
  #repulse=FALSE
  V<-tree
  if(is(V)[1]=="phylo") {V=vcv.phylo(V)}
  s<-dim(V)[1]
  y<-rep(0,s)
  if(repulse){iD<-t(chol(ginv(V)))} else {iD<-t(chol(V))}

  if(random)
  {
    rff<-round(runif(2,min=1,max=s))
    #while(length(unique(rff))<2) {rff<-round(runif(2,min=1,max=s))}
    y[rff]<-1
    rnums<-t(t(runif(s)))
    rand<-runif(1)
    y[(rand<rnums/(1+rnums))]<-1
  } else {
    #rnums<-t(t(rnorm(s)))
    #p<-sc*iD%*%rnums
    #exp(p)/(1+exp(p))
    rnums<-t(t(rnorm(s)))
    p<-sc*iD%*%rnums
    rand<-runif(1)
    y[rand<(exp(p)/(1+exp(p)))]<-1
  }
  return(y)
}

#Make a pa matrix with communities that are phylogenetically structured
pa.struct<-function(num=100,tree=stree.b(),sc=1,repulse=FALSE,random=FALSE)
{
  pa<-replicate(num,sim.phylostruct(tree=tree,sc=sc,repulse=repulse,random=random))
  rownames(pa)<-tree$tip
  colnames(pa)<-1:num
  return(t(pa))
}

#Simulates a community with abundances phylogenetically structured
sim.phylostruct.ab<-function(tree=stree.b(32),sc=1,abund.max=50,zeros.likely=FALSE,trun.scale=2,repulse=FALSE,random=FALSE)
{
  V<-tree
  if(is(V)[1]=="phylo") {V=vcv.phylo(V)}
  s<-dim(V)[1]
  y<-rep(0,s)
  if(repulse){iD<-t(chol(ginv(V)))} else {iD<-t(chol(V))}
  rownames(iD)<-rownames(V)
  colnames(iD)<-colnames(V)
  if(random)
  {
    rff<-round(runif(2,min=1,max=s))
    y[rff]<-1
    rnums<-t(t(runif(s)))
    rand<-runif(1)
    y[(rand<rnums/(1+rnums))]<-1
    if(zeros.likely) {y<-y*round(rnums*abund.max)} else {y<-round(rnums*abund.max)}
    if(length(y[y>0])>0) {y=y-(trun.scale*min(y[y>0])) }
    y[y<1]<-0
    y<-round(y)
  } else {
    rnums<-t(t(rnorm(s)))
    p<-sc*iD%*%rnums
    rand<-runif(1)
    y[rand<(exp(p)/(1+exp(p)))]<-1
    if(zeros.likely) {y<-y*round((exp(p)/(1+exp(p)))*abund.max)} else {y<-round((exp(p)/(1+exp(p)))*abund.max)}
    if(length(y[y>0])>0) {y=y-(trun.scale*min(y[y>0])) }
    y[y<1]<-0
    y<-round(y)
  }
  return(t(y))
}

#Make a pa matrix with communities that are phylogenetically structured
ab.struct<-function(num=10,tree=stree.b(),sc=1,abund.max=50,zeros.likely=TRUE,trun.scale=1,repulse=FALSE,random=FALSE)
{
  pa<-replicate(num,sim.phylostruct.ab(tree=tree,sc=sc,abund.max=abund.max,zeros.likely=zeros.likely,trun.scale=trun.scale,repulse=repulse,random=random),
                    simplify="matrix")
  rownames(pa)<-tree$tip
  colnames(pa)<-1:num
  return(t(pa))
}

#Presence/absence matrix of the feasible set of community structure
# - this is feasible only for small trees, and technically to construct the feasible set a tree is not needed only a species pool size, s
pa.feasible<-function(tree,wrtout=FALSE,file.name="feasible_set.csv")
{
  V=tree

  tate<-function(hold,s)
  {
    ap<-rep(0,s)
    ap[hold]<-1
    ap
  }
  if(is(V)=="phylo") {V=vcv.phylo(V)}
  s<-dim(V)[1]
  if(wrtout)
  {
    hold<-combn(s,2)
    pa<-t(apply(hold,2,tate,s))
    colnames(pa)<-tree$tip
    write.table(pa,file=file.name,sep=",",row.names=FALSE,col.names=TRUE)
    for(i in 3:s)
    {
      hold<-combn(s,i)
      write.table(t(apply(hold,2,tate,s)),file=file.name,append=TRUE,sep=",",row.names=FALSE,col.names=FALSE)
    }
    print(paste("pa feasible set written to file",file.name))

  } else {
    pa<-NULL
    for(i in 2:s)
    {
      hold<-combn(s,i)
      pa<-rbind(pa,t(apply(hold,2,tate,s)))
    }
    colnames(pa)<-tree$tip
    rownames(pa)<-1:dim(pa)[1]
    return(pa)

  }
}

feasible.n<-function(nspp)
{
 k<-1:nspp
 sum(choose(nspp,k))
}

#################
#PSV#############
#################

#Calculates the  psv value for a vector of species names given a phylo covariance matrix
psv.set<-function (index,Cmatrix)
{
  n <- length(index)
  C <- Cmatrix[index, index]
  PSV <- (n * sum(diag(as.matrix(C))) - sum(C))/(n *(n - 1))
  return(c(n,PSV))
}

#PSV of the feasible set

psv.feasible<-function(tree,pa=NULL,pasave=FALSE)
{
  if(is.null(pa)){pa<-pa.feasible(tree,wrtout=FALSE)}
  psv.dat<-psv(pa,tree)
  colnames(psv.dat)[1]<-"PSV"
  psv.dat<-data.frame(psv.dat,PSR=psv.dat[,2]*psv.dat[,1],N=1)
  if(pasave){return(list(aggregate(N~PSV+SR,data=psv.dat,FUN=sum),pa))} else {return(list(aggregate(N~PSV+SR,data=psv.dat,FUN=sum),NULL))}
}

#Aggregated table of the unique PSV values and the number of communties in the feasible set with those values

psv.unique<-function(pa,tree)
{
  psv.dat<-psv(pa,tree)
  colnames(psv.dat)[1]<-"PSV"
  psv.dat<-data.frame(psv.dat,PSR=psv.dat[,2]*psv.dat[,1],N=1)
  psv.dat<-psv.dat[!is.na(psv.dat[,1]),]
  return(aggregate(N~PSV+SR,data=psv.dat,FUN=sum))
}

#################
#METRIC##########
#CLASS WRAPPERS##
#################

#Calculates multiple SHAPE metrics for a community set
# - vecnums chooses the eigenvector to calculate sumvar in Diniz-Filho J.A.F., Cianciaruso M.V., Rangel T.F. & Bini L.M. (2011). Eigenvector estimation of phylogenetic and functional diversity. Functional Ecology, 25, 735-744.
shape.metrics<-function(pa,tree,vecnums=1) 
{
  nspp<-Ntip(tree)
  nsite<-dim(pa)[1]
  SR<-rowSums(pa>0)
  tree.dist<-cophenetic(tree)
  
  pdivs<-psd(pa,tree)[,c(-2:-3,-5:-6)]                      #Metric: PSV (Scaled MPD), PSR, SR
  
  mpd(pa,tree.dist)->MPD                                    #Metric: MPD
  pdivs<-data.frame(pdivs,MPD)

  PD<-pd.(pa,tree)[,1]                                       #Metric: PD
  pdivs<-data.frame(pdivs,PD)

  pdivs<-data.frame(pdivs,regPD=resid(lm(PD~rowSums(pa))))  #Metric: PD standardized by regression for species richness
  
  colless.<-function(pa.vec,tree,nams)                      #Metric: Ic Colless   NOTE THIS IN NOT STANDARDIZED BY YULE ETC.
  {
    if(sum(pa.vec)<3)
    {
      return(NA)
    } else {           
      return(colless(tipsubtree(tree.shape,nams[pa.vec!=0])))
    }
  }
  tree.shape=as.treeshape(tree)                             
  nams<-tree.shape$names
  Ic<-apply(pa,1,colless.,tree.shape,nams)
  pdivs<-data.frame(pdivs,Ic)

  gamma.<-function(pa.vec,tree,nams)                        #Metric: GAMMA
  {
    if(sum(pa.vec)<3)
    {
      return(NA)
    } else {
      return(gammaStat(drop.tip(tree,nams[pa.vec==0])))
    }
  }
  Gamma<-apply(pa,1,gamma.,tree,nams)
  pdivs<-data.frame(pdivs,Gamma)

  kk<-taxondive(pa,tree.dist)                              #Metric: Taxonomic Diversity Index D+ and SRD+
  pdivs<-data.frame(pdivs,TDI_Dplus=kk$Dplus,TDI_SRDplus=kk$Species * kk$Dplus)  #Note this is the same as MPD
  
  x <- PVRdecomp(tree)                                      #Metric: Eigen Vector SUMVAR
  evc<- x@Eigen$vectors
  sumvar<-function(pa.vec,evc,vecnums=1)                        
  {
   pv<-sum(apply(as.matrix(evc[pa.vec>0,vecnums]),2,var))
   if(sum(pa.vec>0)){return(pv)} else {return(NA)}
  }
  
  pdivs<-data.frame(pdivs,PEsumvar=apply(pa,1,sumvar,evc=evc,vecnums=vecnums))
  
  pa.<-pa[rowSums(pa>0)>1,]
  x<-phylo4d(tree,t(pa.))                                 
  pdabund<-data.frame(Eed=eed(x),Hed=hed(x))              #Metric: Eed and Hed #metrics from Cadotte et al 2010
  names(SR)<-paste("X",c(1:nsite),sep="")
  hold<-merge(SR,pdabund,by=0,all=T)
  rownames(hold)<-hold$Row.names
  hold<-hold[names(SR),-1:-2]
  pdivs<-data.frame(pdivs,hold)
  
  return(pdivs)
}

#EVENNESS metrics with only one value for a community set

evenness.metrics.single<-function(abund,tree)
{
  tree.dist<-cophenetic(tree)
  spacodi.calc(sp.plot = as.spacodi(abund), phy = tree)
  return(pdivs)
}

# Calculates multiple EVENNESS metrics for a community set
evenness.metrics<-function(abund,tree)
{
  tree.dist<-cophenetic(tree)
  SR<-rowSums(abund>0)
#  orig.site.names<-rownames(abund)
  nsite<-dim(abund)[1]
  nspp<-Ntip(tree)
  
  PSE<-pse(abund,tree)[,1]                                                             #Metric: PSE
  pdivs<-data.frame(PSE,raoDkk=raoD(abund,tree)$Dkk)                                   #Metric: Raos D Quadradic Entropy/Diversity Dkk for each communitiy 
  
  kk<-taxondive(abund,tree.dist)
  pdivs<-data.frame(pdivs,TDivI_D=kk$D,TDistI_Dstar=kk$Dstar,TDistI_Lambda=kk$Lambda)  #Metric: Taxonomic diversity indicies
  
  pent<-pe(abund,tree)
  pdivs<-data.frame(pdivs,Hp=pent)                                                     #Metric: Phylogenetic Entropy
  
  abund.<-abund[rowSums(abund>0)>1,]
  x<-phylo4d(tree,t(abund.))                                                           #Metric: Phylogenetic abundance evenness  and other metrics from Cadotte et al 2010    
  pdabund<-data.frame(PAE=pae(x),IAC=iac(x),Haed=haed(x),Eaed=eaed(x),Pst=simpson(x,"phylogenetic")) #Metric: Pst simpson diversity from Hardy and Senterre 2007
  names(SR)<-paste("X",c(1:nsite),sep="")
  hold<-merge(pdabund,SR,by=0,all=T)
  rownames(hold)<-hold$Row.names
  hold<-hold[names(SR),-1]
  colnames(hold)[dim(hold)[2]]<-"SR"
  pdivs<-data.frame(pdivs,hold)
  return(pdivs)
}

#Calculates multiple DISPERSION metrics for a community set
dispersion.metrics<-function(pa,tree)
{
  tree.dist<-cophenetic(tree)
  NRI<-ses.mpd(pa,tree.dist)$mpd.obs.z        #NRI
  NTI<-ses.mntd(pa,tree.dist)$mntd.obs.z      #NTI
  sesPD<- ses.pd(pa,tree)$pd.obs.z            #sesPD
  INND<-ses.mpd(pa,1/tree.dist)$mpd.obs.z     #INND
  
  pdivs<-data.frame(NRI,NTI,INND,sesPD)                    
}

#Dissimilarity
dissimilarity.metrics<-function(pa,tree,reps=1000)
{
  tree.dist<-cophenetic(tree)
  pa<-pa[rowSums(pa)>0,]
  UF<-unifrac(pa,tree)                             #UniFrac
  UF<-rowMeans(as.matrix(UF))

  PCD<-pcd.(pa,tree,reps=reps)                             #PCD
  
  PCDf<-rowMeans(as.matrix(as.dist(t(PCD$PCD),upper=TRUE)),na.rm=TRUE)
  PCDc<-rowMeans(as.matrix(as.dist(t(PCD$PCDc),upper=TRUE)),na.rm=TRUE)
  PCDp<-rowMeans(as.matrix(as.dist(t(PCD$PCDp),upper=TRUE)),na.rm=TRUE)

  PSor<-phylosor(pa,tree)
  PSor<-rowMeans(as.matrix(PSor))                  #PhyloSorensen

  pdivs<-data.frame(UF,PCDf,PCDc,PCDp,PSor)
  return(pdivs)                    
}
