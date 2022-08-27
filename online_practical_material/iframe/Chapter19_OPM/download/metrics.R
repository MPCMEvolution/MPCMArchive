#Novel code to calculate community phylogenetic metrics
#Taken from Metrics and Models of Community Phylogenetics; Pearse et al. 2014
#Code mostly by Matthew R. Helmus; some (cosmetic) commenting/restructuring/plotting by Will Pearse
# - greatly appreciated code contributions from Gustavo Carvalho and Marc Cadotte (as described in chapter)
#*Important* - at the time of writing, Helmus and Pearse are finishing an R package called 'pez' that implements all of this and more, and is probably (by now) easier to use. Have a look for it!

#################
#HEADERS#########
#################
require(ape)
require(ade4)
require(ecoPD)

#################
#FUNCTIONS#######
#################
#Phylogenetic Entropy by Allen et al 2009
#code adapted from http://www.cerradoecology.com/codes_files/Phylogenetic%20Entropy.R      ## Gustavo Carvalho (gustavo.bio@gmail.com)
#E.g., pe(abund,tree)
pe <- function(abund, tree)
{
  comm<-abund
  tree.phylo <- tree
  tree.phylog <- newick2phylog(write.tree(tree))

  ## checking the community data
  #if (any(rowSums(comm) == 0)) {
  #  warning("Your community data has sites without individuals")
  #  warning("The phylogenetic entropies for these sites will be 0")
  #}

  if (is.vector(comm)) {
    comm <- t(as.data.frame(comm))
    rownames(comm) <- "site"
  }

  ## stop if there is a mismatch between leaves labels and species in
  ## the community data
  if (!all(colnames(comm) %in% names(tree.phylog$leaves))) {
    stop("Leaves and community species do not match for pe (Phylogenetic Entropy)")
  }

  species <- colnames(comm)

  ## initializing the vector which will contain values for each site
  hp.sites <- c()

  ## beginning of the sites iteration
  for (i in 1:dim(comm)[1])
  {
    ## initilizing the vector which will hold values for the site in
    ## the loop
    hp <- c()

    ## species which occured at the site
    site.species <- species[which(comm[i,] > 0)]

    ## species which NOT occured at the site. They will be removed
    ## from the phylogenetic tree.
    other.species <- species[which(comm[i,] == 0)]

    ## proportions of occurrence of each species
    proportions <- comm[i,site.species] / sum(comm[i,])

    ## god, how I hate these tree conversions...
    ## also, this comparison is very ugly, it has to be a better
    ## way...
    ## TODO: Search for a better way to verify if an object is empty
    if (length(other.species) == 0) {
      partial.tree <- tree.phylog
    } else {
      if(length(site.species) == 1) {other.species<-c(other.species,site.species)}
      partial.tree <- drop.tip(tree.phylo, other.species)
      if (all(partial.tree$edge.length[1] == partial.tree$edge.length) | length(site.species) == 2)
      {
        hp.sites[rownames(comm)[i]] <- abs(sum(proportions * log(proportions) * partial.tree$edge.length))
        next
      }
      partial.tree <- newick2phylog(write.tree(drop.tip(tree.phylo, other.species)))
    }
    ## TODO: Some (I think) partial trees are not rooted. The results
    ## seem to be ok, but the paper states that Hp should be
    ## calculated from a rooted tree. Do not know why though. Check
    ## what can be done to keep partial trees rooted.
    if ("Root" %in% names(partial.tree$nodes))
    {
      partial.branches <-
      partial.tree$nodes[-c(length(partial.tree$nodes))]
    } else {
      partial.branches <- partial.tree$nodes
    }
    ## terminal branches sizes
    partial.leaves <- partial.tree$leaves

    ## first part of the calculations. Here we calculate the index for
    ## each terminal branch
    sum.leaves <- sum(partial.leaves *
                      proportions[names(partial.leaves)] *
                      log(proportions[names(partial.leaves)]))

    ## storing the first part of the calculation
    hp <- c(sum.leaves)

    ## initilizing the list that will hold the descending leaves for
    ## each branch
    descending.leaves <- list()

    ## determining the descending leaves for each branch
    for (j in names(partial.branches)) {
      if (all(partial.tree$parts[[j]] %in% names(partial.leaves))) {
        descending.leaves[[j]] <- partial.tree$parts[[j]]
      } else {
        branches <- partial.tree$parts[[j]][!partial.tree$parts[[j]]
                                            %in% names(partial.leaves)]
        leaves <- partial.tree$parts[[j]][partial.tree$parts[[j]] %in%
                                          names(partial.leaves)]
        for (k in branches) {
          leaves <- c(leaves, descending.leaves[[k]])
        }
        descending.leaves[[j]] <- leaves
      }
    }
    ## calculating the index for each internal branch
    for (j in names(partial.branches)) {
      sum.proportions.desc.leaves <-
        sum(proportions[descending.leaves[[j]]])
      hp <- c(hp, (partial.branches[[j]] * sum.proportions.desc.leaves
                   * log(sum.proportions.desc.leaves)))
    }
    ## putting it all together
    hp.sites[rownames(comm)[i]] <- abs(sum(hp))
  }
  #Make the 0 values NAs
  hp.sites[hp.sites==0]<-NA
  ## the end.
  return(hp.sites)
}

#Phylogenetic Community Dissimilarity from
#Ives A.R. & Helmus M.R. (2010). Phylogenetic metrics of community similarity. The American Naturalist, 176, E128-E142.
# - to be used in preference to the version in the package 'picante'
pcd. <- function(samp, tree, PSVmncd=NULL, PSVpool=NULL, reps=10^4)
{
  SSii<-PSVmncd
  SCii<-PSVpool

  # Make samp matrix a pa matrix
  samp[samp>0]<-1

	# convert trees to VCV format
	if (is(tree)[1] == "phylo")
  {
		if (is.null(tree$edge.length)) {tree <- compute.brlen(tree, 1)}    #If phylo has no given branch lengths
		tree <- prune.sample(samp, tree)
		V <- vcv.phylo(tree, cor = TRUE)
		samp <- samp[, tree$tip.label]
	} else {
		V <- tree
    species<-colnames(samp)
    preval<-colSums(samp)/sum(samp)
    species<-species[preval>0]
    V<-V[species,species]
    samp<-samp[,colnames(V)]
  }

  if(any(rowSums(samp)==0))  #it is possible that some species are not included in the phylogeny and as such their removal creates communities of zero spp
  {
   warning("Some of the communities have zero species.")
  }
  samp.orig<-samp
  samp<-samp[rowSums(samp)>0,]
  m <- dim(samp)[1]
	n <- dim(samp)[2]


  if (!is.null(SSii) & length(SSii)!=max(rowSums(samp))) {
        stop("The length of PSVmncd does not match the community with the highest species richness")
  }

	# m=number of communities; n=number of species; nsr=maximum sr value across all communities

  nsr <-max(rowSums(samp))
  if(is.null(SSii) & is.null(SCii))   #If the user already has calculated the mean conditional PSV values for all levels of SR
  {                                   #and the PSV of the species pool
    SSii <- array(0,nsr)
  	n1 <- 2
  	for (n2 in 1:nsr)
  	{
      temp <- array(0,reps)
  		for (t in 1:reps)
  		{
 			 rp <- sample(n)
		   pick1 <- rp[1:n1]

       rp <- sample(n)
		   pick2 <- rp[1:n2]

  			C11 <- V[pick1,pick1]
  			C22 <- V[pick2,pick2]
  			C12 <- V[pick1,pick2]

  			invC22 <- solve(C22)
  			S11 <- C11 - C12%*%invC22%*%t(C12)
  			SS11 <- (n1*sum(diag(S11))-sum(S11))/(n1*(n1-1))
  			temp[t] <- SS11
  		}
  		SSii[n2] <- mean(temp)
  	}
  	SCii=1-(sum(V)-sum(diag(V)))/(n*(n-1))
  }

	# calculate PCD
	PCD <- array(NA,c(m,m))
	PCDc <- array(NA,c(m,m))
	PCDp <- array(NA,c(m,m))
	for (i in 1:(m-1))
	{
		for (j in (i+1):m)
		{
			pick1 <- (1:n)[samp[i,]==1]
 			pick2 <- (1:n)[samp[j,]==1]

			n1 <- length(pick1)
			n2 <- length(pick2)

			C <- V[c(pick1, pick2),c(pick1, pick2)]

			C11 <- C[1:n1,1:n1]
			C22 <- C[(n1+1):(n1+n2),(n1+1):(n1+n2)]
			C12 <- C[1:n1,(n1+1):(n1+n2)]
			if(is.null(dim(C12)))
      {
        if(is.null(dim(C22))){C12<-as.matrix(C12)} else {C12<-t(as.matrix(C12))}
      }

			invC11 <- solve(C11)
			S22 <- C22 - t(C12)%*%invC11%*%C12

			invC22 <- solve(C22)
			S11 <- C11 - C12%*%invC22%*%t(C12)
      if(n1>1)
      {
			 SC11 <- (n1*sum(diag(C11))-sum(C11))/(n1*(n1-1))
       SS11 <- (n1*sum(diag(S11))-sum(S11))/(n1*(n1-1))
			} else {          #Added to deal with communities of only one species
       SC11 <- 0
       SS11 <- S11
      }
      if(n2>1)
      {
        SC22 <- (n2*sum(diag(C22))-sum(C22))/(n2*(n2-1))
        SS22 <- (n2*sum(diag(S22))-sum(S22))/(n2*(n2-1))
      } else {          #Added to deal with communities of only one species
        SC22 <- 0
        SS22 <- S22
      }

      if((n1+n2)==2){ #both communities have only one species
        D=(n1*SS11 + n2*SS22) #we do not standardize by the unconditional PSVs, since the PSV of a 1 spp community is zero (or undefined)
      } else {
        D=(n1*SS11 + n2*SS22)/(n1*SC11 + n2*SC22)   #if one of the communities is of one species, then the unconditional PSV is 0
      }
			a <- length(unique(c(pick1, pick2)))
			b <- length(pick1)-a
			cc <- length(pick2)-a
			dsor <- 2*a/(2*a+b+cc) - 1

			pred.D <- (n1*SSii[n2]+n2*SSii[n1])/(n1*SCii+n2*SCii)
			pred.dsor <- 1 - 2*n1*n2/((n1+n2)*n)

			PCD[i,j] <- D/pred.D
			PCDc[i,j] <- dsor/pred.dsor
			PCDp[i,j] <- PCD[i,j]/PCDc[i,j]
		}
	}
	colnames(PCD)<-rownames(samp)
  rownames(PCD)<-rownames(samp)
  colnames(PCDc)<-rownames(samp)
  rownames(PCDc)<-rownames(samp)
  colnames(PCDp)<-rownames(samp)
  rownames(PCDp)<-rownames(samp)

  return(list(PCD=PCD,PCDc=PCDc,PCDp=PCDp,PSVmncd=SSii,PSVpool=SCii))
}

#################################################################################################################################################################################################################
# Quick plot of some of the different PD values
caterpillar <- function(data, center=c("mean", "median"), sort=TRUE,
  ...) {

  method <- attr(data, "method")
  center <- match.arg(center)
  if(class(data)=="matrix") {
    data <- as.data.frame(data)
  }
  if(sort) {
    data <- data[order(data$pd.obs),]
  }

  # set up plotting area
  ind <- seq_len(nrow(data))
  xlab <- paste("PD",
    if(!is.null(method)) paste(" (", method, ")", sep="") else NULL,
    sep="")
  plot(data[, "0.025"], ind, type = "n", yaxt = "n", xlab = xlab,
    ylab = NA, xlim = range(c(data[,"0.025"], data[, "0.975"])),
    ...)
  axis(2, at = ind, labels = rownames(data), las = 2)

  # plot 95% bounds
  arrows(data[, "0.025"], ind, data[, "0.975"], ind, code = 3,
    angle = 90, length = 0.01)

  # plot mean or median
  if (center == "mean") {
    points(data[, "mean"], ind, pch = 3)
  } else if (center == "median") {
    points(data[, "median"], ind, pch = 3)
  }

  # plot actual pd values
  points(data[, "pd.obs"], ind, pch = 1, col = "red")

}
#########################################################################################################################################################################################################
##
## Shannon's index
##


#The Cadotte metrics define a pd function, but the picante pd function is still used.
pd.<-function (samp, tree, include.root = TRUE) 
{
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute pd")
    }
    species <- colnames(samp)
    SR <- rowSums(ifelse(samp > 0, 1, 0))
    nlocations = dim(samp)[1]
    nspecies = dim(samp)[2]
    PDs = NULL
    for (i in 1:nlocations) {
        present <- species[samp[i, ] > 0]
        treeabsent <- tree$tip.label[which(!(tree$tip.label %in% 
            present))]
        if (length(present) == 0) {
            PD <- 0
        }
        else if (length(present) == 1) {
            if (!is.rooted(tree) || !include.root) {
                warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species sampunities. Single species sampunity assigned PD value of NA.")
                PD <- NA
            }
            else {
                PD <- node.age(tree)$ages[which(tree$edge[, 2] == 
                  which(tree$tip.label == present))]
            }
        }
        else if (length(treeabsent) == 0) {
            PD <- sum(tree$edge.length)
        }
        else {
            sub.tree <- drop.tip(tree, treeabsent)
            if (include.root) {
                if (!is.rooted(tree)) {
                  stop("Rooted tree required to calculate PD with include.root=TRUE argument")
                }
                sub.tree.depth <- max(node.age(sub.tree)$ages)
                orig.tree.depth <- max(node.age(tree)$ages[which(tree$edge[, 
                  2] %in% which(tree$tip.label %in% present))])
                PD <- sum(sub.tree$edge.length) + (orig.tree.depth - 
                  sub.tree.depth)
            }
            else {
                PD <- sum(sub.tree$edge.length)
            }
        }
        PDs <- c(PDs, PD)
    }
    PDout <- data.frame(PD = PDs, SR = SR)
    rownames(PDout) <- rownames(samp)
    return(PDout)
}

#################
#PGLMMM##########
#################
pglmm.data.<-function (modelflag = 1, sim.dat = NULL, samp = NULL, tree = NULL, traits = NULL, env = NULL, Vcomp = NULL) 
{
    if (!is.null(sim.dat)) {
        tree <- sim.dat$Vphylo
        Vcomp <- sim.dat$Vcomp
        samp <- sim.dat$Y
        traits <- sim.dat$bspp1
        env <- sim.dat$u
    }
    is.empty <- function(x) {
        length(x) == 0
    }
    if (is.empty(samp)) {
        stop("sample matrix (Y) is empty")
    }
    if (is(tree)[1] == "phylo") 
    {
        if (is.null(tree$edge.length)) 
        {
            tree <- compute.brlen(tree, 1)
        }
        tree <- prune.sample(samp, tree)
        samp <- samp[, tree$tip.label]
        V <- vcv.phylo(tree, corr = TRUE)
        species <- colnames(samp)
        preval <- colSums(samp)/sum(sum(samp))
        species <- species[preval > 0]
        V <- V[species, species]
        if(!is.null(Vcomp)){Vcomp <- Vcomp[species, species]}
        samp <- samp[, colnames(V)]
        if(!is.null(traits)){traits <- as.matrix(traits[species, ])}
    } else {
        V <- tree
        species <- colnames(samp)
        preval <- colSums(samp)/sum(sum(samp))
        species <- species[preval > 0]
        V <- V[species, species]
        Vcomp <- Vcomp[species, species]
        samp <- samp[, colnames(V)]
        traits <- as.matrix(traits[species, ])
    }
    Y <- samp
    X <- traits
    nsites <- dim(Y)[1]
    nspp <- dim(Y)[2]
    YY <- t(Y)
    YY <- as.matrix(as.vector(as.matrix(YY)))
    if (modelflag == 1) {
        Vfullspp <- kronecker(diag(nsites), V)
        Vfullsite <- kronecker(diag(nsites), matrix(1, nspp, 
            nspp))
        VV <- list(Vfullspp = Vfullspp, Vfullsite = Vfullsite)
        XX <- kronecker(matrix(1, nsites, 1), diag(nspp))
        return(list(YY = YY, VV = VV, XX = XX))
    }
    if (modelflag == 6) {
        if (is.null(Vcomp)) {
            compscale <- 1
            Vcomp <- solve(V, diag(nspp))
            Vcomp <- Vcomp/max(Vcomp)
            Vcomp <- compscale * Vcomp
        }
        
        Vfullspp <- kronecker(diag(nsites), Vcomp)
        Vfullsite <- kronecker(diag(nsites), matrix(1, nspp, 
            nspp))
        VV <- list(Vfullspp = Vfullspp, Vfullsite = Vfullsite)
        XX <- kronecker(matrix(1, nsites, 1), diag(nspp))
        return(list(YY = YY, VV = VV, XX = XX))
    }

    if (modelflag == 2) {
        u <- scale(U)
        U <- matrix(env, nrow = length(env), ncol = 1)
        U <- kronecker(u, matrix(1, nspp, 1))
        Vfullspp <- kronecker(matrix(1, nsites, nsites), diag(nspp))
        VfullsppV <- kronecker(matrix(1, nsites, nsites), V)
        VfullUCU <- diag(as.vector(U)) %*% Vfullspp %*% diag(as.vector(U))
        VfullUCUV <- diag(as.vector(U)) %*% VfullsppV %*% diag(as.vector(U))
        Vfullsite <- kronecker(diag(nsites), matrix(1, nspp, 
            nspp))
        VV <- list(VfullUCU = VfullUCU, VfullUCUV = VfullUCUV, 
            Vfullsite = Vfullsite)
        XXspp <- kronecker(matrix(1, nsites, 1), diag(nspp))
        XX <- cbind(U, XXspp)
        YY <- as.vector(t(Y))
        return(list(YY = YY, VV = VV, XX = XX))
    }
    if (modelflag == 3) {
        u <- scale(U)
        U <- kronecker(u, matrix(1, nspp, 1))
        Vfullsite <- kronecker(diag(nsites), matrix(1, nspp, 
            nspp))
        if (is.null(Vcomp)) {
            compscale <- 1
            Vcomp <- solve(V, diag(nspp))
            Vcomp <- Vcomp/max(Vcomp)
            Vcomp <- compscale * Vcomp
        }
        Vfullcomp <- kronecker(diag(nsites), Vcomp)
        VV <- list(Vfullcomp = Vfullcomp, Vfullsite = Vfullsite)
        XXspp <- kronecker(matrix(1, nsites, 1), diag(nspp))
        XX <- cbind(XXspp, ((U %*% matrix(1, 1, nspp)) * XXspp))
        YY <- as.vector(t(Y))
        return(list(YY = YY, VV = VV, XX = XX))
    }
    if (modelflag == 4) {
        Vfulltrait <- kronecker(diag(nsites), traits %*% t(traits))
        traitscale4 <- 100
        Vfulltrait <- traitscale4 * Vfulltrait
        Vfullsite <- kronecker(diag(nsites), matrix(1, nspp, 
            nspp))
        VV <- list(Vfulltrait = Vfulltrait, Vfullsite = Vfullsite)
        XXspp <- kronecker(matrix(1, nsites, 1), diag(nspp))
        XX <- XXspp
        YY <- as.vector(t(Y))
        return(list(YY = YY, VV = VV, XX = XX))
    }
    if (modelflag == 5) {
        Vfulltrait <- kronecker(diag(nsites), traits %*% t(traits))
        traitscale5 <- 10
        Vfulltrait <- traitscale5 * Vfulltrait
        Vfullspp <- kronecker(diag(nsites), V)
        Vfullsite <- kronecker(diag(nsites), matrix(1, nspp, 
            nspp))
        VV <- list(Vfulltrait = Vfulltrait, Vfullspp = Vfullspp, 
            Vfullsite = Vfullsite)
        XXspp <- kronecker(matrix(1, nsites, 1), diag(nspp))
        XX <- XXspp
        YY <- as.vector(t(Y))
        return(list(YY = YY, VV = VV, XX = XX))
    }
}
