######insert module
library(geiger)
library(ape)
library(mnormt)
#library(tcltk2)


####construct cato type data
xDataFactor = function(xData, factorName){
  #function which can return the xMatrix which involves the factor matrix
  #factorName is the name of factor variables in the regression analysis
  name = colnames(xData)
  nFactor = length(factorName)
  factorSummary = array(NA, nFactor)
  nDim = dim(xData)[2]
  xDataNew = array(NA, c(dim(xData)[1], 1))
  xDataNew[, 1] = xData[, 1]
  for(i in 2:nDim){
    if(sum(name[i] == factorName) == 0){
    #not sure whether this can deal with the name stuff
      xDataNew = cbind(xDataNew, xData[, i])
    }
  }
  colnames(xDataNew) = colnames(xData)[1:(nDim - length(factorName))]
  indicator = nDim - length(factorName)
  for(i in 1:nFactor){
    dataTemp = unique(xData[, indicator + i])
    factorSummary[i] = length(dataTemp) - 1
    xDataTemp = array(0, c(dim(xData)[1], (length(dataTemp) - 1)))
    for(k in 2:(length(dataTemp))){
      xDataTemp[which(xData[, indicator + i] == dataTemp[k]), k - 1]  = 1
    }
    colnames(xDataTemp) = paste(paste("factor.", factorName[i], sep = ""), dataTemp[2:length(dataTemp)], sep = "")
    xDataNew = cbind(xDataNew, xDataTemp)
  }
  return (list(xDataNew, factorSummary))
}




####missing data cleaning#######
missingDataCleaning = function(xData, yData, treeDataAll, missingList = c()){
  #this function can deal with the data cleaning
  #missing data cleaning should deal with the dataset after the factor name analysis
  #so the xData should be xData after the factor variable process
  #######data initialize#######
  numOfTrees = length(treeDataAll)
  #xDataAll = vector("list", numOfTrees)
  #yDataAll = vector("list", numOfTrees)
  numOfSpecies = dim(xData)[1]
  n = length(missingList)
  nameO = rownames(xData)
  tr = treeDataAll
  trTemp = vector("list", numOfTrees)
  for(i in 1:numOfTrees){
    if (n > 0){
      treeTemp = drop.tip(tr[[i]], missingList)
    }
    if (n == 0){
      treeTemp = treeDataAll[[i]]
    }
    trTemp[[i]] = treeTemp
    mm = getCovarianceMatrixLambda(treeTemp, 0, 1.0)
    names = colnames(mm)
    dataO = array(NA, c(dim(xData)[1] - n, dim(xData)[2]))
    yDataO = array(NA, dim(xData)[1] - n)
    for(k in 1:(numOfSpecies - n)){
      index = which(nameO == names[k])
      dataO[k, ] = xData[index, ]
      yDataO[k] = yData[index]
    }
    colnames(dataO) = colnames(xData)
    rownames(dataO) = names
    xDataAll = dataO
    yDataAll = yDataO   
    #missingList is the name of NA species
  }
  return (list(xDataAll, yDataAll, trTemp))
}

                     
##Some Tree Func#####
getCovarianceMatrixLambda = function(treeData, lambda, lambdaUpperBound = 1.0){
  #this lambda is the logit(lambda), used for MH algorithm
  lambdaScale = exp(lambda)/(1 + exp(lambda)) * lambdaUpperBound
  covariance = vcv.phylo(rescale(treeData, "lambda", lambdaScale))
  return (covariance)
}

getCovarianceMatrixLambdaTwo = function(treeData, lambda){
  #this lambda is true lambda, used for analysis of result
  covariance = vcv.phylo(rescale(treeData, "lambda", lambda))
  return (covariance)
}


#Rho function is not used in this example                
getCovarianceMatrixRho = function(CCMatrix, rho){
  covariance = 1 - (1 - CCMatrix)^rho
  return (covariance)
}

getCovarianceMatrixKappa = function(treeData, kappa, kappaUpperBound = 1.){
  #this kappa is the logit(kappa), used for MH algorithm
  kappaScale = exp(kappa)/(1 + exp(kappa)) * kappaUpperBound
  covariance = vcv.phylo(rescale(treeData, "kappa", kappaScale))
  return (covariance)
}


getCovarianceMatrixKappaTwo = function(treeData, kappa){
  #this kappa is true kappa, used for analysis of result
  covariance = vcv.phylo(rescale(treeData, "kappa", kappa))
  return (covariance)
}



XOXXOY = function(xData, yData, omegaCurrent, gamma){
#add factor variable into the xData, 04/11/10
#For the xData part, here, it should be xData after the data factor process
#added parameter xData, yData, into the fuc, 03/10/10
#generate XOX and XOY, O is the matrix omega.
#omegaCurrent is the inverse of the covariance matrix
#gamma should be the original gamma, 05/07/10
  xDataNew = xData[, which(gamma == 1)]
  XOX = t(xDataNew) %*% omegaCurrent %*% xDataNew
  XOY = t(xDataNew) %*% omegaCurrent %*% yData
  return (list(XOX, XOY))
}


densityFunc = function(xData, yData, omegaCurrent, gamma, sigmaSquareCurrent, lambdakappa){
#returns the log of density function.
  xoxxoy = XOXXOY(xData, yData, omegaCurrent, gamma)
  XOX = xoxxoy[[1]]
  XOY = xoxxoy[[2]][,1]
  YOY = t(yData) %*% omegaCurrent %*% yData
  betaHat = (solve(XOX) %*% XOY)[, 1]
  q = sum(gamma)
  ll = exp(lambdakappa)/(1 + exp(lambdakappa))
  result = log(ll) + log(1 - ll) - (dim(omegaCurrent)[1] - sum(gamma))/2 * log(sigmaSquareCurrent) + 1/2 * log(det(omegaCurrent)) - 1/2 * log(det(XOX)) + (- .5*(YOY - t(XOY) %*% solve(XOX) %*% (XOY))*(1/(0.0 + sigmaSquareCurrent)))[1,1]
  return (result)
}


densityFuncTwo = function(xData, yData, gammaNew, omegaCurrent, sigmaSquareNew, betaNew){
  betaNew[is.na(betaNew)] = 0
  meanBeta = xData %*% (betaNew)
  sigmaX = sigmaSquareNew * solve(omegaCurrent + .0001 * dim(omegaCurrent)[1])
  result = dmvnorm(x = yData, mean = meanBeta, sigma  = sigmaX, log = T)
  return (result)
}

priorGamma = function(gamma){
  #returns the prior distribution of gamma
  return (0)
}

priorSigmaSquare = function(sigmaSquareCurrent){
  #returns the prior distribution of sigma square
  return (0)
}

priorLambda = function(lambda){
  #returns the prior distribution of lambda
  return (0)
}

priorKappa = function(kappa){
  #returns the prior distribution of kappa
  return (0)
}

priorP = function(p){
  #returns the prior distribution of p, which is the model selection variable
  return (0)
}

posteriorFunc = function(xData, yData, omegaCurrent, gamma, sigmaSquareCurrent, lambdakappa){
  #returns the posterior function of model
  density = densityFunc(xData, yData, omegaCurrent, gamma, sigmaSquareCurrent, lambdakappa)
  posterior = density + priorGamma(gamma) + priorSigmaSquare(sigmaSquareCurrent)
  return (posterior)
}


#####preanalysis part######
#pre analysis is the analysis which can find approximate value of parameters which can give good starting value of the analysis
lik = function(xData, yData, lambda, p, treeData){
#xData here should be xData after the factor cleaning part, so the gamma is corresponding vectors, 04/11/10
  if(p == 1){
    omega = solve(getCovarianceMatrixLambdaTwo(treeData, lambda) + .0001 * diag(length(yData)))
  }
  if(p == 0){
    omega = solve(getCovarianceMatrixKappaTwo(treeData, lambda) + .0001 * diag(length(yData)))
  }
  XOX = t(xData) %*% omega %*% xData
  XOY = t(xData) %*% omega %*% yData
  YOY = t(yData) %*% omega %*% yData
  n = dim(omega)[1]
  const = (1/n) * (YOY - t(XOY) %*% solve(XOX) %*% (XOY))
  result = -n/2 * log(const) + 1/2 * log(det(omega))
  return (result[1, 1])
}

estimation = function(xData, yData, coef, p, treeData){
  lambda = coef[1]
  sigmaSquare = coef[2]
  if(p == 1){
    omega = solve(getCovarianceMatrixLambdaTwo(treeData, lambda) + .0001 * diag(length(yData)))
  }
  if(p == 0){
    omega = solve(getCovarianceMatrixKappaTwo(treeData, lambda) + .0001 * diag(length(yData)))
  }
  XOX = t(xData) %*% omega %*% xData
  XOY = t(xData) %*% omega %*% yData
  YOY = t(yData) %*% omega %*% yData
  beta = solve(XOX) %*% (XOY)
  varCoef = diag(sigmaSquare * solve(XOX))
  tStat = beta / sqrt(varCoef)
  return (cbind(beta, varCoef, tStat))
}

preAnalysis = function(xData, yData, treeData){
  N = 200
  lambdaTest = seq(.01, .99, length = N)
  result = array(NA, N)
  for(i in 1:N){
    result[i] = lik(xData, yData, lambdaTest[i], treeData)
  }
  lambda = lambdaTest[which(result == max(result))]
  omega = solve(getCovarianceMatrixLambdaTwo(treeData, lambda) + .0001 * diag(length(yData)))
  XOX = t(xData) %*% omega %*% xData
  XOY = t(xData) %*% omega %*% yData
  YOY = t(yData) %*% omega %*% yData
  n = dim(omega)[1]
  sigmaSquare = (1/n) * (YOY - t(XOY) %*% solve(XOX) %*% (XOY))
  coef = c(lambda, sigmaSquare)
  re = estimation(xData, yData, coef, treeData)
  tStat = re[, 3]
  beta = re[, 1]
  gammaCurrent = array(1, length(tStat))
  #betaCurrent = array(NA, length(gammaCurrent))
  #betaCurrent[which(gammaCurrent == 1)] = beta[which(gammaCurrent == 1)]
  betaCurrent = beta
  lambdaCurrent = kappaCurrent = log(lambda/(1 - lambda))
  sigmaSquareCurrent = sigmaSquare
  pCurrent = 1
  currentValue = list(gammaCurrent, sigmaSquareCurrent, pCurrent, lambdaCurrent, kappaCurrent, betaCurrent, -10)
  return (currentValue)
}


initialValueTest = function(xData, yData, treeData, p, upperBound = 1.){
  N = 200
  lambdaTest = seq(.01, upperBound * .99, length = N)
  result = array(NA, N)
  for(i in 1:N){
    result[i] = lik(xData, yData, lambdaTest[i], p, treeData)
  }
  lambda = lambdaTest[which(result == max(result))]
  if(p == 1){
    omega = solve(getCovarianceMatrixLambdaTwo(treeData, lambda) + .0001 * diag(length(yData)))
  }
  if(p == 0){
    omega = solve(getCovarianceMatrixKappaTwo(treeData, lambda) + .0001 * diag(length(yData)))
  }
  XOX = t(xData) %*% omega %*% xData
  XOY = t(xData) %*% omega %*% yData
  YOY = t(yData) %*% omega %*% yData
  n = dim(omega)[1]
  sigmaSquare = (1/n) * (YOY - t(XOY) %*% solve(XOX) %*% (XOY))
  coef = c(lambda, sigmaSquare)
  re = estimation(xData, yData, coef, p, treeData)
  tStat = re[, 3]
  beta = re[, 1]
  gammaCurrent = array(1, length(tStat))
  #betaCurrent = array(NA, length(gammaCurrent))
  #betaCurrent[which(gammaCurrent == 1)] = beta[which(gammaCurrent == 1)]
  betaCurrent = beta
  lambdaCurrent = kappaCurrent = log(lambda/(upperBound - lambda))
  sigmaSquareCurrent = sigmaSquare
  pCurrent = 1
  currentValue = list(gammaCurrent, sigmaSquareCurrent, pCurrent, lambdaCurrent, kappaCurrent, betaCurrent, -10)
  return (currentValue)
}



#####Bayesian Posterior Draw Part#####
generateBeta = function(xData, yData, omegaCurrent, gammaCurrent, sigmaSquareCurrent){
  #generate the posterior of beta distribution, this distribution is not related to the posterior distributions mentioned above. It is independent
  #the returns of this function will be a vector with dimension equal to gammaCurrent but with some NA values
  #this gamma is the gamma after the gamma selection, which has been incorporated the factor variable analysis, 04/10/10
  betaDim = length(gammaCurrent)
  xoxxoy = XOXXOY(xData, yData, omegaCurrent, gammaCurrent)
  XOX = xoxxoy[[1]]
  XOY = xoxxoy[[2]][,1]
  betaNew = array(NA, betaDim)
  meanBeta = solve(XOX) %*% XOY
  varianceBeta = solve(XOX)
  library(mvtnorm)
  betaNewTemp = rmvnorm(n = 1, mean = meanBeta, sigma = varianceBeta * (sigmaSquareCurrent))
  index = which(gammaCurrent == 1)
  for(i in 1:(length(index))){
    betaNew[index[i]] = betaNewTemp[i]
  }
  return (betaNew)
}

generateGamma = function(xData, yData, gammaCurrent, factorSummary, facrotName, omegaCurrent, sigmaSquareCurrent, lambdakappa, restriction = c(1, array(0, (numberOfVariable - 1)))){
#This function generate gamma vector
#factorSummary required, 05/07/10
    if (length(factorSummary) > 0){
        numberOfVariable = length(gammaCurrent)
        k = length(gammaCurrent)
        gammaNew = array(1, k)
        numVariable = k - sum(factorSummary)
        for(m in 1:numVariable){
            i = restriction[m]
            if (i == 0){
            gammaTemp1 = gammaCurrent
            gammaTemp2 = gammaCurrent
            gammaTemp1[m] = 1
            gammaTemp2[m] = 0
            prob1 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp1, sigmaSquareCurrent, lambdakappa)
            prob2 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp2, sigmaSquareCurrent, lambdakappa)
            gammaNew[m] = ifelse(runif(1) < 1 /(1 + exp(prob1 - prob2)), 0, 1)
            }
            if (i == 1){
            gammaNew[m] = 1
            }
            gammaCurrent[m] = gammaNew[m]
        }
        for(factorIter in 1:length(factorName)){
            i = restriction[factorIter]
            if (i == 0){
                gammaTemp1 = gammaCurrent
                gammaTemp2 = gammaCurrent
                const = ifelse(factorIter == 1, 0, sum(factorSummary[1:(factorIter - 1)]))
                gammaTemp1[numberOfVariable - sum(factorSummary) + 1:factorSummary[factorIter] + const] = 1
                gammaTemp2[numberOfVariable - sum(factorSummary) + 1:factorSummary[factorIter] + const] = 0
                prob1 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp1, sigmaSquareCurrent, lambdakappa)
                prob2 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp2, sigmaSquareCurrent, lambdakappa)
                gammaNew[numberOfVariable - sum(factorSummary) + 1:factorSummary[factorIter] + const] = ifelse(runif(1) < 1 /(1 + exp(prob1 - prob2)), 0, 1)
            }
            if (i == 1){
                const = ifelse(factorIter == 1, 0, sum(factorSummary[1:(factorIter - 1)]))
                gammaNew[numberOfVariable - sum(factorSummary) + 1:factorSummary[factorIter] + const] = 1
            }
            gammaCurrent[numberOfVariable - sum(factorSummary) + 1:factorSummary[factorIter] + const] = gammaNew[numberOfVariable - sum(factorSummary) + 1:factorSummary[factorIter] + const]
        }
    }
    if (length(factorSummary) == 0){
        numberOfVariable = length(gammaCurrent)
        nonFactorVariable = numberOfVariable
        gammaNew = array(NA, numberOfVariable)
        for(m in 1:nonFactorVariable){
            i = restriction[m]
            if (i == 0){
            gammaTemp1 = gammaCurrent
            gammaTemp2 = gammaCurrent
            gammaTemp1[m] = 1
            gammaTemp2[m] = 0
            prob1 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp1, sigmaSquareCurrent, lambdakappa)
            prob2 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp2, sigmaSquareCurrent, lambdakappa)
            gammaNew[m] = ifelse(runif(1) < 1 /(1 + exp(prob1 - prob2)), 0, 1)
            }
            if (i == 1){
            gammaNew[m] = 1
            }
            gammaCurrent[m] = gammaNew[m]
        }
    }
  return (gammaCurrent)
}


generateGammaFactor = function(xData, yData, gammaCurrent, omegaCurrent, sigmaSquareCurrent, lambdakappa, factorSummary, restriction = c(1, array(0, (numberOfVariable - 1)))){
#factor summary is a vector which includes number of levels of each factor name, numeric variable. added 04/10/10
  if (length(factorSummary) > 0){
    numberOfVariable = length(gammaCurrent)
    nonFactorVariable = numberOfVariable - sum(factorSummary)
    gammaNew = array(NA, numberOfVariable)
    for(m in 1:nonFactorVariable){
        i = restriction[m]
        if (i == 0){
          gammaTemp1 = gammaCurrent
          gammaTemp2 = gammaCurrent
          gammaTemp1[m] = 1
          gammaTemp2[m] = 0
          prob1 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp1, sigmaSquareCurrent, lambdakappa)
          prob2 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp2, sigmaSquareCurrent, lambdakappa)
          gammaNew[m] = ifelse(runif(1) < 1 /(1 + exp(prob1 - prob2)), 0, 1)
        }
        if (i == 1){
          gammaNew[m] = 1
        }
        gammaCurrent[m] = gammaNew[m]
    }
    indicator = nonFactorVariable + 1
    for(fact in 1:length(factorSummary)){
      i = restriction[fact]
      if (i == 0){
          gammaTemp1 = gammaCurrent
          gammaTemp2 = gammaCurrent
          gammaTemp1[indicator : (indicator + factorSummary[fact] - 1)] = 1
          gammaTemp2[indicator : (indicator + factorSummary[fact] - 1)] = 0
          prob1 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp1, sigmaSquareCurrent, lambdakappa)
          prob2 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp2, sigmaSquareCurrent, lambdakappa)
          gammaNew[indicator : (indicator + factorSummary[fact] - 1)] = ifelse(runif(1) < 1 /(1 + exp(prob1 - prob2)), rep(0, factorSummary[fact]), rep(1, factorSummary[fact]))
        }
        if (i == 1){
          gammaNew[indicator : (indicator + factorSummary[fact] - 1)] = 1
      }
      gammaCurrent[indicator : (indicator + factorSummary[fact] - 1)] = gammaNew[indicator : (indicator + factorSummary[fact] - 1)]
      indicator = indicator + factorSummary[factor]
    }
  }
  if (length(factorSummary) == 0){
    numberOfVariable = length(gammaCurrent)
    nonFactorVariable = numberOfVariable
    gammaNew = array(NA, numberOfVariable)
    for(m in 1:nonFactorVariable){
        i = restriction[m]
        if (i == 0){
          gammaTemp1 = gammaCurrent
          gammaTemp2 = gammaCurrent
          gammaTemp1[m] = 1
          gammaTemp2[m] = 0
          prob1 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp1, sigmaSquareCurrent, lambdakappa)
          prob2 = posteriorFunc(xData, yData, omegaCurrent, gammaTemp2, sigmaSquareCurrent, lambdakappa)
          gammaNew[m] = ifelse(runif(1) < 1 /(1 + exp(prob1 - prob2)), 0, 1)
        }
        if (i == 1){
          gammaNew[m] = 1
        }
        gammaCurrent[m] = gammaNew[m]
    }
  }
  return (gammaCurrent)
}


generateSigmaSquare = function(xData, yData, gammaCurrent, omegaCurrent){
  #generate the posterior sample of sigma square, it will be a chi-square distribution
  numberOfVariable = length(gammaCurrent)
  XOXXOY1 = XOXXOY(xData, yData, omegaCurrent, gammaCurrent)
  XOX = XOXXOY1[[1]]
  XOY = XOXXOY1[[2]][,1]
  YOY = t(yData) %*% omegaCurrent %*% yData
  betaConst = (.5*(YOY - t(XOY) %*% solve(XOX) %*% (XOY)))[1,1]
  alphaConst = (dim(omegaCurrent)[1] - sum(gammaCurrent) - 1)/2
  gammaVariable = rgamma(1, shape = alphaConst)
  sigmaSquareNew = betaConst/gammaVariable
  return (sigmaSquareNew)
}

generateP = function(xData, yData, treeData, gammaCurrent, sigmaSquareLambdaCurrent, sigmaSquareKappaCurrent, lambdaCurrent, kappaCurrent, lambdaUpperBound, kappaUpperBound) {
  numberOfVariable = length(gammaCurrent)
# P = 1, lambda model, P = 0, kappa model
  varOmega1 = solve(getCovarianceMatrixLambda(treeData, lambdaCurrent, lambdaUpperBound) + .0001 * diag(length(yData)))
  prob1Log = posteriorFunc(xData, yData, varOmega1, gammaCurrent, sigmaSquareLambdaCurrent, lambdaCurrent)
  varOmega2 = solve(getCovarianceMatrixKappa(treeData, kappaCurrent, kappaUpperBound) + .0001 * diag(length(yData)))
  prob2Log = posteriorFunc(xData, yData, varOmega2, gammaCurrent, sigmaSquareKappaCurrent, kappaCurrent)
  PNew = ifelse(runif(1) < 1 /(1 + exp(prob1Log - prob2Log)), 0, 1)
  return (PNew)
}

generateLambda = function(xData, yData, treeData, gammaCurrent, sigmaSquareCurrent, lambdaCurrent, delta = .15, lambdaUpperBound = 1.){
  numberOfVariable = length(gammaCurrent)
  varOmega1 = solve(getCovarianceMatrixLambda(treeData, lambdaCurrent, lambdaUpperBound) + .0001 * diag(length(yData)))
  prob1Log = posteriorFunc(xData, yData, varOmega1, gammaCurrent, sigmaSquareCurrent, lambdaCurrent)
  lambdaNew = rnorm(1, lambdaCurrent, delta)
  varOmega2 = solve(getCovarianceMatrixLambda(treeData, lambdaNew, lambdaUpperBound) + .0001 * diag(length(yData)))
  prob2Log = posteriorFunc(xData, yData, varOmega2, gammaCurrent, sigmaSquareCurrent, lambdaNew)
  lambdaNext = ifelse(runif(1) < exp(prob2Log - prob1Log), lambdaNew, lambdaCurrent)
  return (lambdaNext)
}


generateKappa = function(xData, yData, treeData, gammaCurrent, sigmaSquareCurrent, kappaCurrent, delta = .15, kappaUpperBound = 1.){
  numberOfVariable = length(gammaCurrent)
  varOmega1 = solve(getCovarianceMatrixKappa(treeData, kappaCurrent, kappaUpperBound) + .0001 * diag(length(yData)))
  prob1Log = posteriorFunc(xData, yData, varOmega1, gammaCurrent, sigmaSquareCurrent, kappaCurrent)
  kappaNew = rnorm(1, kappaCurrent, delta)
  varOmega2 = solve(getCovarianceMatrixKappa(treeData, kappaNew, kappaUpperBound) + .0001 * diag(length(yData)))
  prob2Log = posteriorFunc(xData, yData, varOmega2, gammaCurrent, sigmaSquareCurrent, kappaNew)
  kappaNext = ifelse(runif(1) < exp(prob2Log - prob1Log), kappaNew, kappaCurrent)
  return (kappaNext)
}


##########posterior draw part############
posteriorDrawFunc = function(xData, yData, currentValue, treeData, factorSummary, priors, varSelection = "random", lambdaValue = NA, kappaValue = NA, restriction = "no restriction", lambdaUpperBound = 1, kappaUpperBound = 1){
  #varSelection can be "lambda" or "kappa", "random"
  #restriction is a list of col names which includes all the variable you want in.
  gammaCurrent = currentValue[[1]]
  P = currentValue[[2]]
  lambdaCurrent = currentValue[[3]]
  kappaCurrent = currentValue[[4]]
  betaCurrent = currentValue[[5]]
  sigmaSquareLambdaCurrent = currentValue[[6]]
  sigmaSquareKappaCurrent = currentValue[[7]]
  numberOfVariable = length(betaCurrent)
  mu_lambda = priors[1]
  mu_kappa = priors[2]
  mu_sigmaLambda = priors[3]
  mu_sigmaKappa = priors[4]
  if (varSelection == "lambda"){
    pNew = 1
    if (is.na(lambdaValue)){
      lambdaNew = generateLambda(xData, yData, treeData, gammaCurrent, sigmaSquareLambdaCurrent, lambdaCurrent, lambdaUpperBound = lambdaUpperBound)
      kappaNew = rnorm(1, mu_kappa, .1)
    }
    if (!is.na(lambdaValue)){
      lambdaNew = log(lambdaValue/(lambdaUpperBound - lambdaValue))
      kappaNew = rnorm(1, mu_kappa, .1)
    }
    omegaCurrent = solve(getCovarianceMatrixLambda(treeData, lambdaNew, lambdaUpperBound))
    sigmaSquareLambdaNew = generateSigmaSquare(xData, yData, gammaCurrent, omegaCurrent)
    sigmaSquareKappaNew = exp(rnorm(1, mu_sigmaKappa, .05))
    lambdakappa = lambdaNew
    s.temp = sigmaSquareLambdaNew
  }
 
  if (varSelection == "kappa"){
    pNew = 0
    if (is.na(kappaValue)){
      lambdaNew = rnorm(1, mu_lambda, .1)
      kappaNew = generateKappa(xData, yData, treeData, gammaCurrent, sigmaSquareKappaCurrent, kappaCurrent, kappaUpperBound = kappaUpperBound)
    }
    if (!is.na(kappaValue)){
      kappaNew = log(kappaValue/(kappaUpperBound - kappaValue))
      lambdaNew = 0
    }
      omegaCurrent = solve(getCovarianceMatrixKappa(treeData, kappaNew, kappaUpperBound))
      sigmaSquareKappaNew = generateSigmaSquare(xData, yData, gammaCurrent, omegaCurrent)
      sigmaSquareLambdaNew = exp(rnorm(1, mu_sigmaLambda, .05))
      lambdakappa = kappaNew
      s.temp = sigmaSquareKappaNew
  }
 
  if (varSelection == "random"){
    pNew = generateP(xData, yData, treeData, gammaCurrent, sigmaSquareLambdaCurrent, sigmaSquareKappaCurrent, lambdaCurrent, kappaCurrent, lambdaUpperBound, kappaUpperBound)
    lambdaNew = generateLambda(xData, yData, treeData, gammaCurrent, sigmaSquareLambdaCurrent, lambdaCurrent)
    kappaNew = generateKappa(xData, yData, treeData, gammaCurrent, sigmaSquareKappaCurrent, kappaCurrent)
    if(pNew == 1){
       omegaCurrent = solve(getCovarianceMatrixLambda(treeData, lambdaNew, lambdaUpperBound))
       sigmaSquareLambdaNew = generateSigmaSquare(xData, yData, gammaCurrent, omegaCurrent)
       sigmaSquareKappaNew = exp(rnorm(1, mu_sigmaKappa, .05))
       lambdakappa = lambdaNew
       s.temp = sigmaSquareLambdaNew
    }
    if(pNew == 0){
      omegaCurrent = solve(getCovarianceMatrixKappa(treeData, kappaNew, kappaUpperBound))
      sigmaSquareKappaNew = generateSigmaSquare(xData, yData, gammaCurrent, omegaCurrent)
      sigmaSquareLambdaNew = exp(rnorm(1, mu_sigmaLambda, .05))
      lambdakappa = kappaNew
      s.temp = sigmaSquareKappaNew
    }
  }
 
  if (restriction[1] == "no restriction"){
    resCode = c(1, array(0, (numberOfVariable - 1)))
  }
 
  if (restriction[1] != "no restriction"){
    resCode = c(1, array(0, (numberOfVariable - 1)))
    res.length = length(restriction)
    for(i in 1:res.length){
       index = which(colnames(xData) == restriction[i])
       resCode[index] = 1
    }
  }
  gammaNew = generateGamma(xData, yData, gammaCurrent, factorSummary, factorName, omegaCurrent, s.temp, lambdakappa, resCode)
 
 
  if (pNew == 1){
    betaNew = generateBeta(xData, yData, omegaCurrent, gammaNew, sigmaSquareLambdaNew)
    betaNew[is.na(betaNew)] = 0
    m_temp = xData %*% betaNew
    v_temp = getCovarianceMatrixLambda(treeData, lambdaNew, lambdaUpperBound) * sigmaSquareLambdaNew
  }
 
  if (pNew == 0){
    betaNew = generateBeta(xData, yData, omegaCurrent, gammaNew, sigmaSquareKappaNew)
    betaNew[is.na(betaNew)] = 0
    m_temp = xData %*% betaNew
    v_temp = getCovarianceMatrixKappa(treeData, lambdaNew, kappaUpperBound) * sigmaSquareKappaNew
  }
  y_data_transform = array(NA, c(1, length(yData)))
  y_data_transform[1, ] = yData 
  lkhood = dmnorm(y_data_transform, m_temp[, 1] , v_temp, log = T)
  return (list(gammaNew, pNew, lambdaNew, kappaNew, betaNew, sigmaSquareLambdaNew, sigmaSquareKappaNew, lkhood))
}


posteriorSample = function(yDataAll, xDataAll, treeDataAll, factorSummary, nposterior = 30000, burnin = 0, varSelection = "random", lambdaValue = NA, kappaValue = NA, restriction = "no restriction", lambdaUpperBound = 1, kappaUpperBound = 1){
  numberOfVariable = dim(xDataAll)[2]
  N = nposterior
  gammaSample = array(NA, c(N, numberOfVariable))
  betaSample = array(NA, c(N, numberOfVariable))
  pSample = array(NA, N)
  lambdaSample = array(NA, N)
  kappaSample  = array(NA, N)
  sigmaSquareLambdaSample = array(NA, N)
  sigmaSquareKappaSample = array(NA, N)
  lkhoodSample = array(NA, N)
  treeIndex = array(NA, N)
  numOfTrees = length(treeDataAll)
  indi = runif(1, 0, 1)
  initialLambda = initialValueTest(xDataAll, yDataAll, treeDataAll[[1]], 1)
  initialKappa = initialValueTest(xDataAll, yDataAll, treeDataAll[[1]], 0, 3)
  priors = c(initialLambda[[4]], initialKappa[[4]], log(initialLambda[[2]]), log(initialKappa[[2]]))
  if(indi < .5){
    currentValue = list(initialLambda[[1]], initialLambda[[3]], priors[1], priors[2], initialLambda[[6]], exp(priors[3]), exp(priors[4]), initialLambda[[7]], 1)
  }
  if(indi >= .5){
    currentValue = list(initialKappa[[1]], initialKappa[[3]], priors[1], priors[2], initialKappa[[6]], exp(priors[3]), exp(priors[4]), initialKappa[[7]], 1)
  }
  #pb <- winProgressBar(title = "Bayesian Posterior Drawing", min = 0,
  #                  max = N, width = 500)
  quarter_check <- floor(N / 4)
  for(i in 1:N){
    index = floor(runif(1, 0, numOfTrees)) + 1
    treeData = treeDataAll[[index]]
    xData = xDataAll
    yData = yDataAll
    newValue = posteriorDrawFunc(xData, yData, currentValue, treeData, factorSummary, priors, varSelection = varSelection, lambdaValue = lambdaValue, kappaValue = kappaValue, restriction = restriction, lambdaUpperBound = lambdaUpperBound, kappaUpperBound = kappaUpperBound)
    gammaSample[i, ] = newValue[[1]]
    if(newValue[[2]] == 1){
      sigmaSquareLambdaSample[i] = newValue[[6]]
      sigmaSquareKappaSample[i] = NA
      pSample[i] = "Lambda Model"
      lambdaSample[i] = newValue[[3]]
      kappaSample[i] = NA
    }
    if(newValue[[2]] == 0){
      sigmaSquareLambdaSample[i] = NA
      sigmaSquareKappaSample[i] = newValue[[7]]
      pSample[i] = "Kappa Model"
      lambdaSample[i] = NA
      kappaSample[i] = newValue[[4]]
    }
    betaSample[i, ] = newValue[[5]]
    lkhoodSample[i] = newValue[[8]]
    treeIndex[i] = index
    currentValue = newValue
    #setWinProgressBar(pb, i, label=paste( round(i/N*100, 0),
    #                                    "% done"))
	if (floor((i - 1) / quarter_check) != floor((i) / quarter_check)){
		print(paste("regression finished ", floor((i) / quarter_check) * 25, "%", sep = ""))
	}
  }
  #close(pb)
  result = list(gammaSample, pSample, lambdaSample, kappaSample, betaSample, sigmaSquareLambdaSample, sigmaSquareKappaSample, lkhoodSample, treeIndex)
  return (result)
}


#this function can deal with the data cleaning~~
bayesianModel = function(yData, xData, treeDataAll, factorName = c(), missingList = c(), lambdaUpperBound = 1, kappaUpperBound = 1, currentValue = 0, nposterior = 30000, burnin = 0, thin = 30, varSelection = "random", lambdaValue = NA, kappaValue = NA, restriction = "no restriction", path = "default"){
#this is the wrapper function for the bayesian posterior draw, the input will be the basic model input, the output will be the csv files which include the posterior draw. In this function, we will not do any posterior analysis. This is the core part of the function
  #print (factorName)
  if (length(factorName) > 0){
     listData = xDataFactor(xData, factorName)
     xData = listData[[1]]
     factorSummary = listData[[2]]
  }
  if (length(factorName) == 0){
    xData = xData
    factorSummary = c()
  }
  dataAfterCleaning = missingDataCleaning(xData, yData, treeDataAll, missingList)
  xDataAll = dataAfterCleaning[[1]]
  yDataAll = dataAfterCleaning[[2]]
  treeDataAll = dataAfterCleaning[[3]]
  print ("pre-analysis begins...")
  if(!is.na(lambdaValue)){
    lambdaOriginal = lambdaValue
    if(lambdaValue < .001){
      lambdaValue = .001
    }
    if(lambdaValue > .999){
      lambdaValue = .999
    }
  }
  if(!is.na(kappaValue)){
    kappaOriginal = kappaValue
    if(kappaValue < .001){
      kappaValue = .001
    }
    if(kappaValue > .999){
      kappaValue = .999
    }
  }
  numberOfVariable = dim(xData)[2]
  print ("pre-analysis finished, Bayesian posterior draw begins...")
  posteriorSampleList = posteriorSample(yDataAll = yDataAll, xDataAll = xDataAll, treeDataAll = treeDataAll, factorSummary = factorSummary, nposterior, burnin = burnin, varSelection = varSelection, lambdaValue = lambdaValue, kappaValue = kappaValue, restriction = restriction, lambdaUpperBound = lambdaUpperBound, kappaUpperBound = kappaUpperBound)
  print ("Bayesian posterior draw finished, writing files...")
  npos = floor(nposterior / thin)
  start = 1 + floor(burnin / thin)
  gammaSample = posteriorSampleList[[1]][(start:npos) * thin, ]
  betaSample = posteriorSampleList[[5]][(start:npos) * thin, ]
  pSample = posteriorSampleList[[2]][(start:npos) * thin]
  if(!is.na(lambdaValue)){
    lambdaSample = lambdaOriginal
  }
  if(is.na(lambdaValue)){
    lambdaTemp = posteriorSampleList[[3]][(start:npos) * thin]
    lambdaSample = lambdaUpperBound * exp(lambdaTemp) / (1 + exp(lambdaTemp))
  }
  if(!is.na(kappaValue)){
    kappaSample = kappaOriginal
  }
  if(is.na(kappaValue)){
    kappaTemp = posteriorSampleList[[4]][(start:npos) * thin]
    kappaSample = kappaUpperBound * exp(kappaTemp) / (1 + exp(kappaTemp))
  }
  sigmaSquareLambdaSample = posteriorSampleList[[6]][(start:npos) * thin]
  sigmaSquareKappaSample = posteriorSampleList[[7]][(start:npos) * thin]
  lkhoodSample = posteriorSampleList[[8]][(start:npos) * thin]
  treeIndex = posteriorSampleList[[9]][(start:npos) * thin]
  numberOfVariable = dim(gammaSample)[2]
  colnames(betaSample) = paste("coef.", colnames(xData), sep = "")
  result = cbind(gammaSample, pSample, lambdaSample, kappaSample, betaSample, sigmaSquareLambdaSample, sigmaSquareKappaSample, lkhoodSample, treeIndex)
  colnames(result) = c(colnames(xData), "pSample", "lambdaSample", "kappaSample", colnames(betaSample), "sigmaSquareLambdaSample", "sigmaSquareKappaSample", "lkhoodSample", "treeIndex")
  data = cbind(xData, yData)
  colnames(data) = c(colnames(xData), "y")
  dataname = rownames(xData)
  if (path == "default"){
    filename1 = "result.csv"
    filename2 = "data.csv"
    filename3 = "dataname.csv"
    write.csv(result, "result.csv")
    write.csv(data, "data.csv")
    write.csv(dataname, "dataname.csv")
  }
  if (path != "default"){
    filename1 = paste(path, "result.csv")
    filename2 = paste(path, "data.csv")
    filename3 = paste(path, "dataname.csv")
    write.csv(result, filename1)
    write.csv(data, filename2) 
    write.csv(dataname, filename3)
  }
  textPosteriorSample = paste("posterior sample is written in the file", filename1)
  textData = paste("dataset is written in the file", filename2)
  print (textPosteriorSample)
  print (textData)
  print ("files writing completed...")
  return (c(filename1, filename2, filename3))
}

blm = function(formula, data, treeDataAll, factorName = c(), missingList = c(), currentValue = 0, nposterior = 30000, burnin = 0, thin = 30, varSelection = "random", lambdaUpperBound = 1., kappaUpperBound = 1., lambdaValue = NA, kappaValue = NA, restriction = "no restriction", path = "default"){
  N = dim(data)[1]
  initial = strsplit(formula, "~")[[1]]
  initial = gsub("(^ +)|( +$)", "", initial)
  y_index = initial[1]
  yData = array(NA, N)
  for(i in 1:N){
    yData[i] = as.numeric(as.character(data[, which(colnames(data) == y_index)][i]))
  }
  x_index = strsplit(initial[2], "\\+")[[1]]
  x_index = gsub("(^ +)|( +$)", "", x_index)
  xData = array(1, N)
  name = c("intercept")
  n = length(x_index)
  for(i in 1:n){
    if(length(grep("\\*", x_index[i])) > 0){
      x_index_temp = strsplit(x_index[i], "\\*")[[1]]
      x_index_temp = gsub("(^ +)|( +$)", "", x_index_temp)
      m = length(x_index_temp)
      temp = array(1, N)
      for(k in 1:m){
        temp_array = array(NA, N)
        for(mm in 1:N){
          temp_array[mm] = as.numeric(as.character(data[mm, which(colnames(data) == x_index_temp[k])]))
        }
        temp = temp * temp_array
      }
      xData = cbind(xData, temp)
      name = c(name, x_index[i])
    }
    if(length(grep("\\*", x_index[i])) == 0){
      temp_array = array(NA, N)
      for(m in 1:N){
        temp_array[m] = as.numeric(as.character(data[m, which(colnames(data) == x_index[i])]))
      }
      xData = cbind(xData, temp_array)
      name = c(name, x_index[i])
    }
  }
  fan = factorName
  factorName = c()
  for(i in 1:length(fan)){
    if(sum(name == fan[i]) > 0){
		factorName = c(factorName, fan[i])
	}
  }
  colnames(xData) = name
  rownames(xData) = as.character(data[, 1])
  data.list = preDataCleaning(yData, xData, treeDataAll, missingList)
  xData = data.list[[1]]
  yData = data.list[[2]]
  treeDataAll = data.list[[3]]
  missingList = data.list[[4]]
    
  filename = bayesianModel(yData, xData, treeDataAll, factorName = factorName, missingList = missingList, currentValue = currentValue, lambdaUpperBound = lambdaUpperBound, kappaUpperBound = kappaUpperBound, nposterior = nposterior, burnin = burnin, thin = thin, varSelection = varSelection, lambdaValue = lambdaValue, kappaValue = kappaValue, restriction = restriction, path = path)
  return (list(filename, treeDataAll, c(lambdaUpperBound, kappaUpperBound)))
}


preDataCleaning = function(yData, xData, treeDataAll, missingList){
	originalMissingList = missingList
	mm = getCovarianceMatrixLambda(treeDataAll[[1]], 0, 1.)
    name.tree = colnames(mm)
	index.x.nmissing = !is.na(rowMeans(xData))
	index.x.missing = which(is.na(rowMeans(xData)))
	if(sum(index.x.nmissing) < dim(xData)[1]){
		print ("The following species in X has missing data")
		print (rownames(xData)[is.na(rowMeans(xData))])
		print ("Deleted these species")
	}
	xData = xData[index.x.nmissing, ]
	yData = yData[index.x.nmissing]
	name.data = rownames(xData)
	y.missing = name.data[is.na(yData)]
      if(sum(is.na(yData)) > 0){
	   missingList = c(missingList, y.missing)
	   print ("The following species in Y is missing!")
	   print (y.missing)
	   print ("Added these species into missingList for prediction")
      }
	name.final = intersect(name.data, name.tree)
	name.delete = c()
	for(i in 1:length(name.tree)){
		if(sum(name.final == name.tree[i]) == 0){
			name.delete = c(name.delete, name.tree[i])
		}
	}
	if(length(name.delete) > 0){
		print ("The species in trees is not consistent with the species in data")
		print ("The common species are:")
		print (name.final)
		print ("All other species will be deleted")
	}
	missingList = intersect(missingList, name.final)
	if(length(missingList) > 0){
		print ("Those species are deleted from regression")
		print (missingList)
		print ("Those species are in missingList")
		print (originalMissingList)
		
	}
	if(length(missingList) == 0){
		print ("After the adjustment, we do not have any missing list")
	}
	
	xData = xData[name.final, ]
      index.name.final = array(0, length(name.data))
      for(i in 1:length(name.data)){
          if(sum(name.final == name.data[i]) > 0){
            index.name.final[i] = 1
          }
      }
      yData = yData[which(index.name.final == 1)]
	n = length(treeDataAll)
	if(length(name.delete) > 0){
		tr = vector("list", n)
		for(i in 1:n){
			tr[[i]] = drop.tip(treeDataAll[[i]], name.delete)
		}
	}
	else{
		
		tr = treeDataAll
	}
	treeDataAll = tr
	return (list(xData, yData, treeDataAll, missingList))
}

######begin with the model output part#####
###analysis part####
######posterior analysis########
##Find the model################
#df is the number of variables

##initialize the analysis, the file name required
initializeAnalysis = function(filename, lambdaUpperBound, kappaUpperBound){
  filename1 = filename[1]
  filename2 = filename[2]
  filename3 = filename[3]
  sampleO = read.table(filename1, sep = ",", header = T)
  sample = array(NA, dim(sampleO))
  nV = (dim(sampleO)[2] - 8) / 2
  for(i in 1:(dim(sampleO)[2])){
    if (i != (nV + 2)){
      sample[, i] = sampleO[, i]
    }
    if (i == (nV + 2)){
      for(mm in 1:dim(sampleO)[1]){
        sample[mm, i] = ifelse(as.character(sampleO[mm, i]) == "Lambda Model", 1, 0)
      }
    }
  }
  colnames(sample) = colnames(sampleO)
  dataO = read.table(filename2, sep = ",", header = T)
  data = array(NA, dim(dataO))
  for(i in 1:(dim(dataO)[2])){
    data[, i] = dataO[, i]
  }
  colnames(data) = colnames(dataO)
  name = read.table(filename3, sep = ",", head = T)
 
  rownames(data) = name[, 2]
  return (list(sample, data, c(lambdaUpperBound, kappaUpperBound)))
} 
 

initialSample = function(listOfSampleData){
#the input of this function will be the output of the previous function
  sample = listOfSampleData[[1]]
  data = listOfSampleData[[2]]
  lambdaUpperBound = listOfSampleData[[3]][1]
  kappaUpperBound = listOfSampleData[[3]][2]
  
  numberOfVariable = (dim(data)[2] - 2)
  gammaSample = sample[, 2:(1 + numberOfVariable)]
  pSample = sample[, 2 + numberOfVariable]
  lambdaTemp = sample[, 3 + numberOfVariable]
  if(prod(is.na(lambdaTemp)) == 0){
    if(min(lambdaTemp, na.rm = T) < .001){
      lambdaTemp[which(lambdaTemp <= .001)] = .001
    }
    if(max(lambdaTemp, na.rm = T) > .999){
      lambdaTemp[which(lambdaTemp >= .999)] = .999
    }  
  }
  lambdaSample = log(lambdaTemp) - log(lambdaUpperBound - lambdaTemp)
  kappaTemp = sample[, 4 + numberOfVariable]
  #if(prod(is.na(kappaTemp)) == 0){
  #  if(min(kappaTemp, na.rm = T) < .001){
  #    kappaTemp[which(kappaTemp <= .001)] = .001
  #  }
  #  if(max(kappaTemp, na.rm = T) > .999){
  #    kappaTemp[which(kappaTemp >= .999)] = .999
  #  }  
  #}
  kappaSample = log(kappaTemp) - log(kappaUpperBound - kappaTemp)
  betaSample = (sample[, (5 + numberOfVariable):(4 + 2 * numberOfVariable)])
  sigmaSquareLambdaSample = sample[, 5 + 2 * numberOfVariable]
  sigmaSquareKappaSample = sample[, 6 + 2 * numberOfVariable]
  lkhoodSample  = sample[, 7 + 2 * numberOfVariable]
  treeIndex  = sample[, 8 + 2 * numberOfVariable]
  xData = data[, 2:(1 + numberOfVariable)]
  yData = data[, (2 + numberOfVariable)]
  colnames(gammaSample) = colnames(betaSample) = colnames(sample)[2:(numberOfVariable + 1)]
  return (list(numberOfVariable, gammaSample, pSample, lambdaSample, kappaSample, betaSample, sigmaSquareLambdaSample, sigmaSquareKappaSample, lkhoodSample, xData, yData, treeIndex))
}

###
locateBestModel = function(gammaSample, betaSample){
  numberOfVariable = dim(gammaSample)[2]
  index = array(0, 2^numberOfVariable)
  len.gamma = dim(gammaSample)[1]
  gammaSampleIndicator = array(NA, len.gamma)
  binary = 2^((numberOfVariable - 1):0)
  modelCoefficientVec = modelCoefficientVarVec = array(NA, c(2^numberOfVariable, (numberOfVariable + 1)))
  for(i in 1:len.gamma){
    binaryNumber = sum(binary * gammaSample[i, ]) + 1
    gammaSampleIndicator[i] = binaryNumber
    index[binaryNumber] =  index[binaryNumber] + 1
  }
  indexSort = sort(index)
  modelSelectionVec = array(NA, c(2^numberOfVariable, (numberOfVariable + 1)))
  for(i in 1:(2^numberOfVariable)){
    modelIndicator = (i - 1)
    modelRep = array(NA, numberOfVariable)
    for(j in 1:numberOfVariable){
      remainder = modelIndicator %% 2
      modelRep[numberOfVariable + 1 - j] = remainder
      modelIndicator = (modelIndicator - remainder) / 2
    }
    modelSelectionVec[i, (1:numberOfVariable)] = modelRep
    modelSelectionVec[i, (numberOfVariable + 1)] =  modelCoefficientVec[i, (numberOfVariable + 1)] = modelCoefficientVarVec[i, (numberOfVariable + 1)] = index[i]/(dim(gammaSample)[1])  
    if ((length((betaSample[which(gammaSampleIndicator == i), 1]))) > 1){
      betaModel = colMeans(betaSample[which(gammaSampleIndicator == i), ])
      betaModelVar = apply(betaSample[which(gammaSampleIndicator == i), ], 2, var, na.rm = T)
      modelCoefficientVarVec[i, (1:numberOfVariable)] = betaModelVar
    }
   
    if ((length((betaSample[which(gammaSampleIndicator == i), 1]))) == 0){
      betaModel = array(NA, numberOfVariable)
      modelCoefficientVarVec[i, (1:numberOfVariable)] = array(NA, numberOfVariable)
    }

    if ((length((betaSample[which(gammaSampleIndicator == i), 1]))) == 1){
      betaModel = betaSample[which(gammaSampleIndicator == i), ]
      modelCoefficientVarVec[i, (1:numberOfVariable)] = array(NA, numberOfVariable)
    }
    modelCoefficientVec[i, (1:numberOfVariable)] = betaModel
  }
  selectionVariable = modelSelectionVec[order(modelSelectionVec[, numberOfVariable + 1], decreasing = T), ]
  estimateCoefficient = modelCoefficientVec[order(modelSelectionVec[, numberOfVariable + 1], decreasing = T), ]
  estimateVariance = modelCoefficientVarVec[order(modelSelectionVec[, numberOfVariable + 1], decreasing = T), ]
  indexPo = which(selectionVariable[, numberOfVariable + 1] * (dim(gammaSample)[1]) > 1)
  result = list(selectionVariable[indexPo, ], estimateCoefficient[indexPo, ], estimateVariance[indexPo, ])
  return (result)
}


###########main datainput wrapper#####
#currentValue

 
 
 
 
locateModelWrapper = function(gammaSample, betaSample, path = "default"){
  ####model-running####
  #select the model and output the selection as an excel file.
  model = locateBestModel(gammaSample, betaSample)
  ##model is the selection summary
  indicator = model[[1]]
  meanBeta = model[[2]]
  varBeta = model[[3]]
  if(length(dim(indicator)) == 2){
	colnames(indicator) = c(colnames(gammaSample), "prob")
	colnames(meanBeta) = c(colnames(gammaSample), "prob")
	colnames(varBeta) = c(colnames(gammaSample), "prob")
  }
  name = ifelse(path == "default", "modelIndicator.csv", paste(path, "modelIndicator.csv"))
  write.csv(indicator, name)
  model = list(indicator, meanBeta, varBeta)
  return (model)
}

#model = locateModelWrapper(gammaSample, betaSample)
#indicator = model[[1]]
#meanBeta = model[[2]]
#varBeta = model[[3]]



plotModelSelection = function(indicator){
  par(mfrow = c(1, 1))
  if(length(dim(indicator)) == 2){
	numberOfVariable = dim(indicator)[2] - 1
  }
  else{
	numberOfVariable = length(indicator) - 1
  }
  plot(c(0, 0), type = "n", xlim = c(0, 10), ylim = c(0, 13), ylab = "candidate models", xlab = ".", main = "Summary of Model Selection", axes = F)
  axis(1)
  axis(2, at = 1:13, labels = rep("", 13))
  len = 4.2/numberOfVariable
  if(length(dim(indicator)) == 2){
	den = indicator[, numberOfVariable + 1]
  }
  else{
	den = indicator[numberOfVariable + 1]
  }
  upperBound = min(length(den), 6)
  if(length(dim(indicator)) == 2){  
	  for(i in 1:upperBound){
		for(j in 1:numberOfVariable){
		  text((j)*len, (12 - 2*(i - 1)), ifelse(indicator[i, j] == 1, "I", "N"))
		}
		lines(c(5, 5 + den[i] * 5), c((12 - 2*(i - 1)), (12 - 2*(i - 1))), lwd = 12, col = "red")
		text(5 + den[i] * 5 + .5, 14 - 2 * i, floor(den[i]*1000)/1000)
	  }
  }
  else{
	  for(i in 1:upperBound){
		 text((1)*len, (12 - 2*(i - 1)), ifelse(indicator[i] == 1, "I", "N"))
		lines(c(5, 5 + den[i] * 5), c((12 - 2*(i - 1)), (12 - 2*(i - 1))), lwd = 12, col = "red")
		text(5 + den[i] * 5 + .5, 14 - 2 * i, floor(den[i]*1000)/1000)
	  }	
  }
 text(6, 1, "I, N means whether this coefficients is included or not. ")

}

ll = function(gammaSample){
  NN = dim(gammaSample)[1]
  numberOfVariable  = dim(gammaSample)[2]
  index = array(NA, NN)
  for(i in 1:NN){
    index[i] = sum(2^(0:(numberOfVariable - 1)) * gammaSample[i, ])
  }
  po = sum(2^(0:(numberOfVariable - 1)) * indicator[1, (1 : numberOfVariable)])
  indexFind = which(index == po)
  return (indexFind)
}

betaHistogram = function(betaSample, gammaSample, major = F){
  if(major == F){
    betaSampleNew = betaSample
  }
  if (major == T){
    lll = ll(gammaSample)
    betaSampleNew = betaSample[lll, which(indicator[1, 1:(dim(betaSample)[2])] == 1)]
  }
  betaSample = betaSampleNew
  NN = floor(sqrt(dim(betaSample)[2]))
  par(mfrow = c(NN, (NN + 1)))
  for(i in 1:dim(betaSample)[2]){
    if (sum(betaSample[, i], na.rm = T) != 0){
      hist(betaSample[, i], xlab = colnames(betaSample)[i], freq = F, main = paste("histogram for", colnames(betaSample)[i]), col = "blue")
    }
  }
}     

lambdaHistogram = function(pSample, lambdaSample, kappaSample, lambdaUpperBound, kappaUpperBound, major = F){
  if(major == F){
    lambdaTemp = lambdaSample
    kappaTemp = kappaSample
    pTemp = pSample
  }
  if (major == T){
    lll = ll(gammaSample)
    lambdaTemp = lambdaSample[lll]
    kappaTemp = kappaSample[lll]
    pTemp = pSample[lll]
  }
  lambdaSample = lambdaTemp
  kappaSample = kappaTemp
  pSample = pTemp
  par(mfrow = c(2, 1))
  if (sum(pSample) > 0){
    lambda = lambdaSample[which(pSample == 1)]
    hist(lambdaUpperBound * exp(lambda)/(1 + exp(lambda)), xlab = "lambda", freq = F, main = "histogram for lambda", col = "blue")
    text(lambdaUpperBound * max(exp(lambda)/(1 + exp(lambda))) - .1 , 1.5, paste("lambda model", sum(pSample)/(length(pSample))), col = "red")
  }
  if (sum(pSample) != length(lambdaSample)){
    kappa = kappaSample[which(pSample == 0)]
    hist(kappaUpperBound * exp(kappa)/(1 + exp(kappa)), xlab = "kappa", freq = F, main = "histogram for kappa", col = "blue")
    text(kappaUpperBound * max(exp(kappa)/(1 + exp(kappa))) - .1 , 1.5, paste("kappa model", 1 - sum(pSample)/(length(pSample))), col = "red")
  }
}

sigmaSquareHistogram = function(sigmaSquareSample, major = F){
  if(major == F){
    sigmaTemp = sigmaSquareSample
  }
  if (major == T){
    lll = ll(gammaSample)
    sigmaTemp = sigmaSquareSample[lll]
  }
  par(mfrow = c(1, 1))
  sigmaSquareSample = sigmaTemp
  hist(sigmaSquareSample, xlab = "sigmaSquare", freq = F, main = "histogram for sigmasquare", col = "blue")                                                                                         
}




plotCoef = function(coefMatrix, mainX){
  par(mfrow = c(1,1))
  names = rownames(coefMatrix)
  #print (names)
  #print ('----')
  len = dim(coefMatrix)[1]
  line1 = coefMatrix[, 1] + 1.96 * sqrt(coefMatrix[, 2])
  line2 = coefMatrix[, 1] - 1.96 * sqrt(coefMatrix[, 2])
  range = c(min(line2, na.rm = T) - 2, max(line1, na.rm = T))
  plot(c(0, 0), type = "n", xlim = range, ylim = c(0, 2 * len + 2), xlab = "coefficient", ylab = ".", main = mainX, axes = F)
  axis(1)
  for(i in 1:len){
	text(min(line2) - 1, 2 * len + 1 - 2*i, names[i])
    if (is.na(coefMatrix[i, 1]) == F){
      lines(c(line1[i], line2[i]), c(2 * len + 1 - 2*i, 2 * len + 1 - 2*i), col = "blue", lwd = 2)
      points(coefMatrix[i, 1], 2 * len + 1 - 2*i, pch = 19, cex = 2, col = "blue")
    }
  }
  lines(c(0, 0), c(0, 2 * len + 2), col = "red", lwd = 3)
}


summaryStat = function(betaSample, model){
  #plot the coefficient of the model
  par(mfrow = c(1, 1))
#difference the major model and the full model
  indicator = model[[1]]
  #print (indicator)
  #print ('-----')
  meanBeta = model[[2]]
  varBeta = model[[3]]
  if(length(dim(indicator)) == 2){
	  numberOfVariable = dim(indicator)[2] - 1
	  majorModelCoef = cbind(meanBeta[1, 1:numberOfVariable], varBeta[1, 1:numberOfVariable], (meanBeta[1, 1:numberOfVariable]/sqrt(varBeta[1, 1:numberOfVariable])))
	  colnames(majorModelCoef) = c("esti", "var", "tStat")
	  rownames(majorModelCoef) = colnames(betaSample)[1:numberOfVariable]
	#  print ("major Model Coef")
	#  print (majorModelCoef)
	  #plotCoef(majorModelCoef)
	  overallModelCoef = array(NA, c(numberOfVariable, 3))
	  for(i in 1:numberOfVariable){
		meanTemp = sum(meanBeta[, i] * meanBeta[, numberOfVariable + 1], na.rm = T)
		varTemp = sum(varBeta[, i] * varBeta[, numberOfVariable + 1], na.rm = T)
		tStat = meanTemp/sqrt(varTemp)
		overallModelCoef[i, ] = c(meanTemp, varTemp, tStat)
	  }
	  colnames(overallModelCoef) = c("esti", "var", "tStat")
	  rownames(overallModelCoef) = colnames(betaSample)[1:numberOfVariable]
	#  print ("overall Model Coef")
	#  print (overallModelCoef)
	  plotCoef(majorModelCoef, "coefficient for major model")
	  plotCoef(overallModelCoef, "coefficient for adjusted model")
  }
  else{
	  numberOfVariable = length(indicator) - 1
	  majorModelCoef = cbind(meanBeta[1:numberOfVariable], varBeta[1:numberOfVariable], (meanBeta[1:numberOfVariable]/sqrt(varBeta[1:numberOfVariable])))
	#  colnames(majorModelCoef) = c("esti", "var", "tStat")
	  rownames(majorModelCoef) = colnames(betaSample)[1:numberOfVariable]
	#  print ("major Model Coef")
	#  print (majorModelCoef)
	  #plotCoef(majorModelCoef)
	  overallModelCoef = array(NA, c(numberOfVariable, 3))
	  for(i in 1:numberOfVariable){
		meanTemp = sum(meanBeta[i] * meanBeta[numberOfVariable + 1], na.rm = T)
		varTemp = sum(varBeta[i] * varBeta[numberOfVariable + 1], na.rm = T)
		tStat = meanTemp/sqrt(varTemp)
		overallModelCoef[i, ] = c(meanTemp, varTemp, tStat)
	  }
	  colnames(overallModelCoef) = c("esti", "var", "tStat")
	  rownames(overallModelCoef) = colnames(betaSample)[1:numberOfVariable]
	#  print ("overall Model Coef")
	#  print (overallModelCoef)
	  plotCoef(majorModelCoef, "coefficient for major model")
	  plotCoef(overallModelCoef, "coefficient for adjusted model")  
  }
}



plotWholeCoef = function(betaSample, gammaSample){
  #print the overall model, which is the overall mean and variance of all posterior samples
  par(mfrow = c(1, 1))
  names = colnames(gammaSample)
  len = dim(gammaSample)[2]
  numberOfVariable = dim(gammaSample)[2]
  coefMatrix = array(NA, c(numberOfVariable, 2))
  for(i in 1:numberOfVariable){
    coefMatrix[i, 1] = mean(betaSample[, i], na.rm = T)
    coefMatrix[i, 2] = var(betaSample[, i], na.rm = T)
  }
  line1 = coefMatrix[, 1] + 1.96 * sqrt(coefMatrix[, 2])
  line2 = coefMatrix[, 1] - 1.96 * sqrt(coefMatrix[, 2])
  range = c(min(line2, na.rm = T) - 2, max(line1, na.rm = T))
  plot(c(0, 0), type = "n", xlim = range, ylim = c(0, 2 * len + 2), xlab = "coefficients", ylab = "", main = "overall model coefficients", axes = F)
  axis(1)
  for(i in 1:len){
    text(min(line2) - 1, 2 * len + 1 - 2*i, names[i])
    lines(c(line1[i], line2[i]), c(2 * len + 1 - 2*i, 2 * len + 1 - 2*i), col = "blue", lwd = 2)
    points(coefMatrix[i, 1], 2 * len + 1 - 2*i, pch = 19, cex = 1.5, col = "blue")
  }
  lines(c(0, 0), c(0, 2 * len + 2), col = "red", lwd = 3)
}

analysis = function(bmselection, path = "default"){
#analysis wrapper function, the output will be a pdf file
  filename = bmselection[[1]]
  lambdaUpperBound = bmselection[[3]][1]
  kappaUpperBound = bmselection[[3]][2]  
  listOfSampleData = initializeAnalysis(filename, lambdaUpperBound, kappaUpperBound)
  posteriorSampleList = initialSample(listOfSampleData)
  print ("initialize analysis and get posterior sample ends, analysis begins...")
  numberOfVariable = posteriorSampleList[[1]]
  gammaSample = posteriorSampleList[[2]]
  pSample = posteriorSampleList[[3]]
  lambdaSample = posteriorSampleList[[4]]
  kappaSample = posteriorSampleList[[5]]
  betaSample = posteriorSampleList[[6]]
  sigmaSquareLambdaSample = posteriorSampleList[[7]]
  sigmaSquareKappaSample = posteriorSampleList[[8]]
  lkhoodSample = posteriorSampleList[[9]]
  xData = posteriorSampleList[[10]]
  yData = posteriorSampleList[[11]]
  model = locateModelWrapper(gammaSample, betaSample, path)
  nameTemp = "modelIndicator.csv"
  nameModelIndicator = ifelse(path == "default", nameTemp, paste(path, nameTemp))
  textModelIndicator = paste("model indicator written in file", nameModelIndicator)
  print (textModelIndicator)
  indicator = model[[1]]
  meanBeta = model[[2]]
  varBeta = model[[3]]
  print ("initial analysis completed, plot the results...")
  ####Now we can plot the model####
  filenameOutput = ifelse(path == "default", "outputBayesianModel.pdf", paste(path, "outputBayesianModel.pdf"))
  pdf(file = filenameOutput)
  plotModelSelection(indicator)
  betaHistogram(betaSample, gammaSample, major = F)
  lambdaHistogram(pSample, lambdaSample, kappaSample, lambdaUpperBound, kappaUpperBound, major = F)
  if(sum(!is.na(sigmaSquareLambdaSample)) > 0){
    sigmaSquareHistogram(sigmaSquareLambdaSample, major = F)
  }
  if(sum(!is.na(sigmaSquareKappaSample)) > 0){
    sigmaSquareHistogram(sigmaSquareKappaSample, major = F)
  }
  summaryStat(betaSample, model)
  plotWholeCoef(betaSample, gammaSample)
  dev.off()
  aa = paste("plot ploted in the file", filenameOutput)
  print (aa)
  print ("plot completed, analysis finished...")
}

####begins the model checking part####

predictiveDraw = function(xDataAll, yDataAll, treeDataAll, gammaSample, betaSample, pSample, lambdaSample, kappaSample, sigmaSquareLambdaSample, sigmaSquareKappaSample,treeIndexSample, lambdaUpperBound = 1, kappaUpperBound = 1){
  #use the posterior sample in the part 1 to get predictive draw, will be used for model checking
  N = length(lambdaSample)
  numberOfVariable = dim(betaSample)[2]
  predictiveSample = array(NA, c(N, length(yDataAll)))
  betaSampleNaToZero = betaSample
  betaSampleNaToZero[is.na(betaSampleNaToZero)] = 0
  for(i in 1:N){
    xData = xDataAll
    yData = yDataAll
    treeData = treeDataAll[[treeIndexSample[i]]]
    meanY = xData %*% betaSampleNaToZero[i, ]
    if (pSample[i] == 1){
      MM = getCovarianceMatrixLambda(treeData, lambdaSample[i], lambdaUpperBound)
      sigmaSquareSample = sigmaSquareLambdaSample[i]
    }

    if (pSample[i] == 0){
      MM = getCovarianceMatrixKappa(treeData, kappaSample[i], kappaUpperBound)
      sigmaSquareSample = sigmaSquareKappaSample[i]
    }

    MM = MM * sigmaSquareSample
    yData = rmvnorm(n = 1, mean = meanY[,1], sigma = MM)[1, ]
    predictiveSample[i, ] = yData
  }
  return (predictiveSample)
}


plotModelChecking = function(yData, predictiveSample){
  #plot the model checking result based on var, mean, min, max function
  data = array(NA, c(4, dim(predictiveSample)[1]))
  data[1, ] = apply(predictiveSample, 1, var)
  data[2, ] = apply(predictiveSample, 1, mean)
  data[3, ] = apply(predictiveSample, 1, min)
  data[4, ] = apply(predictiveSample, 1, max)
  dataTrue = c(var(yData, na.rm = T), mean(yData, na.rm = T), min(yData, na.rm = T), max(yData, na.rm = T))
  histname = c("model checking for var", "model checking for mean", "model checking for min", "model checking for max")
  for(i in 1:4){
    hist(data[i, ], xlab = "posterior sample", main = histname[i], col = "blue")
    lines(c(dataTrue[i], dataTrue[i]), c(0, 1000), col = "red", lwd = 2)
  }
}


modelChecking = function(bmselection, missingList, path = "default"){
  #wrapper function for model checking
  filename = bmselection[[1]]
  treeDataAll = bmselection[[2]]
  lambdaUpperBound = bmselection[[3]][1]
  kappaUpperBound = bmselection[[3]][2]
  listOfSampleData = initializeAnalysis(filename, lambdaUpperBound, kappaUpperBound)
  posteriorSampleList = initialSample(listOfSampleData)
  numberOfVariable = posteriorSampleList[[1]]
  gammaSample = posteriorSampleList[[2]]
  sigmaSquareLambdaSample = posteriorSampleList[[7]]
  sigmaSquareKappaSample = posteriorSampleList[[8]]
  pSample = posteriorSampleList[[3]]
  lambdaSample = posteriorSampleList[[4]]
  kappaSample = posteriorSampleList[[5]]
  betaSample = posteriorSampleList[[6]]
  lkhoodSample = posteriorSampleList[[9]]
  treeIndexSample = posteriorSampleList[[12]]
  xData = posteriorSampleList[[10]]
  yData = posteriorSampleList[[11]]
  print ("initialize analysis, predictive draw begins...")
  missingList = intersect(missingList, rownames(xData))
  aa = missingDataCleaning(xData, yData, treeDataAll, missingList)
  xDataAll = aa[[1]]
  yDataAll = aa[[2]]
  treeDataAll = aa[[3]]
  predictiveSample = predictiveDraw(xDataAll, yDataAll, treeDataAll, gammaSample, betaSample, pSample, lambdaSample, kappaSample, sigmaSquareLambdaSample, sigmaSquareKappaSample, treeIndexSample)
  print ("predictive draw finished, model checking begins...")
  filenameOutput = ifelse(path == "default", "modelChecking.pdf", paste(path, "modelChecking.pdf"))
  pdf(file = filenameOutput)
  plotModelChecking(yData, predictiveSample)
  dev.off()
  textModelChecking = paste("model checking result finished, written in file", filenameOutput)
  print (textModelChecking)
}


#####end of model checking part######



predictiveDraw2 = function(xDataAllWithMissing, yDataAllWithMissing, treeDataAll, missingList, gammaSample, betaSample, pSample, lambdaSample, kappaSample, sigmaSquareLambdaSample, sigmaSquareKappaSample, treeIndexSample, nameData, lambdaUpperBound = 1., kappaUpperBound = 1.){
  #use the posterior sample in the part 1 to get predictive draw, will be used for model checking
  N = length(lambdaSample)
  numberOfVariable = dim(betaSample)[2]
  predictiveSample = array(NA, c(N, length(missingList)))
  betaSampleNaToZero = betaSample
  betaSampleNaToZero[is.na(betaSampleNaToZero)] = 0
  nameTree =  rownames(getCovarianceMatrixLambda(treeDataAll[[1]], .5))
  numberOfMissingData = length(missingList)
  for(i in 1:N){
    xData = xDataAllWithMissing
    yData = yDataAllWithMissing
    treeData = treeDataAll[[treeIndexSample[i]]]
   
    meanY = xData %*% betaSampleNaToZero[i, ]

    if (pSample[i] == 1){
      MM = getCovarianceMatrixLambda(treeData, lambdaSample[i], lambdaUpperBound)
      sigmaSquareSample = sigmaSquareLambdaSample[i]
    }

    if (pSample[i] == 0){
      MM = getCovarianceMatrixKappa(treeData, kappaSample[i], kappaUpperBound)
      sigmaSquareSample = sigmaSquareKappaSample[i]
    }
    MM = MM * sigmaSquareSample
    nameTemp = colnames(MM)
    indexMissing = c()
    indexNonMissing = c()
    for(iterMissing in 1:length(nameTemp)){
        if(length(which(missingList == nameTemp[iterMissing])) == 0){
            indexNonMissing = c(indexNonMissing, iterMissing)
        }
        if(length(which(missingList == nameTemp[iterMissing])) > 0){
            indexMissing = c(indexMissing, iterMissing)
        }
    }
    MM1 = MM[indexMissing, indexMissing]
    MM2 = MM[indexNonMissing, indexNonMissing]
    MM12 = MM[indexMissing, indexNonMissing]
    mean1 = meanY[indexMissing]
	#check the name match
	#print (c(betaSampleNaToZero[i, ], mean1, xData[indexMissing, ]))
    mean2 = meanY[indexNonMissing]
    munew = mean1 + MM12 %*% solve(MM2) %*% (yData[indexNonMissing] - mean2)
    if(length(indexMissing) == 1){
      sigmanew = MM1 - MM12 %*% solve(MM2) %*% (MM12)
    }
    if(length(indexMissing) > 1){
      sigmanew = MM1 - MM12 %*% solve(MM2) %*% t(MM12)
    }

    yData = rmvnorm(n = 1, mean = munew, sigma = sigmanew)[1, ]
    NN = length(missingList)
    predictiveSample[i, ] = yData
	#print (yData)
  }
  colnames(predictiveSample) = missingList
  return (predictiveSample)
}
predictiveUnknown = function(filename, treeDataAll, missingList, lambdaUpperBound, kappaUpperBound, path = default){
  #wrapper function for model checking
  listOfSampleData = initializeAnalysis(filename, lambdaUpperBound, kappaUpperBound)
  posteriorSampleList = initialSample(listOfSampleData)
  numberOfVariable = posteriorSampleList[[1]]
  gammaSample = posteriorSampleList[[2]]
  sigmaSquareLambdaSample = posteriorSampleList[[7]]
  sigmaSquareKappaSample = posteriorSampleList[[8]]
  pSample = posteriorSampleList[[3]]
  lambdaSample = posteriorSampleList[[4]]
  kappaSample = posteriorSampleList[[5]]
  betaSample = posteriorSampleList[[6]]
  lkhoodSample = posteriorSampleList[[9]]
  treeIndexSample = posteriorSampleList[[12]]
  xData = posteriorSampleList[[10]]
  yData = posteriorSampleList[[11]]
  print ("initialize analysis, predictive draw begins...")
  #aa = missingDataCleaning(xData, yData, treeDataAll, missingList = c())
  xDataAllWithMissing = array(NA, dim(xData))
  yDataAllWithMissing = array(NA, length(yData))
  name.local <- rownames(getCovarianceMatrixLambda(treeDataAll[[1]], .5))
  name.ori <- rownames(xData)
  for(m.iter in 1:dim(xData)[1]){
    index.temp <- which(name.ori == name.local[m.iter])
    xDataAllWithMissing[m.iter, ] <- xData[index.temp, ]
    yDataAllWithMissing[m.iter] <- yData[index.temp]
  }
 
  missingListOrigin = intersect(rownames(xData), missingList)
  missingList = unique(c(missingListOrigin, rownames(xData)[which(is.na(yData))]))
  #treeDataAll = aa[[3]]
  nameData = rownames(xData)
  predictiveSample = predictiveDraw2(xDataAllWithMissing, yDataAllWithMissing, treeDataAll, missingList, gammaSample, betaSample, pSample, lambdaSample, kappaSample, sigmaSquareLambdaSample, sigmaSquareKappaSample, treeIndexSample, nameData)
  print ("predictive draw end, begin for missing data analysis")
  index <- c()
  for(k in 1:length(missingList)){
	if(sum(missingListOrigin == missingList[k]) > 0){
		index <- c(index, k)
	}
  }
  predictiveSampleFinal <- as.matrix(predictiveSample[, index])
  colnames(predictiveSampleFinal) <- missingList[index]
  return (list(predictiveSampleFinal, xData, yData))
}



predictiveSummary = function(predictiveSampleWrapper, missingList, path){
  #give the predictive sample summary of the predictiveSamples
  predictiveSample <- predictiveSampleWrapper[[1]]
  xData <- predictiveSampleWrapper[[2]]
  yData <- predictiveSampleWrapper[[3]]
  n = length(missingList)
  predictive = array(NA, c(n,8))
  name.temp <- rownames(xData)
  if(path == "default"){
    filename <- "histogram_predictive_sample.pdf"
  }
  if(path != "default"){
    filename <- paste(path, "histogram_predictive_sample.pdf", sep = "")
  }
  pdf(filename)
  for(i in 1:n){
    index = which(colnames(predictiveSample) == missingList[i])
    hist(predictiveSample[, i], main = missingList[i], xlab = "predictive samples", col = "blue")
    index.name.check <- which(name.temp == missingList[i])
    if(!is.na(yData[index.name.check])){
        lines(c(yData[index.name.check], yData[index.name.check]), c(0, 100000000), col = "red", lwd = 2)
    }
   
    sample = array(predictiveSample[ ,index])
    sampleSort = sort(sample)
    indexSample = floor(c(0.025, 0.25, 0.75, 0.975) *  length(sample)) + 1
    predictive[i, ] = c(min(sampleSort), sampleSort[indexSample[1]], sampleSort[indexSample[2]], median(sampleSort), mean(sampleSort), sampleSort[indexSample[3]], sampleSort[indexSample[4]], max(sampleSort))
  }
  dev.off()
  colnames(predictive) = c("min", "2.5%q", "25%q", "median", "mean", "75%q", "97.5%q", "max")
  rownames(predictive) = missingList
  print ("missing data analysis end, begin printing out the result")
  return (predictive)
}

predict = function(bmselection, missingList, path = "default"){
#wrapper function for prediction
  filename = bmselection[[1]]
  treeDataAll = bmselection[[2]]
  lambdaUpperBound <- bmselection[[3]][1]
  kappaUpperBound <- bmselection[[3]][2]
  
  name.tree = rownames(getCovarianceMatrixLambda(treeDataAll[[1]], .5))
  missingList = intersect(missingList, name.tree)
  predictiveSampleWrapper = predictiveUnknown(filename, treeDataAll, missingList, lambdaUpperBound, kappaUpperBound)
  outputfilename <- paste(path, "predictions.csv", sep = "")
  write.csv(predictiveSampleWrapper[[1]], outputfilename)
  pre = predictiveSummary(predictiveSampleWrapper, missingList, path)
  #print (pre)
  return (pre)

}

