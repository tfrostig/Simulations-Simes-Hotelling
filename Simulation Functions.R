### Libraries 
### Code taken from the following stack post 
### stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them

list.of.packages <- c("doParallel", "MASS", "highmean", 
                      "highD2pop", "doRNG", "dplyr", "here", 
                      "ggplot2", "stringr", "reshape2", "fda.usc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, library, character.only = TRUE)

### Installing simesHotelling 
if ( !('simesHotelling' %in% installed.packages()[,"Package"]) ) {
  if (  !('devtools' %in% installed.packages()[,"Package"]) ) {
    install.packages('devtools')
  }
  devtools::install_github('tfrostig/simesHotelling', build_vignettes = FALSE)
}

# Data Generation  --------------------------------------------------------


## Create a list of X and Y same covariance matrix 
## Input: n - number of observations, p - number of dimensions, FHR - number of false hypothesis or rate 
##        delta - size of diffence in means (could be vector or scalar) cov.mat - covariance matrix 
## Output: list with matrix X and Y multivariate norm according to the specified requirements 

createData <-  function (non.null.prec, euclid.norm = 2, n.x, n.y = n.x, cov.mat.x, 
            cov.mat.y = cov.mat.y, group.num = 5, signal.type = "Unif", 
            noise.type = "none") {
    p <- ncol(cov.mat.x)
    signal.vec <- CreateSignal(round(non.null.prec * p), euclid.norm, 
                               group.num, signal.type)
    mu.vec <- rep(0, p)
    mu.vec[sample(p, length(signal.vec))] <- signal.vec
    X <- mvrnorm(n.x, mu.vec, cov.mat.x)
    Y <- mvrnorm(n.y, rep(0, p), cov.mat.y)
    if (noise.type == "exp") {
      X <- X + replicate(ncol(X), rexp(nrow(X)))
      Y <- Y + replicate(ncol(X), rexp(nrow(X)))
    }
    if (noise.type == "laplace") {
      X <- X + replicate(ncol(X), rlaplace(nrow(X), 0, sqrt(0.5)))
      Y <- Y + replicate(ncol(X), rlaplace(nrow(X), 0, sqrt(0.5)))
    }
    if (noise.type == "none") {
      X <- X
      Y <- Y
    }
    if (noise.type == "normal") {
      X <- X + replicate(ncol(X), rnorm(nrow(X), 0, 1))
      Y <- Y + replicate(ncol(X), rnorm(nrow(X), 0, 1))
    }
    return(list(X = X, Y = Y))
}



## Creating singal so that signal is of specific norm
CreateSignal <- function(non.null.num, euclid.norm, group.num, signal.type = 'Unif') {
  if (non.null.num %% group.num != 0) {
    group.num <- non.null.num
  }
  group.size  <- non.null.num / group.num
  effect.size <- sqrt(euclid.norm * group.num^2 / (sum((1:group.num)^2 * group.size)))
  if (signal.type == 'Unif') {
    mu.vec <- rep(seq(effect.size / group.num, effect.size, length.out = group.num), group.size)
  }
  return(mu.vec)
}


## Covariance creating 

## Create covariance matrix with exponential decline 
## Input: familly size - dimension size of matrix, first correlation near diagonal
##        variance on diagonal - v.d 
## Output: time series covariance matrix 

corr2 <- function(fam.size , cor, v.d = 1) {
  corr2 <- matrix(0 , ncol = fam.size , nrow = fam.size)
  for (i in 1:(fam.size - 1)) { 
    corr2[i,((i + 1):fam.size)] <- cor ^ (1:(fam.size - i))
  }
  diag(corr2) <- v.d
  corr2[lower.tri(corr2)] <- t(corr2)[lower.tri(corr2)]
  return(corr2)
}


## Full Block 

corr.full <- function(fam.size, rho, v.d = 1){ 
  cov.mat <- matrix(rho, ncol = fam.size, nrow = fam.size)
  diag(cov.mat) <- v.d 
  return(cov.mat)
}

## Covariance matrix with correlation that decrease as function of |i-j|
## Input: familly size - dimensions required, v.d - variance on diagonal 
## Output: Required matrix 


corr1 <- function(fam.size, v.d = 1) {
  corr1 <- matrix(0, ncol = fam.size , nrow = fam.size)
  for (i in 1:fam.size) {
    for (j in i:fam.size) {
      corr1[i,j] <-  (fam.size - j + i) / fam.size 
    }
  }
  diag(corr1) <- v.d
  corr1[lower.tri(corr1)] <- t(corr1)[lower.tri(corr1)]
  return(corr1)
}

## Covariance matrix with correlation that decrease as function of |i-j|
## Input: familly size - dimensions required, roMean of covariance in covariance matrix, 
##        alpha - alpha parameter for beta distribution, lim - mistake allowed for roMean
## Output: Required matrix 

covMat_beta_maker <- function(famSize , roMean , alpha = 10^-6  , lim = 0.05)
{
  covMat_beta <- matrix(-1 , ncol = famSize , nrow = famSize)
  i <- 0
  while ((!is.positive.semi.definite(covMat_beta)) 
         || (mean(covMat_beta[upper.tri(covMat_beta)]) < (roMean - lim))
         || (mean(covMat_beta[upper.tri(covMat_beta)]) > (roMean + lim))) {
    betaMean <- (sqrt(roMean)) 
    beta <- ((1 - betaMean) * alpha) / betaMean 
    vec <- rbeta(famSize , alpha , beta) 
    covMat_beta <- vec %*% t(vec)
    ## Can Be changed , changed to 0.95 from 0.99 
    covMat_beta[which(covMat_beta == 1)] <- 0.95
    diag(covMat_beta) <- 1 
    i <- i + 1
    print(i)
  }
  return(covMat_beta)
}

## Creating Block Covariance Matrix 
## Input: Family size - dimension required, correlation = correlation inside each block, Var - var on diagonal 
## output: required matrix 

corrBlock <- function(famSize , blockSize , correlation = 0.8 , var = 1) {
  if ((famSize < blockSize) | (((famSize / blockSize) %% 1) != 0)) {
    print('Block Size Must Be Smaller Than Familly Size Or The Devision Of Them Isnt a Whole Number')
  }
  corrBlock <-  matrix(0,ncol = famSize , nrow = famSize)
  for (i in 1:(famSize/blockSize)) {
    startPos <- 1 + (i - 1) * blockSize
    corrBlock[startPos:(startPos + blockSize -1),startPos:(startPos + blockSize - 1)] <- correlation 
  }
  diag(corrBlock) <- var
  return(corrBlock)
}

## Covariance model by Cai model 7 
## Input: fam.size - number of dimensions 
## Output: Requested matrix 

CovMakerModel7 <- function(fam.size) { 
  D <- matrix(0, ncol = fam.size, nrow = fam.size)
  diag(D) <- runif(fam.size , 1, 3)
  Sigma <- outer(1:fam.size, 1:fam.size , FUN="-")
  diag(Sigma) <- 1 
  Sigma <- ((abs(Sigma)^(-5)) / 2)  
  diag(Sigma) <- 1
  cov.mat <- sqrt(D) %*% Sigma %*% sqrt(D)
  return(cov.mat)
}

## Covariance model by Cai model 4 
## Input: fam.size - nubmer of dimensions, corr - number of covariance 
## Output: Requested matrix 

CovMakerModel4 <- function(fam.size , corr = 0.6) { 
  Omega <- 0.6^abs(outer(1:fam.size, 1:fam.size, FUN = "-"))
  D <- matrix(0, ncol = fam.size, nrow = fam.size)
  diag(D) <- runif(fam.size , 1, 3)
  return(sqrt(D) %*% solve(Omega) %*% sqrt(D))
}



# Simulation Maker  -----------------------------------------------------

### Iteartion for certain core 
### Input: covariance matrix, number of false hypothesis (or percentage), observations, differnce in mean vector, 
###        number of iterations, functions 

simulation.function <- function(cov.mat , delta, obs, fam.size, iter.num, multi.fun, func.num) {
  p <- ncol(cov.mat) 
  pval.mat <- data.frame(matrix(NA , nrow = iter.num, ncol = func.num))
  for (i in 1:iter.num) { 
    dat.list <- createData(n = obs, p = p, delta = delta, cov.mat = cov.mat)
    pval.mat[i,]  <- multi.fun(X = dat.list$X , Y = dat.list$Y)
  }
  return(pval.mat)
}

## MU creator 
## Creating delta vector according to type ('Ascending' , 'Descending' and 'Uniform')
## Input: non.null - number of false hypothesis, L - largest effect, group.num - number of groups ,type - type of effect 
## Output: Mu vector 

MU.Creator <- function(non.null, L, group.num = 5, type = 'Descend') {
  effect   <- seq(L / group.num, L, length.out = group.num)
  base.num <- round(non.null / sum(1:group.num)) 
  if (type == 'Ascend') { 
    mu.vec <-  unlist(mapply(rep, times = base.num * 1:group.num, x = effect))
    return(mu.vec)
  }
  if (type == 'Descend') {
    mu.vec <- unlist(mapply(rep, times = rev(base.num * 1:group.num), x = effect))
    return(mu.vec)
  }
  if (type == 'Unif') { 
    group.size <- rep(round(non.null / group.num), group.num)
    mu.vec <- as.vector(unlist(mapply(rep, times = group.size, x = effect)))
    return(mu.vec)
  } else { 
    stop('Not known type')
  }
}



### Mean Delta Maker - Samples multiple time mu vectors and return the average value of delta for them  
### Input: mu vector, covariance matrix, number of repeatition 
### Output: Delta 

DeltaMaker <- function(mu.vec, cov.mat, rep.numeric = 250) {
  p <- nrow(cov.mat) 
  k <- length(mu.vec)
  temp.delta <- rep(NA, rep.numeric)
  sov.cov <- solve(cov.mat)
  for (i in 1:rep.numeric) { 
    temp.mu <- rep(0, p)
    temp.mu[sample.int(p, k)] <- mu.vec 
    temp.delta[i] <- t(temp.mu) %*% sov.cov %*% temp.mu
  }
  return(mean(temp.delta))
}

### Finding The power 
PowerHoteling <- function(non.null.prec, cov.mat, delta = 2.5, group.num = 5, type = 'Unif', Mahalnobis = F) { 
  p <- ncol(cov.mat) 
  non.null <- round(non.null.prec * p)
  poss.L <- seq(0.01, 2, 0.001)
  mu.mat <- t(mapply(MU.Creator, L = poss.L, MoreArgs = list('non.null' = non.null, 
                                                             'type' = type,
                                                             'group.num' = group.num)))
  ## removing duplicates 
  mu.mat <- mu.mat[!duplicated(mu.mat), ]
  if (Mahalnobis == T) { 
    avg.delta <- apply(mu.mat, 1, DeltaMaker, cov.mat = cov.mat)
  }
  if (Mahalnobis == F) { 
    avg.delta <- apply(mu.mat, 1, function(x) sum(x^2))
  }
  return(c('L' = poss.L[which.min((avg.delta -  delta)^2)], 
           'Val' = avg.delta[which.min((avg.delta -  delta)^2)]))
}


## Adding Signal 
AddSignal <- function(X, signal) { 
  p <- ncol(X) 
  signal.rand <- rep(0, p) 
  signal.rand[1:length(signal)] <- signal
  signal.rand <- sample(signal.rand)
  return(sweep(X, 2, FUN = '+', signal.rand))
}





## Lopes Test
## Input: X - matrix of multivariate normal, Y - matrix of multivariate norm we want to compare
##        B1 - number of repeatition
## Output: P - value of pvalue of Hotelling test to the mean of repetitions


lopesTest <- function(X ,Y , B1 = 1) {
  n.x <- nrow(X)
  n.y <- nrow(Y)
  n  <- n.x + n.y - 2
  p <- ncol(X)
  k <- floor((n / 2))
  f.statistic <- rep(NA , B1)
  if (p > k) {
    for (i in 1:B1){
      Pk <- matrix(rnorm(k *p), p, k)
      cov.mat <- (cov(X) * (n.x - 1) + cov(Y) * (n.y - 1)) / n
      mean.x <- apply(X , 2 , mean)
      mean.y <- apply(Y , 2 , mean)
      t.statistic <- ((n.x * n.y) / (n.x + n.y)) *
        t(t(Pk) %*% (mean.x - mean.y)) %*%
        solve(t(Pk) %*% cov.mat %*% Pk) %*%
        (t(Pk) %*% (mean.x - mean.y))
      f.statistic[i] <- ((n - k + 1) / (k*n)) * t.statistic
    }
    return(1 - pf(mean(f.statistic), k, n - k))
  } else {
    return(as.numeric(HotellingsT2(X,Y)$p.value))
  }
}


# CLX Test  ---------------------------------------------------------------


## Cai Fucntions
## Center columns using means
center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

### Thresholding the eigen - values in order to recieve positive definite matrix
EigenValMax <- function(x , eps = 10^-5) {
  return(ifelse(x > eps, x , eps))
}

## Extreme Value

## CLX test, based on code sent by Xia Yin.
## Input:  X - matrix of multivariate normal, Y - matrix of multivariate norm we want to compare to
## Output: P-value of test


CLX_test <- function(X , Y ) {
  n1 <- nrow(X)
  n2 <- nrow(Y)
  p  <- ncol(X)
  eX <- array(1:1,dim=c(1,n1))
  eY <- array(1:1,dim=c(1,n2))
  A.tz <- (cov(X) * (n1 - 1) + cov(Y) * (n2 - 1)) / ((n1 + n2) - 2)
  X.bar <- (eX %*% X) / n1;
  Y.bar <- (eY %*% Y) / n2;
  Xcov <- cov(X) * (n1 - 1) / n1;
  Ycov <- cov(Y) * (n2 - 1) / n2;
  XXX <- X - ((t(eX))%*%X.bar);
  YYY <- Y - ((t(eY))%*%Y.bar);
  C1 <- (t(XXX^2))%*%(XXX^2)/n1;
  CC1 <- (t(YYY^2))%*%(YYY^2)/n2;
  C2 <- ((t(XXX))%*%(XXX)/n1)*(Xcov);
  CC2 <- ((t(YYY))%*%(YYY)/n2)*(Ycov);
  C3 <- Xcov^2;
  CC3 <- Ycov^2;
  T1 <- C1-2*C2+C3;
  T2 <- CC1-2*CC2+CC3;
  AA.tz <- sqrt( (2 * (n1 + n2) * (A.tz^2)) /(T1 + T2));
  O1.Hat <- (A.tz * (abs(AA.tz) >= 2 * sqrt(log(p))))
  eigen.decomp <- eigen(O1.Hat)
  est.cov <- eigen.decomp$vectors %*% diag(EigenValMax(eigen.decomp$values)) %*% t(eigen.decomp$vectors)
  O1.Hat <- ginv(est.cov)
  Cov.tz <- (cov((X) %*% O1.Hat) * (n1 - 1) + cov((Y) %*% O1.Hat) * (n2 - 1)) / (n1 + n2 -2)
  Z.tz <- O1.Hat %*% t(X.bar - Y.bar)
  M3est <- (max(abs((Z.tz))/abs(diag(Cov.tz))^{1/2}))^2 * (n1 * n2 / (n1 + n2))
  TSvalue <- M3est - 2 * log(p) + log(log(p))
  p.val  <-  1 - exp(-(1/sqrt(pi) * exp(-TSvalue/2)))
  return('Pval' = p.val)
}


### T.test on each dimensions and combining them using Simes 
t.simes <- function(X, Y, var.equal = FALSE) 
{
  P <- ncol(X)
  pval.vec <- rep(NA, P)
  for (i in 1:P){
    pval.vec[i] <- t.test(X[,i], Y[,i], var.equal = var.equal)$p.value
  } 
  return(simesTest(pval.vec))
}
