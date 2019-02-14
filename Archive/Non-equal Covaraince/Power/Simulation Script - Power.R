
### Creating Simulation
library('doParallel')
library('MASS')
library('highmean')
library('highD2pop')
library('doRNG')
library('dplyr')
library('rmutil')
library('MTSKNN')
library('NonParametricMean')


## Loading Functions 
source('C:/Users/user/Google Drive/Simes Hotelling Paper/Simes Hoteling V2 - No HC/R Script/Full Simulation 21112018/Data Generation.R')
setwd('C:/Users/user/Google Drive/Simes Hotelling Paper/Simes Hoteling V2 - No HC/R Script/Full Simulation 21112018/Non-equal Covaraince/Power')


### Preparing Parallelization 
FuncName <- c('SimesHotelingArma', 'hotelingTest2', 'findMaxHotelingArma', 'simesTest',
              'findCluster', 'subSpaceClusterTest', 'CombineTest', 'PermuteTal', 
              'PermutationTests')
CppFunc  <- c('MahalanobisArma', 'RandomHoteling', 'CombineVectors', 'HotelingFind', 
              'CreateSol', 'Permute', 'RandomClusterSH') 
pacName <- c('MASS', 'NonParametricMean', 'Matrix'  ,'highmean', 'highD2pop', 'rmutil' ,'MTSKNN')


cl = makeCluster(8)
registerDoParallel(cl)
clusterExport(cl, c(CppFunc, FuncName))
  clusterEvalQ(cl, library('MASS'))

## Default Values 
iter.num <- 1000
false.hypo <-  c(0.01, 0.05, 0.15)
obs   <- c(20, 50, 100)
type  <- 'Unif'

## Simulation Functions 
set.seed(999)

## Name Vector 
Name.Vec <- c('Srivastava', 'Chen', 'SH.Ln.LB', 'SH.Sn.LB', 'GCT', 'Simes')
L        <- length(Name.Vec)
## Run tests 
RunTests <- function(X, Y) { 
  p <- ncol(X) 
  n <- nrow(X) + nrow(Y)
  mat.list <- list('X' = X, 'Y' = Y)
  temp.test <-c('Srivastava'  = apval_Sri2008(X, Y)$pval,
                'Chen'        = apval_Chen2010(X, Y, eq.cov = FALSE)$pval,
                'SH.Ln.LB'    = SimesHotelingArma(mat.list, round(n/2) , p * log(p), equal.cov = TRUE), 
                'SH.Sn.LB'    = SimesHotelingArma(mat.list, round(n/4) , p * log(p), equal.cov = TRUE), 
                'GCT'   = GCT.test(X, Y, 25)$pvalue, 
                'Simes' = t.simes(X, Y, var.equal = FALSE)

  )
  return(temp.test)
}

# Simulation 1  -----------------------------------------------------------

## Using Corr2 (exponential decrease in correlation as function of 'Distance')
## Variance is 1 

fam.size <- 600
correlation <- c(0.3, 0.45, 0.6, 0.75, 0.95)
power.simulation.matrix.1 <- data.frame(expand.grid(false.hypo, obs, correlation, type), 
                                        matrix(NA ,ncol = L))
sim.length <- nrow(power.simulation.matrix.1)
colnames(power.simulation.matrix.1) <- c('False Hypotheses' ,
                                         'Obs','Correlation', 'Type', 
                                         Name.Vec)

for (i in 1:nrow(power.simulation.matrix.1)) { 
  temp.false.hypo <- power.simulation.matrix.1[i,'False Hypotheses']
  temp.obs        <- power.simulation.matrix.1[i,'Obs']
  temp.cor        <- power.simulation.matrix.1[i,'Correlation']
  temp.type       <- power.simulation.matrix.1[i,'Type']
  temp.false.num  <- round(temp.false.hypo * fam.size)
  cov.mat         <- corr2(fam.size, temp.cor, v.d = 1)
  temp.res <- foreach(j = 1:iter.num, .combine = rbind, 
                      .packages = pacName, 
                      .options.RNG = 9999, 
                      .noexport = CppFunc) %dorng% { 
                        temp.dat        <-  CreateData(temp.false.hypo, 
                                                       euclid.norm = 2.5 * sqrt(20 / temp.obs), 
                                                       group.num = 5, 
                                                       n.x = temp.obs, 
                                                       n.y = temp.obs, 
                                                       cov.mat.x = cov.mat, 
                                                       cov.mat.y = 2 * cov.mat,
                                                       noise.type = temp.type,
                                                       signal.type = 'Unif')
                        test.pval <- tryCatch(RunTests(temp.dat$X, temp.dat$Y), error = function(e) e) 
                        if (inherits(test.pval, 'error')) {
                          rep(NA, length(Name.Vec)) 
                        }  else { 
                          test.pval   
                        }
                      }
  mean.res <- apply(temp.res, 2, function(x) mean(x < 0.05))
  power.simulation.matrix.1[i,Name.Vec] <- mean.res
  write.csv(power.simulation.matrix.1, 'Power.Simulation.Matrix.1.csv', row.names = FALSE)
  print(paste('Iteration', i, 'Out Of', nrow(power.simulation.matrix.1)))
}




# Simulation 3 - Block covariance matrix  ---------------------------------


fam.size <- 600
block.size <- c(20, 50, 100)
rho <- c(0.5, 0.65, 0.8)


#### Settings 
power.simulation.matrix.3 <- data.frame(expand.grid(false.hypo, obs, block.size, rho, type), 
                                        matrix(NA ,ncol = length(Name.Vec)))
sim.length <- nrow(power.simulation.matrix.3)
colnames(power.simulation.matrix.3) <- c('False Hypotheses' ,
                                         'Obs','Block_Size', 'Rho', 'Type', 
                                         Name.Vec)


for (i in 1:nrow(power.simulation.matrix.3)) { 
  temp.false.hypo <- power.simulation.matrix.3[i,'False Hypotheses']
  temp.obs        <- power.simulation.matrix.3[i,'Obs']
  temp.rho        <- power.simulation.matrix.3[i,'Rho']
  temp.block.size <- power.simulation.matrix.3[i,'Block_Size']
  temp.type       <- power.simulation.matrix.3[i,'Type']
  temp.false.num  <- round(temp.false.hypo * fam.size)
  cov.mat         <- corrBlock(fam.size, temp.block.size, correlation = temp.rho) 
  temp.res <- foreach(j = 1:iter.num, .combine = rbind, 
                      .packages = pacName, 
                      .options.RNG = 9999, 
                      .noexport = CppFunc) %dorng% { 
                        temp.dat        <- CreateData(temp.false.hypo, 
                                                      euclid.norm = 2.5 * sqrt(20 / temp.obs), 
                                                      group.num = 5, 
                                                      n.x = temp.obs, 
                                                      n.y = temp.obs, 
                                                      cov.mat.x = cov.mat, 
                                                      cov.mat.y = 2 * cov.mat,
                                                      noise.type = temp.type,
                                                      signal.type = 'Unif')
                        test.pval <- tryCatch(RunTests(temp.dat$X, temp.dat$Y), error = function(e) e) 
                        if (inherits(test.pval, 'error')) {
                          rep(NA, length(Name.Vec)) 
                        }  else { 
                          test.pval   
                        }
                      }
  mean.res <- apply(temp.res, 2, function(x) mean(x < 0.05))
  power.simulation.matrix.3[i,Name.Vec] <- mean.res
  print(paste('Iteration', i, 'Out Of', nrow(power.simulation.matrix.3)))
}

write.csv(power.simulation.matrix.3, 'Power.Simulation.Matrix.3.csv', row.names = FALSE)


# Simulation 5 - Cai Model 7  ---------------------------------------------


fam.size <- 600

power.simulation.matrix.5 <- data.frame(cbind(expand.grid(false.hypo, obs, type), matrix(NA ,ncol = length(Name.Vec))))
colnames(power.simulation.matrix.5) <- c('False Hypotheses','Number of Obs', 'Type', Name.Vec)
sim.length <- nrow(power.simulation.matrix.5)



## Creating Simulation

for (i in 1:nrow(power.simulation.matrix.5)) {
  false.hypo.temp <- power.simulation.matrix.5[i,1]  
  obs.temp        <- power.simulation.matrix.5[i,2] 
  temp.type       <- power.simulation.matrix.5[i,3]
  cov.mat         <- CovMakerModel7(fam.size)
  ## Creating Simulation 
  temp.res <- foreach(j = 1:iter.num, .combine = rbind, 
                      .packages = pacName, 
                      .options.RNG = 9999, 
                      .noexport = CppFunc) %dorng% { 
                        temp.dat        <- CreateData(temp.false.hypo, 
                                                      euclid.norm = 5 * sqrt(20 / temp.obs),  ## Due to sigma_ii ~ Unif[1,3] E(sigma_ii) = 2
                                                      group.num = 5, 
                                                      n.x = temp.obs, 
                                                      n.y = temp.obs, 
                                                      cov.mat.x = cov.mat, 
                                                      cov.mat.y = 2 * cov.mat,
                                                      noise.type = temp.type,
                                                      signal.type = 'Unif')
                        test.pval <- tryCatch(RunTests(temp.dat$X, temp.dat$Y), error = function(e) e) 
                        if (inherits(test.pval, 'error')) {
                          rep(NA, length(Name.Vec)) 
                        }  else { 
                          test.pval   
                        }
                      }
  mean.res <- apply(temp.res, 2, function(x) mean(x < 0.05))
  power.simulation.matrix.5[i, Name.Vec] <- mean.res
  print(paste(i,'out of',nrow(power.simulation.matrix.5)))
}


write.csv(power.simulation.matrix.5,'power.simulation.matrix.5.csv', row.names = FALSE)


# Simulation 7 - Equicorrelated  ------------------------------------------

## Using Corr2 (exponential decrease in correlation as function of 'Distance')
## Variance is 1 

fam.size <- 600
correlation <- c(0.1, 0.3, 0.5)

power.simulation.matrix.7 <- data.frame(expand.grid(false.hypo, obs, correlation, type), matrix(NA ,ncol = length(Name.Vec)))
colnames(power.simulation.matrix.7) <- c('False Hypotheses','Number of Obs','Correlation','Type', Name.Vec)

## Run Simulation 
for (i in 1:nrow(power.simulation.matrix.7)) {
  false.hypo.temp <- power.simulation.matrix.7[i,1]  
  obs.temp  <- power.simulation.matrix.7[i,2] 
  corr.temp <- power.simulation.matrix.7[i,3] 
  temp.type <- power.simulation.matrix.7[i,4]
  cov.mat   <- corr.full(fam.size, corr.temp)
  temp.res <- foreach(j = 1:iter.num, .combine = rbind, 
                      .packages = pacName, 
                      .options.RNG = 9999, 
                      .noexport = CppFunc) %dorng% { 
                        temp.dat        <- CreateData(false.hypo.temp, 
                                                      euclid.norm = 2.5 * sqrt(20 / obs.temp), 
                                                      group.num = 5, 
                                                      n.x = obs.temp, 
                                                      n.y = obs.temp, 
                                                      cov.mat.x = cov.mat, 
                                                      cov.mat.y = 2 * cov.mat,
                                                      noise.type = temp.type,
                                                      signal.type = 'Unif')
                        test.pval <- tryCatch(RunTests(temp.dat$X, temp.dat$Y), error = function(e) e) 
                        if (inherits(test.pval, 'error')) {
                          rep(NA, length(Name.Vec)) 
                        }  else { 
                          test.pval   
                        }
                      }
  mean.res <- apply(temp.res, 2, function(x) mean(x < 0.05, na.rm = TRUE))
  power.simulation.matrix.7[i, Name.Vec] <- mean.res
  print(paste(i,'out of',nrow(power.simulation.matrix.7)))
}

write.csv(power.simulation.matrix.7,'power.simulation.matrix.7.csv', row.names = FALSE)