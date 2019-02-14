
set.seed(998)
setwd('Non-equal Covaraince/Power')

## Name Vector 
name.vec <- c('SD', 'CQ', 'SH.Ln.LB', 'SH.Sn.LB', 'Simes') 
L        <- length(name.vec)

## Run tests 
RunTests <- function(X, Y) { 
  p <- ncol(X) 
  n <- nrow(X) + nrow(Y)
  temp.test <-c('SD'  = apval_Sri2008(X, Y)$pval,
                'CQ'        = apval_Chen2010(X, Y, eq.cov = FALSE)$pval,
                'SH.Ln.LB'    = SHTest(X, Y, samp.size = round(n / 2) , iterations = p * log(p), equal.cov = FALSE)$`P-value`,
                'SH.Sn.LB'    = SHTest(X, Y, samp.size = round(n / 4) , iterations = p * log(p), equal.cov = FALSE)$`P-value`, 
                'Simes'       = t.simes(X, Y, var.equal = FALSE)
                
  )
  return(temp.test)
}

# Simulation 1 - AR(1) -----------------------------------------------------------

## Using Corr2 (exponential decrease in correlation as function of 'Distance')
## Variance is 1 

fam.size <- 600
correlation <- c(0.3, 0.45, 0.6, 0.75, 0.95)
power.simulation.matrix.1 <- data.frame(expand.grid(false.hypo, obs, correlation, type), 
                                        matrix(NA ,ncol = L))
sim.length <- nrow(power.simulation.matrix.1)
colnames(power.simulation.matrix.1) <- c('False Hypotheses' ,
                                         'Obs','Correlation', 'Type', 
                                         name.vec)
print('Simulation AR(1)')

for (i in 1:nrow(power.simulation.matrix.1)) { 
  temp.false.hypo <- power.simulation.matrix.1[i,'False Hypotheses']
  temp.obs        <- power.simulation.matrix.1[i,'Obs']
  temp.cor        <- power.simulation.matrix.1[i,'Correlation']
  temp.type       <- power.simulation.matrix.1[i,'Type']
  temp.false.num  <- round(temp.false.hypo * fam.size)
  cov.mat         <- corr2(fam.size, temp.cor, v.d = 1)
  temp.res <- foreach(j = 1:iter.num, 
                      .combine = rbind, 
                      .packages = pack.name, 
                      .options.RNG = 9999) %dorng% { 
                        temp.dat        <-  createData(temp.false.hypo, 
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
                          rep(NA, L) 
                        }  else { 
                          test.pval   
                        }
                      }
  mean.res <- apply(temp.res, 2, function(x) mean(x < alpha))
  power.simulation.matrix.1[i,name.vec] <- mean.res
  print(paste('Iteration', i, 'Out Of', nrow(power.simulation.matrix.1)))
}

write.csv(power.simulation.matrix.1, 'nonequal.power.simulation.matrix.1.csv', row.names = FALSE)



# Simulation 3 - Block covariance matrix  ---------------------------------


fam.size <- 600
block.size <- c(20, 50, 100)
rho <- c(0.5, 0.65, 0.8)


#### Settings 
power.simulation.matrix.3 <- data.frame(expand.grid(false.hypo, obs, block.size, rho, type), 
                                        matrix(NA ,ncol = L))
sim.length <- nrow(power.simulation.matrix.3)
colnames(power.simulation.matrix.3) <- c('False Hypotheses' ,
                                         'Obs','Block_Size', 'Rho', 'Type', 
                                         name.vec)

print('Simulation block covariance')

for (i in 1:nrow(power.simulation.matrix.3)) { 
  temp.false.hypo <- power.simulation.matrix.3[i,'False Hypotheses']
  temp.obs        <- power.simulation.matrix.3[i,'Obs']
  temp.rho        <- power.simulation.matrix.3[i,'Rho']
  temp.block.size <- power.simulation.matrix.3[i,'Block_Size']
  temp.type       <- power.simulation.matrix.3[i,'Type']
  temp.false.num  <- round(temp.false.hypo * fam.size)
  cov.mat         <- corrBlock(fam.size, temp.block.size, correlation = temp.rho) 
  temp.res <- foreach(j = 1:iter.num, 
                      .combine = rbind, 
                      .packages = pack.name, 
                      .options.RNG = 9999) %dorng% { 
                        temp.dat        <- createData(temp.false.hypo, 
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
                          rep(NA, L) 
                        }  else { 
                          test.pval   
                        }
                      }
  mean.res <- apply(temp.res, 2, function(x) mean(x < alpha))
  power.simulation.matrix.3[i,name.vec] <- mean.res
  print(paste('Iteration', i, 'Out Of', nrow(power.simulation.matrix.3)))
}

write.csv(power.simulation.matrix.3, 'nonequal.power.simulation.matrix.3.csv', row.names = FALSE)


# Simulation 5 - Cai Model 7  ---------------------------------------------


fam.size <- 600

power.simulation.matrix.5 <- data.frame(expand.grid(false.hypo, obs, type), 
                                              matrix(NA ,ncol = L))
colnames(power.simulation.matrix.5) <- c('False Hypotheses','Obs', 'Type', name.vec)
sim.length <- nrow(power.simulation.matrix.5)



print('Simulation Cai model')

for (i in 1:nrow(power.simulation.matrix.5)) {
  false.hypo.temp <- power.simulation.matrix.5[i,1]  
  obs.temp        <- power.simulation.matrix.5[i,2] 
  temp.type       <- power.simulation.matrix.5[i,3]
  cov.mat         <- CovMakerModel7(fam.size)
  ## Creating Simulation 
  temp.res <- foreach(j = 1:iter.num, 
                      .combine = rbind, 
                      .packages = pack.name, 
                      .options.RNG = 9999) %dorng% { 
                        temp.dat        <- createData(temp.false.hypo, 
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
                          rep(NA, L) 
                        }  else { 
                          test.pval   
                        }
                      }
  mean.res <- apply(temp.res, 2, function(x) mean(x < alpha))
  power.simulation.matrix.5[i, name.vec] <- mean.res
  print(paste(i,'out of',nrow(power.simulation.matrix.5)))
}

 
write.csv(power.simulation.matrix.5, 'nonequal.power.simulation.matrix.5.csv', row.names = FALSE)


# Simulation 7 - Equicorrelated  ------------------------------------------

## Using Corr2 (exponential decrease in correlation as function of 'Distance')
## Variance is 1 

fam.size <- 600
correlation <- c(0.1, 0.3, 0.5)

power.simulation.matrix.7 <- data.frame(expand.grid(false.hypo, obs, correlation, type),
                                        matrix(NA ,ncol = L))
colnames(power.simulation.matrix.7) <- c('False Hypotheses','Obs','Correlation','Type', name.vec)
print('Simulation equicorrelation')

## Run Simulation 
for (i in 1:nrow(power.simulation.matrix.7)) {
  false.hypo.temp <- power.simulation.matrix.7[i, 'False Hypotheses']  
  obs.temp  <- power.simulation.matrix.7[i, 'Obs'] 
  corr.temp <- power.simulation.matrix.7[i, 'Correlation'] 
  temp.type <- power.simulation.matrix.7[i, 'Type']
  cov.mat   <- corr.full(fam.size, corr.temp)
  temp.res <- foreach(j = 1:iter.num, 
                      .combine = rbind, 
                      .packages = pack.name, 
                      .options.RNG = 9999) %dorng% { 
                        temp.dat        <- createData(false.hypo.temp, 
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
                          rep(NA, L) 
                        }  else { 
                          test.pval   
                        }
                      }
  mean.res <- apply(temp.res, 2, function(x) mean(x < alpha, na.rm = TRUE))
  power.simulation.matrix.7[i, name.vec] <- mean.res
  print(paste(i,'out of',nrow(power.simulation.matrix.7)))
}

write.csv(power.simulation.matrix.7,'nonequal.power.simulation.matrix.7.csv', row.names = FALSE)