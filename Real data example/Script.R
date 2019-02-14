## Real Data Example 
library(dplyr)
library(NonParametricMean)
library(highD2pop)
library(highmean)

setwd('C:/Users/Tzviel/Google Drive/Simes Hotelling Paper/R Files/Real Data Example/Cancer Example 02082017')
source('C:/Users/Tzviel/Google Drive/Simes Hotelling Paper/R Files/Simulation 11.07.2017 - Nonequal Covariance matrices/Function Simes Hotelling.R')
cancer.dat <- read.csv('CarcinomaNormaldatasetCancerResearch.csv')
set.seed(999)


## Simes Hotelling Functions 

RunTests <- function(X, Y) { 
  p <- ncol(X) 
  n <- floor(nrow(X) + nrow(Y))
  c('Srivastava' = apval_Sri2008(X, Y)$pval,
    'Chen'  = apval_Chen2010(X, Y)$pval,
    'GCT'   = GCT.test(X, Y, 25)$pvalue,
    'Simes' = t.simes(X, Y))
}

RunRandTest <- function(X, Y, equal.cov = F) {
  p <- ncol(X) 
  n <- floor(nrow(X) + nrow(Y))
  c('Lopes' = lopesTest(X, Y, B1 = 1),
    'SH.Ln.LB'    = SimesHotelingArma(list(X, Y), samp.size = 18, iterations = 60000),
    'SH.Sn.LB'    = SimesHotelingArma(list(X, Y), samp.size = 18, iterations = 60000))
}

### Making Matrices 


## Tumor data SD 
tumor.dat <- t(cancer.dat %>% dplyr::select(contains('Tumor')))
norm.dat  <-  t(cancer.dat %>% dplyr::select(contains('Normal')))




p.val <- RunTests(tumor.dat, norm.dat)
p.val.rand <- replicate(10, RunRandTest(tumor.dat, norm.dat))
avg.pval.rand <- apply(p.val.rand, 1, mean)
sd.pval.rand  <- apply(p.val.rand, 1, sd) 



### Data from 
#Broad patterns of gene expression revealed by 
#clustering of tumor and normal colon tissues probed by oligonucleotide arrays
I2000CSV <- read.csv('I2000CSV.csv', header = FALSE)
ident    <- read.csv('I2000Ident.csv', header = FALSE)

rem.ind  <- which(apply(I2000CSV, 1, function(x) return(all(is.na(x)))))
full.mat <- t(I2000CSV[-rem.ind,])

## Removing Duplicates 

dup.ind <- which(duplicated(t(full.mat)))
norm.dat <- full.mat[which(sign(ident) == 1), -dup.ind]
canc.dat <- full.mat[which(sign(ident) == -1), -dup.ind]

## Calculate 
SH.pval <- replicate(10, SimesHotelingArma(list(norm.dat, canc.dat), 31, 50000))
RunTests(norm.dat, canc.dat)



### Testing type 1 error
samp1 <- sample(40, 20)
samp2 <- setdiff(1:40, samp1)
fake.dat.1 <- canc.dat[samp1,]
fake.dat.2 <- canc.dat[samp2,]


replicate(10, SimesHotelingArma(list(fake.dat.1, fake.dat.2), 20, 50000))
RunTests(fake.dat.1, fake.dat.2)

any(cor(full.mat) == 1 )
duplicated(t(full.mat))


### Calcium Data 
library(fda.usc)
data(MCO)




## Changing for non-equal covariance matrices
RunTests <- function(X, Y, equal.cov = F) { 
  p <- ncol(X) 
  n <- floor(nrow(X) + nrow(Y))
  c('Chen'  = apval_Chen2010(X, Y, eq.cov = F)$pval,
    'GCT'   = GCT.test(X, Y, 25)$pvalue,
    'Simes' = t.simes(X, Y))
}

RunRandTest <- function(X, Y, equal.cov = F) {
  p <- ncol(X) 
  n <- floor(nrow(X) + nrow(Y))
  c('SH.Ln.LB'    = SimesHotellingB1(X, Y, round(n/2) , p * log(p), equal.cov),
    'SH.Sn.LB'    = SimesHotellingB1(X, Y, round(n/5) ,  p * log(p), equal.cov))
}


intact.group.a <- (MCO$intact$data[MCO$classintact == 1, -c(1:18)])
intact.group.b <- (MCO$intact$data[MCO$classintact == 2, -c(1:18)])


p.val.intact <- RunTests(intact.group.a, intact.group.b)
p.val.rand.intact <- replicate(10, RunRandTest(intact.group.a, intact.group.b, equal.cov = F))
avg.pval.rand.intact <- apply(p.val.rand.intact, 1, mean)
sd.pval.rand.intact  <- apply(p.val.rand.intact, 1, sd) 


perma.group.a <- (MCO$permea$data[MCO$classpermea == 1, -c(1:18)])
perma.group.b <- (MCO$permea$data[MCO$classpermea == 2, -c(1:18)])


p.val.perm <- RunTests(perma.group.a, perma.group.b)
p.val.rand.perm <- replicate(10, RunRandTest(perma.group.a, perma.group.b))
avg.pval.rand.perm <- apply(p.val.rand.perm, 1, mean)
sd.pval.rand.perm  <- apply(p.val.rand.perm, 1, sd) 
