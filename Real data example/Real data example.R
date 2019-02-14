## Libraries 
setwd(here::here())
source('Simulation Functions.R')


## Directory
setwd(paste0(here::here(), '/Real data example'))
## Simes Hotelling Functions 

RunTests <- function(X, Y, equal.cov = TRUE) { 
  p <- ncol(X) 
  n <- floor(nrow(X) + nrow(Y))
  c('SD'  = apval_Sri2008(X, Y)$pval,
    'CQ'  = apval_Chen2010(X, Y, eq.cov = equal.cov)$pval,
    'GCT'   = GCT.test(X, Y, 25)$pvalue,
    'Simes' = t.simes(X, Y, var.equal= equal.cov))
}

RunRandTest <- function(X, Y, equal.cov = F) {
  p <- ncol(X) 
  n <- floor(nrow(X) + nrow(Y))
  c('Lopes'       = lopesTest(X, Y, B1 = 1),
    'SH.Ln.LB'    = SHTest(X, Y, samp.size = round((nrow(X) + nrow(Y)) / 2), 
                           iterations = ncol(X) * log(ncol(X)), equal.cov = equal.cov)$`P-value`,
    'SH.Sn.LB'    = SHTest(X, Y, samp.size = round((nrow(X) + nrow(Y)) / 4),
                           iterations = ncol(X) * log(ncol(X)), equal.cov = equal.cov)$`P-value`)
}




# First example  ----------------------------------------------------------
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
set.seed(999)
rand.pval <- replicate(10, RunRandTest(norm.dat, canc.dat))
rand.pval <- apply(rand.pval, 1, mean) 
### Removed GCT (See testing type 1 error)
print('Results for colon data are')
c(RunTests(norm.dat, canc.dat, equal.cov = TRUE)[-3], rand.pval)


### Testing type 1 error
samp1 <- sample(40, 20)
samp2 <- setdiff(1:40, samp1)
fake.dat.1 <- canc.dat[samp1,]
fake.dat.2 <- canc.dat[samp2,]

## Calculate 
set.seed(999)
rand.pval <- replicate(10, RunRandTest(fake.dat.1, fake.dat.2))
rand.pval <- apply(rand.pval, 1, mean) 
### Removed GCT (See testing type 1 error)
print('Results for fake data are')
c(RunTests(fake.dat.1, fake.dat.2), rand.pval)



# Second example  ---------------------------------------------------------
data(MCO)
set.seed(999)

intact.group.a <- (MCO$intact$data[MCO$classintact == 1, -c(1:18)])
intact.group.b <- (MCO$intact$data[MCO$classintact == 2, -c(1:18)])

p.val.intact      <- RunTests(intact.group.a, intact.group.b, equal.cov = FALSE)
p.val.rand.intact <- apply(replicate(10, RunRandTest(intact.group.a, intact.group.b)), 1, mean)
print('Results for intact cells data are')
c(p.val.intact[-1], p.val.rand.intact[-1])

perma.group.a <- (MCO$permea$data[MCO$classpermea == 1, -c(1:18)])
perma.group.b <- (MCO$permea$data[MCO$classpermea == 2, -c(1:18)])


p.val.perm      <- RunTests(perma.group.a, perma.group.b, equal.cov = FALSE)
p.val.rand.perm <- apply(replicate(10, RunRandTest(perma.group.a, perma.group.b, equal.cov = FALSE)), 1, mean)
print('Results for permablized cells data are')

c(p.val.perm[-1], p.val.rand.perm[-1])
