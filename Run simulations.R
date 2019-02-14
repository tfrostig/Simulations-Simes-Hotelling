### Making sure directory is set correctly 
setwd(here::here()) 
## Loading functions 
source('Simulation Functions.R')


##### Running simulations 
### Preparing parallelization 
pack.name <- c('MASS', 'simesHotelling', 'Matrix', 'highmean', 'highD2pop')
cl = makeCluster(parallel::detectCores() - 1) ## Detect cores
registerDoParallel(cl)

## Default values 
iter.num   <- 3
false.hypo <- c(0.01, 0.05, 0.15) 
obs        <- c(20, 50, 100)
type       <- 'Unif'
alpha      <- 0.05


setwd(here::here()) 
false.hypo <- c(0) 
source('Parametric Tests/Type I error/Simulation Script Parametric - Alpha.R')
setwd(here::here()) 
false.hypo <- c(0.01, 0.05, 0.15) 
source('Parametric Tests/Power/Simulation Script Parametric - Power.R')
setwd(here::here()) 
false.hypo <- c(0) 
source('Non-equal Covaraince/Type I error/Simulation Script Nonequal - Alpha.R')
setwd(here::here()) 
false.hypo <- c(0.01, 0.05, 0.15) 
source('Non-equal Covaraince/Power/Simulation Script Nonequal - Power.R')



#### Plotting results 
### Making sure directory is set correctly 
setwd(here::here()) 
## Loading functions 
source('Simulation Functions.R')
### Plotting parametric power results
setwd(here::here())
source('Parametric Tests/Parametric - plotting and princeton.R')
### Plotting parametric power results 
setwd(here::here())
source('Parametric Tests/Parametric - plotting and princeton.R')








