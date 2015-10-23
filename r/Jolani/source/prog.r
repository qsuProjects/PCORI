# required library
require(bindata)
require(mvtnorm)
require(lme4)
require(mice)
require(MASS)

# loading setup values
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/Jolani/source")
load("setup.RData")
source("subfunctions.r")

rep <- 500
seed <- 2781903
fmi <- matrix(NA, nrow = rep, ncol = 3)
result <- array(NA, dim = c(length(setup$fixef),(4+3),5,rep))                
dimnames(result) <- list(c(attr(setup$fixef, "names")),
                         c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "r.int", "r.malign", "r.calf"),
                         c("full", "fca", "tmi", "mifix", "mlmi"), 1:rep)
data <- array(NA, dim = c(3000,4,2,rep))

# This is only for one systematically missing predictor (scenarios 1-4). 
# Different scenarios can be obtained by changing 'stnum' and 'perc'.
for (i in 1:rep){
  temp <- annals(setup = setup, stnum = "small", perc = "low", pred = "one", MAX = 1, M=5)
  result[,,"full",i] <- cbind(temp$f.full, temp$r.full)
  result[,,"fca",i] <- cbind(temp$f.fca, temp$r.fca)
  result[,,"tmi",i] <- cbind(temp$f.tmi, temp$r.tmi)
  result[,,"mifix",i] <- cbind(temp$f.mifix, temp$r.mifix)  
  result[,,"mlmi",i] <- cbind(temp$f.mlmi, temp$r.mlmi)
  fmi[i, ] <- temp$fmi.mlmi
  data[,,1,i] <- as.matrix(temp$fdata)
  data[,,2,i] <- as.matrix(temp$mdata)
}
sysmis <- temp$sysmis
tstmis <- temp$tstmis

final <- list(data = data, result=result, fmi = fmi, sysmis=sysmis, tstmis=tstmis, seed=seed)


# For scenarios 5 and 6 the following code must be added in the loop.
# Scenario 6 can be obtained by changing 'perc = "high"'.
temp <- annals(setup = setup, stnum = "large", perc = "low", pred = "two", MAX = 10, M=5)


