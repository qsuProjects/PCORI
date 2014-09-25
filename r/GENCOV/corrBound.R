
########################### SETUP ###########################

#desired parameters for HIV data
pmeans = c(46,25,308.1,105337,126.85,76.7,99.1,39.02,180.87,0.12,0.08,0.09,0.16,0.11,0.07,0.40,0.12,0.08,0.07,0.05,0.04,0.17,0.22,0.24)
pps = c(0.275, 0.222, 0.198, 0.50, 0.413, 0.206, 0.688, 0.253, 0.24, 0.184, 0.18, 0.103, 0.332, 0.526, 0.544, 0.97)



########################### FIND-BOUND FUNCTIONS ###########################

##### Function: given parameters for 2 binary RVs, return the minimum and maximum correlation.

BB.rBound = function(p1, p2) {
  q1 = 1-p1; q2 = 1-p2
  
  #compute upper bound
  hiBound = min( sqrt( (p1*q2) / (q1*p2) ) ,
                  sqrt( (q1*p2) / (p1*q2) ) )  
  
  #compute lower bound
  loBound = max( -sqrt( (p1*p2) / (q1*q2) ) ,
                 -sqrt( (q1*q2) / (p1*p2) ) )  
  
  return( c( round(loBound, 2),  round(hiBound, 2) ) )
}


##### Function: given parameters for 1 binary and 1 normal RV, return the maximum correlation.
#for binary-normal correlations, the min and max correlations are the same (in abs value)

BN.rBound = function(p) {
  q = 1-p
  
  #compute upper bound
  hiBound = dnorm( qnorm(p) ) / sqrt(p*q)
  return( round(hiBound, 2) )
}



########################### CORRELATION BOUND MATRIX FUNCTION ###########################

##### Function: given parameters for normal and binary variables, make matrix with bound correlation values
#pmeans: vector of population means for normal variables
#pps: vector of population proportions for binary variables

corrBound = function(pmeans, pps) {
  
  #number of binary and normal variables
  nBin = length(pps)
  nNorm = length(pmeans)
  
  #initialize max-correlation matrix
  maxCorr = matrix(nrow = length(pmeans) + length(pps), ncol = length(pmeans) + length(pps))
  
  #step 1: compute bounds for the binary-binary correlations
  for(i in 1:nBin) {
    for (j in 1:nBin) {
      #fill only upper-diagonal elements to avoid redundancy
      if (i < j) maxCorr[i,j] = paste("[", BB.rBound(pps[i], pps[j])[1], ", " , BB.rBound(pps[i], pps[j])[2], "]", sep="")
    }
  }
  
  #step 2: compute bounds for the binary-normal correlations
  for(i in 1:nBin) {
    for (j in (nBin+1):(nBin+nNorm) ) {
      #fill only upper-diagonal elements to avoid redundancy
      if (i < j) maxCorr[i,j] = paste("[-", BN.rBound(pps[i]), ", " , BN.rBound(pps[i]), "]", sep="")
    }
  }
  
  #step 3: put placeholders for the normal-normal correlations
  #since they have no special bounds
  for(i in (nBin+1):(nBin+nNorm) ) {
    for (j in (nBin+1):(nBin+nNorm) ) {
      #fill only upper-diagonal elements to avoid redundancy
      if (j > i) maxCorr[i,j] = "nn"
      #if (j <= i) maxCorr[i,j] = "-"
    }
  }
    
  print(maxCorr)
}


#apply to desired HIV parameters
( table = corrBound(pmeans, pps) )

#export to csv
setwd("~/Desktop")
write.csv(table, file="Correlation bounds matrix.csv")





################################################################################################
########################### TEST CORRBOUND AGAINST DEMIRTAS FUNCTION ###########################

##### Test Case: Binary-Normal #####
pmeans = c(25.5); var=c(5)
pps = c(.11)

#my bounds
corrBound(pmeans, pps)

#try putting in borderline correlations
r = .18
jointly.generate.binary.normal(no.rows=10, no.bin=1, no.nor=1, prop.vec.bin=pps, mean.vec.nor=pmeans, var.nor=var, corr.vec=r)
#my function gives correct bounds



##### Test Case: Binary-Binary #####

#(had to include some normal vars due to function constraints)
pmeans = c(100); var=c(300)
pps = c(.001, .001)

#my bounds
corrBound(pmeans, pps)

#try putting in some borderline correlations
r = c(.99, .05, .05)
jointly.generate.binary.normal(no.rows=10, no.bin=2, no.nor=1, prop.vec.bin=pps, mean.vec.nor=pmeans, var.nor=var, corr.vec=r)
#my function gives correct bounds


