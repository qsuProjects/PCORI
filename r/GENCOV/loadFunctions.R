#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Project: PCORI Missing Data
# 
# -This program creates covariates for simulated data based on existing VA data 
#    (summary stats provided by Vilija at VA, see Q:\Datasets\PCORI\data simulation\real data estimates)
#
# -Using jointly_generate_binary_normal (binNor package) function by Demirtas, simulate clinical characteristics and drug exposure for
#  a large number of patients over time.
#
# -This file loads all necessary functions. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

############################## LOAD PACKAGES ##############################

# use the program by Demirtas (modified by MM)
source("/Volumes/QSU/Datasets/PCORI/data simulation/generate covariates/Current simulation files/jointly_generate_binary_normal_modified_v2.R")

# load packages
library(MBESS)
library(mvtnorm)
library(ICC)
library(miscTools)
library(car)


############################## FUNCTION: OVERRIDE DRUG PROBS ##############################

# For subjects who were never on a given drug, changes the corresponding drug propensity to zero.
# overrides drug probabilities if indicator = 0

# mus0: matrix of subject means
# n.Drugs: number of drugs

overrideDrugProbs = function(mus0, n.Drugs) {
  # number of non-drug variables
  n.NonDrugVars = dim(mus0)[2] - n.Drugs
  
  for (m in 1:n.Drugs) {
    # replace each entry in the part of mus0 corresponding to the drug propensities with 0...
    # ...if the part of mus0 corresponding to the drug indicator is 0
    # otherwise leave the drug propensity alone
    mus0[, m + n.NonDrugVars] = ifelse( mus0[,m] == 0, zero, mus0[, m + n.NonDrugVars])
  }
  
  return(mus0)
}


############################## FUNCTION: UPPER TRI VEC ##############################

# turns matrix into vector of upper-triangular elements
# matrix elements are arranged by row

upperTriVec = function(m) {
  v1 = as.vector( t(m) )
  keepElement = as.vector( t(upper.tri(m) ) ) #use transpose to avoid going by columns
  v2 = as.numeric( v1[keepElement] )
  return(v2)
}


############################## FUNCTION: PROPORTIONIZE ##############################

# Changes a number into a valid proportion. 
# Numbers >= 1 become one.
# Numbers <= 0 become zero.

# x: the number to proportionize
# zero: a very small number that is almost 0
# one: a number that is almost 1

proportionize = function(x, zero, one) {
  y = ifelse(x <= 0, zero, x)
  y = ifelse(y >= 1, one, y)
  return(y)
}


############################## FUNCTION: GENERATE RACE ##############################

#datc: the dataset without race
add_race <- function(d3, n, obs) {
  
  # create dataframe with just the first observation from each patient
  first.dat <- d3[!duplicated(d3$id),]
      
  # estimate logits for black and other-race
  logit.race.b <- bb0 + bb1*first.dat$ever_abac + bb2*first.dat$ever_ataz + bb3*first.dat$ever_dida + bb4*first.dat$ever_efav + bb5*first.dat$ever_emtr +
    bb6*first.dat$ever_indi + bb7*first.dat$ever_lami + bb8*first.dat$ever_lopi + bb9*first.dat$ever_nelf + bb10*first.dat$ever_nevi +
    bb11*first.dat$ever_rito + bb12*first.dat$ever_saqu + bb13*first.dat$ever_stav + bb14*first.dat$ever_teno + bb15*first.dat$ever_zido +
    bb16*first.dat$male + 
    bb17*(first.dat$age_cat=="b.35to45") + bb18*(first.dat$age_cat=="c.45to55") + bb19*(first.dat$age_cat=="d.55to65") + bb20*(first.dat$age_cat=="e.Over65") +  
    bb21*(first.dat$bmi_cat=="b.Under20") + bb22*(first.dat$bmi_cat=="c.25to35") + bb23*(first.dat$bmi_cat=="d.Over30") + 
    bb33*first.dat$bpd + bb34*first.dat$bps + bb35*first.dat$hdl + bb36*first.dat$ldl + bb37*first.dat$trig
  
  logit.race.o <- bo0 + bo1*first.dat$ever_abac + bo2*first.dat$ever_ataz + bo3*first.dat$ever_dida + bo4*first.dat$ever_efav + bo5*first.dat$ever_emtr +
    bo6*first.dat$ever_indi + bo7*first.dat$ever_lami + bo8*first.dat$ever_lopi + bo9*first.dat$ever_nelf + bo10*first.dat$ever_nevi +
    bo11*first.dat$ever_rito + bo12*first.dat$ever_saqu + bo13*first.dat$ever_stav + bo14*first.dat$ever_teno + bo15*first.dat$ever_zido +
    bo16*first.dat$male + 
    bo17*(first.dat$age_cat=="b.35to45") + bo18*(first.dat$age_cat=="c.45to55") + bo19*(first.dat$age_cat=="d.55to65") + bo20*(first.dat$age_cat=="e.Over65") +  
    bo21*(first.dat$bmi_cat=="b.Under20") + bo22*(first.dat$bmi_cat=="c.25to35") + bo23*(first.dat$bmi_cat=="d.Over30") + 
    bo33*first.dat$bpd + bo34*first.dat$bps + bo35*first.dat$hdl + bo36*first.dat$ldl + bo37*first.dat$trig
  
  # probability of being each race
  # each is a vector where the jth entry is the probability for the jth person
  p.race1 <- exp(logit.race.b)/(1+exp(logit.race.b)+exp(logit.race.o))
  p.race2 <- exp(logit.race.o)/(1+exp(logit.race.b)+exp(logit.race.o))
  p.race0 <- 1/(1+exp(logit.race.b)+exp(logit.race.o))
  #p.race0 <- ifelse(p.race0 < 0 ,0, p.race0)
  
  # create vector of lists for each person's race
  race <- vector("list", n)
  
  # generate each person's race multinomially
  for (j in 1:n)  race[[j]] <- rmultinom(1, 1, c(p.race0[j], p.race1[j], p.race2[j]))
  
  race <- t(do.call("cbind", race))
  
  # recode race
  # if the third col is 1, race=2
  # if the second col is 1, race=1
  # if the first col is 1, race=0
  # the i,jth entry of race now corresponds to whether person i was of race j
  racecat <- ifelse(race[,3] == 1, 2, ifelse(race[,2] == 1, 1, 0))
  
  ### expand race
  racecatexp = rep(NA, n*obs)
  # for each subject, fill in the relevant race entries of the vector, repeating for all the observations
  for (i in 1:n ) {
    # number of entries already used: number of previous subjects * obs per subject
    entriesUsed = (i-1)*obs 
    racecatexp[ (entriesUsed + 1) : (entriesUsed + obs)] = racecat[i]
  }
  
  
  ### step 3.5 - combine race
  d4 <- cbind(d3, racecatexp)
  
  return(d4)
}



############################## FUNCTION: EXPAND SUBJECTS ##############################

#subjectMus: a matrix of subject means for each variable
#n.OtherNorms: number of non-drug normals
#wcorin: within-subject correlation matrix
#obs: number of observations per subject to generate

expandSubjects = function(mus3, n.OtherNorms, n.OtherBins, n.Drugs, wcorin, obs) {

  #number of subjects
  n = dim(mus3)[1]
  
  #total number of variables
  n.Vars = n.OtherNorms + n.OtherBins + n.Drugs
  
  dat <- vector("list", n)
  
  for (s in 1:n) {
    staticBins <- mus3[ s, 1:n.OtherBins ] #subset the matrix to get static (non-time-varying) binaries for subject s
    drugProbs <- mus3[ s, (n.OtherBins + n.OtherNorms + 1):n.Vars ] # subset the matrix to get binaries for subject s
    normMeans <- mus3[s, c( (n.OtherBins + 1) : (n.OtherBins + n.OtherNorms) )] # subset the matrix to get the non-drug normals for subject s
    
    zerodrugs <- which(drugProbs==zero) # which drugs have p = 0?
    
    # within-subject correlation matrix
    wcorin2 = as.matrix(wcorin)
    
    # change the correlations to 0s where pdrugs = zero
    for (r in zerodrugs){ 
      wcorin2[,r] <- 0
      wcorin2[r,] <- 0
    }
    
    # convert correlation matrix to vector of just upper-tri elements
    newwcorvec = upperTriVec(wcorin2) 
    
    #create a list with a dataset (of length obs) for each subject
    dat[[s]] <- mod.jointly.generate.binary.normal(no.rows=obs, no.bin=length(drugProbs), no.nor=length(normMeans),
                                                   prop.vec.bin=drugProbs, mean.vec.nor=normMeans,
                                                   var.nor=wvar, corr.vec=newwcorvec, adjust.corrs=T)
    
    
    # put non-time-varying binaries back in
    dat[[s]] = cbind( "gender" = staticBins, dat[[s]] )
    
  }
  
  dat <- do.call("rbind", dat)
  return(dat)
}



############################## FUNCTION: GENERATE DATA ##############################


### function to create time-varying covariates
#will need to take as arguments:
#wcorin
#obs
#n.Drugs

#n: number of subjects to generate
#obs: number of observations to generate per subject
#pps: vector of population proportions for binary variables, starting with drug indicators and ending with any non-drug binaries (gender)
#pmeans: vector of population means for normal variables, starting with non-drug normals and ending with drug propensities
#n.Drugs: number of drugs
#pcor: population correlation vector
#wcorin: population correlation matrix
# pvars: vector of variances for normal variables

gendata <- function(n, obs, pps, pmeans, n.Drugs, pcor, wcorin, pvars) {  
  
  ### step 0 - number of different types of variables
  n.BinVars = length(pps) # number binary variables
  n.NormVars = length(pmeans) # number normal variables
  n.OtherNorms = n.NormVars - n.Drugs # number of non-drug normal variables
  n.OtherBins = n.BinVars - n.Drugs # number of non-drug binary variables
  n.Vars = n.OtherBins + n.OtherNorms + n.Drugs #total number of variables in study
  
  ### step 1 - generate mu for each person
  mus0 <- mod.jointly.generate.binary.normal( no.rows=n, no.bin=n.BinVars, no.nor=n.NormVars,
                                              prop.vec.bin=pps, mean.vec.nor=pmeans,
                                              var.nor=pvars, corr.vec=pcor )
  
  ### step 1.1 - if drug indicator is 0, then convert probability of receiving drug to 0
  mus1 = overrideDrugProbs(mus0, n.Drugs)
  
  
  ### step 1.2 - set aside ever-use indicators and expand them
  everUser <- mus1[, ( 1:n.Drugs ) ]
  everUserExp <- expandMatrix(everUser, obs)
  
  ### step 1.25 - temporarily remove drug indicators from matrix
  mus2 <- mus1[, -c( 1:n.Drugs ) ]
  
  
  ### step 1.3 - "proportionize" normal drug variables (force them to be strictly between 0 and 1)
  bins = mus2[ , (n.OtherBins + n.OtherNorms + 1):dim(mus2)[2] ]  # just binaries
  bins.prop = proportionize(bins, zero, one)  # proportionized version of the binaries
  mus3 = cbind( mus2[, 1:(n.OtherBins + n.OtherNorms)], bins.prop )  # new matrix with the proportionized binaries
  
  ### step 2 - generate time-varying data for each person
  d1 <- expandSubjects(mus3, n.OtherNorms, n.OtherBins, n.Drugs, wcorin, obs)
 
  ### step 2.5 - add patient id and ever-use indicators
  id <- rep(1:n, each=obs)
  d2 <- as.data.frame( cbind(id, d1, everUserExp) ) 
  
  ### step 2.6 - dummy-code variables for race model
  d3 = add_dummy_vars(d2)
  
  ### step 3 - add race ###
  d4 = add_race(d3, n, obs)

  sim <- list("data" = d4, "everUser" = everUser)
  return(sim)
}



##### Function: longitudinally expand a matrix of single observations by subject
# repeat each subject's entry in each row for obs number of times

expandMatrix <- function(matrix, obs) {
  
  n <- nrow(matrix)
  expanded <- matrix(c(NA), nrow = n*obs, ncol = ncol(matrix) )
  
  for ( i in 1:nrow(expanded) ) {
    for (j in 1:ncol(expanded) ) {
      id <- ceiling(i/obs)  # which subject id corresponds to the ith row in the expanded matrix?
      expanded[i,j] <- matrix[id,j]
    }
  }
  print(expanded)
}

# example
( mat <- matrix( seq(1:10), nrow=2, byrow=F) )
expandMatrix(mat, 4)


######################### FUNCTION: ADD DUMMY VARIABLES FOR USE IN RACE MODEL #########################

add_dummy_vars <- function(d2) {
  
  names(d2) <- varnames
  
  library(car)
  
  # dummy code age 
  age_cat <- recode(d2$age, "0:35='a.Under35'; 35:45='b.35to45'; 45:55='c.45to55'; 55:65='d.55to65'; 65:120='e.Over65'")
  
  # dummy code bmi
  bmi_cat <- recode(d2$bmi, "0:20='b.Under20'; 20:25='a.20to25'; 25:30='c.25to30'; 30:100='d.Over30'")
  
  # dummy code cd4
  # >500 is the reference category
  #cd4_lt50 <- as.numeric(d$cd4 < 50)
  #cd4_50to100 <- as.numeric(d$cd4 >= 50 & d$cd4 < 100)
  #cd4_100to200 <- as.numeric(d$cd4 >= 100 & d$cd4 < 200)
  #cd4_200to350 <- as.numeric(d$cd4 >= 200 & d$cd4 < 350)
  #cd4_350to500 <- as.numeric(d$cd4 >= 200 & d$cd4 < 350)
  cd4_cat <- recode(d2$cd4, "0:50='a.Under50'" )  # temporary!
  
  # dummy vln
  vln_cat <- recode(d2$vln, "0:400='a.Under400'; 400:3500='b.400to3500'; 3500:10000='c.3500to10K';
                    10000:50000='d.10Kto50K'; 50000:300000='e.Over50K'")
  
  d3 = cbind(d2, age_cat, bmi_cat, vln_cat)
  
  return(d3)
}




######################### FUNCTION: CREATE PROPORTION-OF-TIME-ON-DRUG DATAFRAME #########################

# d: dataset
# n: number of subjects in dataset
# n.Drugs: number of drugs
# n.OtherBins: number of non-drug binary variables
# ever_use: a matrix generated from makeEverUse function.
#           Has n rows indicating whether each subject was ever on given drug.

make_prop_drug_time = function(sim, n, n.Drugs, n.OtherBins) {
 
  # extract dataframes
  d <- sim$data
  
  # proportion of time on drug among ever-users (pmeans)
  # initialize matrix
  prop_drug_time <- sim$everUser
  
  # iterate through subjects and drugs, filling in prop_drug_time matrix
  for (j in 1:n.Drugs) { # i: rows; j: columns
    for (i in 1:n) {
      
      # for subjects who never used the drug, set relevant entry to NA
      if (sim$everUser[i,j] == 0) {
        prop_drug_time[i,j] <- NA
      }
      
      else {
        # pull out drug indicator column for drug of interest
        drugIndicator <- d[, n.OtherBins + 1 + j]
        # fill in proportion of time subject i was on drug j (mean of their indicators)
        prop_drug_time[i,j] <- mean( drugIndicator[d$id == i] )
        
      }
    } 
  }
  return(prop_drug_time)
}



######################### FUNCTION: COMPUTE PERFORMANCE STATISTICS FOR A SINGLE DATASET ######################### 

# d: dataset
# n: number of subjects in dataset
# obs: number of observations per subject
# n.Drugs: number of drugs
# n.OtherBins: number of non-drug binary variables
# n.OtherBins: number of non-drug normal variables
# pmeans: vector of population means for normal variables, starting with non-drug normals and ending with drug propensities
# pmeans.target: vector of desired population means; default is same as pmeans (this lets you adjust for bias)
# pps: vector of population proportions for binary variables, starting with drug indicators and ending with any non-drug binaries (gender)
# pps.target: vector of desired population proportions; default is same as pps (this lets you adjust for bias)
# pvars: vector of variances for normal variables

gendata_performance = function(sim, n, obs, n.Drugs, n.OtherBins, n.OtherNorms, pmeans.target=pmeans, pps.target=pps, pvars) {
  
  #### Set Up ######

  # make proportion-of-time-on-drug dataframe
  prop_drug_time = make_prop_drug_time(sim, n, n.Drugs, n.OtherBins)
 
  # get first observations for each subject
  first.dat <- sim$data[!duplicated(sim$data$id),]
  
  # get subsets of matrix #CHANGED
  OtherNorms = first.dat[ , ( (2+n.OtherBins+n.Drugs) : (1+n.OtherBins+n.Drugs+n.OtherNorms) ) ]  # just the non-drug normals; use only first observation per subject
  #Drugs = d[ , (2+n.OtherBins) : (1+n.OtherBins+n.Drugs) ]  # just the drugs
  
  # set up list of values to return
  performance = make_result_list()
  
  
  ###### Performance Statistics: Proportion Male #####
  
  target = pps.target[n.Drugs + 1]  # population proportion
  n.Males = length(which(first.dat$male==1))  # number of males
  ours = n.Males / n  # sample proportion
  
  performance$prop.male$abs.bias = ours - target
  performance$prop.male$std.bias = (ours - target) / sqrt(target*(1-target)/n)
  
  CIbounds = prop.test(n.Males, n)$conf.int
  performance$prop.male$coverage = target < CIbounds[2] & target > CIbounds[1]  # coverage is 1 if target is between the CIbounds
  
  
  ###### Performance Statistics: Normal Variables ######
  
  target = pmeans.target[1:n.OtherNorms]
  ours = apply(OtherNorms, 2, mean)
  trueSEs = sqrt(pvars[1:n.OtherNorms])
  
  performance$nor.vars$std.bias = (ours-target)/trueSEs
  performance$nor.vars$abs.bias = (ours-target)
  
  # coverage
  CIbounds = apply( OtherNorms, 2, function(x) t.test(x)$conf.int )
  performance$nor.vars$coverage = target < CIbounds[2,] & target > CIbounds[1,] # coverage is 1 if target is between the CIbounds 
  
  
  ###### Performance Statistics: Drug Ever-Use Variables ######
 
  # proportion of ever-users for each drug (pps)
  target = pps.target[1:n.Drugs]
  ours = apply(sim$everUser, 2, mean)
  trueSEs = sqrt( ( target * (1 - target) ) / n ) #CHANGED
  
  performance$drug.evers$abs.bias = ours - target
  performance$drug.evers$std.bias = (ours - target) / trueSEs
  
  n.EverUsers = apply(sim$everUser, 2, sum)  # number of people ever on each drug
  CIbounds = vapply( n.EverUsers, function(x) prop.test(x, n)$conf.int, c(0,0) )  # get CI for proportion ever-users
  performance$drug.evers$coverage = ( target < CIbounds[2,] ) & ( target > CIbounds[1,] )  # coverage is TRUE if target is between the CIbounds
  
  
  
  ###### Performance Statistics: Proportion Drug Time Variables ######
  
  # proportion of time on drug among ever-users (pmeans)
  target = pmeans.target[ (n.OtherNorms + 1) : length(pmeans.target)]
  ours = apply(prop_drug_time, 2, function(x) mean(x, na.rm=T) ) # compute mean of proportion drug time among ever-users
  sampleSDs = apply(prop_drug_time, 2, function(x) sd(x, na.rm=T) )
  
  # compute sample SEs of mean
  sampleSEs = sampleSDs / sqrt(n)  #CHANGED
  
  performance$prop.drug.time$abs.bias = ours-target
  performance$prop.drug.time$std.bias = (ours-target) / sampleSEs
  
  # coverage
  halfWidth = qt(0.975, df=n-1) * sampleSEs
  
  upperBounds = ours + halfWidth
  lowerBounds = ours - halfWidth
  performance$prop.drug.time$coverage = (target < upperBounds) & (target > lowerBounds) # coverage is 1 if target is between the CIbounds 
  
  
  ###### Performance Statistics: Race Summary Statistics ######
  
  performance$race$prop.table = prop.table( table( first.dat$racecatexp ) )
  
  
  ###### Performance Statistics: Population Correlation Matrix ######
  
  performance$race$prop.table = prop.table( table( first.dat$racecatexp ) )
  
  
  
  
  ###### Performance Statistics: Race Model Coefficients ######
  
  # fit multinomial models for races 1 and 2 vs. 0
  #library(glm2)
  
  # this will need to be updated! 
  # mod.race1 <- glm( (racecatexp==1) ~  ever_abac + ever_ataz + ever_dida + ever_efav + ever_emtr + ever_indi + ever_lami + ever_lopi + ever_nelf + ever_nevi +
  #       ever_rito + ever_saqu + ever_stav + ever_teno + ever_zido + male + age_cat + bmi_cat + bpd + bps + hdl + ldl + trig, 
  #     data=first.dat, family=binomial(link="logit"))
  
  # mod.race2 <- glm( (racecatexp==2) ~  ever_abac + ever_ataz + ever_dida + ever_efav + ever_emtr + ever_indi + ever_lami + ever_lopi + ever_nelf + ever_nevi +
  #                    ever_rito + ever_saqu + ever_stav + ever_teno + ever_zido + male + age_cat + bmi_cat + bpd + bps + hdl + ldl + trig, 
  #                  data=first.dat, family=binomial(link="logit"))
  
  
  return(performance)
  
}



######################### FUNCTION: RUN SEVERAL SIMULATIONS AND RETURN RESULTS LIST ######################### 

###### This is the function that the user would call. ######


#n: number of subjects to generate
#obs: number of observations to generate per subject
#pps: vector of population proportions for binary variables, starting with drug indicators and ending with any non-drug binaries (gender)
# pps.target: vector of desired population proportions; default is same as pps (this lets you adjust for bias)
#pmeans: vector of population means for normal variables, starting with non-drug normals and ending with drug propensities
# pmeans.target: vector of desired population means; default is same as pmeans (this lets you adjust for bias)
#n.Drugs: number of drugs
#pcor: population correlation vector
#wcorin: population correlation matrix
# pvars: vector of variances for normal variables
# n.Reps: number of datasets to generate
# varnames: vector of variable names
# race.names: names of races
# write.data: should R write all the generated datasets to csv files?


repeat_sim <- function(n, obs, pps, pps.target=pps, pmeans, pmeans.target=pmeans, n.Drugs, pcor, wcorin, pvars, n.Reps, varnames, race.names, write.data=FALSE) {
  
  ##### compute numbers of different types of variables #####
  n.BinVars = length(pps) # number binary variables
  n.NormVars = length(pmeans) # number normal variables
  n.OtherNorms = n.NormVars - n.Drugs # number of non-drug normal variables
  n.OtherBins = n.BinVars - n.Drugs # number of non-drug binary variables
  
  ##### initialize results list #####
  results = make_result_list()
  
  ##### simulate data n.Reps times, adding each entry to results list #####

  for (i in 1:n.Reps) {
    sim <- gendata(n, obs, pps, pmeans, n.Drugs, pcor, wcorin, pvars)
    newEntry <- gendata_performance(sim, n, obs, n.Drugs, n.OtherBins, n.OtherNorms, pmeans=pmeans, pmeans.target=pmeans.target, pps=pps, pps.target=pps.target, pvars=pvars)
    
    # add the new entry as a new "row" in the results list
    results <- Map( function(x,y) Map(rbind, x, y), results, newEntry )
    
    # optionally, write the dataset to a csv file in current working directory
    file.name = paste(Sys.Date(), "_dataset_", i, sep="" )
    write.csv( dput(sim), file.name )
  }
  
  
  ##### put in variable names  #####
  normal.names <- varnames[ (n.OtherBins + n.Drugs + 2) : (n.OtherBins + n.Drugs + 1 + n.OtherNorms) ]
  other.bin.names <- varnames[ 2 : (1 + n.OtherBins) ]
  drug.names <- varnames[ (n.OtherBins + 2) : (n.OtherBins + 1 + n.Drugs) ]
  
  names(results$prop.male$abs.bias) = other.bin.names
  names(results$prop.male$std.bias) = other.bin.names
  names(results$prop.male$coverage) = other.bin.names
  
  names(results$nor.vars$abs.bias) = normal.names
  names(results$nor.vars$std.bias) = normal.names
  names(results$nor.vars$coverage) = normal.names
  
  names(results$drug.evers$abs.bias) = drug.names
  names(results$drug.evers$std.bias) = drug.names
  names(results$drug.evers$coverage) = drug.names
  
  names(results$prop.drug.time$abs.bias) = drug.names
  names(results$prop.drug.time$std.bias) = drug.names
  names(results$prop.drug.time$coverage) = drug.names
  
  names(results$race$prop.table) = race.names
  
  return(results)
}



######################### FUNCTION: COMPUTE MEAN PERFORMANCE STATS ######################### 

# results: a list of performance results returned by repeat_sim

mean_performance = function(results) {
  
  # initialize mean results list
  mean_results = make_result_list()
                      
  for ( i in 1:length(results) ) {
    
    for (j in 1:length(results[[i]] )) {
      
      avg = apply(results[[i]][[j]], 2, mean)  # compute the average across all simulation reps
      mean_results[[i]][[j]] = avg  # put the averages in the mean_results list
      
    }  
  }
  return(mean_results)
}



######################### FUNCTION: ADJUST POPULATION PARAMETERS GIVEN ESTIMATED BIAS #########################

# only adjust drug variables, not normals or gender
# adjust population parameters by the estimated (observed) bias

adjust_parameters = function(pmeans, pps, mean_results, n.OtherNorms, n.OtherBins) {
  
  # adjust drug propensity variables
  bias = as.vector( mean_results$prop.drug.time$abs.bias )
  adjustment = c( rep(0, n.OtherNorms), bias)  # no adjustment (0) for normal variables
  pmeans.adj = pmeans - adjustment
  pmeans.adj[pmeans.adj <= 0] <- zero  # since these are proportions, replace any negative values with zero
  
  # adjust drug ever-use variables
  bias = as.vector( mean_results$drug.evers$abs.bias )
  adjustment = c( bias, rep(0, n.OtherBins ))  # no adjustment (0) for non-drug binaries
  pps.adj = pps - adjustment
  pps.adj[pps.adj <= 0] <- zero  # since these are proportions, replace any negative values with zero
  
  return( list("pmeans.adj" = pmeans.adj, "pps.adj" = pps.adj) )
}




######################### FUNCTIONS TO PLOT PERFORMANCE STATS WITH 95% CIs #########################


plotBias = function(bias, x.names, yaxp, ylim, pch=20, xlab="", ylab="Bias", main="", refline=0, col="red") {
  
  # were multiple simulations run over which to average? 
  multiple = !is.null( dim(bias) )
  
  # compute mean bias across all simulations
  # if there was only 1 simulation, just use that
  values = if (multiple) apply(bias, 2, mean) else bias
  
  # make plot
  plot(values, pch=pch, xlab=xlab, ylab=ylab, main=main, yaxp=yaxp, ylim=ylim, xaxt="n", type="b", col=col)
  axis(1, at=1:length(bias), labels=x.names)
  abline(h=refline, lty=2)
  
  # empirical 95% CI based on quantiles of sample
  # applies only if n.Reps > 1, aka bias has non-null dimensions
  if (multiple) {
    loBound <- apply(bias, 2, function(x) quantile(x, 0.025))
    hiBound <- apply(bias, 2, function(x) quantile(x, 0.975))
    arrows(x0 = seq(1:length(values)), y0=loBound, y1=hiBound, angle=90, code=3, length=0.05)
  }
}

plotCoverage = function(coverage, n.Reps, x.names, pch=20, xlab="", ylab="Coverage", main="", refline=0.95, col="red") {
  
  # were multiple simulations run over which to average? 
  multiple = n.Reps > 1
  
  # compute mean coverage across all simulations
  # if there was only 1 simulation, just use that
  values = if (multiple) apply(coverage, 2, mean) else coverage
  
  # make plot
  plot(values, ylim=c(0,1), yaxp=c(0,1,10), pch=pch, xlab=xlab, ylab=ylab, main=main, xaxt="n", type="b", col=col)
  axis(1, at=1:length(coverage), labels=x.names)
  abline(h=refline, lty=2)
  
  # 95% CI bars based on binomial distribution and sample coverage proportion
  # applies only if n.Reps > 1, aka bias has non-null dimensions
  if (multiple) {
    SE <- sqrt( (values * (1-values) ) / n.Reps )
    loBound <- values - ( SE * qnorm(0.975) )  # CI bounds
    hiBound <- values + ( SE * qnorm(0.975) ) 
    # error bars
    arrows(x0 = seq(1:length(values)), y0=loBound, y1=hiBound, angle=90, code=3, length=0.05)
  }
}



######################### FUNCTION: INITIALIZE EMPTY RESULTS LIST #########################

make_result_list = function() {
  empty.list <- list( "prop.male" = list("abs.bias"=data.frame(), "std.bias"=data.frame(), "coverage"=data.frame()),
              "nor.vars" = list("abs.bias"=data.frame(), "std.bias"=data.frame(), "coverage"=data.frame()), 
              "drug.evers" = list("abs.bias"=data.frame(), "std.bias"=data.frame(), "coverage"=data.frame()),
              "prop.drug.time" = list("abs.bias"=data.frame(), "std.bias"=data.frame(), "coverage"=data.frame()),
              "race" = list("prop.table"=data.frame())
  ) 
  return(empty.list)
}

