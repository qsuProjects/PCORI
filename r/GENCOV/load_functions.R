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

# load packages
library(MBESS)
library(mvtnorm)
library(ICC)
library(miscTools)
library(car)

# load jointly generate function
# source("jointly_generate_binary_normal_modified_v2.R")


############################## FUNCTION: OVERRIDE DRUG PROBS ##############################

# For subjects who were never on a given drug, changes the corresponding drug propensity to zero.
# overrides drug probabilities if indicator = 0

# mus0: matrix of subject means
# n.Drugs: number of drugs

override_drug_probs = function(mus0, n.Drugs) {
  # number of non-drug variables
  n.NonDrugVars = dim(mus0)[2] - n.Drugs
  
  for (m in 1:n.Drugs) {
    # replace each entry in the part of mus0 corresponding to the drug propensities with 0...
    # ...if the part of mus0 corresponding to the drug indicator is 0
    # otherwise leave the drug propensity alone
    mus0[, n.Drugs + n.OtherBins + m] = ifelse( mus0[,m] == 0, zero, mus0[, n.Drugs + n.OtherBins + m])
  }
  
  return(mus0)
}


############################## FUNCTION: UPPER TRI VEC ##############################

# turns matrix into vector of upper-triangular elements
# matrix elements are arranged by row

upper_tri_vec = function(m) {
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
add_race = function(d3, n, obs) {
  
  # create dataframe with just the first observation from each patient
  first.dat = d3[!duplicated(d3$id),]
      
  # estimate logits for black and other-race
  logit.race.b = bb0 + bb1*first.dat$ever_abac + bb2*first.dat$ever_ataz + bb3*first.dat$ever_dida + bb4*first.dat$ever_efav + bb5*first.dat$ever_emtr +
    bb6*first.dat$ever_indi + bb7*first.dat$ever_lami + bb8*first.dat$ever_lopi + bb9*first.dat$ever_nelf + bb10*first.dat$ever_nevi +
    bb11*first.dat$ever_rito + bb12*first.dat$ever_saqu + bb13*first.dat$ever_stav + bb14*first.dat$ever_teno + bb15*first.dat$ever_zido +
    bb16*first.dat$male + 
    bb17*(first.dat$age_cat=="b.35to45") + bb18*(first.dat$age_cat=="c.45to55") + bb19*(first.dat$age_cat=="d.55to65") + bb20*(first.dat$age_cat=="e.Over65") +  
    bb21*(first.dat$bmi_cat=="b.Under20") + bb22*(first.dat$bmi_cat=="c.25to35") + bb23*(first.dat$bmi_cat=="d.Over30") + 
    bb33*first.dat$bpd + bb34*first.dat$bps + bb35*first.dat$hdl + bb36*first.dat$ldl + bb37*first.dat$trig
  
  logit.race.o = bo0 + bo1*first.dat$ever_abac + bo2*first.dat$ever_ataz + bo3*first.dat$ever_dida + bo4*first.dat$ever_efav + bo5*first.dat$ever_emtr +
    bo6*first.dat$ever_indi + bo7*first.dat$ever_lami + bo8*first.dat$ever_lopi + bo9*first.dat$ever_nelf + bo10*first.dat$ever_nevi +
    bo11*first.dat$ever_rito + bo12*first.dat$ever_saqu + bo13*first.dat$ever_stav + bo14*first.dat$ever_teno + bo15*first.dat$ever_zido +
    bo16*first.dat$male + 
    bo17*(first.dat$age_cat=="b.35to45") + bo18*(first.dat$age_cat=="c.45to55") + bo19*(first.dat$age_cat=="d.55to65") + bo20*(first.dat$age_cat=="e.Over65") +  
    bo21*(first.dat$bmi_cat=="b.Under20") + bo22*(first.dat$bmi_cat=="c.25to35") + bo23*(first.dat$bmi_cat=="d.Over30") + 
    bo33*first.dat$bpd + bo34*first.dat$bps + bo35*first.dat$hdl + bo36*first.dat$ldl + bo37*first.dat$trig
  
  # probability of being each race
  # each is a vector where the jth entry is the probability for the jth person
  p.race1 = exp(logit.race.b)/(1+exp(logit.race.b)+exp(logit.race.o))
  p.race2 = exp(logit.race.o)/(1+exp(logit.race.b)+exp(logit.race.o))
  p.race0 = 1/(1+exp(logit.race.b)+exp(logit.race.o))
  #p.race0 = ifelse(p.race0 < 0 ,0, p.race0)
  
  # create vector of lists for each person's race
  race = vector("list", n)
  
  # generate each person's race multinomially
  for (j in 1:n)  race[[j]] = rmultinom(1, 1, c(p.race0[j], p.race1[j], p.race2[j]))
  
  race = t(do.call("cbind", race))
  
  # recode race
  # if the third col is 1, race=2
  # if the second col is 1, race=1
  # if the first col is 1, race=0
  # the i,jth entry of race now corresponds to whether person i was of race j
  racecat = ifelse(race[,3] == 1, 2, ifelse(race[,2] == 1, 1, 0))
  
  ### expand race
  racecatexp = rep(NA, n*obs)
  # for each subject, fill in the relevant race entries of the vector, repeating for all the observations
  for (i in 1:n ) {
    # number of entries already used: number of previous subjects * obs per subject
    entriesUsed = (i-1)*obs 
    racecatexp[ (entriesUsed + 1) : (entriesUsed + obs)] = racecat[i]
  }
  
  
  ### step 3.5 - combine race
  d4 = cbind(d3, racecatexp)
  
  return(d4)
}



############################## FUNCTION: EXPAND SUBJECTS ##############################

# subjectMus: a matrix of subject means for each variable
# n.OtherNorms: number of non-drug normals
# wcorin: within-subject correlation matrix
# obs: number of observations per subject to generate

expand_subjects = function(mus3, n.OtherNorms, n.OtherBins, n.Drugs, wcorin, obs) {
  
  # number of subjects
  n = dim(mus3)[1]
  
  # total number of variables
  n.Vars = n.OtherNorms + n.OtherBins + n.Drugs
  
  dat = vector("list", n)
  
  for (s in 1:n) {
    staticBins = mus3[ s, 1:n.OtherBins ] #subset the matrix to get static (non-time-varying) binaries for subject s
    drugProbs = mus3[ s, (n.OtherBins + 1):(n.OtherBins + n.Drugs) ] # subset the matrix to get binaries for subject s
    
    ### FOR DEBUGGING ONLY ###
    if ( ( min(drugProbs)<=0 ) | ( max(drugProbs)>=1 ) ) {
      string = paste("SUBJECT", s, "HAD ILLEGAL DRUG PROBS:")
       print(string)
      print(drugProbs)
    }
    
    
    normMeans = mus3[ s, (n.OtherBins + n.Drugs + 1) : ncol(mus3) ] # subset the matrix to get the non-drug normals for subject s
    
    zerodrugs = which(drugProbs==zero) # which drugs have p = 0?
    
    # within-subject correlation matrix
    wcorin2 = as.matrix(wcorin)
    
    # change the correlations to 0s where pdrugs = zero
    # SHOULD WE ALSO DO THIS WHEN PDRUGS = 1? DOES THIS EVEN HAPPEN?
    for (r in zerodrugs){ 
      wcorin2[,r] = 0
      wcorin2[r,] = 0
    }
    
    # convert correlation matrix to vector of just upper-tri elements
    newwcorvec = upper_tri_vec(wcorin2) 
    
    # create a list with a dataset (of length obs) for each subject
    dat[[s]] = mod.jointly.generate.binary.normal(no.rows=obs, no.bin=length(drugProbs), no.nor=length(normMeans),
                                                   prop.vec.bin=drugProbs, mean.vec.nor=normMeans,
                                                   var.nor=parameters$within.var[parameters$type == "normal.other" | 
                                                                                   parameters$type == "time.function"], 
                                                   corr.vec=newwcorvec, adjust.corrs=T)
 
    # put non-time-varying binaries back in
    dat[[s]] = cbind( staticBins, dat[[s]] )
  }
  
  dat = do.call("rbind", dat)
  return(dat)
}



############################## FUNCTION: GENERATE DATA ##############################

# n: number of subjects to generate
# obs: number of observations to generate per subject
# parameters: parameters dataframe
# n.Drugs: number of drugs
# pcor: across-subject correlation matrix
# wcorin: within-subject correlation matrix

make_one_dataset = function(n, obs, parameters, n.Drugs, pcor, wcorin) {  

  ### step 0.1 - extract parameter vectors from given dataframe
  bin.props = parameters$prop[ parameters$type == "bin.other" | parameters$type == "bin.drug" ]  # = bin.props
  nor.means = parameters$across.mean[ parameters$type %in% c("normal.drug", "normal.other", "time.function") ]  # = nor.means
  across.vars = parameters$across.var[ parameters$type %in% c("normal.drug", "normal.other", "time.function") ]  # = across.vars

  ### step 0.2 - number of different types of variables
  n.BinVars = length( bin.props )  # number binary variables (16)
  n.NormVars = length( nor.means )  # number normal variables
  n.OtherNorms = n.NormVars - n.Drugs  # number of non-drug normal variables
  n.OtherBins = n.BinVars - n.Drugs  # number of non-drug binary variables
  n.Vars = n.OtherBins + n.OtherNorms + n.Drugs  # total number of variables in study (not double-counting drugs)
  
  ### step 1.0 - generate mu for each person
  mus0 = mod.jointly.generate.binary.normal( no.rows = n, no.bin = n.BinVars, no.nor = n.NormVars,
                                              prop.vec.bin = bin.props, 
                                              mean.vec.nor = nor.means,
                                              var.nor = across.vars, corr.vec = pcor )
  
  ### step 1.1 - if drug indicator is 0, then convert probability of receiving drug to 0
  mus1 = override_drug_probs(mus0, n.Drugs)
  
  ### step 1.2 - set aside ever-use indicators and expand them
  everUser = mus1[, ( 1:n.Drugs ) ]
  everUserExp = expand_matrix(everUser, obs)
  
  ### step 1.3 - temporarily remove drug indicators from matrix
  mus2 = mus1[, -c( 1:n.Drugs ) ]
  
  
  ### step 1.4 - "proportionize" normal drug variables (force them to be strictly between 0 and 1)
  bins = mus2[ , (n.OtherBins + 1):(n.OtherBins + n.Drugs) ]  # just binaries
  bins.prop = proportionize(bins, zero, one)  # proportionized version of the binaries
  mus3 = mus2
  mus3[ , (n.OtherBins + 1):(n.OtherBins + n.Drugs) ] = bins.prop
  
  
  ### step 2 - generate time-varying data for each person
  d1 = expand_subjects(mus3, n.OtherNorms, n.OtherBins, n.Drugs, wcorin, obs)
 
  ### step 3 - add patient id, ever-use indicators, and variable names
  id = rep(1:n, each=obs)
  d2 = as.data.frame( cbind(id, d1, everUserExp) ) 
  
  names = c( "id", as.character(parameters$name[parameters$type=="bin.other"]), 
             as.character(parameters$name[parameters$type=="normal.drug"]), 
            as.character(parameters$name[parameters$type %in% c("normal.other", "time.function")]), 
            as.character(parameters$name[parameters$type=="bin.drug"]))
  names(d2) = names
  
  ### step 3.1 - dummy-code variables for race model
  d3 = add_dummy_vars(d2)

  ### step 3.2 - add race ###
  d4 = add_race(d3, n, obs)
  
  ### step 3.3 - add time-function variables ###
  d5 = add_time_function_vars(d4, obs, parameters)
    
  sim = list("data" = d5, "everUser" = everUser)
  return(sim)
}



##### Function: longitudinally expand a matrix of single observations by subject
# repeat each subject's entry in each row for obs number of times

expand_matrix = function(matrix, obs) {
  
  n = nrow(matrix)
  expanded = matrix(c(NA), nrow = n*obs, ncol = ncol(matrix) )
  
  for ( i in 1:nrow(expanded) ) {
    for (j in 1:ncol(expanded) ) {
      id = ceiling(i/obs)  # which subject id corresponds to the ith row in the expanded matrix?
      expanded[i,j] = matrix[id,j]
    }
  }
  print(expanded)
}

# example
#( mat = matrix( seq(1:10), nrow=2, byrow=F) )
#expand_matrix(mat, 4)


######################### FUNCTION: ADD DUMMY VARIABLES FOR USE IN RACE MODEL #########################

add_dummy_vars = function(d2) {
  
  # dummy code age 
  age_cat = recode(d2$age, "0:35='a.Under35'; 35:45='b.35to45'; 45:55='c.45to55'; 55:65='d.55to65'; 65:120='e.Over65'")
  
  # dummy code bmi
  bmi_cat = recode(d2$bmi, "0:20='b.Under20'; 20:25='a.20to25'; 25:30='c.25to30'; 30:100='d.Over30'")
  
  # dummy code cd4
  # >500 is the reference category
  #cd4_lt50 = as.numeric(d$cd4 < 50)
  #cd4_50to100 = as.numeric(d$cd4 >= 50 & d$cd4 < 100)
  #cd4_100to200 = as.numeric(d$cd4 >= 100 & d$cd4 < 200)
  #cd4_200to350 = as.numeric(d$cd4 >= 200 & d$cd4 < 350)
  #cd4_350to500 = as.numeric(d$cd4 >= 200 & d$cd4 < 350)
  cd4_cat = recode(d2$log_cd4, "0:50='a.Under50'" )  # temporary!
  
  # dummy vln
  vln_cat = recode(d2$log_vln, "0:400='a.Under400'; 400:3500='b.400to3500'; 3500:10000='c.3500to10K';
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
  d = sim$data
  
  # proportion of time on drug among ever-users (nor.means)
  # initialize matrix
  prop_drug_time = sim$everUser
  
  # iterate through subjects and drugs, filling in prop_drug_time matrix
  for (j in 1:n.Drugs) { # i: rows; j: columns
    for (i in 1:n) {
      
      # for subjects who never used the drug, set relevant entry to NA
      if (sim$everUser[i,j] == 0) {
        prop_drug_time[i,j] = NA
      }
      
      else {
        # pull out drug indicator column for drug of interest
        drugIndicator = d[, n.OtherBins + 1 + j]
        # fill in proportion of time subject i was on drug j (mean of their indicators)
        prop_drug_time[i,j] = mean( drugIndicator[d$id == i] )
        
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
# mean target: vector of desired population means; default is same as mean.target (this lets you adjust for bias)
# bin.props: vector of population proportions for binary variables, starting with drug indicators and ending with any non-drug binaries (gender)
# across.var.target: vector of desired population proportions; default is same as across.vars (this lets you adjust for bias)
# across.vars: vector of variances for normal variables

dataset_performance = function(sim, n, obs, n.Drugs, n.OtherBins, n.OtherNorms, 
                               mean.target, prop.target, bin.props, nor.means, across.vars, var.names) {
  
  #### Set Up ######
  # if proportion and mean performance targets aren't set, use the parameters around which we're actually generating
  if ( is.null(prop.target) ) prop.target=bin.props
  if ( is.null(mean.target) ) mean.target=nor.means
  
  # make proportion-of-time-on-drug dataframe
  prop_drug_time = make_prop_drug_time(sim, n, n.Drugs, n.OtherBins)
 
  # get first observations for each subject
  first.dat = sim$data[!duplicated(sim$data$id),]
  
  # get subsets of matrix
  OtherNorms = first.dat[ , ( (2+n.OtherBins+n.Drugs) : (n.OtherBins+n.Drugs+n.OtherNorms+1) ) ]  # just the non-drug normals; use only first observation per subject
  #Drugs = d[ , (2+n.OtherBins) : (1+n.OtherBins+n.Drugs) ]  # just the drugs
  
  # set up list of values to return
  performance = make_result_list()
  
   
  
  ###### Performance Statistics: Proportion Male #####
  target = prop.target[n.Drugs + 1]  # extract proportion males
  n.Males = length(which(first.dat$male==1))  # number of males
  ours = n.Males / n  # sample proportion
  
  performance$prop.male$abs.bias = ours - target
  performance$prop.male$std.bias = (ours - target) / sqrt(target*(1-target)/n)
  
  CIbounds = prop.test(n.Males, n)$conf.int
  performance$prop.male$coverage = ( target < CIbounds[2] ) & ( target > CIbounds[1] )  # coverage is 1 if target is between the CIbounds
  
  
  ###### Performance Statistics: Normal Variables ######
  target = mean.target[ (n.Drugs + 1) : length(mean.target) ]
  ours = apply(OtherNorms, 2, mean)
  trueSEs = sqrt(across.vars[1:n.OtherNorms])
  
  performance$other.normals$std.bias = (ours-target)/trueSEs
  performance$other.normals$abs.bias = (ours-target)
  
  # coverage
  CIbounds = apply( OtherNorms, 2, function(x) t.test(x)$conf.int )
  performance$other.normals$coverage = target < CIbounds[2,] & target > CIbounds[1,] # coverage is 1 if target is between the CIbounds 
  
  
  ###### Performance Statistics: Drug Ever-Use Variables ######
  # proportion of ever-users for each drug
  target = prop.target[1:n.Drugs]
  ours = apply(sim$everUser, 2, mean)
  trueSEs = sqrt( ( target * (1 - target) ) / n )
  
  performance$drug.evers$abs.bias = ours - target
  performance$drug.evers$std.bias = (ours - target) / trueSEs
  
  n.EverUsers = apply(sim$everUser, 2, sum)  # number of people ever on each drug
  CIbounds = vapply( n.EverUsers, function(x) prop.test(x, n)$conf.int, c(0,0) )  # get CI for proportion ever-users
  performance$drug.evers$coverage = ( target < CIbounds[2,] ) & ( target > CIbounds[1,] )  # coverage is TRUE if target is between the CIbounds
  
  
  
  ###### Performance Statistics: Proportion Drug Time Variables ######
  
  # proportion of time on drug among ever-users (nor.means)
  target = mean.target[1:n.Drugs]
  ours = apply(prop_drug_time, 2, function(x) mean(x, na.rm=T) ) # compute mean of proportion drug time among ever-users
  sampleSDs = apply(prop_drug_time, 2, function(x) sd(x, na.rm=T) )
  
  # compute sample SEs of mean
  sampleSEs = sampleSDs / sqrt(n) 
  
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
  # mod.race1 = glm( (racecatexp==1) ~  ever_abac + ever_ataz + ever_dida + ever_efav + ever_emtr + ever_indi + ever_lami + ever_lopi + ever_nelf + ever_nevi +
  #       ever_rito + ever_saqu + ever_stav + ever_teno + ever_zido + male + age_cat + bmi_cat + bpd + bps + hdl + ldl + trig, 
  #     data=first.dat, family=binomial(link="logit"))
  
  # mod.race2 = glm( (racecatexp==2) ~  ever_abac + ever_ataz + ever_dida + ever_efav + ever_emtr + ever_indi + ever_lami + ever_lopi + ever_nelf + ever_nevi +
  #                    ever_rito + ever_saqu + ever_stav + ever_teno + ever_zido + male + age_cat + bmi_cat + bpd + bps + hdl + ldl + trig, 
  #                  data=first.dat, family=binomial(link="logit"))
  
  return(performance)
}



######################### WRAPPER FUNCTION: RUN SEVERAL SIMULATIONS AND RETURN RESULTS LIST ######################### 

###### This is the function that the user would call. ######


#n: number of subjects to generate
#obs: number of observations to generate per subject
# props: vector of population proportions for binary variables, starting with drug indicators and ending with any non-drug binaries (gender)
# prop.target: vector of desired population proportions; default is same as bin.props (this lets you adjust for bias)
# means: vector of population means for normal variables, starting with non-drug normals and ending with drug propensities
# mean.target: vector of desired population means; default is same as nor.means (this lets you adjust for bias)
# across.vars: vector of variances for normal variables

#n.Drugs: number of drugs
#pcor: population correlation vector
#wcorin: population correlation matrix

# n.Reps: number of datasets to generate
# varnames: vector of variable names
# race.names: names of races
# write.data: should R write all the generated datasets to csv files?


repeat_sim = function(n, obs, parameters, prop.target = NULL, mean.target = NULL, n.Drugs, 
                       pcor, wcorin, n.Reps, race.names, write.data=FALSE, name_prefix) {
  
  ##### extract parameters from parameter matrix #####
  bin.props = parameters$prop[parameters$type == "bin.other" | parameters$type == "bin.drug"]  # = bin.props
  nor.means = parameters$across.mean[ parameters$type %in% c("normal.drug", "normal.other", "time.function") ]  # = nor.means
  across.vars = parameters$across.var[ parameters$type %in% c("normal.drug", "normal.other", "time.function") ]  # = across.vars
  
  #### extract variable names from parameter matrix ####
  other.bin.names = parameters$name[parameters$type=="bin.other"]
  normal.names = parameters$name[ parameters$type %in% c("normal.other", "time.function") ]
  drug.ever.names = parameters$name[parameters$type=="bin.drug"]
  drug.prop.names = parameters$name[parameters$type=="normal.drug"]

  ### if proportion and mean performance targets aren't set, use the parameters around which we're actually generating
  if ( is.null(prop.target) ) prop.target=bin.props
  if ( is.null(mean.target) ) mean.target=nor.means
  
  ##### compute numbers of different types of variables #####
  n.BinVars = length(bin.props) # number binary variables
  n.NormVars = length(nor.means) # number normal variables
  n.OtherNorms = n.NormVars - n.Drugs # number of non-drug normal variables
  n.OtherBins = n.BinVars - n.Drugs # number of non-drug binary variables
  
  ##### initialize results list #####
  results = make_result_list()

  ##### simulate data n.Reps times, adding each entry to results list #####
  for (i in 1:n.Reps) {
    sim = make_one_dataset(n, obs, parameters, n.Drugs, pcor, wcorin)

    newEntry = dataset_performance(sim, n, obs, n.Drugs, n.OtherBins, n.OtherNorms, 
                                   mean.target, prop.target, bin.props, nor.means, across.vars, var.names)
    
    # add the new entry as a new "row" in the results list
    results = Map( function(x,y) Map(rbind, x, y), results, newEntry )
    
    # optionally, write the dataset to a csv file in current working directory
    if(write.data) {
      file.name = paste(Sys.Date(), name_prefix, "dataset", i, sep="_" )
      write.csv( sim$data, file.name )
    }
  }
  
  
  # 
  
  ##### put in variable names  #####
  names(results$prop.male$abs.bias) = other.bin.names
  names(results$prop.male$std.bias) = other.bin.names
  names(results$prop.male$coverage) = other.bin.names
  
  names(results$other.normals$abs.bias) = normal.names
  names(results$other.normals$std.bias) = normal.names
  names(results$other.normals$coverage) = normal.names
  
  names(results$drug.evers$abs.bias) = drug.ever.names
  names(results$drug.evers$std.bias) = drug.ever.names
  names(results$drug.evers$coverage) = drug.ever.names
  
  names(results$prop.drug.time$abs.bias) = drug.prop.names
  names(results$prop.drug.time$std.bias) = drug.prop.names
  names(results$prop.drug.time$coverage) = drug.prop.names
  
  # NOTE: THIS LINE GIVES DIMENSION-MATCH ERROR IF WE DON'T GENERATE PEOPLE OF ALL 3 RACES
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

adjust_parameters = function(nor.means, bin.props, mean_results, n.OtherNorms, n.OtherBins) {
  
  # adjust drug propensity variables
  bias = as.vector( mean_results$prop.drug.time$abs.bias )
  adjustment = c( rep(0, n.OtherNorms), bias)  # no adjustment (0) for normal variables
  nor.means.adj = nor.means - adjustment
  nor.means.adj[nor.means.adj <= 0] = zero  # since these are proportions, replace any negative values with zero
  
  # adjust drug ever-use variables
  bias = as.vector( mean_results$drug.evers$abs.bias )
  adjustment = c( bias, rep(0, n.OtherBins ))  # no adjustment (0) for non-drug binaries
  bin.props.adj = bin.props - adjustment
  bin.props.adj[bin.props.adj <= 0] = zero  # since these are proportions, replace any negative values with zero
  
  return( list("nor.means.adj" = nor.means.adj, "bin.props.adj" = bin.props.adj) )
}




######################### FUNCTIONS TO PLOT PERFORMANCE STATS WITH 95% CIs #########################


plot_bias = function(bias, x.names, yaxp, ylim, pch=20, xlab="", ylab="Bias", main="", refline=0, col="red") {
  
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
    loBound = apply(bias, 2, function(x) quantile(x, 0.025))
    hiBound = apply(bias, 2, function(x) quantile(x, 0.975))
    arrows(x0 = seq(1:length(values)), y0=loBound, y1=hiBound, angle=90, code=3, length=0.05)
  }
}

plot_coverage = function(coverage, n.Reps, x.names, pch=20, xlab="", ylab="Coverage", main="", refline=0.95, col="red") {
  
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
    SE = sqrt( (values * (1-values) ) / n.Reps )
    loBound = values - ( SE * qnorm(0.975) )  # CI bounds
    hiBound = values + ( SE * qnorm(0.975) ) 
    # error bars
    arrows(x0 = seq(1:length(values)), y0=loBound, y1=hiBound, angle=90, code=3, length=0.05)
  }
}



######################### FUNCTION: INITIALIZE EMPTY RESULTS LIST #########################

make_result_list = function() {
  empty.list = list( "prop.male" = list("abs.bias"=data.frame(), "std.bias"=data.frame(), "coverage"=data.frame()),
              "other.normals" = list("abs.bias"=data.frame(), "std.bias"=data.frame(), "coverage"=data.frame()), 
              "drug.evers" = list("abs.bias"=data.frame(), "std.bias"=data.frame(), "coverage"=data.frame()),
              "prop.drug.time" = list("abs.bias"=data.frame(), "std.bias"=data.frame(), "coverage"=data.frame()),
              "race" = list("prop.table"=data.frame())
  ) 
  return(empty.list)
}



######################### FUNCTION: COMPLETE PARAMETER DATAFRAME GIVEN N ######################### 

complete_parameters = function(parameters, n) {
  
  # calculate SDs for proportions based on sample size
  parameters$across.SD[parameters$type=="normal.drug"] = sqrt( parameters$across.mean[parameters$type=="normal.drug"]
                                                               * (1-parameters$across.mean[parameters$type=="normal.drug"]) / n )
  # calculate vars based on SDs
  parameters$across.var = parameters$across.SD ^ 2
  
  # ARBITRARILY SET WITHIN-S VARIANCE TO 1/3 OF ACROSS-S VARIANCE
  var.index = parameters$type %in% c("normal.other", "time.function")  # index of variables to consider
  within.var.vector = parameters$within.var[var.index]
  across.var.vector = parameters$across.var[var.index]
  within.var.vector[ is.na(within.var.vector) ] = across.var.vector[ is.na(within.var.vector) ] / 3 
  parameters$within.var[var.index] = within.var.vector

  return(parameters)
}


######################### FUNCTION: MAKE A COVARIATE AS A FUNCTION OF TIME #########################

# given slope and intercept for this subject, increase the variable as a function of time
# and add random error ~ N(0, error.SD) 

# change parameters matrix so that it has a normal.time option for type and an error.SD column for these

add_time_function_vars = function(d4, obs, parameters) { 
  
  # if no time-function variable, return dataframe unchanged
  if ( length( parameters$type[parameters$type == "time.function"] ) == 0 ) return(d4)

  # extract time-function variables
  time.list = parameters$name[parameters$type == "time.function"]
  first.dat = d4[!duplicated(d4$id),]
  time.vars = as.data.frame( first.dat[ , names(first.dat) %in% time.list ] )
  names(time.vars) = parameters$name[parameters$type=="time.function"]  # needed in case there's only 1 time function variable, in which case the subset above is an unnamed vector
  
  # extract parameters for time-function variables
  error.SD = parameters$error.SD[parameters$type == "time.function"]

  d5 = d4
  
  for ( s in unique(d4$id) ) {  # for each subject
    for ( i in 1:ncol(time.vars) ) {  # for each time-function variable
      
      # intercept is the first observation generated for this subject
      intercept = time.vars[s, i]
      
      # pull slope from first.dat
      # what variable name are we looking for?
      # assumes the slope variables are named XXX_slope where XXX is corresponding normal variable
      var.name = paste( names(time.vars)[i], "_slope", sep="" )
      # pull the generated slopes for subject s and take the first one
      subj.slope = d4[[var.name]][d4$id == s][1] 
      
      # generate error vector for all obs
      # error ~ N(0, error.SD)
      errors = rnorm(n=obs, mean = 0, sd = error.SD[i])
      
      # compute vector of observations for this subject as a linear function of time
      v = intercept + subj.slope*seq( 0, (obs-1) ) + errors
      
      # put in d5
      var.index = which(names(d5)==names(time.vars)[i])  # column index for the variable we're doing
      d5[d5$id==s, var.index] = v
    }  
  }
  return(d5)
}






