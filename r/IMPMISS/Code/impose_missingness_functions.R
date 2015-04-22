#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Project: PCORI Missing Data
# Author: Maya Mathur 
#
# This algorithm:
#   1.) Creates auxiliary variables related to cardio, HIV, and etiology variables.
#       These variables affect missingness. 
#   2.) Creates missingness indicators as a function of the covariates and auxiliary variables.
#   3.) Imposes missingness according to the missingess indicators.
#
# Usage notes:
#   1.) Outcome must be coded as 0/1
#   2.) Everyone must have event (pre-censoring)
#   3.) Aux parameter matrix must use var name "time2event" if it uses that as parameter
#   4.) Variables in missingness matrix must use name format "miss.myVar" where myVar is the name of
#       the variable to be made missing.
#   5.) For both parameter matrices, each variable must have a predictor variable called "intercept".
#   6.) Miss parameter matrix must use var name "subj.rand.int" if it uses subject random intercepts.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

########################## WRAPPER FUNCTION: IMPOSE MISSINGNESS ##########################

##### Arguments #####

# data: original dataset
# outcome.name: quoted name of the binary outcome variable
# start.time.name: quoted name of the time defining beginning of each time interval
# stop.time.name: quoted name of the time defining end of each time interval
# id.var.name: quoted name of subject id variable

# aux matrix: a matrix defining the parameters for the auxiliary variables
#   it should have the following columns:
#   aux.var: lists names of auxiliary variables (each will have multiple rows)
#   parameter: lists names of existing variables in dataset that are used to compute the aux var (one must be called "intercept")
#   beta: the coefficient for the variable
#   error.SD: the error SD for the aux var; only first entry will be used

# miss.matrix: a matrix defining the parameters for the missingness indicators
#   it should have the following columns:
#   miss.var: name of missingness indicator; must follow the form "miss.myVar" where myVar is the exact name of the variable to be made missing
#   parameter: lists names of existing variables in dataset that are used to compute the aux var (one must be called "intercept")
#   beta: the coefficient for the variable

# make.relateds.missing: TRUE/FALSE to indicate whether any variable containing the name of the var to be made missing
#   should also be made missing (this takes care of categorical versions of the var)


impose_missingness = function( data, outcome.name, start.time.name, stop.time.name, 
                               id.var.name, aux.matrix, miss.matrix, make.relateds.missing=FALSE,
                               rand.int.SD=NA ) {
  
  outcome = data[[outcome.name]]
  start.time = data[[start.time.name]]
  stop.time = data[[stop.time.name]]
  id = data[[id.var.name]]
 
  ##### Step 1: Create Time-to-Event Variable #####
  d2 = create_time_vars(data=data, id=id, outcome=outcome, start.time=start.time, stop.time=stop.time)
  
  ##### Step 2: Create Subject Random Effect Variable #####
  # only do this if an SD that isn't NA is specified
  if (!is.na(rand.int.SD)) {d.rand = create_random_effects(data=d2, id.var.name=id.var.name,
                                                      rand.int.SD=rand.int.SD)}

  ##### Step 2: Create All Auxiliary Variables ####
  d3 = make_aux_vars(data=d.rand, aux.matrix=aux.matrix)
  
  ##### Step 3: Create All Missingness Indicators ####
  d4 = make_missing_indics(data=d3, miss.matrix=miss.matrix)
  
  ##### Step 4: Override Observations using Missingness Indicator ####
  d5 = missingify(d4, make.relateds.missing)
  
  ##### Step 5: Print Baseline Proportion Missing for a Few Vars ####
  first.dat = d5[ d5$baseline, ]
  cat("\n")

  for (i in c("ldl", "log_vln", "cd4")) {
    prop.missing.baseline = sum( is.na(first.dat[[i]]) ) / length( is.na(first.dat[[i]]) )
    cat( "\nProportion of baseline observations made missing for", i, ":", prop.missing.baseline )
  }
  
  ##### Step 5: Print Info on Global Missingness Targets ####
  proportion_dropped(data=d5, miss.matrix=miss.matrix)
  prop_subj_dropped(data=d5, id.var.name=id.var.name, miss.matrix=miss.matrix)
  
  return(d5)
}


########################## STEP 1 FUNCTION: CREATE TIME-TO-EVENT VARIABLE ##########################

create_time_vars = function( data, id, outcome, start.time, stop.time ) {
  # since we assume everyone has the event, time to event is just the largest observed time for each subject
  maxes = aggregate( stop.time ~ id, FUN=max)[,2]
  
  # number observations for each subject
  reps = as.vector(table(id))
  
  # put in dataframe
  data$time2event = rep(maxes, times=reps)
  
  # make indicator for whether it's the baseline observation
  data$baseline = (start.time == 0)
  
  return(data)
}

# test - WORKS :)
#create_time2event_var( data=data, id=data$id, outcome=data$d, stop.time=data$t)



########################## STEP 2 FUNCTION: CREATE RANDOM INTERCEPT VARIABLE ##########################

#fake = rnorm(mean=0, sd=rand.int.SD, n=length(unique(data[[id]])))
#length(unique(data[[id]]))

create_random_effects = function( data, id.var.name, rand.int.SD ) {
  
  # generate random intercepts ~ N(0, rand.int.SD)
  rand.intercepts = rnorm(mean=0, sd=rand.int.SD, n=length(unique(data[[id.var.name]])))
  
  # number of observations for each subject (vector)
  obs.per.subj = as.numeric(table(data[[id.var.name]]))
  
  # put in data
  fake = data
  fake$subj.rand.int = rep(rand.intercepts, obs.per.subj)
  return(fake)
}



########################## STEP 2 FUNCTION: CREATE ALL AUX VARIABLES ##########################

# given the dataset, apply function below to create all auxiliary variables
#  and append to dataset

make_aux_vars = function( data, aux.matrix ) {
  # extract names of auxiliary variables from matrix
  aux.names = unique( aux.matrix$aux.var )
  
  # start with passed dataset
  # iterate through auxiliary variables, making and appending them one at a time
  current.data = data
  for (i in aux.names) {
    current.data = make_one_aux_var( data=current.data, aux.var.name=i, aux.matrix=aux.matrix )
  }
  
  return(current.data) 
}


# test it - WORKS! :)
#data=d2
#fake = make_aux_vars(d2, aux.matrix)


########################## STEP 2 FUNCTION: CREATE A SINGLE AUX VARIABLE VECTOR ##########################

# given the input matrix, name of aux variable, and number observations, apply betas
#  and add noise to get vector for the aux variable with length n.obs

make_one_aux_var = function( data, aux.var.name, aux.matrix ) {

  n.obs = nrow(data)
  
  # pull out relevant parameters except intercept
  if ( !any(miss.matrix$parameter == "intercept") ) stop("Intercept is missing in missing variable matrix")
  matrix2 = aux.matrix[ aux.matrix$aux.var==aux.var.name & aux.matrix$parameter!="intercept", ]
  
  # pull out intercept
  intercept = aux.matrix$beta[ aux.matrix$aux.var==aux.var.name & aux.matrix$parameter=="intercept" ]
  
  # pull out model type (linear or logistic)
  model = as.character( aux.matrix$model[ aux.matrix$aux.var==aux.var.name ][1] )
  
  # calculate linear predictor, not including any error term
  linear.pred = intercept
  for (p in matrix2$parameter) {
    beta = matrix2$beta[ matrix2$parameter==p ]  # beta is length 1; recycled
    linear.pred = linear.pred + ( data[[p]] * beta )  # linear predictor (no error term yet)
  }
  
  # if auxiliary variable is built according to a linear model
  if (model == "linear") {
    # random error vector
    error.SD = aux.matrix$error.SD[aux.matrix$aux.var==aux.var.name][1]  # use just first entry since constant within an aux variable
    #error.SD = matrix2$error.SD[1]  # use just first entry since constant within an aux variable
    errors = rnorm(n=n.obs, mean=0, sd=error.SD )
    
    # observed values
    values = linear.pred + errors
  } 
  
  # if auxiliary variable is built according to a logistic model 
  if (model == "logistic") {
    # calculate probability of missingness, pi
    pi = exp(linear.pred) / (1 + exp(linear.pred))
    
    if ( any(is.na(linear.pred)) ) warning( cat("\n\nLinear predictor to construct", aux.var.name,
                                                "has missing values.\nLinear predictor may have been too large, leading to Inf values when calculating predicted probability.\n") )
    
    # observed values (Bernouilli trials)
    values = rbinom(n=n.obs, size=1, p=pi)
  }
  
  # put new variable in dataframe
  data[[aux.var.name]] = values
  return(data)
}


# test it - WORKS! :)
#fake = make_one_aux_var(data=data, aux.var.name="aux.cardio", aux.matrix=aux.matrix)
#fake = make_one_aux_var(data=data, aux.var.name="aux.hiv", aux.matrix=aux.matrix)

#d2 = create_time2event_var(data=data, id=id, outcome=outcome, stop.time=stop.time)
#fake = make_one_aux_var(data=d2, aux.var.name="aux.etiol", aux.matrix=aux.matrix)



########################## STEP 3 FUNCTION: CREATE A SINGLE MISSINGNESS INDICATOR ##########################

make_one_missing_indic = function( data, miss.var.name, miss.matrix ) {
  
  n.obs = nrow(data)
  
  # pull out relevant parameters except intercept
  if ( !any(miss.matrix$parameter == "intercept") ) stop("Intercept is missing in missing variable matrix")
  matrix2 = miss.matrix[ miss.matrix$miss.var==miss.var.name &
                           miss.matrix$parameter!="intercept", ]
  
  # pull out intercept
  intercept = miss.matrix$beta[ miss.matrix$miss.var==miss.var.name &
                                  miss.matrix$parameter=="intercept" ]
  
  # start with just intercept
  linear.pred = intercept
  
  # for each parameter, multiply it by corresponding beta and add to values vector
  # linear.pred becomes vector with same length as dataset
  for (p in matrix2$parameter) {
    beta = matrix2$beta[ matrix2$parameter==p ]  # beta is length 1; recycled
    linear.pred = linear.pred + ( data[[p]] * beta )  
  }

  # calculate probability of missingness, pi
  pi = exp(linear.pred) / (1 + exp(linear.pred))
  
  # if linear predictor so large that it couldn't be exponentiated, set pi to 1 to avoid missing values
  if ( any( exp(linear.pred) == Inf ) ) warning( cat("\n\nLinear predictor to construct", miss.var.name,
                                              "had enormous values; setting pi to 1 for each of these.\n\n") )
  pi[ exp(linear.pred) == Inf ] = 1
  
  # simulate Bernouilli trials based on this
  values = rbinom(n=n.obs, size=1, p=pi)
  
  # print proportion missing and non-missing
  prop.missing = length(values[values==1]) / length(values)
  cat( "\nProportion of overall observations made missing for", miss.var.name, ":", prop.missing )
  
  # put in dataframe
  data[[miss.var.name]] = values
  return(data)
}

# test - WORKS :)
#data=d3
#miss.var.name="miss.bmi"
#fake = make_one_missing_indic( data=d3, miss.var.name="miss.bmi", miss.matrix=miss.matrix )


########################## STEP 3 FUNCTION: CREATE ALL MISSINGNESS INDICATORS ##########################

make_missing_indics = function( data, miss.matrix ) {
  # extract names of auxiliary variables from matrix
  
  miss.names = unique( miss.matrix$miss.var )
  
  # start with passed dataset
  # iterate through auxiliary variables, making and appending them one at a time
  current.data = data
  for (i in miss.names) {
    current.data = make_one_missing_indic( data=current.data,
                                           miss.var.name=i, miss.matrix=miss.matrix )
  }
  
  return(current.data) 
}

# test
#data = d3




########################## STEP 4 FUNCTION: OVERRIDE MISSING ONES ##########################

# assumes that missing indicators are titled "miss.myVar"

missingify = function( data, make.relateds.missing ) { 
  # extract missing indicators
  miss.indic.names = names(data)[ substr( names(data), 0, 4 )=="miss" ]

  # for each missing indicator, find the corresponding variable and make it missing
  for (i in miss.indic.names) {
    orig.var.name = substr( i, 6, nchar(i) )  # extract name of original variable that will be missingified
    is.missing = (data[[i]]==1)  # vector of missing indicators
    data[[orig.var.name]][is.missing] = NA
    
    if (make.relateds.missing) {
      # find names of "related" variables
      # where a "related" variable is one whose name contains the original variable name
      related.names = names(data)[ grep( orig.var.name, names(data) ) ]
      
      cat("\n\nVariables made missing:\n")
      cat(related.names)

      # impose missingness
      data[ is.missing, names(data) %in% related.names ] = NA
    }
  }
  return(data)
}

# test - WORKS :)
#data = d4


########################## STEP ## FUNCTION: COMPUTE PROPORTION OF OBS DROPPED IN CC ANALYSIS ##########################

# want 84%
proportion_dropped = function(data, miss.matrix) {

  # keep only vars of interest
  vars.of.interest = unique(miss.matrix$miss.var)
  temp = data[ , names(data) %in% vars.of.interest ]

  # compute proportion dropped in complete-case analysis
  prop.dropped = sum( complete.cases(temp) ) / length( complete.cases(temp) )
  cat( "\n\nProportion of all observations dropped in CC analysis:", prop.dropped )
  return(prop.dropped)
}


########################## STEP ## FUNCTION: COMPUTE PROPORTION OF SUBJECTS DROPPED IN CC ANALYSIS ##########################

# subject gets dropped if they have a missing value somewhere in each row
# or: subject is kept if they have at least one row that's 100% complete

prop_subj_dropped = function(data, id.var.name, miss.matrix) {

  # keep only vars of interest
  vars.of.interest = unique(miss.matrix$miss.var)
  temp = data[ , names(data) %in% vars.of.interest ]
  
  # dataset of only complete rows
  d.complete = data[complete.cases(temp),]  
  
  # keep only S who are represented in this complete dataset
  n.kept = length( unique( d.complete[[id.var.name]] ) )
  prop.dropped = n.kept / length( unique( data[[id.var.name]] ) )
  
  cat( "\n\nProportion of all subjects dropped in CC analysis:", prop.dropped )
  return(prop.dropped)
}

