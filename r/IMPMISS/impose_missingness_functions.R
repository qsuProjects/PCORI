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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TEST: 4:10 


########################## READ IN DATA ##########################

# TEMPORARY - EVENTUALLY MOVE TO RUN FILE

# t: time at end of interval
# t0: time at beginning of interval
# d: whether event occurred at that interval

setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/IMPMISS/Data from KK")
data = read.csv("SURV_2015-02-01_job_10_dataset_1.csv")

setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/IMPMISS")
aux.matrix = read.csv("aux_var_parameters_matrix_v2.csv")
miss.matrix = read.csv("missing_var_parameters_matrix_v4.csv")


########################## TEST PERFORMANCE ##########################

#data=data; outcome.name="d"; start.time.name="t0";
#stop.time.name="t"; id.var.name="id"; aux.matrix=aux.matrix;
#miss.matrix=miss.matrix; make.relateds.missing=TRUE

d.miss = impose_missingness( data=data, outcome.name="d", start.time.name="t0",
                             stop.time.name="t", id.var.name="id", aux.matrix=aux.matrix,
                             miss.matrix=miss.matrix, make.relateds.missing=TRUE )

# refit models that created missingness to see how we're doing
# looks good!
#summary( glm( is.na(bmi) ~ aux.etiol + time2event + d, data=d.miss, family=binomial(link="logit") ) ) 


##### Did We Hit Individual Variable Missingness Targets? #####

# compare probability of missingness for VLN at baseline vs. overall
first.dat = d.miss[!duplicated(d.miss$id),]
prop.table( table( is.na(first.dat$log_vln) ) )  # at baseline

# compare probability of missingness for CD4 at baseline vs. overall
prop.table( table( is.na(first.dat$cd4) ) )  # at baseline

# compare probability of missingness for LDL at baseline vs. overall
prop.table( table( is.na(first.dat$ldl) ) )  # at baseline


##### Did We Hit Global Missingness Targets? #####

# proportion of observations missing at least 1 variable
# want 84%
d.miss$row.has.missing = ( is.na(d.miss$bmi) | is.na(d.miss$log_vln) | is.na(d.miss$ldl) | 
                          is.na(d.miss$hdl) | is.na(d.miss$cd4) )
prop.table( table(d.miss$row.has.missing) )

# proportion of subjects where each observation is missing at least 1 variable 
# want 18%
gets.dropped = c()
for (i in unique(d.miss$id) ) {
  gets.dropped[i] = all( d.miss$row.has.missing[d.miss$id==i] )
}
prop.table(table(gets.dropped))

# proportion of subjects who are systematically missing each var
dt = data.table(d.miss)
dt[, prop.miss.bmi := mean( is.na(miss.bmi) ), by=id ]
dt[, prop.miss.ldl := mean( is.na(miss.ldl) ), by=id ]
dt[, prop.miss.hdl := mean( is.na(miss.hdl) ), by=id ]
dt[, prop.miss.log_vln := mean( is.na(miss.log_vln) ), by=id ]
dt[, prop.miss.cd4 := mean( is.na(miss.cd4) ), by=id ]

first.dat = dt[ !duplicated(dt$id), ]
prop.table( table(first.dat$prop.miss.bmi == 1) )


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
                               id.var.name, aux.matrix, miss.matrix, make.relateds.missing=FALSE ) {
  
  outcome = data[[outcome.name]]
  start.time = data[[start.time.name]]
  stop.time = data[[stop.time.name]]
  id = data[[id.var.name]]
  
  ##### Step 1: Create Time-to-Event Variable #####
  d2 = create_time2event_var(data=data, id=id, outcome=outcome, stop.time=stop.time)
  
  ##### Step 2: Create All Auxiliary Variables ####
  d3 = make_aux_vars(data=d2, aux.matrix=aux.matrix)
  
  ##### Step 3: Create All Missingness Indicators ####
  d4 = make_missing_indics(data=d3, miss.matrix=miss.matrix)
  
  ##### Step 4: Override Observations using Missingness Indicator ####
  d5 = missingify(d4, make.relateds.missing)
}


########################## STEP 1 FUNCTION: CREATE TIME-TO-EVENT VARIABLE ##########################

create_time2event_var = function( data, id, outcome, stop.time ) {
  # since we assume everyone has the event, time to event is just the largest observed time for each subject
  maxes = aggregate( stop.time ~ id, FUN=max)[,2]
  
  # number observations for each subject
  reps = as.vector(table(id))
  
  # put in dataframe
  data$time2event = rep(maxes, times=reps)
  return(data)
}

# test - WORKS :)
#create_time2event_var( data=data, id=data$id, outcome=data$d, stop.time=data$t)



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
    error.SD = matrix2$error.SD[1]  # use just first entry since constant within an aux variable
    errors = rnorm(n=n.obs, mean=0, sd=error.SD )
    
    # observed values
    values = linear.pred + errors
  } 
  
  # if auxiliary variable is built according to a logistic model 
  if (model == "logistic") {
    # calculate probability of missingness, pi
    pi = exp(linear.pred) / (1 + exp(linear.pred))
    
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
  
  # simulate Bernouilli trials based on this
  values = rbinom(n=n.obs, size=1, p=pi)
  
  # print proportion missing and non-missing
  prop.missing = length(values[values==1]) / length(values)
  cat( "\nProportion of observations made missing for", miss.var.name, ":", prop.missing )
  
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
  
  # browser()
  
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
