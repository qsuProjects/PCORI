 
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
#
# USAGE NOTES
#  1.) If the type field in parameters matrix has "static" in its name, it will be overridden using S' first observation.
#  2.) Must put in "ref" as the beta for one entry in categorical parameters matrix.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


############################## PACKAGES ##############################

library(survival)
library(mice)
library(coxme)
library(data.table)


############################## LOCAL TEST ##############################
# # 
# write.path="~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/local-test"
# file.path = "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/local-test/SURV_2015-08-01_simple_covs_dataset_1.csv"
# miss.matrix.hi.path = "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/for-sherlock/missing_var_parameters_matrix_high_Y.csv"
# miss.matrix.lo.path = "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/for-sherlock/missing_var_parameters_matrix_low_Y.csv"
# file.name="fake_file_name_2_@.csv"
# 
# d = read.csv(file.path)
# 
# # create variables for log-time interactions
# d$time.X = log(d$t) * d$X
# d$time.Z = log(d$t) * d$Z
# 
# ( miss.matrix.hi = read.csv(miss.matrix.hi.path) )
# ( miss.matrix.lo = read.csv(miss.matrix.lo.path) )
# 
# name.prefix = "dataset_2"
# 
# time.name="t"
# event.name="d"
# cluster.name="id"
# cox.predictors = c("X")
# na.methods = c("complete.case", "complete.case.by.subj",
#                "naive", "frailty", "log-t", "full", "no.time")  # misnomer because some of these don't involve an NA estimator (e.g., CC)
# 
# impute.with = c("id", "X", "Z", "time.X", "time.Z", "Z", cox.predictors)
# 
# #~~~~~~~~~ CHANGED 12/12
# #make.miss.if.contains="X"
# 
# 
# do_one_dataset(.d=d, .source.file.name=file.name,
#                .miss.matrix.hi=miss.matrix.hi, 
#                .miss.matrix.lo=miss.matrix.lo,
#                .aux.matrix=aux.matrix,
#                .time.name=time.name,
#                .event.name=event.name, .cluster.name=cluster.name,
#                .cox.predictors=cox.predictors, .name.prefix=name.prefix,
#                .na.methods=na.methods, .write.path=write.path,
#                .impute.with = impute.with,
#                .make.miss.if.contains=make.miss.if.contains
#               )
# 
# 
# stitch_files( "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/local-test",
#               "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/local-test",
#               .name.prefix="results",
#               .stitch.file.name="stitched.csv"
#               )
# 

############################## FUNCTION: DO ONE DATASET ##############################

do_one_dataset = function(.d, .source.file.name, .miss.matrix.hi, .miss.matrix.lo, .aux.matrix,
                          .time.name, .event.name, .cluster.name,
                          .cox.predictors, .name.prefix, .na.methods, .write.path,
                          .impute.with, .make.miss.if.contains=NULL, .p.censor) {

  ##### Overwrite Z to make it more correlated with X (r=0.86) #####
  .d$Z = .d$X + rnorm(dim(.d)[1], mean=0, sd=.6)

  ##### Impose Missingness (MAR) ####
      # survival time by subject
      library(data.table)
      dt = data.table(.d)
      dt[, surv := max(t), by=id ]  # get survival time for each subject (last observed time point)
      d2 = data.frame(dt)

  ##### V3: SPLIT ON X AND MAKE MISSINGNESS CONDITIONAL ON X #####
  # split roughly on median X
  hi = d2[d2$X > 0,] 
  lo = d2[d2$X <= 0,] 
  
  # impose missingness by X group
  hi = make_missing_indics(data=hi, miss.matrix=.miss.matrix.hi)
  lo = make_missing_indics(data=lo, miss.matrix=.miss.matrix.lo)

#   JUST FOR TEST PURPOSES
#   print( table(hi$miss.ind_X) )  # how much missingness among high X?
#   print( table(lo$miss.ind_X) )  # and low X?
# 
#   print( median(hi$surv[hi$miss.ind_X==0]) )  # mean of observed high-X
#   print( median(hi$surv) )  # mean of all high-X
#   print( median(lo$surv[lo$miss.ind_X==0]) )  # mean of observed low-X
#   print( median(lo$surv) )  # mean of all low-X

  # glue together
  d2 = rbind(hi, lo)
  
  # make missing based on Frankenstein dataset
  d3 = missingify(d2, make.relateds.missing=FALSE, make.miss.if.contains=.make.miss.if.contains)
  
  for (i in .na.methods) {  # do once for each NA estimator method
    
    cat( "\n\nBeginning loop for method: ", i)
    
    d4 = d3
    
    ##### Get NA Estimator ####
        if (i == "naive") { d4$estim = calc_NA(.time.name=.time.name,
                                              .event.name=.event.name, .type="naive",
                                              .data=d4) }
  
        if (i == "frailty") { d4$estim = calc_NA(.time.name=.time.name,
                                               .event.name=.event.name, 
                                               .cluster.name="id",
                                               .type="frailty",
                                              .data=d4) }
    
        if (i == "log-t") d4$estim = log( d4[[.time.name]] )
        
        # log
        cat( "\n\nFinished making estimator for method: ", i)

    ##### Impute Using NA Estimator #####
        #if (!i %in% c("complete.case", "full") ) {  # unless doing a non-imputation approach
    
        if (i %in% c("naive", "frailty", "log-t", "no.time") ) {  # if using an imputation-based approach
          imp = impute(.data=d4, .method="pmm", .cluster.name="id", .na.name="estim", .impute.with)
          cat( "\nFinished imputing")
        }
 
    ##### Fit Cox Frailty Model #####
        # formula is same regardless of method
        .coxme_formula <- paste0("Surv(t0, t, d) ~ ",
                                 paste0(.cox.predictors, collapse = " + "), " + (1|id) ")   
        
        # different model call depending on whether we need to pool over imputations or not
        if (i == "complete.case") {
          rs1 = coxme( as.formula(.coxme_formula), data=d3 )
          coefs = coxme_coefs(.coxme_object=rs1)
        }

        if (i == "complete.case.by.subj") {  # new as of 10/9
          dt = data.table(d3)
          # indicator for whether a subject has any missingness
          dt[, has.miss := mean( is.na(d3[ , 1:dim(d3)[2] ] ) ) > 0, by=id ]
          rs1 = coxme( as.formula(.coxme_formula), data=d3 )
          coefs = coxme_coefs(.coxme_object=rs1)
        }
    
        if (i == "full") {
          rs1 = coxme( as.formula(.coxme_formula), data=d2 )  # use d2, pre-missingness, instead of d3
          coefs = coxme_coefs(.coxme_object=rs1)
        }
        
        # all imputation-based approaches get treated the same
        if (i %in% c("naive", "frailty", "log-t", "no.time")) {
          rs1 = with( imp, coxme( as.formula(.coxme_formula) ) )
          rs2 = pool(rs1); print(rs2)
        coefs = pooled_coxme_coefs(.coxme_object=rs2)
      }
    
    cat( "\n\nFinished imputing for method: ", i)
    
    ##### Export Useful Info About Dataset ####
    coefs$source.file = .source.file.name
    coefs$method = i
    coefs$prop.rows.missing = sum(is.na(d4$X)) / length(d4$X)  # proportion of observations (NOT subjects) with missing X
    coefs$XZ.corr = cor(.d$X, .d$Z)  # X-Z correlation in full data
    
    cat( "\nFinished fitting Cox model")

    # write single-line file of results
    setwd(.write.path)
    write.csv( coefs, paste(Sys.Date(), .name.prefix, "NA", i, "results.csv", sep="_") ) 
  }
}



############################## FUNCTION: RANDOM CENSORING ##############################

##~~~~~~~~ BOOKMARK: CHANGE THIS. CURRENTLY P.CENSOR IS ACTUALLY P OF EVENT
# .d = dataset
# .p.censor = probability of censoring
# random_censor = function(.d, .id.var.name, .p.censor) {
#   dt = data.table(.d)
#   dt[, event := rbinom(n=1, size=1, prob=.p.censor), by=.id.var.name ]
#   return( as.data.frame(dt) )
# }

# test
#fake = random_censor(d, "id", 0.999)

################################# MAKE NELSON-AALEN ESTIMATORS #################################

# create naive (unclustered) Nelson-Aalen estimator (naive; ignoring study clustering)
# or frailty NA estimator
# this function is modified from SAS & R blog (http://sas-and-r.blogspot.com/2010/05/example-739-nelson-aalen-estimate-of.html?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+SASandR+%28SAS+and+R%29)

calc_NA = function(.time.name, .event.name, .type, .cluster.name=NA, .data) {
  
  # Cox formula
  if (.type=="naive") RHS = "1"
  if (.type=="frailty") RHS = paste( "(1|", .cluster.name, ")" )
  
  # for censored data, simply use the binary event indicator for .event.name
  form = paste("Surv(", .time.name, ",", .event.name, ") ~", RHS, sep="")
  
  # get NA estimate from coxph
  na.fit = survfit( coxph( as.formula(form), data=.data ), type="aalen" )
  
  # times at which risk set changes
  time = .data[[.time.name]]
  jumps = c( 0, na.fit$time, max(time) ) 
  
  # cumulative survival at each jump time
  # need to be careful at the beginning and end
  surv = c(1, na.fit$surv, na.fit$surv[length(na.fit$surv)])
  
  # apply appropriate transformation
  neglogsurv = -log(surv)   
  
  # create placeholder of correct length
  naest = numeric(length(time))  
  
  # for each time, apply the appropriate NA estimate
  for (i in 2:length(jumps)) {
    naest[which(time>=jumps[i-1] & time<=jumps[i])] = 
      neglogsurv[i-1]   # snag the appropriate value
  }
  return(naest)
}



################################# FUNCTION: IMPUTE #################################

impute = function( .data, .method, .cluster.name, .na.name, .impute.with ) {
  # remove the variables that aren't needed for imputation
  # except keep the ones we need for modeling
  
  keep.in.data = c( .impute.with, "t", "t0", "d", "estim" ) 
  #keep.in.data = c( .impute.with, "t", "t0", "d", "estim" )  # FOR UNCENSORED VERSION
  .data = .data[ , names(.data) %in% keep.in.data ]
  
  imp=mice(.data)
  
  # print out pred matrix
  cat("\nPred matrix:\n")
  
  print(imp$pred)
  
  return(imp)
}



# PREVIOUS VERSION OF IMPUTE FUNCTION THAT MICROMANAGES THE PRED MATRIX. 
# data: raw dataset
# method: "2l.norm.me" for Resche-Rigon's MICE-RE, "2l.norm" for MICE's native 2l.norm, or "pmm" for MICE default (primary analysis)
# 
# impute = function( .data, .method, .cluster.name, .na.name, .impute.with ) {
#   
#   browser()
#   # remove the variables that aren't needed for imputation
#   # except keep the ones we need for modeling
#   keep.in.data = c( .impute.with, "t", "t0", "d", "estim" )
#   .data = .data[ , names(.data) %in% keep.in.data ]
#   
#   cat("\nVariables in imputation data (not all used for imputation):")
#   print(names(.data))
#   
#   # if using 2l.norm, put in the required constant term
#   if (.method == "2l.norm") .data$const=1
#   
#   # first fit normal MICE to get predictor matrix
#   ini = mice(.data, maxit = 0)
#   pred = ini$predictorMatrix
#   
#   # don't impute with these ones
#   # (most are already removed except for those needed for Cox model)
#   pred[, "t0"] = 0
#   pred[, "t"] = 0
#   
#   #LOG
#   cat("\nEntered impute function; finished dry run")
#   
#   # if not using default PMM, modify the predictor matrix to specify multilevel model
#   if (.method != "pmm") {
#     # treat trial as the cluster term
#     pred[ , .cluster.name] = -2
#     
#     # fixed effects for N-A estimator
#     pred[ , .na.name ] = 1
#  
#     # if using 2l.norm, code constant as random effect
#     if (.method=="2l.norm") pred[ , "const" ] = 2
#   }
#   
#   # change method based on user specification
#   method = ini$method; method[method == "pmm"] = .method
#   
#   # impute with specified method
#   imp = mice(.data, maxit = 0, pred = pred, method = method)
#   
#   #LOG
#   cat(imp$pred)
#   
#   return(imp)
# }


################################# FUNCTIONS: GET COEFFICIENTS FROM COX MODELS #################################

pooled_coxme_coefs = function(.coxme_object) {
  s = summary(.coxme_object)
  
  coefs = data.frame( s[,c("est", "se", "lo 95", "hi 95")] )
  names(coefs) = row.names(s)  # use actual variable names
  return(coefs)
}


# modified from Kris' function
coxme_coefs <- function (.coxme_object) 
{
  .tmp <- ""
  .beta <- .coxme_object$coefficients
  .nvar <- length(.beta)
  .nfrail <- nrow(.coxme_object$var) - .nvar
  .omit <- .coxme_object$na.action
  if (.nvar > 0) {
    .se <- sqrt(bdsmatrix::diag(.coxme_object$var)[.nfrail + 1:.nvar])
    .tmp <- cbind(.beta, 
                  .se, 
                  .beta + qnorm(.025)*.se,
                  .beta + qnorm(.975)*.se)
    dimnames(.tmp) <- list(names(.beta), c("est", 
                                           "se", 
                                           "lo 95",
                                           "hi 95"))
  }
  coefs = data.frame( t( .tmp ) )
  return(coefs)
}

