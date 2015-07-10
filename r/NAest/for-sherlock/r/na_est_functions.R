 
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


############################## LOCAL TEST ##############################
# 
# write.path="~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/local-test"
# file.path = "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/local-test/SURV_2015-02-01_job_10_dataset_1.csv"
# miss.matrix.path = "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/local-test/missing_var_parameters_matrix.csv"
# file.name="fake_file_name_2_@.csv"
# 
# # my impose missingness code
# #source("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/IMPMISS/Code/impose_missingness_functions.R")
# 
# d = read.csv(file.path)
# d = d[1:2000,]  # small one just for debugging
# miss.matrix = read.csv(miss.matrix.path)
# 
# name.prefix = "dataset_2"
# 
# 
# time.name="t"
# event.name="d"
# cluster.name="id"
# cox.predictors = c("ind_cd4_50_100", "ind_cd4_350_500", "ind_cd4_200_350", "ind_cd4_100_200", "d_dida")
# na.methods = c("naive", "frailty", "log-t")
# 
# # RUN BY MANISHA!
# impute.with = c("id", "d", cox.predictors)
# 
# make.miss.if.contains="cd4"
# 
# 
# 
# do_one_dataset(.d=d, .source.file.name=file.name, .miss.matrix=miss.matrix, .time.name=time.name,
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

do_one_dataset = function(.d, .source.file.name, .miss.matrix, .time.name, .event.name, .cluster.name,
                          .cox.predictors, .name.prefix, .na.methods, .write.path,
                          .impute.with, .make.miss.if.contains=NULL) {

  ##### Impose Missingness (MAR) ####
  d2 = make_missing_indics(data=.d, miss.matrix=.miss.matrix)
  d3 = missingify(d2, make.relateds.missing=FALSE, make.miss.if.contains=.make.miss.if.contains)
  
  for (i in .na.methods) {  # do once for each NA estimator method
    d4 = d3
    
    ##### Get NA Estimator ####
    if (i == "naive") { d4$estim = calc_NA(.time.name=.time.name,
                                          .event.name=.event.name, .type="naive",
                                          .data=d4) }
    if (i == "strat") { d4$estim = calc_NA_strat(time=d4[[.time.name]],
                                                 event=d4[[.event.name]], cluster=d4[[.cluster.name]], data=d4) }
    
    cat( "\nAbout to do frailty NA")
    
    if (i == "frailty") { d4$estim = calc_NA(.time.name=.time.name,
                                           .event.name=.event.name, 
                                           .cluster.name="id",
                                           .type="frailty",
                                          .data=d4) }
    cat( "\nDid frailty NA")
    
    if (i == "log-t") d4$estim = log( d4[[.time.name]] )
    
    # LOG
    cat( "\nFinished making estimator")
    
    ##### Impute Using NA Estimator ####
    imp = impute(.data=d4, .method="2l.norm", .cluster.name="id", .na.name="estim", .impute.with)
    
    # LOG
    cat( "\nFinished imputing")
    
    #print( head(imp$data) )
    
    ##### Fit Cox Frailty Model ####
    .coxme_formula <- paste0("Surv(t0, t, d) ~ ",
                             paste0(.cox.predictors, collapse = " + "), " + (1|id) ")   
    rs1 = with( imp, coxme( as.formula(.coxme_formula) ) )
    
    print(rs1)
    
    rs2 = pool(rs1)
    
    print(rs2)

    coefs = pooled_coxme_coefs(.coxme_object=rs2)
    coefs$source.file = .source.file.name
    coefs$method = i
    # SPECIFIC TO THIS SIM: PROPORTION OF MISSING VALUES
    coefs$prop.missing = sum(is.na(d4$ind_cd4_50_100)) / length(d4$ind_cd4_50_100)
    
    # LOG
    cat( "\nFinished fitting Cox model")

    # write single-line file of results
    setwd(.write.path)
    write.csv( coefs, paste(Sys.Date(), .name.prefix, "NA", i, "results.csv", sep="_") ) 
  }
}






################################# MAKE NELSON-AALEN ESTIMATORS #################################

# create naive (unclustered) Nelson-Aalen estimator (naive; ignoring study clustering)
# or frailty NA estimator
# this function is modified from SAS & R blog (http://sas-and-r.blogspot.com/2010/05/example-739-nelson-aalen-estimate-of.html?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+SASandR+%28SAS+and+R%29)

calc_NA = function(.time.name, .event.name, .type, .cluster.name=NA, .data) {
  
  # Cox formula
  if (.type=="naive") RHS = "1"
  if (.type=="frailty") RHS = paste( "(1|", .cluster.name, ")" )
  form = paste("Surv(", .time.name, ",", .event.name, ") ~", RHS, sep="")
  #form = paste("Surv(time, event) ~", RHS)
  
  # get NA estimate from coxph
  #na.fit = survfit( coxph( Surv(time, event) ~ 1 ), type="aalen" )
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


# create stratified Nelson-Aalen estimator (computed separately within each study)
# this function is modified from SAS & R blog (http://sas-and-r.blogspot.com/2010/05/example-739-nelson-aalen-estimate-of.html?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+SASandR+%28SAS+and+R%29)
# calc_NA_strat = function(time, event, cluster, data) {
# 
#   final.na.est = c()
#   
#   for ( j in unique(cluster) ) {
#     temp = data[cluster==j, ]  # subsets with only the desired cluster
#     time.temp = time[cluster==j]
#     event.temp = event[cluster==j]
#     
#     na.fit = survfit( coxph( Surv(time.temp, event.temp) ~ 1 ), type="aalen" )
#     
#     jumps = c( 0, na.fit$time, max(time.temp) ) 
#     # need to be careful at the beginning and end
#     surv = c(1, na.fit$surv, na.fit$surv[length(na.fit$surv)])
#     
#     # apply appropriate transformation
#     neglogsurv = -log(surv)   
#     
#     # create placeholder of correct length
#     naest = numeric(length(time.temp))  
#     for (i in 2:length(jumps)) {
#       naest[ which(time.temp>=jumps[i-1] & time.temp<=jumps[i]) ] = 
#         neglogsurv[i-1]   # snag the appropriate value
#     }
#     final.na.est[cluster==j] = naest
#   }
#   return(final.na.est)
# }


################################# FUNCTION: IMPUTE #################################

# data: raw dataset
# method: "2l.norm.me" for Resche-Rigon's MICE-RE, "2l.norm" for MICE's native 2l.norm, or "pmm" for MICE default (primary analysis)

impute = function( .data, .method, .cluster.name, .na.name, .impute.with ) {
  
  #cat("\nEntered impute function")

  # remove the variables that aren't needed for imputation
  # except keep the ones we need for modeling
  keep.in.data = c( .impute.with, "t", "t0", "d", "estim" )
  .data = .data[ , names(.data) %in% keep.in.data ]
  
  cat("\nVariables in imputation data (not all used for imputation):")
  print(names(.data))
  
  # if using 2l.norm, put in the required constant term
  if (.method == "2l.norm") .data$const=1
  
  # first fit normal MICE to get predictor matrix
  ini = mice(.data, maxit = 0)
  pred = ini$predictorMatrix
  
  # don't impute with these ones
  # (most are already removed except for those needed for Cox model)
  # NEEDS TO ALSO BE t
  pred[, "t0"] = 0
  pred[, "t"] = 0
  
  #LOG
  cat("\nEntered impute function; finished dry run")
  
  # if not using default PMM, modify the predictor matrix to specify multilevel model
  if (.method != "pmm") {
    #pred[ pred==1 ] = 2  # give random effects to all variables used for prediction
    
    # treat trial as the cluster term
    col = pred[ , .cluster.name]; col[col==2] = -2; pred[, .cluster.name] = col
    
    # fixed effects for N-A estimator
    col = pred[ , .na.name ]; col[col==2] = 1; pred[ , .na.name ] = col
 
    # if using 2l.norm, code constant as random effect
    if (.method=="2l.norm") pred[ , "const" ] = 2
  }
  
  #LOG
  #cat(pred)
  
  # change method based on user specification
  method = ini$method; method[method == "pmm"] = .method
  
  # convert class column to integer (required)
  #.data[[.cluster]] = as.integer(d$trial2)
  
  #LOG
  #cat("\nAbout to impute")
  
  # impute with specified method
  imp = mice(.data, maxit = 0, pred = pred, method = method)
  
  return(imp)
}




pooled_coxme_coefs = function(.coxme_object) {
  s = summary(.coxme_object)
  coefs = data.frame( t( s[,c("est", "se", "lo 95", "hi 95")] ) )
  #coefs$var = row.names(s)
  return(coefs)
}



