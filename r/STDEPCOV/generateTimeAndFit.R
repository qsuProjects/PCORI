if(FALSE) {
  rm(list = ls())
  source("~/.Rprofile")
}

###########################################
#this script loads simulated covariate data,
#uses those data to generate survival times,
#merges those survival times with covariate data,
#then fits various models on results

#It is written to be sourced from a script 
##executed on a multi-core system.

#Sourcing script should initialize:

#number_cores: Integer, number of cores to use
#the_seed: Integer, the seed
#covariate_data_path: String, path and filename for covariate data
#results_write_directory: String, where to write results
#beta: Data frame, one row, each column has a corresponding column
##in the covariate file. Names column names should be identical.
##values are the values to use to generate linear predictor
##for use in survival time generation algorithm
#nreps: Integer, number of replications
#sim_results_name: String, file name for simulation results


require(msm)
require(survival)
require(plyr)

setwd("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/")

#load data
covariates <- read.table(file = covariate_data_path, sep = ",", header = T, stringsAsFactors = F)

#generate dummy variables for race
covariates$race_black <- 0
covariates$race_other <- 0

covariates$race_black[covariates$racecatexp == 1] <- 1
covariates$race_other[covariates$racecatexp == 2] <- 1

#determine how many subjects are in the data
n_subjects <- length(unique(covariates$id))

#set up parallelization
require(doSNOW)
cl<-makeCluster(number_cores)
registerDoSNOW(cl)

n_clusters
getDoParWorkers()

clusterApply(cl, seq(along=cl), function(.id, .results_write_directory, .sim_results_name) {
  WORKER.ID <<- paste0("core_", .id)
  set.seed(the_seed*(.id^8))
  cat("coef,se,p,var,type,rep\n", file = paste0(.results_write_directory, WORKER.ID, "_", .sim_results_name, ".csv"), append = F)
}, results_write_directory, sim_results_name)

# 
# ##define beta
# beta <- data.frame(male = 0.2,
#                    d_abac = 0.7, 
#                    d_ataz = 0.5,
#                    d_dida = 0.3,
#                    d_efav = 0,
#                    d_emtr = 0,
#                    d_indi = 0,
#                    d_lami = 0,
#                    d_lopi = 0,
#                    d_nelf = 0,
#                    d_nevi = 0,
#                    d_rito = 0,
#                    d_saqu = 0,
#                    d_stav = 0,
#                    d_teno = 0,
#                    d_zido = 0,
#                    age = 0.02,
#                    bmi = 0.02,
#                    cd4 = 0.02,
#                    vln = 0.0002,
#                    bps = 0.02,
#                    bpd = 0.02,
#                    ldl = 0.02,
#                    hdl = 0.02,
#                    trig = 0.02,
#                    race_black = -0.2,
#                    race_other = 0)


l_ply(1:nreps, .parallel = T, function(.rep, .covariates, .n_subjects, .beta, .results_write_directory, .sim_results_name) {
  
  #load required packages
  require(msm)
  require(survival)
  require(plyr)
  require(MBESS)
  
  
  #generate time-dependent lambda
  .covariates <- .covariates[ , c(names(.beta), "id")]
  .covariates$linpred <- as.matrix(.covariates[ ,names(.beta)]) %*% t(.beta)
  .covariates$xB <-  exp(.covariates$linpred)
  
  #put data frame in list form by id
  .covars_list <- dlply(.covariates, .(id))
  
  #define simulation parameters
  .ZBbar <- mean(.covariates$linpred)
  .nu <- 3
  .mediangoal <- 200
  
  .lambda <- (log(2)/exp(.ZBbar))*.mediangoal^(-.nu)
  
  # the g function is defined as the inverse of the baseline cummulative hazard from
  ## a Weibull with shape nu and scale lambda defined above
  .g <- function(x){
    ((1/.lambda)*x)^(1/.nu)
  }
  .g_inv <- function(x){
    lambda*(x^.nu)
  }  
  
  .t <- 1:350
  .t_diff <- (t[-1] - t[1:(length(.t) - 1)])[-(length(.t) - 1)]
  .g_inv_t <- .g_inv(.t)
  .g_inv_t_diff <- (.g_inv(.t[-1]) - .g_inv(.t[1:(length(.t) - 1)]))[-(length(.t) - 1)]
  
  #CREATING THE BOUNDS OF TRUNCATION
  .t_max <- 350
  .t_min <- 20
  
  .g_inv_t_max <- .g_inv(t_max)
  .g_inv_t_min <- .g_inv(t_min)
  
  
  #K function applies ACCEPT-REJECT algorithm
  .k <- function(..x, ..m, ..M, ..rates, ..t){
    ifelse(..x <= ..m | ..x >= ..M, 0, dpexp(..x, ..rates, ..t))
  }
  
  #define survival time generation function
  .gen_y <- function(.x, .g_inv_t,  .g_inv_t_min, .g_inv_t_max) {
    .x1 <- .x$xB
    #   print(length(.x$xB))
    #   print(length(.g_inv_t))
    .d <- ppexp(.g_inv_t_max, .x1, .g_inv_t) - ppexp(.g_inv_t_min, .x1, .g_inv_t)
    .M <- 1 / .d
    .r <- 60
    .count<-0
    #counter of times repeat is run
    while (.count<1000) {
      .count <- .count+1
      .y <- rpexp(.r, .x1, .g_inv_t)
      .u <- runif(.r)
      .t <- .M * (k(.y, .g_inv_t_min, .g_inv_t_max, .x1, .g_inv_t) / .d / dpexp(.y, .x1, .g_inv_t))
      .y <- .y[.u <= .t][1]
      if (!is.na(.y)) {break}
    }
    .y
  }
  
  .survival_times <- ldply(.covars_list, .gen_y, .g_inv_t,  .g_inv_t_min, .g_inv_t_max)
  
  .survival_times$g.y <- g(.survival_times[ ,2])
  
  if (sum(is.na(.survival_times$V1)) == 0) {
    
    #create uncensored model-ready dataset
    .data_none <-  ldply(1:.n_subjects, function(..subject, ..survival_times, ..covariates) {
        ..survival_time <- ceiling(..survival_times$g_y[..subject])
        ..to_return <- ..covariates[[..subject]]
        ..to_return$id <- ..subject
        ..to_return$t <- c(1:..survival_time)
        ..to_return$t0 <- 0:(..survival_time - 1)
        ..to_return$d <- c( rep(0, ..survival_time - 1), 1)
        ..to_return$proportion_censored = 0
        
        return(..to_return)
        
    }, .survival_times, .covariates)
    
    #create randomly censored model-ready dataset
    .data_censored <-  llply(1:3, function(..censor, ..survival_times, ..covariates, ..n_subjects) {
      ..proportion_censored <- c(0.2, 0.5, 0.8)[..censor]
      ..events <- rbinom(.n_subjects, size = 1, p = (1 - ..proportion_censored) )
      
      ldply(1:..n_subjects, function(...subject, ...survival_times, ...covariates, ...events, ...proportion_censored) {
        ...survival_time <- ceiling(...survival_times$g_y[...subject])
        ...to_return <- ...covariates[[...subject]]
        ...to_return$id <- ...subject
        ...to_return$t <- c(1:...survival_time)
        ...to_return$t0 <- 0:(...survival_time - 1)
        ...to_return$d <- ...events[...subject]
        ...to_return$proportion_censored = ...proportion_censored
        return(...to_return, ..proportion_censored)
        
      }, ..survival_times, ..covariates, ..events)
    }, .survival_times, .covariates, .n_subjects)
      
    
    #initialize model container list
    .uncensored_results <- list()
    #fit uncensored models
    
    .cox_formula <- as.formula(paste0("Surv(t0, t, d_none) ~ ", paste0(names(.beta), collapse = " + ")))
    .unpooled_cox_results <- data.frame(summary(coxph(.cox_formula, data = .data_none))$coef[ ,c(1,3,5)])
    names(.unpooled_cox_results) <- c("coef", "se", "p")
    .unpooled_cox_results$var <- row.names(.unpooled_cox_results)
    .unpooled_cox_results$type <- "Unpooled Cox"
    .unpooled_cox_results$proportion_censored <- 0
    
    .uncensored_results[[1]] <- .unpooled_cox_results
    
    .pooled_cox_formula <- as.formula(paste0("Surv(t0, t, d_none) ~ ", paste0(names(.beta), collapse = " + "), " + cluster(id) "))
    .pooled_cox_results <- data.frame(summary(coxph(.pooled_cox_formula, data = .data_none), )$coef[ ,c(1,4,6)])
    names(.pooled_cox_results) <- c("coef", "se", "p")
    .pooled_cox_results$var <- row.names(.pooled_cox_results)
    .pooled_cox_results$type <- "Pooled Cox"
    .pooled_cox_results$proportion_censored <- 0
    
    .uncensored_results[[1]] <- .pooled_cox_formula
    
    .uncensored_results <- ldply(.uncensored_results)
    
    .censored_results <- ldply(.data_censored, function(..data_censored, ..unpooled_cox_formula, ..pooled_cox_formula) {
      
      ..proportion_censored <- ..data_censored$proportion_censored[1]
      
      ..censored_results <- list()
      
      ..unpooled_cox_results <- data.frame(summary(coxph(..unpooled_cox_formula, data = ..data_censored))$coef[ ,c(1,3,5)])
      names(..unpooled_cox_results) <- c("coef", "se", "p")
      ..unpooled_cox_results$var <- row.names(..unpooled_cox_results)
      ..unpooled_cox_results$type <- "Unpooled Cox"
      ..unpooled_cox_results$proportion_censored <- ..proportion_censored
      
      ..censored_results[[1]] <- ..unpooled_cox_results
      
      ..pooled_cox_results <- data.frame(summary(coxph(..pooled_cox_formula, data = ..data_censored), )$coef[ ,c(1,4,6)])
      names(..pooled_cox_results) <- c("coef", "se", "p")
      ..pooled_cox_results$var <- row.names(.pooled_cox_results)
      ..pooled_cox_results$type <- "Pooled Cox"
      ..pooled_cox_results$proportion_censored <- ..proportion_censored
      
      ..censored_results[[1]] <- ..pooled_cox_results
      
      ldply(..censored_results)
    })
    to_write <- rbind(..uncensored_results, ..censored_results)
    
    write.table(to_write, file = paste0(.results_write_directory, WORKER.ID, "_", .sim_results_name, ".csv"), 
                append = T, row.names = F, col.names = F, sep = ",")
  } else {
    print("This won't print, but that's cool, right?")
  }
  
}, covariates, n_subjects, beta, results_write_directory, sim_results_name)

#setwd("/Volumes/QSU/Datasets/PCORI/data simulation/data/Sim 9/")

all_results <- read.table(file = "bootstrapResults.csv", stringsAsFactors = F, header = T, sep = ",")
all_results$stat <- rep(c("coef", "se", "p"), n_subjects)
all_results[1,]




summarized_results <- ldply(names(beta), function(.var, .all_results, .beta) {
  .estimates <- .all_results[.all_results$stat == "coef", ]
  .ses <- .all_results[.all_results$stat == "se", ]
  .ps <- .all_results[.all_results$stat == "p", ]
  
  return(data.frame(var = .var,
                    true = round(.beta[[.var]], 4),
                    mean_estimate = mean(.estimates[[.var]], na.rm = T),
                    n_converge = sum(!is.na(.estimates[[.var]])),
                    bias = mean(.estimates[[.var]], na.rm = T) - .beta[[.var]],
                    std_bias = 100 * (mean(.estimates[[.var]], na.rm = T) - .beta[[.var]]) / mean(.ses[[.var]], na.rm = T),
                    mse = mean( (mean(.estimates[[.var]], na.rm = T) - .beta[[.var]])^2),
                    reject_prob = mean(.ps[[.var]] < 0.05, na.rm = T),
                    stringsAsFactors = F)
  )}, all_results, beta)
summarized_results



names(beta)[13]



summary(coxph(Surv(t0, t, d_none) ~  z1 + z2, data = subject_data))$coef[ ,c(1,3,5)]

data_none[1,]
names(data_none) <- c('id','z1','z2','t','t0','d_none')


source_url('https://raw.github.com/kikapp/GeneralRScripts/master/MiscellaneousRFunctions.R')


