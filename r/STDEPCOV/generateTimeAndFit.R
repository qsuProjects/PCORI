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
#range_low: Integer, first rep to do
#range_high: Integer, last rep to do
#sim_results_name: String, file name for simulation results


require(msm)
require(survival)
require(plyr)

#load data
covariates <- read.table(file = covariate_data_path, sep = ",", header = T, stringsAsFactors = F)

#generate dummy variables for race
covariates$race_black <- 0
covariates$race_other <- 0

covariates$race_black[covariates$racecatexp == 1] <- 1
covariates$race_other[covariates$racecatexp == 2] <- 1

#make age constant for each subject
covariates <- ddply(covariates, .(id), function(.df) {
  .df$age <- .df$age[1]
  return(.df)
})

#determine how many subjects are in the data
n_subjects <- length(unique(covariates$id))

#set up parallelization
require(doSNOW)
cl<-makeCluster(number_cores)
registerDoSNOW(cl)

number_cores
getDoParWorkers()

clusterApply(cl, seq(along=cl), function(.id, .results_write_directory, .sim_results_name, .the_seed) {
  WORKER.ID <<- paste0("core_", .id)
  set.seed(.the_seed*(2^.id))
  cat("coef,se,p,var,type,proportion_censored,rep\n", file = paste0(.results_write_directory, "/", WORKER.ID, "_", .sim_results_name, ".csv"), append = F)
}, results_write_directory, sim_results_name, the_seed)

l_ply(range_low:range_high, .parallel = T, function(.rep, .covariates, .n_subjects, .beta, .results_write_directory, .sim_results_name, .log_write_directory) {
  
  #load required packages
  require(msm)
  require(survival)
  require(plyr)
  
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
    .lambda*(x^.nu)
  }  
  
  .t <- 1:350
  .t_diff <- (.t[-1] - .t[1:(length(.t) - 1)])[-(length(.t) - 1)]
  .g_inv_t <- .g_inv(.t)
  .g_inv_t_diff <- (.g_inv(.t[-1]) - .g_inv(.t[1:(length(.t) - 1)]))[-(length(.t) - 1)]
  
  #CREATING THE BOUNDS OF TRUNCATION
  .t_max <- 350
  .t_min <- 20
  
  .g_inv_t_max <- .g_inv(.t_max)
  .g_inv_t_min <- .g_inv(.t_min)
  
  
  #K function applies ACCEPT-REJECT algorithm
  .k <- function(..x, ..m, ..M, ..rates, ..t){
    ifelse(..x <= ..m | ..x >= ..M, 0, dpexp(..x, ..rates, ..t))
  }
  
  #define survival time generation function
  .gen_y <- function(.x, .g_inv_t,  .g_inv_t_min, .g_inv_t_max) {
    .x1 <- .x$xB
    .M <- 1 / .d
    .r <- 60
    .count<-0
    #counter of times repeat is run
    while (.count<1000) {
      .count <- .count+1
      .y <- rpexp(.r, .x1, .g_inv_t)
      .u <- runif(.r)
      .t <- .M * (.k(.y, .g_inv_t_min, .g_inv_t_max, .x1, .g_inv_t) / .d / dpexp(.y, .x1, .g_inv_t))
      .y <- .y[.u <= .t][1]
      if (!is.na(.y)) {break}
    }
    .y
  }
  
  cat(paste0("Covariates processed for rep ", .rep, " at ", Sys.time(), "\n"), file = paste0(.log_write_directory, "/", WORKER.ID, ".txt"), append = T)
  .survival_times <- ldply(.covars_list, .gen_y, .g_inv_t,  .g_inv_t_min, .g_inv_t_max)
  cat(paste0("Survival times generated for rep ", .rep, " at ", Sys.time(), "\n"), file = paste0(.log_write_directory, "/", WORKER.ID, ".txt"), append = T)
  
  .survival_times$g_y <- .g(.survival_times[ ,2])
  
  if (sum(is.na(.survival_times$V1)) == 0) {
    
    #create uncensored model-ready dataset
    .data_none <-  ldply(1:.n_subjects, function(..subject, ..survival_times, ..covars_list) {
      ..survival_time <- ceiling(..survival_times$g_y[..subject])
      ..to_return <- ..covars_list[[..subject]][1:..survival_time, ]
      ..to_return$id <- ..subject
      ..to_return$t <- c(1:..survival_time)
      ..to_return$t0 <- 0:(..survival_time - 1)
      ..to_return$d <- c( rep(0, ..survival_time - 1), 1)
      ..to_return$proportion_censored = 0
      
      return(..to_return)
      
    }, .survival_times, .covars_list)
    
    cat(paste0("Uncensored data set generated for rep ", .rep, " at ", Sys.time(), "\n"), file = paste0(.log_write_directory, "/", WORKER.ID, ".txt"), append = T)
    
    #create randomly censored model-ready dataset
    .data_censored <-  llply(1:3, function(..censor, ..survival_times, ..covars_list, ..n_subjects) {
      ..proportion_censored <- c(0.2, 0.5, 0.8)[..censor]
      ..events <- rbinom(.n_subjects, size = 1, p = (1 - ..proportion_censored) )
      
      ldply(1:..n_subjects, function(...subject, ...survival_times, ...covars_list, ...events, ...proportion_censored) {
        ...survival_time <- ceiling(...survival_times$g_y[...subject])
        ...to_return <- ...covars_list[[...subject]][1:...survival_time, ]
        ...to_return$id <- ...subject
        ...to_return$t <- c(1:...survival_time)
        ...to_return$t0 <- 0:(...survival_time - 1)
        ...to_return$d <- ...events[...subject]
        ...to_return$proportion_censored = ...proportion_censored
        return(...to_return)
        
      }, ..survival_times, ..covars_list, ..events, ..proportion_censored)
    }, .survival_times, .covars_list, .n_subjects)
    
    cat(paste0("Censored data sets generated for rep ", .rep, " at ", Sys.time(), "\n"), file = paste0(.log_write_directory, "/", WORKER.ID, ".txt"), append = T)
    
    #initialize model container list
    .uncensored_results <- list()
    
    #fit uncensored models
    .unpooled_cox_formula <- as.formula(paste0("Surv(t0, t, d) ~ ", paste0(names(.beta), collapse = " + ")))
    .unpooled_cox_results <- data.frame(summary(coxph(.unpooled_cox_formula, data = .data_none))$coef[ ,c(1,3,5)])
    names(.unpooled_cox_results) <- c("coef", "se", "p")
    .unpooled_cox_results$var <- row.names(.unpooled_cox_results)
    .unpooled_cox_results$type <- "Unpooled Cox"
    .unpooled_cox_results$proportion_censored <- 0
    
    .uncensored_results[[1]] <- .unpooled_cox_results
    
    .pooled_cox_formula <- as.formula(paste0("Surv(t0, t, d) ~ ", paste0(names(.beta), collapse = " + "), " + cluster(id) "))
    .pooled_cox_results <- data.frame(summary(coxph(.pooled_cox_formula, data = .data_none), )$coef[ ,c(1,4,6)])
    names(.pooled_cox_results) <- c("coef", "se", "p")
    .pooled_cox_results$var <- row.names(.pooled_cox_results)
    .pooled_cox_results$type <- "Pooled Cox"
    .pooled_cox_results$proportion_censored <- 0
    
    .uncensored_results[[2]] <- .pooled_cox_results
    
    .uncensored_results <- ldply(.uncensored_results)
    
    cat(paste0("Uncensored models fit for rep ", .rep, " at ", Sys.time(), "\n"), file = paste0(.log_write_directory, "/", WORKER.ID, ".txt"), append = T)
    
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
      
      ..censored_results[[2]] <- ..pooled_cox_results
      
      ldply(..censored_results)
    }, .unpooled_cox_formula, .pooled_cox_formula)
    
    cat(paste0("Censored models fit for rep ", .rep, " at ", Sys.time(), "\n"), file = paste0(.log_write_directory, "/", WORKER.ID, ".txt"), append = T)
    
    .to_write <- rbind(.uncensored_results, .censored_results)
    .to_write$rep <- .rep
    
    write.table(.to_write, file = paste0(.results_write_directory, "/", WORKER.ID, "_", .sim_results_name), 
                append = T, row.names = F, col.names = F, sep = ",")
    
    cat(paste0("Results written for rep ", .rep, " at ", Sys.time(), "\n"), file = paste0(.log_write_directory, "/", WORKER.ID, ".txt"), append = T)
    
  } else {
    cat(paste0("Convergence fail for rep ", .rep, " at ", Sys.time(), "\n"), file = paste0(.log_write_directory, "/", WORKER.ID, ".txt"), append = T)
  }
  
}, covariates, n_subjects, beta, results_write_directory, sim_results_name, log_write_directory)
