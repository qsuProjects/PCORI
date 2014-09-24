if (FALSE) {
  rm(list = ls())
  source("~/.Rprofile")
}
args <- commandArgs(trailingOnly = TRUE)
print(args)

input_location <- args[1]
n_clusters <-  as.numeric(args[2])
n_segs <- as.numeric(args[3])
seg <- as.numeric(args[4])
output_name_stem <- args[5]
log_name_stem <- args[6]

if(is.na(input_location)) {input_location <- "/scratch/users/kikapp/dump/"}
if(is.na(n_clusters)) {n_clusters <- 16}
if(is.na(n_segs)) {n_segs <- 1}
if(is.na(seg)) {seg <- 1}


library(plyr)
require(doSNOW)
cl<-makeCluster(n_clusters)
registerDoSNOW(cl)

n_clusters
getDoParWorkers()

clusterApply(cl, seq(along=cl), function(id, scenario) WORKER.ID <<- paste0("job", scenario, "_core", id), seg)



  # input_location <- paste0("/scratch/users/kikapp/PCORI/data/lim1_large_betas")
  
  ##get list of available data files;
  data_files <- list.files(paste0(input_location, "/data/"))
  data_files <- data_files[grepl("data", data_files)]
  # data_files <- paste0(input_location, "/", data_files)
  length(data_files)
  print(data_files[1])
  
  total_files <- length(data_files)
  
  segments <- ceiling(quantile(0:(total_files-1), seq(0, 1, length.out = (n_segs+1))))
  
  #define specified range of file indices
  analysis_range <- c(segments[seg],segments[seg+1]-1)
  #if segment is rightmost, add last file to end
  if (seg == n_segs) { analysis_range[2] <- analysis_range[2] + 1}
  cat(paste0("Log file for seg: ", seg, "\n",
             "Working on analysis range: ", analysis_range[1], "-", analysis_range[2], "\n"),
      file = paste0(input_location, "/model_fit/log/", seg, "_log.txt"),
      append = F)
  
  
  
  system.time({
    l_ply(c(analysis_range[1]:analysis_range[2]), .parallel=T, function(.item, .data_files, .n_clusters, .input_location, .segnum, .nsegs) {
      #     print(paste0(.item))
      .data_file <- .data_files[.item+1]
      library(msm)
      library(survival)
      library(ggplot2)
      library(reshape2)
      library(geepack)
      library(gee)
      file_location <-  paste0(.input_location, "/data/")
      output_dir <- paste0(.input_location, "/model_fit/results/")
      output_file <- paste0(output_dir, .data_file)
      progress_dir <- paste0(.input_location, "/model_fit/log/")
      progress_file <- paste0(progress_dir, "/progress_", WORKER.ID)
#       progress_file <- paste0(progress_dir, "/progress_", .data_file)
      
      cat(paste0("seg ", .segnum, " file ", .item, ": ", paste0(file_location, "/", .data_file), " started at ", Sys.time(), "\n"), file = progress_file, append = T)
      
      #write header
      #       if ((.item %% .n_clusters) == .item) {
      cat("full_model,estimate,se,p,var,data_file,read_time,write_time,filenum,segnum\n", file = output_file, append = F)
      #       }
      #     cat( paste0(.input_location, "/", .data_file), file = output_file, append = T)
      
      #     temp_results <- list()
      subject_data <- read.table(file = paste0(file_location, "/", .data_file), sep = ",", header = T, stringsAsFactors = F)
      read_time  <- Sys.time()
      cat(paste0("seg ", .segnum, " file ", .item, " data read at ", Sys.time(), "\n"), file = progress_file, append = T)
      
      #parse scenario details
      #     if (SCENARIO_1) {
      #     scenario <- strsplit(gsub("data|_", "", .data_file), split = "")[[1]]
      #     g <- scenario[1]
      #     lim <-scenario[2]
      #     w = scenario[3]
      #     md = scenario[4]
      #     rv <- scenario[5]
      #     set <- gsub("data[:0-9:]{5}_|_[:a-z:]{1,}[:0-9:]{0,1}|.csv", "", .data_file)
      #     }
      #     if (SCENARIO_2) {
      #       scenario <- strsplit(gsub("data|_", "", .data_file), split = "")[[1]]
      #       g <- scenario[1]
      #       lim <-scenario[2]
      #       w = scenario[3]
      #       md = scenario[4]
      #       rv <- scenario[5]
      #       set <- gsub("data[:0-9:]{5}_|_[:a-z:]{1,}[:0-9:]{0,1}|.csv", "", .data_file)
      #     }
      #     if (SCENARIO_3) {
      #       scenario <- strsplit(gsub("data|_", "", .data_file), split = "")[[1]]
      #       g <- scenario[1]
      #       lim <-scenario[2]
      #       w = scenario[3]
      #       md = scenario[4]
      #       rv <- scenario[5]
      #       set <- gsub("data[:0-9:]{5}_|_[:a-z:]{1,}[:0-9:]{0,1}|.csv", "", .data_file)
      #     }
      result_list <- list()
      
      #if data are 'random' censoring
      if ( !grepl("trad|admin", .data_file) ) {    
        #fit models
        cat(paste0("seg ", .segnum, " file ", .item, " fitting models at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["cox_model_nocensor"]] <-    summary(coxph(Surv(t0, t, d_none) ~  z1 + z2, data = subject_data))$coef[ ,c(1,3,5)]
        result_list[["cox_model_random1"]] <-    summary(coxph(Surv(t0, t, d_random1) ~  z1 + z2, data = subject_data))$coef[ ,c(1,3,5)]
        result_list[["cox_model_random2"]] <-    summary(coxph(Surv(t0, t, d_random2) ~  z1 + z2, data = subject_data))$coef[ ,c(1,3,5)]
        result_list[["cox_model_random3"]] <-    summary(coxph(Surv(t0, t, d_random3) ~  z1 + z2, data = subject_data))$coef[ ,c(1,3,5)]
        
        cat(paste0("seg ", .segnum, " file ", .item, " cox models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["pooled_cox_model_nocensor"]] <-    summary(coxph(Surv(t0, t, d_none) ~  z1 + z2 + cluster(id), data = subject_data), )$coef[ ,c(1,4,6)]
        result_list[["pooled_cox_model_random1"]] <-    summary(coxph(Surv(t0, t, d_random1) ~  z1 + z2 + cluster(id), data = subject_data))$coef[ ,c(1,4,6)]
        result_list[["pooled_cox_model_random2"]] <-    summary(coxph(Surv(t0, t, d_random2) ~  z1 + z2 + cluster(id), data = subject_data))$coef[ ,c(1,4,6)]
        result_list[["pooled_cox_model_random3"]] <-    summary(coxph(Surv(t0, t, d_random3) ~  z1 + z2 + cluster(id), data = subject_data))$coef[ ,c(1,4,6)]
        
        cat(paste0("seg ", .segnum, " file ", .item, " pooled cox models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["poisson_model_nocensor"]] <-    summary(glm(data = subject_data, d_none ~ z1 + z2, family = "poisson"))$coef[-1,c(1,2,4)]
        result_list[["poisson_model_random1"]] <-    summary(glm(data = subject_data, d_random1 ~ z1 + z2, family = "poisson"))$coef[-1,c(1,2,4)]
        result_list[["poisson_model_random2"]] <-   summary(glm(data = subject_data, d_random2 ~ z1 + z2, family = "poisson"))$coef[-1,c(1,2,4)]
        result_list[["poisson_model_random3"]] <-    summary(glm(data = subject_data, d_random3 ~ z1 + z2, family = "poisson"))$coef[-1,c(1,2,4)]
        
        cat(paste0("seg ", .segnum, " file ", .item, " poisson models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["pooled_poisson_model_nocensor"]] <-    summary(geeglm(data = subject_data, d_none ~ z1 + z2, id = id, corstr = "independence", family = "poisson"))$coef[-1,c(1,2,4)]
        result_list[["pooled_poisson_model_random1"]] <-   summary(geeglm(data = subject_data, d_random1 ~ z1 + z2, id = id, corstr = "independence", family = "poisson"))$coef[-1,c(1,2,4)]
        result_list[["pooled_poisson_model_random2"]] <-   summary(geeglm(data = subject_data, d_random2 ~ z1 + z2, id = id, corstr = "independence", family = "poisson"))$coef[-1,c(1,2,4)]
        result_list[["pooled_poisson_model_random3"]] <-    summary(geeglm(data = subject_data, d_random3 ~ z1 + z2, id = id, corstr = "independence", family = "poisson"))$coef[-1,c(1,2,4)]    
        
        cat(paste0("seg ", .segnum, " file ", .item, " pooled poisson models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["logistic_model_nocensor"]] <-    summary(glm(data = subject_data, d_none ~ z1 + z2, family = "binomial"))$coef[-1,c(1,2,4)]
        result_list[["logistic_model_random1"]] <-    summary(glm(data = subject_data, d_random1 ~ z1 + z2, family = "binomial"))$coef[-1,c(1,2,4)]
        result_list[["logistic_model_random2"]] <-   summary(glm(data = subject_data, d_random2 ~ z1 + z2, family = "binomial"))$coef[-1,c(1,2,4)]
        result_list[["logistic_model_random3"]] <-    summary(glm(data = subject_data, d_random3 ~ z1 + z2, family = "binomial"))$coef[-1,c(1,2,4)]
        
        cat(paste0("seg ", .segnum, " file ", .item, " logistic models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["pooled_logistic_model_nocensor"]] <-    summary(geeglm(data = subject_data, d_none ~ z1 + z2, id = id, corstr = "independence", family = "binomial"))$coef[-1,c(1,2,4)]
        result_list[["pooled_logistic_model_random1"]] <-    summary(geeglm(data = subject_data, d_random1 ~ z1 + z2, id = id, corstr = "independence", family = "binomial"))$coef[-1,c(1,2,4)]
        result_list[["pooled_logistic_model_random2"]] <-   summary(geeglm(data = subject_data, d_random2 ~ z1 + z2, id = id, corstr = "independence", family = "binomial"))$coef[-1,c(1,2,4)]
        result_list[["pooled_logistic_model_random3"]] <-    summary(geeglm(data = subject_data, d_random3 ~ z1 + z2, id = id, corstr = "independence", family = "binomial"))$coef[-1,c(1,2,4)]
        
        cat(paste0("seg ", .segnum, " file ", .item, " models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        #parse and join results
        results_df <- ldply(result_list, function(.df) {      
          .df <- data.frame(.df)
          names(.df) <- c("estimate", "se", "p")
          .df$var <- rownames(.df)
          return(.df)
        })
        #       results_df$g <- g
        #       results_df$lim <- lim
        #       results_df$w <- w
        #       results_df$md <- md
        #       results_df$rv <- rv
        #       results_df$set <- set
        #       results_df$pc <- as.numeric( gsub("_|[:a-z:]", "0", results_df$.id))
        #       results_df$censor_type <-  gsub("^[:a-z:]{1,}_[:a-z:]{1,}_|model_|[:0-9:]", "", results_df$.id) 
        #       results_df$model <- gsub("_|model|random|nocensor|[:0-9:]", "", results_df$.id)   
        results_df$data_file <- .data_file
        write_time  <- Sys.time()
        results_df$read_time <- read_time
        results_df$write_time <- write_time
        results_df$filenum <- .item
        results_df$segnum <- .segnum
        cat(paste0("seg ", .segnum, " file ", .item, " results parsed at ", Sys.time(), "\n"), file = progress_file, append = T)
        cat(paste0("seg ", .segnum, " file ", .item, " writing results to ", output_file, " at ", Sys.time(), "\n"), file = progress_file, append = T)
        write.table(results_df, file = output_file, sep = ",", append=T, row.names = F, col.names = F)
        cat(paste0("seg ", .segnum, " file ", .item, " results written at ", Sys.time(), "\n"), file = progress_file, append = T)
        #       temp_results[[.data_file]] <- results_df
      }
      
      #if data are admin or traditional
      if ( grepl("trad|admin", .data_file) ) {
        
        cat(paste0("seg ", .segnum, " file ", .item,  " fitting models at ", Sys.time(), "\n"), file = progress_file, append = T)
        #fit models
        result_list[["cox_model"]] <-    summary(coxph(Surv(t0, t, d) ~  z1 + z2, data = subject_data))$coef[ ,c(1,3,5)]
        cat(paste0("seg ", .segnum, " file ", .item, " cox models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["pooled_cox_model"]] <-    summary(coxph(Surv(t0, t, d) ~  z1 + z2 + cluster(factor(id)), data = subject_data))$coef[ ,c(1,4,6)]
        cat(paste0("seg ", .segnum, " file ", .item, " pooled cox models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["poisson_model"]] <-    summary(glm(data = subject_data, d ~ z1 + z2, family = "poisson"))$coef[-1,c(1,2,4)]
        cat(paste0("seg ", .segnum, " file ", .item, " poisson models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["pooled_poisson_model"]] <-    summary(geeglm(data = subject_data, d ~ z1 + z2, id = id, corstr = "independence", family = "poisson"))$coef[-1,c(1,2,4)]
        cat(paste0("seg ", .segnum, " file ", .item, " pooled poisson models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["logistic_model"]] <-    summary(glm(data = subject_data, d ~ z1 + z2, family = "binomial"))$coef[-1,c(1,2,4)]
        cat(paste0("seg ", .segnum, " file ", .item, " logistic models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        result_list[["pooled_logistic_model"]] <-    summary(geeglm(data = subject_data, d ~ z1 + z2, id = id, corstr = "independence", family = "binomial"))$coef[-1,c(1,2,4)]
        cat(paste0("seg ", .segnum, " file ", .item, " models fit at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        #     summary(gee(data = subject_data, d ~ z1 + z2, id = id, corstr = "independence", family = "binomial"))$coef
        #     summary(geeglm(data = subject_data, d ~ z1 + z2, id = id, corstr = "independence", family = "binomial"))
        #     summary(glm(data = subject_data, d ~ z1 + z2, family = "binomial"))
        
        #parse and join results
        results_df <- ldply(result_list, function(.df) {      
          .df <- data.frame(.df)
          names(.df) <- c("estimate", "se", "p")
          .df$var <- rownames(.df)
          return(.df)
        })
        #       results_df$g <- g
        #       results_df$lim <- lim
        #       results_df$w <- w
        #       results_df$md <- md
        #       results_df$rv <- rv
        #       results_df$set <- set
        #       results_df$pc <- as.numeric( gsub("[:a-z:]{1,}", "", gsub("data[:0-9:]{1,}_[:0-9:]{1,}_|.csv", "", .data_file)))
        #       results_df$censor_type <- gsub("[:0-9:]{1,}", "", gsub("data[:0-9:]{1,}_[:0-9:]{1,}_|.csv", "", .data_file))
        #       results_df$model <- gsub("_model", "", results_df$.id)    
        results_df$data_file <- .data_file
        write_time  <- Sys.time()
        results_df$read_time <- read_time
        results_df$write_time <- write_time        
        results_df$filenum <- .item
        results_df$segnum <- .segnum
        
        cat(paste0("seg ", .segnum, " file ", .item, " results parsed at ", Sys.time(), "\n"), file = progress_file, append = T)
        cat(paste0("seg ", .segnum, " file ", .item, " writing results to ", output_file, " at ", Sys.time(), "\n"), file = progress_file, append = T)
        
        write.table(results_df, file = output_file, sep = ",", append=T, row.names = F, col.names = F)    
        
        cat(paste0("seg ", .segnum, " file ", .item, " written to  ", output_file, " at ", Sys.time(), "\n"), file = progress_file, append = T)
      }
      
    }, data_files, n_clusters, input_location, seg, n_segs)
  })
