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

if(is.na(input_location)) {input_location <- "/scratch/users/kikapp/PCORI/data/small_beta_corr/"}
if(is.na(n_clusters)) {n_clusters <- 16}
if(is.na(n_segs)) {n_segs <- 5}
if(is.na(seg)) {seg <- 1}


library(plyr)
require(doSNOW)
cl<-makeCluster(n_clusters)
registerDoSNOW(cl)

n_clusters
getDoParWorkers()

clusterApply(cl, seq(along=cl), function(id, scenario) WORKER.ID <<- paste0("job_", scenario, "_core_", id), seg)



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
# cat(paste0("Log file for seg: ", seg, "\n",
#            "Working on analysis range: ", analysis_range[1], "-", analysis_range[2], "\n"),
#     file = paste0(input_location, "/model_fit/log/", seg, "_log.txt"),
#     append = F)



system.time({
  l_ply(c(analysis_range[1]:analysis_range[2]), .parallel=T, function(.item, .data_files, .n_clusters, .input_location, .segnum, .nsegs, .first) {
    #     print(paste0(.item))
    .data_file <- .data_files[.item+1]
    #     WORKER.ID <- "TESTWORKERID"
    file_location <-  paste0(.input_location, "/data/")
    output_dir <- paste0(.input_location, "/runtime/output/")
    output_file <- paste0(output_dir, "/runtimes_", WORKER.ID, ".csv")
    progress_dir <- paste0(.input_location, "/runtime/log/")
    progress_file <- paste0(progress_dir, "/log_", WORKER.ID, ".log")
    
#     if (.item == .first) {
#       cat("scenario,runtime_S,runtime_none,runtime_random,runtime_trad,runtime_admin,n_runtimes,data_file,read_time,write_time,filenum,segnum\n", file = output_file, append = F)
#     }
    subject_data <- read.table(file = paste0(file_location, "/", .data_file), sep = ",", header = T, stringsAsFactors = F)
    read_time  <- Sys.time()
    cat(paste0("seg ", .segnum, " file ", .item, " data (", .data_file, ") read at ", Sys.time(), "\n"), file = progress_file, append = T)
    
    #if data are traditional censoring
    if ( grepl("trad", .data_file) ) {    
      #fit models
      write_time  <- Sys.time()
      results_df <- data.frame(scenario = paste0(unlist(strsplit(.data_file, split = ""))[5:10], collapse = ""), stringsAsFactors = F)
      results_df$runtime_S <- NA
      results_df$runtime_none <- NA
      results_df$runtime_random <- NA
      results_df$runtime_trad <- subject_data$runtime[1]
      results_df$runtime_admin <- NA
      results_df$n_runtimes <- length(unique(subject_data$runtime))
      results_df$data_file <- .data_file
      results_df$read_time <- read_time
      results_df$write_time <- write_time
      results_df$filenum <- .item
      results_df$segnum <- .segnum
      
      cat(paste0("seg ", .segnum, " file ", .item, " writing results to ", output_file, " at ", Sys.time(), "\n"), file = progress_file, append = T)
      write.table(results_df, file = output_file, sep = ",", append=T, row.names = F, col.names = F)
      cat(paste0("seg ", .segnum, " file ", .item, " results written at ", Sys.time(), "\n"), file = progress_file, append = T)
      #       temp_results[[.data_file]] <- results_df
    }
    #if data are admin censoring
    if ( grepl("admin", .data_file) ) {    
      #fit models
      write_time  <- Sys.time()
      results_df <- data.frame(scenario = paste0(unlist(strsplit(.data_file, split = ""))[5:10], collapse = ""), stringsAsFactors = F)
      results_df$runtime_S <- NA
      results_df$runtime_none <- NA
      results_df$runtime_random <- NA
      results_df$runtime_trad <- NA
      results_df$runtime_admin <- subject_data$runtime[1]
      results_df$n_runtimes <- length(unique(subject_data$runtime))
      results_df$data_file <- .data_file
      results_df$read_time <- read_time
      results_df$write_time <- write_time
      results_df$filenum <- .item
      results_df$segnum <- .segnum
      
      cat(paste0("seg ", .segnum, " file ", .item, " writing results to ", output_file, " at ", Sys.time(), "\n"), file = progress_file, append = T)
      write.table(results_df, file = output_file, sep = ",", append=T, row.names = F, col.names = F)
      cat(paste0("seg ", .segnum, " file ", .item, " results written at ", Sys.time(), "\n"), file = progress_file, append = T)
      #       temp_results[[.data_file]] <- results_df
    }
    #if data are admin or traditional
    if ( !grepl("trad|admin", .data_file) ) {
      write_time  <- Sys.time()
      results_df <- data.frame(scenario = paste0(unlist(strsplit(.data_file, split = ""))[5:10], collapse = ""), stringsAsFactors = F)
      results_df$runtime_S <-  subject_data$runtime_S[1]
      results_df$runtime_none <-  subject_data$runtime_none[1]
      results_df$runtime_random <-  subject_data$runtime_random[1]
      results_df$runtime_trad <- NA
      results_df$runtime_admin <-NA
      results_df$n_runtimes <- length(unique(subject_data$runtime_S))
      results_df$data_file <- .data_file
      results_df$read_time <- read_time
      results_df$write_time <- write_time
      results_df$filenum <- .item
      results_df$segnum <- .segnum
      
      cat(paste0("seg ", .segnum, " file ", .item, " writing results to ", output_file, " at ", Sys.time(), "\n"), file = progress_file, append = T)
      write.table(results_df, file = output_file, sep = ",", append=T, row.names = F, col.names = F)
      cat(paste0("seg ", .segnum, " file ", .item, " results written at ", Sys.time(), "\n"), file = progress_file, append = T)
    }
    
  }, data_files, n_clusters, input_location, seg, n_segs, analysis_range[1])
})
