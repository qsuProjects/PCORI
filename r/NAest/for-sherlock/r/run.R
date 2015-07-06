
########################### LOAD COMMAND-LINE ARGUMENTS ########################### 

# load command line arguments
args = commandArgs(trailingOnly = TRUE)
print(args)

data.path = args[1]
name.prefix = args[2]


############################## SET LESS COMMON ARGUMENTS ##############################

time.name="t"
event.name="d"
cluster.name="id"
cox.predictors = c("d_abac")
na.methods = c("naive", "frailty", "log-t")
write.path = "/share/PI/manishad/naEst/output"
miss.matrix = read.csv("/share/PI/manishad/naEst/missing_var_parameters_matrix.csv")

# RUN BY MANISHA!
dont.impute.with = c("X", "bmi", "bmi_slope", "cd4", "cd4_slope", "log_vln_slope",
                     "ldl_slope", "hdl_slope", "log_vln_5", "log_vln_4",
                     "log_vln_3", "log_vln_2", "cd4_cuts", "racecatexp",
                     "linpred", "frailty", "xB", "t", "t0", "proportion_censored",
                     "source_file", "vln_cat", "bmi_cuts", "log_vln",
                     
                     "ind_cd4_350_500", "ind_cd4_200_350", "ind_cd4_100_200",
                     "ind_cd4_50_100", "ind_bmi_gt_30", 
                     
                     "ind_bmi_25_30", "ind_bmi_lt_20"
                     )


# code
source("/share/PI/manishad/naEst/r/na_est_functions.R")
source("/share/PI/manishad/naEst/r/impose_missingness_functions.R")

# read in complete dataset
d = read.csv(data.path)

do_one_dataset(.d=d, .miss.matrix=miss.matrix, .time.name=time.name,
                          .event.name=event.name, .cluster.name=cluster.name,
                          .cox.predictors=cox.predictors, .name.prefix=name.prefix,
                          .na.methods=na.methods, .write.path=write.path,
                          .dont.impute.with = dont.impute.with
               )
  