
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
cox.predictors = c("ind_cd4_50_100", "ind_cd4_350_500", "ind_cd4_200_350", "ind_cd4_100_200", "d_dida")
na.methods = c("complete.case", "naive", "frailty", "log-t")
write.path = "/share/PI/manishad/naEst/output"
miss.matrix = read.csv("/share/PI/manishad/naEst/missing_var_parameters_matrix.csv")

# code
source("/share/PI/manishad/naEst/r/na_est_functions.R")
source("/share/PI/manishad/naEst/r/impose_missingness_functions.R")

# read in complete dataset
d = read.csv(data.path)
file.name = substring(data.path, 54)  # start at 54th character to erase annoying beginning of path

# RUN BY MANISHA!
# don't impute with anything that's not in the Cox model
impute.with = c("id", "d", cox.predictors)

make.miss.if.contains="cd4"

do_one_dataset(.d=d, .source.file.name=file.name, .miss.matrix=miss.matrix, .time.name=time.name,
               .event.name=event.name, .cluster.name=cluster.name,
               .cox.predictors=cox.predictors, .name.prefix=name.prefix,
               .na.methods=na.methods, .write.path=write.path,
               .impute.with = impute.with,
               .make.miss.if.contains=make.miss.if.contains
)



  