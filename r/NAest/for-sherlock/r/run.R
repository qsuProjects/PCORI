
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
cox.predictors = c("X")
na.methods = c("complete.case", "naive", "frailty", "log-t", "full")
write.path = "/share/PI/manishad/naEst/output"
miss.matrix.hi = read.csv("/share/PI/manishad/naEst/missing_var_parameters_matrix_good_survivors.csv")
miss.matrix.lo = read.csv("/share/PI/manishad/naEst/missing_var_parameters_matrix_bad_survivors.csv")
aux.matrix = read.csv("/share/PI/manishad/naEst/aux_var_parameters_matrix_Z.csv")
impute.with = c("id", "d", "Z", cox.predictors)
make.miss.if.contains="X"


# code
source("/share/PI/manishad/naEst/r/na_est_functions.R")
source("/share/PI/manishad/naEst/r/impose_missingness_functions.R")

# read in complete dataset
d = read.csv(data.path)
file.name = substring(data.path, 58)  # start at 58th character to erase annoying beginning of path

do_one_dataset(.d=d, .source.file.name=file.name, 
               .miss.matrix.hi=miss.matrix.hi, 
               .miss.matrix.lo=miss.matrix.lo,
               .aux.matrix=aux.matrix,
               .time.name=time.name,
               .event.name=event.name, .cluster.name=cluster.name,
               .cox.predictors=cox.predictors, .name.prefix=name.prefix,
               .na.methods=na.methods, .write.path=write.path,
               .impute.with = impute.with,
               .make.miss.if.contains=make.miss.if.contains
)



  