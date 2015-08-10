
################################# MAKE SIMPLE COVARIATES #################################

# generate a random intercept for each subject ~ N(0, something)
# X ~ rand.int + N(0, something)
# Z ~ X + N(0, something)

# load code for creating random effects
setwd("/share/PI/manishad/naEst/r")
source("impose_missingness_functions_general.R")

# load parameters matrices
setwd("/share/PI/manishad/naEst")
( mat.X = read.csv("aux_var_parameters_matrix_main_X.csv") )
( mat.Z = read.csv("aux_var_parameters_matrix_Z.csv") )

# function to generate one simple dataset with X and Z
make_one_dataset = function(.n, .obs) {
  d1 = data.frame( id=rep(1:.n, each=.obs) )  # initialize dataframe
  
  # gnerate subject random intercepts
  d1 = create_random_effects(data=d1, id.var.name="id",
                             rand.int.SD=1)
  
  # generate clustered X as a function of the random intercept
  d2 = make_aux_vars(data=d1, aux.matrix=mat.X)
  
  # generate auxiliary Z as a function of X
  d2 = make_aux_vars(data=d2, aux.matrix=mat.Z)
  return(d2)
}

# set simulation parameters
n = 1000
obs = 200
reps = 100
name.prefix = paste(Sys.Date(), "simple_covs_dataset", sep="_")
write.path = "/share/PI/manishad/naEst/output/datasets"

# simulate multiple times
for (i in 1:reps) {
  name = paste(name.prefix, i, sep="_")
  
  d = make_one_dataset(n, obs)
  write.csv( d, paste(write.path, name, sep="/") )
}


################################# LOCAL TEST #################################
# 
# # load code for creating random effects
# setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/IMPMISS/Code")
# source("impose_missingness_functions_general.R")
# 
# # load other stuff
# setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/for-sherlock")
# ( mat.X = read.csv("aux_var_parameters_matrix_main_X.csv") )
# ( mat.Z = read.csv("aux_var_parameters_matrix_Z.csv") )
# write.path = "~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/local-test"
# name = "test_data.csv"
# 
# n = 1000
# obs = 200
# 
# # generate 1 dataset
# d = make_one_dataset(n, obs)
# write.csv( d, paste(write.path, name, sep="/") )
# 
# ##### Check That Data Generation Worked #####
# 
# # fit data generation model to make sure
# library(lme4)
# rs = lmer(X ~ (1|id), data=d2); summary(rs)
# rs = lm(Z ~ X, data=d2); summary(rs)
# # looks pretty good
# 
# #rs = glmer(y ~ x + (1|id), data=d3, family="binomial"(link="logit") ); summary(rs)
# # works if no random slope; does not work if includes random slopes
# 
# # look at clustering by subject
# boxplot(X ~ as.factor(id), data=d2, xlab="id", ylab="y")
# 
# # how correlated are X and Z?
# plot(X~Z, d2)
# cor(d2$X, d2$Z)






