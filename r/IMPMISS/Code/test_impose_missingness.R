
########################## READ IN DATA ##########################

# t: time at end of interval
# t0: time at beginning of interval
# d: whether event occurred at that interval

# load impose missingness code
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/IMPMISS/Code")
source("impose_missingness_functions.R")

setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/IMPMISS/Data from KK")
data = read.csv("SURV_2015-02-01_job_10_dataset_1.csv")

setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/IMPMISS")
aux.matrix = read.csv("aux_var_parameters_matrix.csv")
miss.matrix = read.csv("missing_var_parameters_matrix.csv")

#### TEST ONLY!!
#miss.matrix = read.csv("missing_var_parameters_matrix_test.csv")
miss.matrix$beta[miss.matrix$parameter=="subj.rand.int"] = 5; miss.matrix  # TEST ONLY!
#miss.matrix$beta[miss.matrix$parameter=="intercept"] = 0; miss.matrix

rand.int.SD = as.numeric( read.table("rand_intercepts_sd.txt") )

# load required libraries
library(data.table)


########################## TEST PERFORMANCE ##########################

#data=data; outcome.name="d"; start.time.name="t0";
#stop.time.name="t"; id.var.name="id"; aux.matrix=aux.matrix;
#miss.matrix=miss.matrix; make.relateds.missing=TRUE

d.miss = impose_missingness( data=data, outcome.name="d", start.time.name="t0",
                             stop.time.name="t", id.var.name="id", aux.matrix=aux.matrix,
                             miss.matrix=miss.matrix, make.relateds.missing=TRUE,
                             rand.int.SD = rand.int.SD )

# refit models that created missingness to see how we're doing
# looks good!
summary( glm( is.na(bmi) ~ aux.etiol + time2event + d + subj.rand.int,
              data=d.miss, family=binomial(link="logit") ) ) 

summary( glm( is.na(ldl) ~ t0 + aux.etiol + subj.rand.int,
              data=d.miss, family=binomial(link="logit") ) ) 

summary( glm( is.na(log_vln) ~ t0 + aux.etiol,
              data=d.miss, family=binomial(link="logit") ) ) 



##### Try Varying Coefficients for Random Intercepts #####

coefs = c(0, 2, 5, 10, 20)
rs = as.data.frame( matrix(nrow=length(coefs), ncol=3) )
names(rs) = c("coef", "prop.obs.drop", "prop.s.drop")
rs$coef = coefs

for ( i in 1:length(rs$coef) ) {
  
  # change missingness matrix
  miss.temp = miss.matrix
  miss.matrix$beta[miss.matrix$parameter=="subj.rand.int"] = rs$coef[i]
  
  # impose missingness using updated missingness matrix
  d.temp = impose_missingness( data=data, outcome.name="d", start.time.name="t0",
                               stop.time.name="t", id.var.name="id", aux.matrix=aux.matrix,
                               miss.matrix=miss.matrix, make.relateds.missing=TRUE,
                               rand.int.SD = rand.int.SD )
  
  rs$perc.obs.drop[i] = proportion_dropped(data=d.temp, miss.matrix=miss.matrix) * 100
  rs$perc.s.drop[i] = prop_subj_dropped(data=d.temp, id.var.name="id", miss.matrix=miss.matrix) * 100
}



