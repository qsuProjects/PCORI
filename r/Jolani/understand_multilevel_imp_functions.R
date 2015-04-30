
# This file tests 3 multilevel imputation approaches (Van Buuren's 2l.norm,
#  Resche-Rigon's 2l.norm.me, and Jolani's 2l.bin) on binary data with either
#  sporadic-only or systematic-only missing data. It also explains how to 
#  use each function. 
###############################################################################

# load Jolani code
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/Jolani/source")
source("subfunctions.r")

# load IMPISS code
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/IMPMISS/Code")
source("impose_missingness_functions.R")

# load parameters for IMPMISS
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/Jolani")
#aux.matrix = read.csv("aux_var_parameters_matrix.csv")
aux.matrix = read.csv("aux_var_parameters_matrix_slopes.csv")  # version with random slopes by subject
#aux.matrix$beta[aux.matrix$parameter=="slope.by.x"]=1

#rand.int.SD = 8

# load imputation functions
library(mice)
source("~/Dropbox/QSU/Mathur/PCORI/2015-01-28 problems with 2l.norm/mice.impute.2l.norm.me-MM_commented.R")

library(lme4)
library(data.table)


################################# MAKE FAKE DATA: SYSTEMATIC MISSINGNESS ONLY #################################

n = 100
obs = 20
d1 = data.frame( id=rep(1:n, each=obs) )  # initialize dataframe
d1$const = 1  # this is required by 2l.norm

# generate X ~ N(50, 10) unclustered by subject
d1$x = rnorm(n*obs, 50, 10)

# gnerate subject random intercepts
d2 = create_random_effects(data=d1, id.var.name="id",
                               rand.int.SD=0.4)
# generate subject random slopes 
d2.2 = create_random_effects(data=d1, id.var.name="id",
                             rand.int.SD=0.1)
d2$subj.rand.slope = d2.2$subj.rand.int

# put in the interaction variable
d2$slope.by.x = d2$x * d2$subj.rand.slope

d3 = make_aux_vars(data=d2, aux.matrix=aux.matrix)


##### Check That Data Generation Worked #####

# see distribution of y's
table(d3$y)
boxplot(d3$x ~ d3$y)

# fit data generation model to make sure
rs = glmer(y ~ x + (x|id), data=d3, family="binomial"(link="logit") ); summary(rs)
# looks pretty good

#rs = glmer(y ~ x + (1|id), data=d3, family="binomial"(link="logit") ); summary(rs)
# works if no random slope; does not work if includes random slopes

# look at clustering by subject
boxplot(y ~ as.factor(id), data=d3, xlab="id", ylab="y")


###### Prep for Imputation #####

# remove random subject variables prior to imputation
d3 = d3[ , !names(d3) %in% c("subj.rand.int", "subj.rand.slope", "slope.by.x") ]

# impose systematic missingness only for binary variable
d.sys = d3
p.miss = 0.4
losers = sample( d.sys$id, replace=FALSE, size=p.miss*n ); length(losers)  # choose which subjects will be systematically missing
d.sys$y[ d.sys$id %in% losers ] = NA




################################# MAKE FAKE DATA: SPORADIC MISSINGNESS ONLY #################################

d.spor = d3

# impose sporadic heavy MCAR missingness on missing variables
p.missing = 0.4
#for (c in 2:ncol(d.spor)) {  # don't put missingness in y column
for (c in 4) {  # only put missingness in bin1
  set.seed(c); missing = rbinom(n=n*obs, size=1, prob=p.missing)
  d.spor[missing==1,c] = NA
}


# check that there's no systematic missingness by chance
dt = data.table(d.spor)
dt[, num.missing := sum(is.na(y)), by=id ]
max(dt$num.missing) == obs  # is anyone missing all obs?


################################# IMPUTATION FUNCTION #################################

multilevel_impute = function(.data, .method.name, .pred) {
  
  # first fit normal MICE to get predictor matrix
  set.seed(1); ini = mice(.data, maxit=0)
  pred = ini$predictorMatrix
  
  # change method
  method = ini$method; method[method == "pmm"] = .method.name
  
  # treat id as the cluster term
  col = pred[ , "id"]; col[col==1] = -2; pred[, "id"] = col
  pred["y",] = c(-2, 0, 2, 0)
  
  if (.method.name == "2l.norm") pred["y",] = c(-2, 2, 2, 0)  # 2l.norm wants a constant term specified as random effect
  
  # impute with given method
  set.seed(1)
  print(pred)
  return( mice(.data, pred = pred, method = method) )
}



################################# HOW TO MAKE PREDICTOR MATRICES MANUALLY #################################

# first, do a dry run just to grab the predictor matrix
#ini = mice(d.spor, maxit=0)
#ini$pred  # this is the predictor matrix

# the pred matrix works the same way for all 3 multilevel methods:
# the i,jth entry says whether variable j is used to impute variable i
# 2 = random slope by cluster (and fixed effect)
# 1 = fixed effect
# -2 = the cluster variable

# make predictor matrix for 2l.norm
#pred.vb = ini$pred
#pred.vb[pred.vb == 1] = 2  # first change all 1s to 2s
# treat id as the cluster term
#col = pred.vb[ , "id"]; col[col==2] = -2; pred.vb[, "id"] = col
# annoying required constant term; see VB's JStatSoft paper
#col = pred.vb[ , "const"]; col = 2; pred.vb[, "const"] = col
# don't worry about the extra entries for vars that don't get imputed
# these will just get set to 0

# make predictor matrix for the other two methods
#pred2 = pred.vb
#pred2[,"const"] = 0
# same as for VB except don't use const


################################# IMPUTE SPORADIC DATA WITH 3 DIFFERENT METHODS #################################

spor.rr = multilevel_impute(.data=d.spor, .method.name="2l.norm.me")  # Resche-Rigon
# not intended case b/c sporadic
# runs without error

spor.vb = multilevel_impute(.data=d.spor, .method.name="2l.norm")  # Van-Buuren
# not intended case b/c binaries?
# runs without error

spor.j = multilevel_impute(.data=d.spor, .method.name="2l.bin")  # Jolani
# not intended case b/c sporadic
# runs without error


################################# IMPUTE SYSTEMATIC DATA WITH 3 DIFFERENT METHODS #################################

syst.rr = multilevel_impute(.data=d.sys, .method.name="2l.norm.me")  # Resche-Rigon
# not intended case b/c binaries
# runs without error

syst.vb = multilevel_impute(.data=d.sys, .method.name="2l.norm")  # Van-Buuren
# not intended case b/c systematic
# breaks; error about Cholesky decomp

syst.j = multilevel_impute(.data=d.sys, .method.name="2l.bin")  # Jolani
# INTENDED CASE
# runs without error


################################# DIAGNOSTICS #################################

diagnose = function(imput) {
  # diagnostics
 print( head(imput$imp$y) )
  bwplot(imput, y~.imp|id)
 
 # regression using imputed datasets
 # result is a mira object
 #rs = pool( with( data = spor.rr, exp = glmer( y ~ x + (1|id), family=binomial(link="logit") ) ) )
 #print(rs)
 # can't do this because y aren't imputed as binary :(
 
 # pool the results
 # result is a mipo object
 #mipo = pool(mira)
 #summary(mipo)
  # proportion imputed as a 1 for each subject
 # print( rowMeans(imput$imp$y) )
}

diagnose(spor.rr)  # yes, looks clustered
diagnose(spor.vb)
diagnose(spor.j)

diagnose(syst.rr)
diagnose(syst.j)



################################# REFIT MODEL #################################

##### Sporadic Missingness #####

rs.spor.j = pool( with( data = spor.j, exp = glmer( y ~ x + (x|id),
                                                    family=binomial(link="logit") ) ) )

fake1 = complete(spor.j)
rs.spor = with( data = fake1, exp = glmer( y ~ x + (x|id), family=binomial(link="logit") ) ); summary(rs.syst.j)






##### Systematic Missingness #####


rs.syst.j = pool( with( data = syst.j, exp = glmer( y ~ x + (x|id), family=binomial(link="logit") ) ) )
# see if it works with 1 completed dataset

fake2 = complete(syst.j)
rs.syst.j = with( data = fake2, exp = glmer( y ~ x + (x|id), family=binomial(link="logit") ) ); summary(rs.syst.j)


