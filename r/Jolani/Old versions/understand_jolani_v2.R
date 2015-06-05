

# load Jolani code
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/Jolani/source")
source("subfunctions-mm_edits.r")

# load IMPISS code
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/IMPMISS/Code")
source("impose_missingness_functions.R")

# load parameters for IMPMISS
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/Jolani")
aux.matrix = read.csv("aux_var_parameters_matrix.csv")
#aux.matrix = read.csv("aux_var_parameters_matrix_slopes.csv")  # version with random slopes by subject
rand.int.SD = 8

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
                               rand.int.SD=rand.int.SD)
# generate subject random slopes 
d2.2 = create_random_effects(data=d1, id.var.name="id",
                             rand.int.SD=1)
d2$subj.rand.slope = d2.2$subj.rand.int

d3 = make_aux_vars(data=d2, aux.matrix=aux.matrix)

# remove random intercept
d3 = d3[ , names(d3) != "subj.rand.int" ]

# see distribution of y's
table(d3$y)


# fit data generation to make sure
rs = glmer(y ~ x + (1|id), data=d3, family="binomial"(link="logit") ); summary(rs)
# yes, works :)

# look at clustering by subject
boxplot(y ~ as.factor(id), data=d3, xlab="id", ylab="y")

# impose systematic missingness only for binary variable
d.sys = d3
p.miss = 0.4
losers = sample( d.sys$id, replace=FALSE, size=p.miss*n ); length(losers)  # choose which subjects will be systematically missing
d.sys$y[ d.sys$id %in% losers ] = NA


################################# MAKE FAKE DATA: SPORADIC MISSINGNESS ONLY #################################

d.spor = d3

# impose sporadic heavy MCAR missingness on missing variables
p.missing = 0.4
#for (c in 2:ncol(d.spor)) {  # don't put missingness in id column
for (c in 4) {  # only put missingness in bin1
  set.seed(c); missing = rbinom(n=n*obs, size=1, prob=p.missing)
  d.spor[missing==1,c] = NA
}


# check that there's no systematic missingness by chance
dt = data.table(d.spor)
dt[, num.missing := sum(is.na(y)), by=id ]
max(dt$num.missing) == obs  # is anyone missing all obs?


################################# IMPUTATION FUNCTION #################################

multilevel_impute = function(.data, .method.name) {
  
  # first fit normal MICE to get predictor matrix
  set.seed(1); ini = mice(.data, maxit=0)
  pred = ini$predictorMatrix
  
  # change method
  method = ini$method; method[method == "pmm"] = .method.name
  
  # treat id as the cluster term
  col = pred[ , "id"]; col[col==1] = -2; pred[, "id"] = col
  pred["y",] = c(-2, 0, 1, 0)
  
  if (.method.name == "2l.norm") pred["y",] = c(-2, 2, 1, 0)  # 2l.norm wants a constant term specified as random effect
  
  # impute with given method
  set.seed(1)
  print(pred)
  return( mice(.data, pred = pred, method = method) )
}



################################# IMPUTE SPOR & SYST DATA WITH 3 DIFFERENT METHODS #################################

##### Sporadic Data #####

spor.rr = multilevel_impute(.data=d.spor, .method.name="2l.norm.me")  # Resche-Rigon
# not intended case b/c sporadic
# runs without error

spor.vb = multilevel_impute(.data=d.spor, .method.name="2l.norm")  # Van-Buuren
# not intended case b/c binaries?
# runs without error

spor.j = multilevel_impute(.data=d.spor, .method.name="2l.bin")  # Jolani
# not intended case b/c sporadic
# runs, but complains that "glmer cannot be run"
# BOOKMARK: Was running through this with debug(mice.impute.2l.bin)
# Although the glmer fits apparently with no problem, gets hung up later on eigenvalue problem



##### Systematic Data #####

syst.rr = multilevel_impute(.data=d.sys, .method.name="2l.norm.me")  # Resche-Rigon
# not intended case b/c binaries
# runs without error

syst.vb = multilevel_impute(.data=d.sys, .method.name="2l.norm")  # Van-Buuren
# not intended case b/c systematic
# breaks; error about Cholesky decomp??

syst.j = multilevel_impute(.data=d.sys, .method.name="2l.bin")  # Jolani
# INTENDED CASE
# runs, but complains that "glmer cannot be run"
# SHOULD WORK??


################################# DIAGNOSTICS #################################

diagnose = function(imput) {
  # diagnostics
 print( head(imput$imp$y) )
  bwplot(imput, y~.imp|id)
  
  # proportion imputed as a 1 for each subject
 # print( rowMeans(imput$imp$y) )
}

diagnose(spor.rr)
diagnose(spor.vb)
diagnose(spor.j)

diagnose(syst.rr)
diagnose(syst.j)



################################# REFIT MODEL #################################

##### Sporadic Missingness #####

rs.spor.rr = pool( with( data = spor.rr, exp = glmer( y ~ x + (1|id),
                                                 family=binomial(link="logit") ) ) )
# error about y being outside [0,1]  
  
rs.spor.vb = pool( with( data = spor.vb, exp = glmer( y ~ x + (1|id),
                                                    family=binomial(link="logit") ) ) )
# error about y being outside [0,1]    

rs.spor.j = pool( with( data = spor.j, exp = glmer( y ~ x + (1|id),
                                                    family=binomial(link="logit") ) ) )



##### Systematic Missingness #####

rs.syst.rr = pool( with( data = syst.rr, exp = glmer( y ~ x + (1|id),
                                                      family=binomial(link="logit") ) ) )
# error about y being outside [0,1]  
# happens also with just a single completed dataset

rs.syst.j = pool( with( data = syst.j, exp = glmer( y ~ x + (1|id), family=binomial(link="logit") ) ) )
# see if it works with 1 completed dataset
fake = complete(syst.j)
rs.syst.j = with( data = fake, exp = glmer( y ~ x + (1|id), family=binomial(link="logit") ) ); summary(rs.syst.j)
# completely wrong

