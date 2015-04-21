

# loading setup values
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/Jolani/source")
source("subfunctions.r")


################################# MAKE FAKE DATA: SYSTEMATIC MISSINGNESS ONLY #################################

n = 8
obs = 20
set.seed(1); sex = sample(c(0,1), n*obs, replace=TRUE)
set.seed(1); norm = rnorm(n*obs, mean=50, sd = 8)

b0 = 0  # fixed intercept
set.seed(1); g0 = rnorm(n, mean=0, sd=10)  # subject random intercept
b1 = 20  # coef for norm
set.seed(1); g1 = rnorm(n, mean=0, sd=3)  # subject random slope
b2 = 8  # coef for sex
set.seed(1); g2 = rnorm(n, mean=0, sd=3)
set.seed(1); e = rnorm(n*obs, mean=0, sd=10)

y = b0 + rep(g0, each=obs) + b1*norm + rep(g1, each=obs)*norm + 
  b2*sex + rep(g2, each=obs)*sex + e

d.sys = data.frame(id=rep(1:n, each=obs), y, sex, norm)

# look at clustering by subject
boxplot(y ~ as.factor(id), data=d.sys, xlab="id", ylab="y")
boxplot(norm ~ as.factor(id), data=d.sys, xlab="id", ylab="y")
boxplot(sex ~ as.factor(id), data=d.sys, xlab="id", ylab="y")

# impose sporadic heavy MCAR missingness
# only impose missingness in binary variable, as Jolani did in paper
d.sys$sex[d.sys$id==1] = NA  # binary completely missing within subject 1
#d.sys$y[d.sys$id==4] = NA  # normal completely missing within subject 1




################################# SYSTEMATIC WITH JOLANI METHOD #################################

# first fit normal MICE to get predictor matrix
set.seed(1); ini = mice(d.sys, maxit=0)
pred = ini$predictorMatrix

# change method to 2l.norm.me
method = ini$method; method[method == "pmm"] = "2l.bin"

# treat trial as the cluster term
col = pred[ , "id"]; col[col==1] = -2; pred[, "id"] = col
#pred["y",] = c(-2, 0, 2, 2)
pred["sex",] = c(-2, 2, 0, 2)
pred["norm",] = c(-2, 2, 2, 0)

# impute with 2l.bin
set.seed(1)
imp.j = mice(d.sys, pred = pred, method = method)

head(imp.j$imp$sex)
# sex is imputed as binary

bwplot(imp.j, sex~.imp|id)

# APPEARS TO WORK? I ONLY USED A BINARY VARIABLE.




################################# MAKE FAKE DATA: SPORADIC MISSINGNESS ONLY #################################

n = 8
obs = 20
set.seed(1); sex = sample(c(0,1), n*obs, replace=TRUE)
set.seed(1); norm = rnorm(n*obs, mean=50, sd = 8)

b0 = 0
set.seed(1); g0 = rnorm(n, mean=0, sd=10)
b1 = 20
set.seed(1); g1 = rnorm(n, mean=0, sd=3)
b2 = 8
set.seed(1); g2 = rnorm(n, mean=0, sd=3)
set.seed(1); e = rnorm(n*obs, mean=0, sd=10)

y = b0 + rep(g0, each=obs) + b1*norm + rep(g0, each=obs)*norm + 
  b2*sex + rep(g2, each=obs)*sex + e

d.spor = data.frame(id=rep(1:n, each=obs), y, sex, norm)

# look at clustering by subject
boxplot(y ~ as.factor(id), data=d.spor, xlab="id", ylab="y")

# impose sporadic heavy MCAR missingness on all 3 variables
p.missing = 0.4
#for (c in 2:ncol(d.spor)) {  # don't put missingness in id column
for (c in 3) {  # only put missingness in sex variable
  set.seed(c); missing = rbinom(n=n*obs, size=1, prob=p.missing)
  d.spor[missing==1,c] = NA
}



################################# SPORADIC WITH 2L.NORM #################################

# first fit normal MICE to get predictor matrix
set.seed(1); ini = mice(d.spor, maxit=0)
pred = ini$predictorMatrix

# change method to 2l.norm
method = ini$method; method[method == "pmm"] = "2l.bin"

# treat trial as the cluster term
col = pred[ , "id"]; col[col==1] = -2; pred[, "id"] = col
#pred["y",] = c(-2, 0, 2, 2)
pred["sex",] = c(-2, 2, 0, 2)
#pred["norm",] = c(-2, 2, 2, 0)

# impute with 2l.norm.bin
set.seed(1); imp.rr = mice(d.spor, pred = pred, method = method)

head(imp.rr$imp$sex)
# appears to work

#bwplot(imp.vb, y~.imp|id)
bwplot(imp.rr, sex~.imp|id)


