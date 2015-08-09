
########################### PLOT STITCHED RESULTS ###########################

setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/stitched")

# reshape the wide-format ones


nf = read.csv("NA_frailty_stitched.csv")
nn = read.csv("NA_naive_stitched.csv")
lt = read.csv("NA_log-t_stitched.csv")


nf.w = reshape(nf, varying=c("est", "se", "lo.95", "hi.95"), 
               v.names="X",
               times = c("est", "se", "lo.95", "hi.95"),
               timevar="stat",
               idvar="source.file",
               direction="long"); dim(nf.w)

nn.w = reshape(nn, varying=c("est", "se", "lo.95", "hi.95"), 
               v.names="X",
               times = c("est", "se", "lo.95", "hi.95"),
               timevar="stat",
               idvar="source.file",
               direction="long"); dim(nf.w)

lt.w = reshape(lt, varying=c("est", "se", "lo.95", "hi.95"), 
               v.names="X",
               times = c("est", "se", "lo.95", "hi.95"),
               timevar="stat",
               idvar="source.file",
               direction="long"); dim(nf.w)

cc = read.csv("CC_stitched.csv"); cc = cc[,-1]
f = read.csv("full_stitched.csv"); f = f[,-1]

r = rbind(nf.w, nn.w, lt.w, cc, f)


# r = rbind( read.csv("NA_frailty_stitched.csv"), 
#            read.csv("NA_naive_stitched.csv"),
#            read.csv("NA_log-t_stitched.csv"),
#            
#            read.csv("CC_stitched.csv"),
#            read.csv("full_stitched.csv")
#            )
r = r[,-1]; dim(r)

# how many datasets were run?
( n = length(unique(r$source.file)) )

# true betas
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest")
b = read.csv("true_betas.csv")


library(reshape2)
library(ggplot2)
library(data.table)



########################### DESCRIPTIVE ###########################

# mean % missing in the datasets

temp = r[!duplicated(r$source.file),]
mean(temp$prop.missing)


########################### RESHAPE RESULTS ###########################

# names of variables that should go into long form
varying.names = names(r)[!names(r) %in% c("stat", "source.file", "method") ]

# long form
l = reshape(r, varying=varying.names, 
            v.names="value",
            timevar="var",
            times = varying.names,
            direction="long"); dim(l)

w = reshape(l, timevar="stat",
            idvar=c("method", "var", "source.file"),
            direction="wide")

# remove dumb columns
w = w[ ,!names(w) %in% c("id.est", "id.se", "id.lo 95", "id.hi 95")]
names(w)[4:7] = c("est", "se", "lo.95", "hi.95")

# merge in betas
w = merge(w, b)
w$bias = w$beta - w$est
w$std.bias = (w$bias/w$se) * 100


########################### MAKE DATASET OF AVERAGES ###########################

w$covers = (w$beta <= w$hi.95) & (w$beta >= w$lo.95)

# col for mean est, mean se, 
m = data.table(w)
m[, est.mean := mean(est), by=list(method, var) ]
m[, est.lo.95 := quantile(est, .025), by=list(method, var) ]
m[, est.hi.95 := quantile(est, .975), by=list(method, var) ]

m[, se.mean := mean(se), by=list(method, var) ]
m[, se.lo.95 := quantile(se, .025), by=list(method, var) ]
m[, se.hi.95 := quantile(se, .975), by=list(method, var) ]

m[, bias.mean := mean(bias), by=list(method, var) ]
m[, bias.lo.95 := quantile(bias, .025), by=list(method, var) ]
m[, bias.hi.95 := quantile(bias, .975), by=list(method, var) ]

m[, std.bias.mean := mean(std.bias), by=list(method, var) ]
m[, std.bias.lo.95 := quantile(std.bias, .025), by=list(method, var) ]
m[, std.bias.hi.95 := quantile(std.bias, .975), by=list(method, var) ]

m[, cov.mean := mean(covers), by=list(method, var) ]
m = data.frame(m)

# make a dumb variable that identifies unique method-variable combos
# and subset down to just 1 row per method-variable combo
m$delete = paste(m$method, m$var)
m2 = m[ !duplicated(m$delete), ]; dim(m2)

# remove unneeded variables
toss = c("source.file", "est", "se", "lo.95", "hi.95", "delete", "std.bias", "covers", "bias")
m2 = m2[, !names(m2) %in% toss]


########################### PLOT RESULTS ###########################

point.size = 4
xlab="Variable"
pd = position_dodge(width=0.3)
lwd=1.2
error.width = .2

p = list()

ylab="Std Bias (95% empirical CI)"
p[[1]] = ggplot(m2, aes(x=var, y=std.bias.mean, color=method) ) +
  theme_bw() + xlab(xlab) + ylab(ylab) +
  coord_flip() +
  geom_hline(yintercept=0, lty=2, lwd=lwd) +
  geom_hline(yintercept=-40, lty=2, col="red") +
  geom_hline(yintercept=40, lty=2, col="red") +
  geom_point(size=point.size, position=pd) +
  geom_errorbar( aes(ymin=std.bias.lo.95, ymax=std.bias.hi.95), width=error.width, lwd=lwd, position=pd )
  #theme(axis.title = element_text(size=22),
  #      legend.text = element_text(size=18), legend.title = element_text(size=18) )

ylab="Bias (95% empirical CI)"
p[[2]] =ggplot(m2, aes(x=var, y=bias.mean, color=method) ) + geom_point(size=point.size, position=pd) +
  theme_bw() + xlab(xlab) + ylab(ylab) +
  geom_errorbar( aes(ymin=bias.lo.95, ymax=bias.hi.95), width=error.width, lwd=lwd, position=pd ) +
  coord_flip() +
  geom_hline(yintercept=0, lty=2, lwd=lwd)

ylab="SE (95% empirical CI)"
p[[3]] =ggplot(m2, aes(x=var, y=se.mean, color=method) ) + geom_point(size=point.size, position=pd) +
  theme_bw() + xlab(xlab) + ylab(ylab) +
  geom_errorbar( aes(ymin=se.lo.95, ymax=se.hi.95), width=error.width, lwd=lwd, position=pd ) +
  coord_flip()

ylab="Coverage (95% theoretical CI)"
p[[4]] =ggplot(m2, aes(x=var, y=cov.mean, color=method) ) + geom_point(size=point.size, position=pd) +
  theme_bw() + xlab(xlab) + ylab(ylab) +
  geom_errorbar( aes(ymin = cov.mean - sqrt( (cov.mean*(1-cov.mean)) / n),
                     ymax = cov.mean + sqrt( (cov.mean*(1-cov.mean)) / n) ),
                 width=error.width, lwd=lwd, position=pd ) +
  coord_flip() +
  geom_hline(yintercept=0.95, lty=2, lwd=lwd)


setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/results")
library(gridExtra)
ggsave( paste(Sys.Date(), "cd4_and_teno_results.pdf", sep="_"),
        do.call("marrangeGrob", c(p,nrow=1,ncol=1) ), width=9, height=9 )
