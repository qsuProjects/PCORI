
########################### PLOT STITCHED RESULTS ###########################

setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/stitched")

l = read.csv("stitched.csv")
l= l[,-1]; dim(l)  # should be 2800 (4 lines * 7 methods * 100 datasets)

# how many datasets were run?
( n = length(unique(l$source.file)) )

# true betas
setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest")
( b = read.csv("true_betas.csv") )


library(reshape2)
library(ggplot2)
library(data.table)



########################### DESCRIPTIVE ###########################

# mean % missing in the datasets

temp = l[!duplicated(l$source.file),]
mean(temp$prop.rows.missing)

# mean X-Z correlation
mean(temp$XZ.corr)


########################### RESHAPE RESULTS ###########################

# names of variables that should go into wide form
varying.names = names(l)[!names(l) %in% c("stat", "source.file", "method", "prop.rows.missing", "XZ.corr") ]

w = reshape(l, timevar="stat",
            idvar=c("method", "source.file"),
            direction="wide")

# remove dumb columns
w = w[ ,!names(w) %in% c("prop.rows.missing.se", "prop.rows.missing.lo.95",
"prop.rows.missing.hi.95", "prop.rows.missing.lo 95", "prop.rows.missing.hi 95",
"XZ.corr.se", "XZ.corr.lo 95",
"XZ.corr.hi 95")]
names(w)[c(3,6:8)] = c("est", "se", "lo.95", "hi.95")

# merge in betas
w = merge(w, b)
w$bias = w$beta - w$est
w$std.bias = (w$bias/w$se) * 100


########################### MAKE DATASET OF AVERAGES ###########################

# compute coverage (0/1) for each sim
w$covers = (w$beta <= w$hi.95) & (w$beta >= w$lo.95)

# compute MSE for each sim
w$mse = (w$se)^2 + (w$bias)^2

# beta hat
m = data.table(w)
m[, est.mean := mean(est), by=list(method, var) ]
m[, est.lo.95 := quantile(est, .025), by=list(method, var) ]
m[, est.hi.95 := quantile(est, .975), by=list(method, var) ]

# SE
m[, se.mean := mean(se), by=list(method, var) ]
m[, se.lo.95 := quantile(se, .025), by=list(method, var) ]
m[, se.hi.95 := quantile(se, .975), by=list(method, var) ]

# bias
m[, bias.mean := mean(bias), by=list(method, var) ]
#m[, bias.lo.95 := quantile(bias, .025), by=list(method, var) ]  # NOT USING EMPIRICAL BIAS QUANTILES
#m[, bias.hi.95 := quantile(bias, .975), by=list(method, var) ]

# std bias
m[, std.bias.mean := mean(std.bias), by=list(method, var) ]
m[, std.bias.lo.95 := quantile(std.bias, .025), by=list(method, var) ]
m[, std.bias.hi.95 := quantile(std.bias, .975), by=list(method, var) ]

# MSE
m[, mse.mean := mean(mse), by=list(method, var) ]
m[, mse.lo.95 := quantile(mse.mean, .025), by=list(method, var) ]
m[, mse.hi.95 := quantile(mse.mean, .975), by=list(method, var) ]

# coverage
m[, cov.mean := mean(covers), by=list(method, var) ]
m = data.frame(m)

# make a dumb variable that identifies unique method-variable combos
# and subset down to just 1 row per method-variable combo
m$delete = paste(m$method, m$var)
m2 = m[ !duplicated(m$delete), ]; dim(m2)

# remove unneeded variables
toss = c("source.file", "est", "se", "lo.95", "hi.95", "delete", "std.bias", "covers", "bias")
m2 = m2[, !names(m2) %in% toss]


########################### NEW VARIABLES FOR PLOTTING PURPOSES ###########################

# variable to group MI and CC methods together
m2$model.type = ifelse( m2$method %in% c("complete.case", "complete.case.by.subj"), "B.CC", "C.MI" )
m2$model.type[m2$method == "full"] = "A.FM"

# use alphabetization to get nice ordering
library(car)
m2$method2 = recode(m2$method, " 'complete.case' = 'F.complete.case'; 
                    'complete.case.by.subj' = 'E.complete.case.by.subj';
                    'frailty' = 'C.frailty';
                    'log-t' = 'B.log-t';
                    'naive' = 'D.naive';
                    'full' = 'G.full';
                    'no.time' = 'A.no.time' ")


########################### PLOT RESULTS ###########################
point.size = 6
xlab="Variable"
pd = position_dodge(width=0.3)
lwd=1.2
error.width = .2
colors = c("darkcyan", "deepskyblue", "cornflowerblue", "blue4", "coral4", "firebrick2", "black", "green")
shapes = c(19,19,19,19,17,17,18,18)
  
labels1 = rev( c("FM", "CC-R", "CC-S", "MI-N", "MI-F", "MI-LT", "MI-NT") )
p = list()

ylab="Std Bias (95% empirical CI)"
p[[1]] = ggplot(m2, aes(x=method2, y=std.bias.mean, color=method2, shape=method2) ) +
  theme_bw() + xlab(xlab) + ylab(ylab) +
  coord_flip() +
  geom_hline(yintercept=0, lty=2, lwd=lwd) +
  geom_hline(yintercept=-40, lty=2, col="red") +
  geom_hline(yintercept=40, lty=2, col="red") +
  geom_point(size=point.size, position=pd) +
  geom_errorbar( aes(ymin=std.bias.lo.95, ymax=std.bias.hi.95), width=error.width, lwd=lwd, position=pd ) +
  scale_color_manual(values = colors, guide = guide_legend(title="Model", reverse=TRUE), labels=labels1) +
  scale_shape_manual(values = shapes, guide = guide_legend(title="Model", reverse=TRUE), labels=labels1) +
  xlab("") +
  scale_x_discrete(labels=labels1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #theme(axis.title = element_text(size=22),
   #     legend.text = element_text(size=18), legend.title = element_text(size=18) )


ylab="Bias (95% theoretical CI)"
p[[2]] = ggplot(m2, aes(x=method2, y=bias.mean, color=method2, shape=method2) ) +
  theme_bw() + xlab(xlab) + ylab(ylab) +
  coord_flip() +
  geom_hline(yintercept=0, lty=2, lwd=lwd) +
  geom_point(size=point.size, position=pd) +
  geom_errorbar( aes(ymin = bias.mean - qnorm(.025)*(se.mean/sqrt(n)),
                     ymax = bias.mean + qnorm(.025)*(se.mean/sqrt(n))),
                 width=error.width, lwd=lwd, position=pd ) +
  scale_color_manual(values = colors, guide = guide_legend(title="Model", reverse=TRUE), labels=labels1) +
  scale_shape_manual(values = shapes, guide = guide_legend(title="Model", reverse=TRUE), labels=labels1) +
  xlab("") +
  scale_x_discrete(labels=labels1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#theme(axis.title = element_text(size=22),
#     legend.text = element_text(size=18), legend.title = element_text(size=18) )


ylab="SE (95% empirical CI)"
p[[3]] = ggplot(m2, aes(x=method2, y=se.mean, color=method2, shape=method2) ) +
  theme_bw() + xlab(xlab) + ylab(ylab) +
  coord_flip() +
  geom_point(size=point.size, position=pd) +
  geom_errorbar( aes(ymin=se.lo.95, ymax=se.hi.95), width=error.width, lwd=lwd, position=pd ) +
  scale_color_manual(values = colors, guide = guide_legend(title="Model", reverse=TRUE), labels=labels1) +
  scale_shape_manual(values = shapes, guide = guide_legend(title="Model", reverse=TRUE), labels=labels1) +
  xlab("") +
  scale_x_discrete(labels=labels1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#theme(axis.title = element_text(size=22),
#     legend.text = element_text(size=18), legend.title = element_text(size=18) )

ylab="Coverage (95% theoretical CI)"
p[[4]] = ggplot(m2, aes(x=method2, y=cov.mean, color=method2, shape=method2) ) +
  theme_bw() + xlab(xlab) + ylab(ylab) +
  coord_flip() +
  geom_hline(yintercept=0.95, lty=2, lwd=lwd) +
  geom_point(size=point.size, position=pd) +
  geom_errorbar( aes(ymin = cov.mean - sqrt( (cov.mean*(1-cov.mean)) / n),
                     ymax = cov.mean + sqrt( (cov.mean*(1-cov.mean)) / n) ),
                 width=error.width, lwd=lwd, position=pd ) +
  scale_color_manual(values = colors, guide = guide_legend(title="Model", reverse=TRUE), labels=labels1) +
  scale_shape_manual(values = shapes, guide = guide_legend(title="Model", reverse=TRUE), labels=labels1) +
  xlab("") +
  scale_x_discrete(labels=labels1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#theme(axis.title = element_text(size=22),
#     legend.text = element_text(size=18), legend.title = element_text(size=18) )


ylab="MSE (95% empirical CI)"
p[[5]] = ggplot(m2, aes(x=method2, y=mse.mean, color=method2, shape=method2) ) +
  theme_bw() + xlab(xlab) + ylab(ylab) +
  coord_flip() +
  geom_point(size=point.size, position=pd) +
  geom_errorbar( aes(ymin=mse.lo.95, ymax=mse.hi.95), width=error.width, lwd=lwd, position=pd ) +
  scale_color_manual(values = colors, guide = guide_legend(title="Model", reverse=TRUE), labels=labels1) +
  scale_shape_manual(values = shapes, guide = guide_legend(title="Model", reverse=TRUE), labels=labels1) +
  xlab("") +
  scale_x_discrete(labels=labels1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#theme(axis.title = element_text(size=22),
#     legend.text = element_text(size=18), legend.title = element_text(size=18) )


setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/results")
library(gridExtra)
ggsave( paste(Sys.Date(), "simple_cov_results_version11.pdf", sep="_"),
        do.call("marrangeGrob", c(p,nrow=1,ncol=1) ), width=9, height=9 )

