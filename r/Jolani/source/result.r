#summarizing simulation results
# first the results of a scenario must be uploaded. The corresponding array is called 'final'.
# second the setup value is uploaded in 'setup'.
load("scenario1.RData")
load("setup.RData")

#parameter estimates (fixed effects)
est <- apply(final$result[,c(1:2,5:7),,],1:3,mean)
fix <- est[,1,]
se.fix <- est[,2,]

# relative bias
bias <- fix - setup$fixef
rb <- round(((fix - setup$fixef)/setup$fixef),3)
# rmse
var <- apply(final$result[,c(1:2,5:7),,],1:3,var)[,1,]
rmse1 <- sqrt((fix - setup$fixef)^2 + var)
#coverage
tempcr <- ifelse(setup$fixef > final$result[,1,,] - 1.96*final$result[,2,,] & 
                   setup$fixef <= final$result[,1,,] + 1.96*final$result[,2,,], 1, 0)
cr <- apply(tempcr, 1:2, mean)

# Fixed effect table 
round(t(rbind(fix[1,], rb[1,], rmse[1,], cr[1,],
              fix[2,], rb[2,], rmse[2,], cr[2,],
              fix[3,], rb[3,], rmse[3,], cr[3,])),3)


# parameter estimates (random effects)
est.ran <- apply(sqrt(abs(final$result[,5:7,,])),1:3,mean)
random <- (cbind(diag(est.ran[,,1]), diag(est.ran[,,2]), diag(est.ran[,,3]), 
                 diag(est.ran[,,4]), diag(est.ran[,,5])))
colnames(random) <- dimnames(est.ran)[[3]]
rownames(random) <- dimnames(est.ran)[[1]]

# random effect table
random <- t(random)

# Graphs
require(lattice)
load("scenario2.RData")
finh <- final$result
load("scenario4.RData")
finl <- final$result

scenario <- gl(2,(500*5*2), labels = c("2", "4"))
para <- rep(gl(2, (500*5), labels = c("tau0","tau1")),2)
meth <- rep(rep(gl(5, 500, labels = c("REF","CCA","TMI","SMI","MLMI")),2),2)
outc <- sqrt(c(c(t(finl[1,5,,])), c(t(finl[2,6,,])),  
               c(t(finh[1,5,,])), c(t(finh[2,6,,]))))
bwplot(outc~meth|para + scenario,
       ylab = c("Scenario 4", "Scenario 2"),
       xlab.top = c(expression(tau[0]), expression(tau[1])),
       col.line = "darkgrey", col.symbol = "black", strip = FALSE,
       between = list(x = .5, y = .5),
      #ylim = list(c(0,1.45),c(0,1.45),c(0,2.2),c(0,2.2)),
       scales = list(#y = list(relation = "free"), 
                     x = list(alternating = 1, rot = 45)))

