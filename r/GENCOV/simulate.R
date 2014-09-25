#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Project: PCORI Missing Data
# 
# -This program creates covariates for simulated data based on existing VA data 
#    (summary stats provided by Vilija at VA, see Q:\Datasets\PCORI\data simulation\real data estimates)
#
# -Using jointly_generate_binary_normal (binNor package) function by Demirtas, simulate clinical characterisitics and drug exposure for
#  a large number of patients over time.
#
# -This file runs the simulation and looks at performance.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~s


######################### SIMULATION TEMPLATE #########################


n=
obs=
n.Drugs=15
n.Reps=

t1 = Sys.time()
results = repeat_sim(n, obs, pps=pps, pps.target=pps, pmeans=pmeans, pmeans.target=pmeans, n.Drugs, pcor, wcorin, pvars, n.Reps, varnames, write.data=FALSE)
t2 = Sys.time()

# time to complete simulation
t2 - t1

mean_results = mean_performance(results)



######################### PLOT PERFORMANCE #########################

# extract variable names
normal.names <- varnames[ (n.OtherBins + n.Drugs + 2) : (n.OtherBins + n.Drugs + 1 + n.OtherNorms) ]
other.bin.names <- varnames[ 2 : (1 + n.OtherBins) ]
drug.names <- varnames[ (n.OtherBins + 2) : (n.OtherBins + 1 + n.Drugs) ]


# make performance plots
par(mfrow=c(1,3))

plotBias(results$nor.vars$abs.bias, normal.names, ylim=c(-15,15), yaxp=c(-15,15,10), ylab="Absolute bias", main="Normal Variables")
plotBias(results$nor.vars$std.bias, normal.names, ylim = c(-3,3), yaxp=c(-3,3,6), ylab="Bias/SE", main = paste("n=", n, ", obs=", obs, ", reps=", n.Reps) )
plotCoverage(results$nor.vars$coverage, n.Reps=50, normal.names, main=Sys.Date() )

plotBias(results$drug.evers$abs.bias, drug.names, ylim=c(-1,1), yaxp=c(-1,1,4), ylab="Absolute bias", main="Drug Ever-Use")
plotBias(results$drug.evers$std.bias, drug.names, ylim=c(-6,6), yaxp=c(-6,6,6), ylab="Bias/SE", main=paste("n=", n, ", obs=", obs, ", reps=", n.Reps) )
plotCoverage(results$drug.evers$coverage, n.Reps=50, drug.names, main=Sys.Date() )

plotBias(results$prop.drug.time$abs.bias, drug.names, ylim=c(-.1,.1), yaxp=c(-.1,.1,8), ylab="Absolute bias", main="Proportion Drug Time")
plotBias(results$prop.drug.time$std.bias, drug.names, ylim=c(-10,10), yaxp=c(-10,10,10), ylab="Bias/SE", main=paste("n=", n, ", obs=", obs, ", reps=", n.Reps) )
plotCoverage(results$prop.drug.time$coverage, n.Reps=50, drug.names, main=Sys.Date() )

plotBias(results$prop.male$abs.bias, other.bin.names, ylab="Absolute bias", main="Male")
plotBias(results$prop.male$std.bias, other.bin.names, ylim=c(-3,3), ylab="Bias/SE", main="Male")
plotCoverage(results$prop.male$coverage, n.Reps=50, other.bin.names)


par(mfrow=c(1,1))

