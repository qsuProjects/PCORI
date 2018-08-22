

######################### LOAD PACKAGES ########################

library(metafor)
library(mvtnorm)
library(ICC)
library(miscTools)
library(car)
library(plyr)
library(corpcor)
library(psych)
library(ggplot2)
library(SimTimeVar)


######################### TIME-VARYING ONLY; NO CATEGORICALS ########################

# cut down the parameters matrices provided with the package
static.names = c("height", "weight", "male")
pcor1 = pcor[ , !names(pcor) %in% static.names]
wcor1 = wcor[ , !names(wcor) %in% static.names]
params1 = params[ !params$name %in% static.names, ]
cat.params1 = cat.params[ !cat.params$parameter == "male", ]


data1 = make_one_dataset(n=10,
                        obs=10,
                        n.TBins=2,
                        pcor=pcor1,
                        wcor=wcor1, 
                        parameters=complete_parameters(params1, n=10) )$data



######################### TIME-VARYING AND STATIC; NO CATEGORICALS ########################

# cut down the parameters matrices provided with the package
pcor2 = pcor[ , names(pcor) != "cat.static" ]
wcor2 = wcor[ , names(wcor) != "cat.static"]
params2 = params[ params$type != "cat.static", ]
cat.params2 = cat.params[ !cat.params$parameter == "male", ]

data2 = make_one_dataset(n=10,
                        obs=10,
                        n.TBins=2,
                        pcor=pcor2,
                        wcor=wcor2, 
                        parameters=complete_parameters(params2, n=10) )$data



######################### ALL VARIABLES ########################

# this is now the same as the example in ?make_one_dataset
# increase sample size as well
# store the entire returned object this time for upcoming plots
sim = make_one_dataset(n=1000,
                       obs=10,
                       n.TBins=2,
                       pcor=pcor,
                       wcor=wcor, 
                       parameters=complete_parameters(params, n=10),
                       cat.parameters=cat.params)



######################### QUICK TOUR THROUGH THE SIMULATED DATA ########################

##### Example of Conditional Probabilty of Binary = 1 Given Ever-Use #####
# these should be similar:
d = sim$data
mean( d$drug1[ d$drug1_s == 1] ); parameters$across.mean[ parameters$name == "drug1" ]
mean( d$drug2[ d$drug2_s == 1] ); parameters$across.mean[ parameters$name == "drug2" ]


##### Example of Categorical Variable #####
# refit the models that generated race from cat.params
coef( glm(black ~ male, data=d, family=binomial(link="logit")) )[["male"]]
coef( glm(other ~ male, data=d, family=binomial(link="logit")) )[["male"]]
# pretty good


##### Example of Time-Function Variable #####
# take random sample of subjects
n = 20; keepers = sample( unique(d$id), size=n )
temp = d[ d$id %in% keepers, ]

# add a time variable for plotting purposes
temp$t = rep(1:10, n)

# plot growth curves (time trajectories) for blood pressure
ggplot( data=temp, aes(x=t, y=bp) ) + geom_point() + geom_line() +
  facet_wrap(~id) + theme_bw() + xlab("Time") + ylab("Blood pressure")


##### Example of Normal Variable Clustered within a Subject #####
temp = d[ d$id %in% 1:6, ]  # look at first 6 subjects

# boxplots of weight by subject
ggplot( data=temp, aes(x=id, y=weight, group=id) ) + geom_boxplot() +
  theme_bw() + xlab("Subject ID") + ylab("Weight")

# it is nicely normal
ggplot(data=d, aes(x=weight)) + geom_histogram() + theme_bw() + xlab("Weight")

# it is also highly correlated with height, as we specified
ggplot(data=d, aes(x=weight, y=height)) + geom_point() + theme_bw() + xlab("Weight") +
  ylab("Height")
# note how there are little horizontal lines (multiple weights associated with the same height)
# this is because height is specified as static, so one person has a static height but fluctuating weight

