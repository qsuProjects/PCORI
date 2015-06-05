
######################### SET SIMULATION PARAMETERS #########################

library(ggplot2)


# name prefix for all datasets
name_prefix = "ex1"

# n (number of subjects)
n.Subj = 1000

# obs (number of observations per subject)
obs = 7

# n.Reps (number of datasets that each worker should generate)
n.Reps = 1

# n.Drugs (number of drug variables)
n.Drugs = 2

setwd("/Users/mmathur/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV/vignettes")

# read in and complete parameters dataframe
parameters = complete_parameters( read.csv("ex1_parameters.csv"), n.Subj )

cat.parameters = read.csv("ex1_categorical_parameters.csv")


# within-subject correlation matrix
wcor = read.csv("ex1_wcor.csv", header=FALSE)[-1,-1]

# population correlation matrix
pcor.old = read.csv("ex1_pcor.csv", header=FALSE)[-1,-1]
pcor = read.csv("ex1_pcor.csv", header=TRUE)[,-1]

setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV")
source("jointly_generate_binary_normal_modified_v2.R", local=TRUE)
source("load_functions.R", local=TRUE)


######################### SIMULATE ########################

# simulate results
sim = repeat_sim(n=n.Subj, obs=obs, parameters=parameters, prop.target=NULL,
                       mean.target=NULL, n.Drugs=n.Drugs, 
                       pcor=pcor, wcor=wcor, n.Reps=n.Reps,
                       write.data=TRUE,
                       #name_prefix= paste( .name_prefix, WORKER.ID, sep="_" ),  # used with Sherlock
                       name_prefix=name_prefix,
                       cat.parameters=cat.parameters )

# extract dataset
d = sim$data
head(d)


######################### HOW'D WE DO? ########################

##### Example of Categorical Variable #####
# refit the models that generated race
summary( glm(black ~ male, data=d, family=binomial(link="logit")) )
summary( glm(other ~ male, data=d, family=binomial(link="logit")) )
# pretty good


##### Example of Time-Function Variable #####
# take random sample of subjects
n = 20; keepers = sample( unique(d$id), size=n )
temp = d[ d$id %in% keepers, ]

# add a time variable for plotting purposes
temp$t = rep(1:obs, n)

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


##### Correlation Matrix Across Subjects #####
# extract observed across-subjects correlation matrix
sim$corr

# how biased were the correlations?
sim$corr.bias


