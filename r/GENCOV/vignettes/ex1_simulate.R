
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
wcorin = read.csv("ex1_wcor.csv", header=FALSE)[-1,-1]
wcor = as.numeric(t(wcorin)[lower.tri(wcorin)]) #it gets read in by rows

# population correlation matrix
pcorin = read.csv("ex1_pcor.csv", header=FALSE)[-1,-1]
pcor = as.numeric(t(pcorin)[lower.tri(pcorin, diag=F)]) # need to transpose and read in the lower half to convert the a matrix into vector by row

setwd("~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV")
source("jointly_generate_binary_normal_modified_v2.R", local=TRUE)
source("load_functions.R", local=TRUE)
#source("init_variables.R", local=TRUE)  # this is PCORI-specific


######################### SIMULATE ########################

# simulate results
d = repeat_sim(n=n.Subj, obs=obs, parameters=parameters, prop.target=NULL,
                       mean.target=NULL, n.Drugs=n.Drugs, 
                       pcor=pcor, wcorin=wcorin, n.Reps=n.Reps,
                       write.data=TRUE,
                       #name_prefix= paste( .name_prefix, WORKER.ID, sep="_" ),  # used with Sherlock
                       name_prefix=name_prefix,
                       cat.parameters=cat.parameters )


######################### HOW'D WE DO? ########################

##### Refit Race Model #####
summary( glm(black ~ male, data=d, family=binomial(link="logit")) )
summary( glm(other ~ male, data=d, family=binomial(link="logit")) )
# pretty good


##### Example of Time-Function Variable #####
# look at time-trajectories of blood pressure
# take random sample of subjects
n = 20; keepers = sample( unique(d$id), size=n )
temp = d[ d$id %in% keepers, ]

# add a time variable for plotting purposes
temp$t = rep(1:obs, n.Subj)

# plot growth curves
ggplot( data=temp, aes(x=t, y=bp) ) + geom_point() + geom_line() +
  facet_wrap(~id) + theme_bw() + xlab("Time") + ylab("Blood pressure")


##### Example of Normal Variable Clustered within a Subject #####
temp = d[ d$id %in% 1:6, ]  # look at first 6 subjects
ggplot( data=temp, aes(x=id, y=weight, group=id) ) + geom_boxplot() +
  theme_bw() + xlab("Subject ID") + ylab("Weight")


##### Look at Correlation Matrix Across Subjects #####
# put in Kris' code

