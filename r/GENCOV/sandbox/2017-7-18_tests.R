


######################### LOAD FUNCTIONS #########################



######################### SET SIMULATION PARAMETERS #########################

# name prefix for all datasets
name_prefix = "ex1"

# n (number of subjects)
n.Subj = 10

# obs (number of observations per subject)
obs = 7

# n.Reps (number of datasets to generate)
n.Reps = 1

# n.Drugs (number of drug variables)
n.Drugs = 2


# use Vignette parameters
setwd("~/Dropbox/QSU/Mathur/PCORI/Git/PCORI/r/GENCOV/vignettes")
parameters = complete_parameters( read.csv("ex1_parameters.csv"), n.Subj )
cat.parameters = read.csv("ex1_categorical_parameters.csv")


# within-subject correlation matrix
wcor = read.csv("ex1_wcor.csv", header=FALSE)[-1,-1]

# population correlation matrix
pcor = read.csv("ex1_pcor.csv", header=TRUE)[,-1]



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