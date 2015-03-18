
######################### RUN IN PARALLEL #########################

# load command line arguments
args = commandArgs(trailingOnly = TRUE)
print(args)

# the number of clusters to use
n_clusters = as.numeric(args[1])

# name prefix for datasets
name_prefix = args[2]

# n (number of subjects)
( n.Subj  = as.numeric(args[3]) )

# obs (number of observations per subject)
obs = as.numeric(args[4])

# n.Reps (number of datasets that each worker should generate)
n.Reps = as.numeric(args[5])

# n.Drugs (number of drug variables)
n.Drugs = as.numeric(args[6])



# load required libraries
library(plyr)
library(doSNOW)

#set up parallelization backend
cl = makeCluster(n_clusters)
registerDoSNOW(cl)

#these should be equal
n_clusters
getDoParWorkers()

#give each worker an ID 
#this is optional, but occasionally is useful 
#in situations where you need workers to write data
#every worker writing to the same file will result in 
#uncaught errors and partially overwritten data
#worker ids enable each worker to write to its own set of 
#files
clusterApply(cl, seq(along=cl), function(id) WORKER.ID <<- paste0("worker_", id))


#example:
#each worker here runs a t test on simulated data and
#writes to a worker-specific file

time = system.time({

# EDITED TO PASS THE SBATCH ARGUMENTS ALONG WITH .ITEM TO FUNCTION
# AFTER DEBUGGING, REMOVE THE INFORM=T ARGUMENT! SLOWS THINGS DOWN!
l_ply( c( 1:getDoParWorkers() ), .parallel=T, .inform=T, function(.item, .n.Subj, .obs, .n.Reps, .n.Drugs) {  
  #l_ply iterates through each element of its first argument, here c(1:1000), and passes each element to
  #function() as .item    
  #instead of c(1:1000), you could pass l_ply() a list where each element is a set of subject specific means
  #then each worker/iteration would expand those means into a full set of data for that subject
  #avoid workers writing to the same file
  #if each worker/iteration generates data for one subject, you will also need to write code that stitches
  #all subject data together into one data file.
  
  #everything in here will be performed by workers in parallel    
  #packages and source functions used by workers should be loaded within the this block
  
  # FOR SOME REASON THIS NEEDS TO BE SET GLOBALLY - WHY?
	n=2000

  setwd("/share/PI/manishad/genCov")
	source("jointly_generate_binary_normal_modified_v2.R", local=TRUE)
  source("load_functions.R", local=TRUE)
  source("init_variables.R", local=TRUE)

	# DELETE THIS - testing only
	write.csv(.n.Subj, "n.Subj.csv")
  
  # simulate results
  setwd("/share/PI/manishad/genCov/datasets")
  results = repeat_sim(n=.n.Subj, obs=.obs, parameters=parameters, prop.target=NULL, mean.target=NULL, n.Drugs=.n.Drugs, 
                       pcor=pcor, wcorin=wcorin, n.Reps=.n.Reps, race.names=race.names, write.data=TRUE, WORKER.ID) 

}, n.Subj, obs, n.Reps, n.Drugs)

})

# write a file about the simulation parameters
write(
x = paste( "Simulation ", name_prefix, " ended on ", Sys.Date(),
           "\nThere were ", getDoParWorkers(), " workers",
          "\nDatasets were generated with n=", n.Subj, ", obs=", obs, ", n.Reps=",
          n.Reps, ", n.Drugs=", n.Drugs,
          "\nThe entire process took ", time[3]/(60*60), " hours",
          sep="" )
, file="simulation_info.txt"
)


