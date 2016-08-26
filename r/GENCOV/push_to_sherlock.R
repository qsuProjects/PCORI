# this pushes the files from local GitHub repo to Sherlock

# 
cd /share/PI/manishad/genCov/sbatch_files
rm rm*
sbatch -p manishad 8.sbatch -p manishad
squeue -p manishad

# push all the individual files
system( "scp -r ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV/* mmathur@sherlock:/share/PI/manishad/genCov" )

# push just R files
scp -r ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV/*.R mmathur@sherlock:/share/PI/manishad/genCov

# push the sbatch folder
#system( "scp ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV/sbatch_files mmathur@sherlock:/share/PI/manishad/genCov" )

# see the sbatch files
# cd /share/PI/manishad/genCov/sbatch_files

# run one of them
# sbatch -p manishad /share/PI/manishad/genCov/sbatch_files/6.sbatch


# see the datasets
# cd /share/PI/manishad/genCov/datasets

# move one dataset to local machine
system( "scp mmathur@sherlock:/share/PI/manishad/genCov/datasets/2016-08-25_job_1_worker_1_dataset_1 ~/Desktop" )

# move one dataset with Kris' survival times
#system( "scp mmathur@sherlock:/scratch/PI/manishad/PCORI/bigSim/covariatesS6/SURV_2015-04-01_job_10_worker_1_dataset_1 ~/Desktop" )


# clean up the directory
rm /share/PI/manishad/genCov/datasets/*
  rm /share/PI/manishad/genCov/sbatch_files/rm*
  
  
rm /share/PI/manishad/genCov/sbatch_files/*


############################ LOOK AT EXAMPLE DATASET ############################ 

# look at example dataset
d = read.csv("~/Desktop/2015-08-25_job_6_worker_1_dataset_1")

dim(d)

# number subjects
length(unique(d$id))

# look at vln
hist(d$log_vln)
summary(d$log_vln)

