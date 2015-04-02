# this pushes the files from local GitHub repo to Sherlock

# push all the individual files
system( "scp -r ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV/* mmathur@sherlock:/share/PI/manishad/genCov" )

# push the sbatch folder
#system( "scp ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/GENCOV/sbatch_files mmathur@sherlock:/share/PI/manishad/genCov" )

# see the sbatch files
# cd /share/PI/manishad/genCov/sbatch_files

# run them
# sbatch /share/PI/manishad/genCov/sbatch_files/manual_test.sbatch -p manishad
# sbatch 2.sbatch -p manishad


# see the datasets
# cd /share/PI/manishad/genCov/datasets

# move one dataset to local machine
system( "scp mmathur@sherlock:/share/PI/manishad/genCov/datasets/2015-04-01_job_10_worker_1_dataset_1 ~/Desktop" )




############################ LOOK AT EXAMPLE DATASET ############################ 

# look at example dataset
d = read.csv("~/Desktop/2015-04-01_job_10_worker_1_dataset_1")

dim(d)

# number subjects
length(unique(d$id))

# look at vln
hist(d$log_vln)
summary(d$log_vln)

