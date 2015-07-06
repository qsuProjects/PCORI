#!/bin/bash
#################
#set partition 
#SBATCH -p manishad
#set a job name  
#SBATCH --job-name=na0001
#################  
#a file for job output, you can check job progress
#SBATCH --output=na0001.out
#################
# a file for errors from the job
#SBATCH --error=na0001.err
#################
#time you think you need; default is one hour
#SBATCH --time=0:30:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB
#SBATCH --mem=16000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=mmathur@stanford.edu
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks=1
#################
#SBATCH --cpus-per-task=8
#now run normal batch commands

ml load R

srun R -f /share/PI/manishad/naEst/r/run.R --args /share/PI/manishad/naEst/test-data.csv test_job

