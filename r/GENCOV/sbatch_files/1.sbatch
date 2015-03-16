#!/bin/bash
#################
#set a job name  
#SBATCH --job-name=job_1
#################  
#a file for job output, you can check job progress
#SBATCH --output=rm_1.out
#################
# a file for errors from the job
#SBATCH --error=rm_1.err
#################
#time you think you need; default is one hour
#SBATCH --time=48:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB
#SBATCH --mem=32000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=mmathur@stanford.edu
#################
#task to run per node; each node has 16 cores
# CHANGE TO 16 TO USE ALL CORES
#SBATCH --ntasks=2
#################
#SBATCH --cpus-per-task=1
#now run normal batch commands

ml load R

# arguments: number of  workers, job title prefix, n, obs, n.Reps, n.Drugs
srun R -f /share/PI/manishad/sim_28_sher/simulate.R --args 2 job_1 2000 200 1 15 /share/PI/manishad/sim_28_sher
