
########### NA Estimator ############

# this pushes the files from local GitHub repo to Sherlock
# note: will give lots of permissions errors, but it does work

# push just R files
system( "scp -r ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/for-sherlock/r/* mmathur@sherlock:/share/PI/manishad/naEst/r" )


# push all files (SLOW)
system( "scp -r ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/for-sherlock/* mmathur@sherlock:/share/PI/manishad/naEst" )


# run test file
sbatch /share/PI/manishad/naEst/test-sbatch.sbatch -p manishad

# pull them back locally
scp mmathur@sherlock:/share/PI/manishad/naEst/output/* ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/unstitched

# delete all output
rm /share/PI/manishad/naEst/output/*
rm /share/PI/manishad/naEst/sbatch/slurm*
rm /share/PI/manishad/naEst/slurm*
 
rm /share/PI/manishad/naEst/na00*
  
# remove unstitched and stitched local files
rm ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/stitched/*
rm ~/Dropbox/QSU/Mathur/PCORI/PCORI_git/r/NAest/from-sherlock/unstitched/*

  
cd /share/PI/manishad/naEst/sbatch